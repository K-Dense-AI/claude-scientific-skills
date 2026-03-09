#!/usr/bin/env python3
"""
WebSocket 실시간 상태 서버
team-orchestrator 팀 상태를 실시간으로 브로드캐스트
"""
import asyncio
import json
import sys
import time
import logging
from pathlib import Path
from typing import Set
from datetime import datetime

from fastapi import FastAPI, WebSocket, WebSocketDisconnect, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse, Response
import uvicorn

# sys.path에 부모 디렉토리 추가 (state, team_manager 임포트)
sys.path.insert(0, str(Path(__file__).parent))

import state as st
import team_manager
import event_bus as eb

# 구조화된 로깅 설정
logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

app = FastAPI(title="Team Orchestrator WebSocket")

# HTML 대시보드 로드
_DASHBOARD_HTML = ""
_dashboard_path = Path(__file__).parent / "notion_dashboard.html"
if _dashboard_path.exists():
    _DASHBOARD_HTML = _dashboard_path.read_text(encoding="utf-8")

# CORS 설정 (Notion 페이지 등에서 WebSocket 연결 가능)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# 연결된 클라이언트 관리
class ConnectionManager:
    def __init__(self):
        self.active_connections: Set[WebSocket] = set()
        self.heartbeat_tasks: dict = {}  # websocket -> asyncio.Task mapping

    async def connect(self, websocket: WebSocket):
        await websocket.accept()
        self.active_connections.add(websocket)
        logger.info(
            f"WebSocket_connect",
            extra={
                "timestamp": datetime.utcnow().isoformat(),
                "total_connections": len(self.active_connections)
            }
        )

    def disconnect(self, websocket: WebSocket):
        self.active_connections.discard(websocket)

        # heartbeat 태스크 취소
        if websocket in self.heartbeat_tasks:
            task = self.heartbeat_tasks.pop(websocket)
            if not task.done():
                task.cancel()

        logger.info(
            f"WebSocket_disconnect",
            extra={
                "timestamp": datetime.utcnow().isoformat(),
                "remaining_connections": len(self.active_connections)
            }
        )

    async def _send_heartbeat(self, websocket: WebSocket):
        """
        30초마다 ping 프레임 전송 (WebSocket keep-alive).
        클라이언트가 응답하지 않으면 자동으로 연결 종료.
        """
        heartbeat_interval = 30  # 초
        attempt = 0

        try:
            while True:
                await asyncio.sleep(heartbeat_interval)
                try:
                    await websocket.send_text(
                        json.dumps({
                            "type": "heartbeat",
                            "timestamp": time.time()
                        })
                    )
                    logger.debug(
                        f"Heartbeat_sent",
                        extra={
                            "websocket_id": id(websocket),
                            "timestamp": datetime.utcnow().isoformat()
                        }
                    )
                    attempt = 0  # 성공 시 리셋
                except Exception as e:
                    attempt += 1
                    logger.warning(
                        f"Heartbeat_failed",
                        extra={
                            "websocket_id": id(websocket),
                            "attempt": attempt,
                            "error_type": type(e).__name__,
                            "error_details": str(e),
                            "timestamp": datetime.utcnow().isoformat()
                        }
                    )
                    if attempt >= 3:
                        self.disconnect(websocket)
                        break
        except asyncio.CancelledError:
            logger.debug(f"Heartbeat_cancelled for {id(websocket)}")
        except Exception as e:
            logger.error(
                f"Heartbeat_error",
                extra={
                    "websocket_id": id(websocket),
                    "error_type": type(e).__name__,
                    "error_details": str(e),
                    "timestamp": datetime.utcnow().isoformat()
                }
            )

    async def _exponential_backoff(self, attempt: int) -> float:
        """
        지수 백오프 계산 (1s -> 2s -> 4s -> ... -> 30s).
        네트워크 불안정 시 재시도 대기 시간.
        """
        base_delay = 1.0
        max_delay = 30.0
        delay = min(base_delay * (2 ** attempt), max_delay)
        jitter = (delay * 0.1) * (asyncio.get_event_loop().time() % 1.0)
        return delay + jitter

    async def broadcast(self, data: dict):
        """모든 연결된 클라이언트에게 상태 브로드캐스트"""
        disconnected = set()
        for connection in self.active_connections:
            try:
                await connection.send_text(
                    json.dumps(data, ensure_ascii=False)
                )
            except Exception as e:
                logger.warning(
                    f"Broadcast_failed",
                    extra={
                        "websocket_id": id(connection),
                        "error_type": type(e).__name__,
                        "error_details": str(e),
                        "timestamp": datetime.utcnow().isoformat()
                    }
                )
                disconnected.add(connection)

        # 연결이 끊긴 클라이언트 제거
        for conn in disconnected:
            self.disconnect(conn)


manager = ConnectionManager()


def _build_status_message(project_id: str = None) -> dict:
    """현재 팀 상태를 JSON으로 빌드"""
    teams = team_manager.get_status_all(project_id)
    current_time = time.time()

    # 상태별 카운트
    status_counts = {"running": 0, "pending": 0, "done": 0, "failed": 0, "paused": 0}
    for t in teams:
        status_counts[t["status"]] = status_counts.get(t["status"], 0) + 1

    # 팀 상세 정보
    teams_data = []
    for t in teams:
        team = st.get_team(t["team_id"])
        elapsed = 0.0
        if t["status"] == "running":
            if team and team.started_at > 0:
                elapsed = current_time - team.started_at
        else:
            if team:
                elapsed = team.elapsed_seconds or 0.0

        total_tokens = ((team.input_tokens or 0) + (team.output_tokens or 0)) if team else 0

        teams_data.append({
            "team_id": t["team_id"],
            "name": t["name"],
            "type": t["type"],
            "status": t["status"],
            "elapsed_seconds": round(elapsed, 1),
            "task": "",  # 손상된 UTF-8 데이터 제외
            "model": team.model_used if team and team.model_used else "",
            "input_tokens": team.input_tokens if team else 0,
            "output_tokens": team.output_tokens if team else 0,
            "total_tokens": total_tokens,
            "total_cost_usd": round(team.total_cost_usd or 0, 4) if team else 0,
            "tok_per_sec": round(total_tokens / max(elapsed, 1), 1) if team else 0,
        })

    return {
        "timestamp": time.time(),
        "project_id": project_id,
        "status_counts": status_counts,
        "teams": teams_data,
    }


@app.websocket("/ws")
async def websocket_endpoint(
    websocket: WebSocket,
    project_id: str = Query(None, description="선택: 특정 프로젝트 ID만 필터링"),
):
    """
    WebSocket 엔드포인트 (UTF-8; charset=utf-8)

    쿼리 파라미터:
    - project_id (선택): 특정 프로젝트만 모니터링

    예:
    - ws://localhost:8765/ws
    - ws://localhost:8765/ws?project_id=abc123

    기능:
    - 초기 상태 전송 직후 heartbeat 태스크 시작
    - 100ms 간격 상태 브로드캐스트
    - 30초마다 heartbeat ping 전송
    - 네트워크 오류 시 자동 재연결 지원 (클라이언트 측)
    """
    await manager.connect(websocket)
    logger.info(
        f"WebSocket_endpoint_connected",
        extra={
            "project_id": project_id,
            "websocket_id": id(websocket),
            "timestamp": datetime.utcnow().isoformat()
        }
    )

    try:
        # 초기 상태 전송
        message = _build_status_message(project_id)
        await websocket.send_text(json.dumps(message, ensure_ascii=False))

        # Heartbeat 태스크 생성 및 저장
        heartbeat_task = asyncio.create_task(manager._send_heartbeat(websocket))
        manager.heartbeat_tasks[websocket] = heartbeat_task

        # 100ms 간격으로 상태 브로드캐스트
        last_broadcast = 0.0
        last_event_id = eb.get_latest_id()
        while True:
            await asyncio.sleep(0.1)  # CPU 부하 방지

            # 새 이벤트 감지 시 즉시 브로드캐스트
            new_events = eb.get_events_since(last_event_id)
            if new_events:
                last_event_id = new_events[-1]["id"]
                last_broadcast = 0.0  # 즉시 상태 전송 트리거

            current_time = time.time()
            if current_time - last_broadcast >= 0.1:  # 100ms마다
                message = _build_status_message(project_id)

                # 해당 클라이언트에만 전송 (직접 연결이므로)
                try:
                    await websocket.send_text(json.dumps(message, ensure_ascii=False))
                    last_broadcast = current_time
                except Exception as e:
                    logger.error(
                        f"WebSocket_send_failed",
                        extra={
                            "project_id": project_id,
                            "websocket_id": id(websocket),
                            "error_type": type(e).__name__,
                            "error_details": str(e),
                            "timestamp": datetime.utcnow().isoformat()
                        }
                    )
                    break

    except WebSocketDisconnect:
        manager.disconnect(websocket)
        logger.info(
            f"WebSocket_disconnected",
            extra={
                "project_id": project_id,
                "websocket_id": id(websocket),
                "reason": "client_disconnect",
                "timestamp": datetime.utcnow().isoformat()
            }
        )
    except Exception as e:
        manager.disconnect(websocket)
        logger.error(
            f"WebSocket_error",
            extra={
                "project_id": project_id,
                "websocket_id": id(websocket),
                "error_type": type(e).__name__,
                "error_details": str(e),
                "timestamp": datetime.utcnow().isoformat()
            }
        )


@app.get("/", response_class=HTMLResponse)
async def dashboard():
    """대시보드 HTML 페이지 (Notion Embed / 브라우저 직접 접속)"""
    if _dashboard_path.exists():
        return HTMLResponse(content=_dashboard_path.read_text(encoding="utf-8"))
    return HTMLResponse(content="<h1>Dashboard HTML not found</h1>", status_code=404)


@app.get("/health")
async def health():
    """헬스 체크 (UTF-8; charset=utf-8)"""
    response_data = {
        "status": "running",
        "clients": len(manager.active_connections),
        "timestamp": time.time(),
    }
    return Response(
        content=json.dumps(response_data, ensure_ascii=False),
        media_type="application/json; charset=utf-8"
    )


def _clean_surrogates(text: str) -> str:
    """Orphaned surrogate 문자 제거"""
    if not text:
        return ""
    try:
        # Python 3.1+의 'surrogatepass'를 사용한 올바른 인코딩
        # UTF-8로 인코딩하되, surrogate는 무시
        return text.encode('utf-8', errors='replace').decode('utf-8')
    except:
        return ''.join(c for c in text if ord(c) < 0xdc00 or ord(c) > 0xdfff)

@app.get("/status")
async def status(project_id: str = Query(None)):
    """HTTP 폴링용 상태 조회 (WebSocket이 없을 때)"""
    data = _build_status_message(project_id)

    # 모든 task 필드에 대해 surrogate 정제
    for team in data.get("teams", []):
        if team.get("task"):
            team["task"] = _clean_surrogates(team["task"])

    json_str = json.dumps(data, ensure_ascii=False)
    return Response(content=json_str, media_type="application/json; charset=utf-8")


if __name__ == "__main__":
    port = 8765
    print(f"WebSocket 서버 시작: ws://localhost:{port}/ws", flush=True)
    print(f"헬스 체크: http://localhost:{port}/health", flush=True)
    print(f"HTTP 폴링: http://localhost:{port}/status", flush=True)

    st.init_db()

    uvicorn.run(
        app,
        host="127.0.0.1",
        port=port,
        log_level="warning",
    )
