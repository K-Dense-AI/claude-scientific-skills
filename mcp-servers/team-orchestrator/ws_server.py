#!/usr/bin/env python3
"""
WebSocket 실시간 상태 서버
team-orchestrator 팀 상태를 실시간으로 브로드캐스트
"""
import asyncio
import json
import sys
import time
from pathlib import Path
from typing import Set

from fastapi import FastAPI, WebSocket, WebSocketDisconnect, Query
from fastapi.middleware.cors import CORSMiddleware
import uvicorn

# sys.path에 부모 디렉토리 추가 (state, team_manager 임포트)
sys.path.insert(0, str(Path(__file__).parent))

import state as st
import team_manager

app = FastAPI(title="Team Orchestrator WebSocket")

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

    async def connect(self, websocket: WebSocket):
        await websocket.accept()
        self.active_connections.add(websocket)

    def disconnect(self, websocket: WebSocket):
        self.active_connections.discard(websocket)

    async def broadcast(self, data: dict):
        """모든 연결된 클라이언트에게 상태 브로드캐스트"""
        disconnected = set()
        for connection in self.active_connections:
            try:
                await connection.send_json(data)
            except Exception:
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
        elapsed = 0.0
        if t["status"] == "running":
            team = st.get_team(t["team_id"])
            if team and team.started_at > 0:
                elapsed = current_time - team.started_at
        else:
            team = st.get_team(t["team_id"])
            if team:
                elapsed = team.elapsed_seconds or 0.0

        teams_data.append({
            "team_id": t["team_id"],
            "name": t["name"],
            "type": t["type"],
            "status": t["status"],
            "elapsed_seconds": round(elapsed, 1),
            "task": t["task"][:100] if t["task"] else "",
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
    WebSocket 엔드포인트

    쿼리 파라미터:
    - project_id (선택): 특정 프로젝트만 모니터링

    예:
    - ws://localhost:8765/ws
    - ws://localhost:8765/ws?project_id=abc123
    """
    await manager.connect(websocket)
    print(f"[WebSocket] 클라이언트 연결 (project_id={project_id})", flush=True)

    try:
        # 초기 상태 전송
        message = _build_status_message(project_id)
        await websocket.send_json(message)

        # 100ms 간격으로 상태 브로드캐스트
        last_broadcast = 0.0
        while True:
            await asyncio.sleep(0.1)  # CPU 부하 방지

            current_time = time.time()
            if current_time - last_broadcast >= 0.1:  # 100ms마다
                message = _build_status_message(project_id)

                # 해당 클라이언트에만 전송 (직접 연결이므로)
                try:
                    await websocket.send_json(message)
                    last_broadcast = current_time
                except Exception as e:
                    print(f"[WebSocket] 전송 실패: {e}", flush=True)
                    break

    except WebSocketDisconnect:
        manager.disconnect(websocket)
        print(f"[WebSocket] 클라이언트 연결 해제 (project_id={project_id})", flush=True)
    except Exception as e:
        manager.disconnect(websocket)
        print(f"[WebSocket] 에러: {e}", flush=True)


@app.get("/health")
async def health():
    """헬스 체크"""
    return {
        "status": "running",
        "clients": len(manager.active_connections),
        "timestamp": time.time(),
    }


@app.get("/status")
async def status(project_id: str = Query(None)):
    """HTTP 폴링용 상태 조회 (WebSocket이 없을 때)"""
    return _build_status_message(project_id)


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
