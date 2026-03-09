#!/usr/bin/env python3
"""
Team Orchestrator MCP Server
CEO 역할의 Claude Code에서 여러 팀을 실시간으로 관리하는 MCP 서버
"""
# [DIAGNOSTIC] 가장 첫 번째 로그 — import 전에 실행되는지 확인
import os as _os_diag, time as _time_diag
_DIAG_LOG = _os_diag.path.join(_os_diag.path.dirname(_os_diag.path.abspath(__file__)), "startup_diag.log")
with open(_DIAG_LOG, "a", encoding="utf-8") as _f:
    _f.write(f"[{_time_diag.strftime('%H:%M:%S')}] server.py 시작 | PID={_os_diag.getpid()} | cwd={_os_diag.getcwd()} | USERPROFILE={_os_diag.environ.get('USERPROFILE','없음')}\n")
del _os_diag, _time_diag, _DIAG_LOG, _f

import sys
from pathlib import Path
# 서버 디렉토리를 sys.path에 추가 (VSCode 확장에서 cwd가 달라도 import 가능)
sys.path.insert(0, str(Path(__file__).parent))

import asyncio
import atexit
import json
import logging
import os
import subprocess
import threading
import time
import traceback
import uuid

# 파일 로거 (MCP stdio 오염 없이 Claude Code 연결 시 디버그 가능)
_LOG_FILE = Path(__file__).parent / "server_debug.log"
logging.basicConfig(
    filename=str(_LOG_FILE),
    level=logging.DEBUG,
    format="%(asctime)s %(levelname)s %(message)s",
    encoding="utf-8",
)
_logger = logging.getLogger("team-orchestrator")
_logger.info("=== server.py 시작 ===")

def _global_exception_hook(exc_type, exc_value, exc_tb):
    _logger.critical("Uncaught exception: %s", "".join(traceback.format_exception(exc_type, exc_value, exc_tb)))
    sys.__excepthook__(exc_type, exc_value, exc_tb)

sys.excepthook = _global_exception_hook

# FastMCP 사용
from mcp.server.fastmcp import FastMCP

import state as st
import task_router
import team_manager
import shutdown_handler
import notion_logger as nl
import telegram_notify as tg

SKILLS_BASE = Path.home() / "claude-scientific-skills" / "scientific-skills"

# ── PID 기록 (디버깅용) ───────────────────────────────────────────────────────
_PID_FILE = Path.home() / ".claude" / "team-orchestrator.pid"
_TG_LISTENER_PID_FILE = Path(__file__).parent / "tg_listener.pid"
_PID_FILE.write_text(str(os.getpid()))


def _cleanup_on_exit():
    """MCP 서버 종료 시 자식 프로세스(ws_server, tg_listener) 정리"""
    _PID_FILE.unlink(missing_ok=True)
    _TG_LISTENER_PID_FILE.unlink(missing_ok=True)
    for proc in _child_procs:
        try:
            proc.terminate()
        except Exception:
            pass


atexit.register(_cleanup_on_exit)


def _inject_skill_hints(task: str, skills: list[str]) -> str:
    """태스크 텍스트에 스킬 SKILL.md 경로 힌트를 추가"""
    if not skills:
        return task
    lines = [task, "", "[Skill hint -- use Read tool to check the relevant SKILL.md for API/library usage]"]
    for skill in skills:
        path = SKILLS_BASE / skill / "SKILL.md"
        lines.append(f"- {skill}: {path}")
    return "\n".join(lines)

# DB 초기화
st.init_db()

# 서버 시작 시 2시간 이상 running/pending 상태인 stale 에이전트 자동 정리
_stale = st.cleanup_stale_agents(max_age_hours=2.0)
if _stale:
    print(f"[startup] stale 에이전트 {len(_stale)}개 정리: {_stale}", file=sys.stderr, flush=True)

# Notion 자동 동기화 비활성화 — 실시간 대시보드(Cloudflare Tunnel)로 대체
# nl.start_auto_sync(interval=30.0)
# print("[startup] Notion 자동 동기화 시작 (30s 간격)", file=sys.stderr, flush=True)

# ── 대시보드 (WS서버 + cloudflared) 백그라운드 자동 시작 ─────────────────────
_DASHBOARD_PUBLIC_URL = None  # 외부에서 조회 가능하도록 모듈 레벨 저장
_CLOUDFLARED_PATH = str(Path.home() / "cloudflared.exe")

_child_procs: list[subprocess.Popen] = []  # 자식 프로세스 추적 (atexit 정리용)


def _is_port_in_use(port: int) -> bool:
    """포트가 이미 사용 중인지 확인"""
    import socket
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(("127.0.0.1", port)) == 0


def _start_dashboard_background():
    """WS 서버 + Cloudflare Tunnel을 백그라운드로 시작"""
    global _DASHBOARD_PUBLIC_URL
    script_dir = str(Path(__file__).parent)

    # 0) 환경변수에서 고정 URL 우선 사용 (ngrok 유료 고정 도메인 등)
    static_url = os.environ.get("NGROK_STATIC_URL", "").strip()

    # 1) WebSocket 서버 (ws_server.py) 서브프로세스 — 이미 실행 중이면 스킵
    if _is_port_in_use(8765):
        print("[startup] WebSocket 서버 이미 실행 중 (포트 8765) — 스킵", file=sys.stderr, flush=True)
        _DASHBOARD_PUBLIC_URL = static_url or "http://localhost:8765"
    else:
        try:
            _cflags = 0x08000000 if sys.platform == "win32" else 0  # CREATE_NO_WINDOW
            proc = subprocess.Popen(
                [sys.executable, os.path.join(script_dir, "ws_server.py")],
                cwd=script_dir,
                stdin=subprocess.DEVNULL,   # MCP stdin 파이프 상속 차단
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                creationflags=_cflags,
            )
            _child_procs.append(proc)
            print("[startup] WebSocket 대시보드 서버 시작 (포트 8765)", file=sys.stderr, flush=True)
            # 로컬 호스트는 항상 사용 가능
            _DASHBOARD_PUBLIC_URL = static_url or "http://localhost:8765"
        except Exception as e:
            print(f"[startup] WebSocket 서버 시작 실패: {e}", file=sys.stderr, flush=True)
            return

    # 고정 URL이 설정되어 있으면 Cloudflare 터널 스킵
    if static_url:
        print(f"[startup] 고정 URL 사용: {static_url} — Cloudflare 터널 스킵", file=sys.stderr, flush=True)
        _update_notion_embed(static_url)
        return

    # 2) Cloudflare Tunnel (무료, 인터스티셜 없음) — 선택사항
    time.sleep(2)
    if not os.path.exists(_CLOUDFLARED_PATH):
        print(f"[startup] cloudflared 미설치 ({_CLOUDFLARED_PATH}) — 로컬 포트만 사용 (http://localhost:8765)", file=sys.stderr, flush=True)
        return
    try:
        import re
        _cflags2 = 0x08000000 if sys.platform == "win32" else 0  # CREATE_NO_WINDOW
        proc = subprocess.Popen(
            [_CLOUDFLARED_PATH, "tunnel", "--url", "http://localhost:8765"],
            stdin=subprocess.DEVNULL,    # MCP stdin 파이프 상속 차단
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            creationflags=_cflags2,
        )
        # cloudflared 출력에서 URL 추출 (최대 30초 대기)
        for _ in range(60):
            line = proc.stdout.readline()
            if not line:
                break
            m = re.search(r"(https://[a-z0-9-]+\.trycloudflare\.com)", line)
            if m:
                _DASHBOARD_PUBLIC_URL = m.group(1)
                print(f"[startup] Cloudflare 터널 활성: {_DASHBOARD_PUBLIC_URL}", file=sys.stderr, flush=True)
                _update_notion_embed(_DASHBOARD_PUBLIC_URL)
                break
        if _DASHBOARD_PUBLIC_URL == "http://localhost:8765":
            print("[startup] Cloudflare 터널 URL 추출 실패 — 로컬 포트만 사용", file=sys.stderr, flush=True)
    except Exception as e:
        print(f"[startup] cloudflared 시작 실패: {e} — 로컬 포트만 사용", file=sys.stderr, flush=True)


# Notion 대시보드 페이지 Embed URL 자동 갱신
_NOTION_DASHBOARD_PAGE_ID = "31cf91aca96f8107819bef1d4b05900a"

def _update_notion_embed(new_url: str):
    """Notion 대시보드 페이지의 embed 블록 URL을 자동 갱신"""
    try:
        import requests as _req
        secrets_path = Path.home() / ".claude" / "secrets.json"
        notion_token = json.loads(secrets_path.read_text()).get("NOTION_TOKEN", "")
        if not notion_token:
            print("[startup] NOTION_TOKEN 없음 — Notion embed 갱신 생략", file=sys.stderr, flush=True)
            return

        headers = {
            "Authorization": f"Bearer {notion_token}",
            "Notion-Version": "2022-06-28",
        }

        # 1) 페이지의 블록 children 조회
        resp = _req.get(
            f"https://api.notion.com/v1/blocks/{_NOTION_DASHBOARD_PAGE_ID}/children?page_size=100",
            headers=headers,
        )
        if resp.status_code != 200:
            print(f"[startup] Notion 블록 조회 실패: {resp.status_code}", file=sys.stderr, flush=True)
            return

        # 2) embed 블록 찾기
        for block in resp.json().get("results", []):
            if block.get("type") == "embed":
                block_id = block["id"]
                # 3) embed URL 업데이트
                update_resp = _req.patch(
                    f"https://api.notion.com/v1/blocks/{block_id}",
                    headers={**headers, "Content-Type": "application/json"},
                    json={"embed": {"url": new_url}},
                )
                if update_resp.status_code == 200:
                    print(f"[startup] Notion embed URL 자동 갱신 완료: {new_url}", file=sys.stderr, flush=True)
                else:
                    print(f"[startup] Notion embed 갱신 실패: {update_resp.status_code}", file=sys.stderr, flush=True)
                return

        print("[startup] Notion 페이지에 embed 블록 없음 — 갱신 생략", file=sys.stderr, flush=True)
    except Exception as e:
        print(f"[startup] Notion embed 갱신 오류: {e}", file=sys.stderr, flush=True)

def _is_tg_listener_running() -> bool:
    """tg_listener.py가 이미 실행 중인지 PID 파일로 확인"""
    if not _TG_LISTENER_PID_FILE.exists():
        return False
    try:
        pid = int(_TG_LISTENER_PID_FILE.read_text().strip())
        # Windows: 프로세스 존재 여부 확인
        import ctypes
        handle = ctypes.windll.kernel32.OpenProcess(0x1000, False, pid)  # PROCESS_QUERY_LIMITED_INFORMATION
        if handle:
            ctypes.windll.kernel32.CloseHandle(handle)
            return True
    except Exception:
        pass
    _TG_LISTENER_PID_FILE.unlink(missing_ok=True)
    return False


# ── Telegram 리스너 자동 시작 (답장 수신용) ──────────────────────────────────
def _start_tg_listener():
    """tg_listener.py를 서브프로세스로 시작 — 이미 실행 중이면 스킵"""
    if _is_tg_listener_running():
        print("[startup] Telegram 리스너 이미 실행 중 — 스킵", file=sys.stderr, flush=True)
        return
    try:
        script = str(Path(__file__).parent / "tg_listener.py")
        _cflags3 = 0x08000000 if sys.platform == "win32" else 0  # CREATE_NO_WINDOW
        proc = subprocess.Popen(
            [sys.executable, script],
            cwd=str(Path(__file__).parent),
            stdin=subprocess.DEVNULL,    # MCP stdin 파이프 상속 차단
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            creationflags=_cflags3,
        )
        _TG_LISTENER_PID_FILE.write_text(str(proc.pid))
        _child_procs.append(proc)
        print("[startup] Telegram 리스너 시작", file=sys.stderr, flush=True)
    except Exception as e:
        print(f"[startup] Telegram 리스너 시작 실패: {e}", file=sys.stderr, flush=True)

mcp = FastMCP("team-orchestrator")


# ── 프로젝트 시작 ────────────────────────────────────────────────────────────

@mcp.tool()
def start_project(description: str, feedback_timeout: int = 10) -> str:
    """
    프로젝트를 분석해 팀을 자동 배정하고 병렬 실행을 시작합니다.

    Args:
        description: 프로젝트 설명 (무엇을 해야 하는지 자연어로)
        feedback_timeout: 팀 완료 후 의존성 팀 자동 시작 전 대기 시간(초). 기본 10초.
                          0으로 설정 시 즉시 시작.

    Returns:
        생성된 팀 목록과 project_id
    """
    project_id = str(uuid.uuid4())[:8]
    st.create_project(project_id, description)
    nl.log_project(project_id, description)
    team_manager.set_project_feedback_timeout(project_id, feedback_timeout)

    # Claude로 태스크 분석 및 팀 배정
    plan = task_router.analyze_and_split(description)

    teams_created = []
    for team_def in plan.get("teams", []):
        team_id = f"{project_id}-{len(teams_created)+1}"
        name = team_def["name"]
        team_type = team_def["type"]
        task = _inject_skill_hints(team_def["task"], team_def.get("skills", []))
        depends_on = team_def.get("depends_on", [])

        # 의존성 없는 팀은 즉시 시작, 의존성 있는 팀은 pending 상태로 등록
        if not depends_on:
            team_manager.start_team(team_id, name, team_type, task, project_id,
                                    feedback_timeout=feedback_timeout)
        else:
            # pending 상태로만 등록 (의존성 해결 시 자동 시작됨)
            team = st.Team(
                team_id=team_id, name=name, team_type=team_type,
                status="pending", task=task, project_id=project_id,
                created_at=time.time(), depends_on=json.dumps(depends_on)
            )
            st.upsert_team(team)

        teams_created.append({
            "team_id": team_id,
            "name": name,
            "type": team_type,
            "status": "running" if not depends_on else "pending (depends on: " + ", ".join(depends_on) + ")",
        })

    result = {
        "project_id": project_id,
        "project_summary": plan.get("project_summary", ""),
        "teams": teams_created,
        "dashboard_url": _DASHBOARD_PUBLIC_URL or "N/A",
        "message": (
            f"{len(teams_created)} agents created. Use get_status() to check progress.\n"
            f"[CEO 지시] 모든 팀리드 완료 보고 수신 후 → shutdown(project_id='{project_id}') 호출하여 Notion 보고 에이전트를 실행할 것."
        )
    }
    return json.dumps(result, ensure_ascii=False, indent=2)


# ── CEO 피드백 ────────────────────────────────────────────────────────────────

@mcp.tool()
def send_feedback(team_id: str, message: str) -> str:
    """
    실행 중인 팀에 CEO 추가 지시나 피드백을 전송합니다.

    Args:
        team_id: 팀 ID (get_status()에서 확인)
        message: CEO의 추가 지시 또는 피드백

    Returns:
        전송 결과
    """
    try:
        team_manager.send_feedback(team_id, message)
        team = st.get_team(team_id)
        return f"[Sent] Feedback delivered to agent '{team.name}'. Will be applied on next response turn."
    except ValueError as e:
        return f"[Error] {e}"


# ── 상태 조회 ─────────────────────────────────────────────────────────────────

@mcp.tool()
def get_status(project_id: str = None) -> str:
    """
    전체 팀 또는 특정 프로젝트의 팀 상태를 조회합니다.

    Args:
        project_id: (선택) 특정 프로젝트 ID. 없으면 전체 조회.

    Returns:
        팀별 상태 요약 + 타이머
    """
    teams = team_manager.get_status_all(project_id)
    current_time = time.time()

    # 자동 보고 트리거 (모든 팀 완료 시)
    if project_id:
        _trigger_auto_reporting(project_id)
    else:
        # 전체 조회 시 각 프로젝트별로 체크
        all_projects = set(t["team_id"].split("-")[0] for t in teams)
        for proj_id in all_projects:
            _trigger_auto_reporting(proj_id)

    # 상태별 카운트
    status_counts = {"running": 0, "pending": 0, "done": 0, "failed": 0, "paused": 0}
    for t in teams:
        status_counts[t["status"]] = status_counts.get(t["status"], 0) + 1

    # 읽지 않은 완료 알림 먼저 표시
    notifications = st.get_pending_notifications()
    notif_header = ""
    if notifications:
        notif_lines = ["## 🔔 완료 알림\n"]
        for n in notifications:
            notif_lines.append(f"- {n['message']}")
        notif_lines.append("")
        notif_header = "\n".join(notif_lines) + "\n"

    if not teams:
        return notif_header + "No active agents." if notif_header else "No active agents."

    # 요약 헤더
    summary_parts = []
    if status_counts["running"] > 0:
        summary_parts.append(f"🔄 running: {status_counts['running']}")
    if status_counts["pending"] > 0:
        summary_parts.append(f"⏳ pending: {status_counts['pending']}")
    if status_counts["done"] > 0:
        summary_parts.append(f"✅ done: {status_counts['done']}")
    if status_counts["failed"] > 0:
        summary_parts.append(f"❌ 실패: {status_counts['failed']}")
    summary = " | ".join(summary_parts) if summary_parts else "상태 없음"

    lines = [f"## Agent Status — {summary}\n"]
    done_leads = []
    for t in teams:
        status_icon = {
            "running": "🔄", "done": "✅", "failed": "❌",
            "pending": "⏳", "paused": "⏸️"
        }.get(t["status"], "❓")

        # 타이머: running 상태이면 경과 시간 표시
        timer_info = ""
        elapsed = 0.0
        if t["status"] == "running":
            team = st.get_team(t["team_id"])
            if team and team.started_at > 0:
                elapsed = current_time - team.started_at
                minutes = int(elapsed // 60)
                seconds = int(elapsed % 60)
                if minutes > 0:
                    timer_info = f" (🕐 {minutes}m {seconds}s 경과)"
                else:
                    timer_info = f" (🕐 {seconds}s 경과)"
                # Notion 실시간 동기화
                nl.update_agent_status(t["team_id"], t["status"], elapsed)
        else:
            # running 아닌 상태도 최종 경과 시간 동기화
            team = st.get_team(t["team_id"])
            if team:
                elapsed = team.elapsed_seconds or 0.0
                nl.update_agent_status(t["team_id"], t["status"], elapsed)

        lines.append(
            f"{status_icon} **{t['name']}** (ID: `{t['team_id']}`){timer_info}\n"
            f"   role: {t['type']} | status: {t['status']}\n"
            f"   task: {t['task']}\n"
        )

        # 리드 에이전트(team_id에 -w 없음)가 done이면 output 자동 포함
        is_lead = "-w" not in t["team_id"]
        if t["status"] == "done" and is_lead:
            team = st.get_team(t["team_id"])
            if team and team.output.strip():
                done_leads.append((t["name"], t["team_id"], team.output.strip()))

    result = notif_header + "\n".join(lines)
    for name, team_id, output in done_leads:
        notion_hint = "\n\n---\n> 📋 To update Notion, please add this result to the Notion page."
        result += f"\n\n---\n## 📋 {name} — 최종 결과\n\n{output}{notion_hint}"

    return result


# ── 팀 산출물 조회 ────────────────────────────────────────────────────────────

@mcp.tool()
def get_output(team_id: str) -> str:
    """
    팀의 현재까지 산출물을 가져옵니다.

    Args:
        team_id: 팀 ID

    Returns:
        팀 산출물 전체 텍스트
    """
    team = st.get_team(team_id)
    if not team:
        return f"Agent {team_id} not found."

    output = team.output.strip()
    if not output:
        return f"Agent '{team.name}' has no output yet. (status: {team.status})"

    return f"## {team.name} Output\n\n{output}"


# ── 팀 일시 정지 ──────────────────────────────────────────────────────────────

@mcp.tool()
def pause_team(team_id: str) -> str:
    """
    실행 중인 팀을 일시 정지합니다.

    Args:
        team_id: 팀 ID

    Returns:
        정지 결과
    """
    team = st.get_team(team_id)
    if not team:
        return f"Agent {team_id} not found."
    team_manager.pause_team(team_id)
    return f"Agent '{team.name}' paused."


# ── 종료 및 Notion 동기화 ────────────────────────────────────────────────────

@mcp.tool()
def shutdown(project_id: str = None) -> str:
    """
    모든 팀(또는 특정 프로젝트)을 종료하고 Notion에 스프린트 리포트를 생성합니다.
    - 팀 산출물 수집
    - Notion 리포트 에이전트 실행

    Args:
        project_id: (선택) 특정 프로젝트만 종료. 없으면 전체.

    Returns:
        종료 결과 요약
    """
    teams = st.list_teams(project_id)

    # 실행 중 팀 정지
    for team in teams:
        if team.status == "running":
            team_manager.pause_team(team.team_id)
            st.update_team_status(team.team_id, "done")

    # 산출물 수집
    team_outputs = [
        {
            "name": t.name,
            "type": t.team_type,
            "task": t.task,
            "output": t.output
        }
        for t in teams
    ]

    # Notion 스프린트 리포트 — haiku 특무 에이전트 비동기 실행
    notion_msg = ""
    try:
        outputs_summary = "\n\n".join(
            f"### {t['name']} ({t['type']})\n**태스크:** {t['task'][:200]}\n\n{t['output'][:1000]}"
            for t in team_outputs if t.get("output", "").strip()
        )
        notion_task = (
            f"[Report]\n## Sprint Result Report\n\n{outputs_summary}\n\n"
            f"[태스크 컨텍스트]\n팀명: sprint-shutdown | 에이전트 수: {len(team_outputs)}개"
        )
        nr = _spawn_solo_agent(
            task=notion_task,
            team_type="report-notion",
            name="NotionReport",
        )
        notion_agent_id = nr["agent_id"]
        # Notion 에이전트 완료까지 최대 120초 대기
        import time as _t
        deadline = _t.time() + 120
        while _t.time() < deadline:
            _t.sleep(5)
            agent = st.get_team(notion_agent_id)
            if agent and agent.status in ("done", "failed"):
                break
        agent = st.get_team(notion_agent_id)
        status = agent.status if agent else "unknown"
        notion_msg = f"Notion 리포트 완료 [{status}] (agent_id: {notion_agent_id})"
    except Exception as e:
        notion_msg = f"Notion 리포트 생략: {e}"

    result = (
        f"## Shutdown Complete\n\n"
        f"**Notion:** {notion_msg}"
    )

    # Notion 보고 완료 후 MCP 서버 자동 종료
    def _exit_after_delay():
        import time as _t
        _t.sleep(5)
        print("[shutdown] MCP 서버 자동 종료.", file=sys.stderr, flush=True)
        import os as _os
        _os.exit(0)

    threading.Thread(target=_exit_after_delay, daemon=True).start()
    return result


# ── 팀 수동 시작 (이미 pending인 팀) ─────────────────────────────────────────

@mcp.tool()
def start_team(team_id: str) -> str:
    """
    pending 상태인 팀을 수동으로 시작합니다 (의존성 해결 후).

    Args:
        team_id: 팀 ID

    Returns:
        시작 결과
    """
    team = st.get_team(team_id)
    if not team:
        return f"Agent {team_id} not found."
    if team.status != "pending":
        return f"Agent '{team.name}' is currently {team.status}. Only pending agents can be started."

    team_manager.start_team(team.team_id, team.name, team.team_type, team.task, team.project_id)
    return f"Agent '{team.name}' started."


# ── 팀 간 대화 조회 ────────────────────────────────────────────────────────────

@mcp.tool()
def get_board_messages(project_id: str = None, limit: int = 30) -> str:
    """
    팀 간 공유 메시지 보드의 최근 메시지를 조회합니다.
    에이전트들이 서로 주고받은 질문, 답변, 중간 결과를 확인할 수 있습니다.

    Args:
        project_id: (선택) 특정 프로젝트 팀만 필터. 없으면 전체.
        limit: 최근 N개 메시지 (기본 30)
    """
    import sqlite3
    from pathlib import Path

    board_db = Path(__file__).parent / "board.db"
    if not board_db.exists():
        return "Board DB not found. No agents have started yet."

    conn = sqlite3.connect(str(board_db))
    conn.row_factory = sqlite3.Row

    if project_id:
        # 해당 프로젝트의 팀 ID 목록 가져오기
        orch_db = Path(__file__).parent / "orchestrator.db"
        orch_conn = sqlite3.connect(str(orch_db))
        team_ids = [r[0] for r in orch_conn.execute(
            "SELECT team_id FROM teams WHERE project_id=?", (project_id,)
        ).fetchall()]
        orch_conn.close()
        if not team_ids:
            conn.close()
            return f"No agents found for project {project_id}."
        placeholders = ",".join("?" * len(team_ids))
        rows = conn.execute(f"""
            SELECT from_agent, to_agent, content, timestamp
            FROM board_messages
            WHERE from_agent IN ({placeholders}) OR to_agent IN ({placeholders})
            ORDER BY timestamp DESC LIMIT ?
        """, team_ids + team_ids + [limit]).fetchall()
    else:
        rows = conn.execute("""
            SELECT from_agent, to_agent, content, timestamp
            FROM board_messages
            ORDER BY timestamp DESC LIMIT ?
        """, (limit,)).fetchall()

    conn.close()

    if not rows:
        return "No board messages yet."

    import time as _time
    lines = [f"## Agent Message Board (recent {len(rows)})\n"]
    for r in reversed(rows):
        ts = _time.strftime("%H:%M:%S", _time.localtime(r["timestamp"]))
        target = f"→ {r['to_agent']}" if r["to_agent"] else "→ all"
        lines.append(f"[{ts}] **{r['from_agent']}** {target}")
        lines.append(f"   {r['content'][:200]}")
        lines.append("")
    return "\n".join(lines)


# ── 자동 보고 트리거 (모든 팀 완료 시 Notion 업데이트) ──────────────────────

def _trigger_auto_reporting(project_id: str):
    """
    프로젝트의 모든 팀이 완료되었을 때 Notion 보고 에이전트를 자동 실행.
    get_status()에서 호출됨.
    """
    teams = st.list_teams(project_id)
    if not teams:
        return

    # 모든 팀이 완료/실패 상태인지 확인
    all_done = all(t.status in ("done", "failed") for t in teams)
    if not all_done:
        return

    # 이미 보고 에이전트가 실행 중이거나 완료된지 확인
    report_teams = [t for t in teams if t.team_type == "report-notion"]
    if report_teams:
        return  # 이미 보고 에이전트 존재

    # 산출물 수집
    team_outputs = [
        {
            "name": t.name,
            "type": t.team_type,
            "task": t.task,
            "output": t.output
        }
        for t in teams if t.output.strip()
    ]

    if not team_outputs:
        return

    # Telegram 프로젝트 완료 알림
    try:
        total_cost = sum(t.total_cost_usd or 0.0 for t in teams)
        total_elapsed = max((t.elapsed_seconds or 0.0 for t in teams), default=0.0)
        done_count = sum(1 for t in teams if t.status == "done")
        fail_count = sum(1 for t in teams if t.status == "failed")
        team_names = ", ".join(t.name for t in teams if t.team_type != "report-notion")
        tg.send(
            f"📋 프로젝트 완료: {project_id}\n"
            f"팀: {team_names}\n"
            f"결과: ✅{done_count} ❌{fail_count} | ${total_cost:.4f} | {total_elapsed:.0f}초"
        )
    except Exception:
        pass

    # Notion 보고 에이전트 spawn
    try:
        outputs_summary = "\n\n".join(
            f"### {t['name']} ({t['type']})\n**태스크:** {t['task'][:200]}\n\n{t['output'][:1000]}"
            for t in team_outputs
        )
        notion_task = (
            f"[보고서]\n## 프로젝트 완료 리포트\n\n{outputs_summary}\n\n"
            f"[태스크 컨텍스트]\n프로젝트ID: {project_id} | 팀 수: {len(team_outputs)}개"
        )
        result = _spawn_solo_agent(
            task=notion_task,
            team_type="report-notion",
            name=f"NotionReport-{project_id[:6]}",
            project_id=project_id,
        )
        notion_agent_id = result.get("agent_id") or result.get("team_id")
        notion_proj_id = result.get("project_id")
    except Exception as e:
        print(f"[warning] Notion 보고 자동 실행 실패: {e}", file=sys.stderr, flush=True)


# ── 내부 헬퍼: solo 에이전트 subprocess 실행 ──────────────────────────────────

def _spawn_solo_agent(task: str, team_type: str = "general", name: str = None,
                      skills: list[str] = None, project_id: str = None) -> dict:
    """solo 에이전트를 detached subprocess로 실행하고 ID 정보 반환 (내부 헬퍼)"""
    import os
    import subprocess
    import tempfile

    DETACHED_PROCESS = 0x00000008
    CREATE_NEW_PROCESS_GROUP = 0x00000200
    CREATE_NO_WINDOW = 0x08000000
    RUN_AGENT_SCRIPT = str(Path(__file__).parent / "run_agent.py")

    proj_id = project_id or f"quick-{str(uuid.uuid4())[:6]}"
    agent_id = f"{proj_id}-1"
    agent_name = name or f"Quick-{team_type}"

    st.create_project(proj_id, f"[특무] {task[:80]}")
    st.upsert_team(st.Team(
        team_id=agent_id, name=agent_name, team_type=team_type,
        status="pending", task=task, project_id=proj_id,
        created_at=time.time(),
    ))

    full_task = _inject_skill_hints(task, skills or [])
    tmp = tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", suffix=".txt",
        prefix=f"task_{agent_id}_", delete=False
    )
    tmp.write(full_task)
    tmp.close()

    clean_env = {k: v for k, v in os.environ.items() if k != "CLAUDECODE"}
    subprocess.Popen(
        [sys.executable, RUN_AGENT_SCRIPT,
         agent_id, team_type, proj_id, str(Path.home()), tmp.name, "solo"],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        env=clean_env,
        creationflags=DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP | CREATE_NO_WINDOW,
        close_fds=True,
    )
    return {"agent_id": agent_id, "name": agent_name, "project_id": proj_id}


# ── 특무 에이전트 (단독 실행) ──────────────────────────────────────────────────

@mcp.tool()
def run_quick_agent(task: str, team_type: str = "general", name: str = None,
                    skills: list[str] = None) -> str:
    """
    특무 에이전트 — 단일 에이전트를 board 없이 즉시 실행합니다.
    팀 리드/워커 계층 없음. 단순·빠른 작업에 최적.

    start_project vs run_quick_agent 선택 기준:
    - run_quick_agent: 단일 태스크, 서브태스크 분해 불필요, 에이전트 간 소통 불필요
    - start_project: 여러 병렬 태스크, 리드+워커 계층 필요, 팀 간 의존성 있음

    Args:
        task: 실행할 태스크 (구체적으로 작성)
        team_type: 에이전트 유형. general|research|code|ops|bioinformatics|data-analysis|
                   writing|lab-protocol|literature|db|seq|struct|omics|stats|viz|ml|
                   bioeng|primer|git|asana|env (기본: general)
        name: 에이전트 이름 (기본: Quick-{team_type})
        skills: 참고할 스킬 이름 목록 (선택)

    Returns:
        agent_id — get_output(agent_id)로 결과 조회
    """
    result = _spawn_solo_agent(task=task, team_type=team_type, name=name, skills=skills)
    nl.log_project(result["project_id"], task)
    return json.dumps({
        **result,
        "team_type": team_type,
        "message": f"특무 에이전트 시작. get_output('{result['agent_id']}')로 결과 확인."
    }, ensure_ascii=False, indent=2)


# ── 스마트 라우팅 ──────────────────────────────────────────────────────────────

@mcp.tool()
def start_smart(description: str, feedback_timeout: int = 10) -> str:
    """
    태스크 복잡도를 분석해 자동으로 실행 방식을 선택합니다.
    - 단순 태스크 (1팀, 의존성 없음) → run_quick_agent (solo, 빠름, 토큰 절약)
    - 복잡 태스크 (다팀, 의존성 있음) → start_project (리드+워커, 병렬 처리)

    Args:
        description: 태스크 설명 (자연어)
        feedback_timeout: start_project 라우팅 시 피드백 대기 시간(초)

    Returns:
        라우팅 결과 + agent_id 또는 project_id
    """
    plan = task_router.analyze_and_split(description)
    teams = plan.get("teams", [])
    has_deps = any(t.get("depends_on") for t in teams)

    if len(teams) <= 1 and not has_deps:
        # 단순 → 특무 에이전트
        team = teams[0] if teams else {"type": "general", "task": description, "skills": []}
        result = _spawn_solo_agent(
            task=_inject_skill_hints(team["task"], team.get("skills", [])),
            team_type=team["type"],
            name=team.get("name", f"Smart-{team['type']}"),
        )
        return json.dumps({
            **result,
            "mode": "solo",
            "team_type": team["type"],
            "message": (
                f"단순 태스크 → 특무 에이전트 실행.\n"
                f"get_output('{result['agent_id']}')로 결과 확인."
            )
        }, ensure_ascii=False, indent=2)
    else:
        # 복잡 → 프로젝트 팀 구성
        return start_project(description, feedback_timeout)


# ── 출력 실시간 조회 ───────────────────────────────────────────────────────────

@mcp.tool()
def tail_output(team_id: str, since_offset: int = 0) -> str:
    """
    에이전트 출력을 실시간으로 조회합니다 (폴링용).
    since_offset을 이용해 마지막으로 읽은 위치 이후 내용만 가져올 수 있습니다.

    Args:
        team_id: 에이전트 ID
        since_offset: 마지막으로 읽은 출력 길이 (0이면 전체 반환)

    Returns:
        JSON: {output_chunk, total_offset, status, done, elapsed_seconds}
        - output_chunk: since_offset 이후 새 내용
        - total_offset: 현재 전체 출력 길이 (다음 호출 시 since_offset으로 사용)
        - elapsed_seconds: 경과 시간 (running 상태일 때만 실시간 계산)
        - done: True이면 에이전트 완료 (더 이상 폴링 불필요)
    """
    team = st.get_team(team_id)
    if not team:
        return json.dumps({"error": f"Agent {team_id} not found."})

    output = team.output
    total_len = len(output)
    chunk = output[since_offset:] if since_offset < total_len else ""

    # 실시간 경과 시간 계산
    current_time = time.time()
    elapsed = 0.0
    if team.status == "running" and team.started_at > 0:
        elapsed = current_time - team.started_at
    else:
        elapsed = team.elapsed_seconds or 0.0

    return json.dumps({
        "team_id": team_id,
        "name": team.name,
        "status": team.status,
        "output_chunk": chunk,
        "total_offset": total_len,
        "elapsed_seconds": round(elapsed, 1),
        "done": team.status in ("done", "failed", "stuck"),
    }, ensure_ascii=False)


# ── 사용량 리포트 ──────────────────────────────────────────────────────────────

@mcp.tool()
def get_usage_report(project_id: str = None) -> str:
    """
    에이전트별 토큰 사용량·모델·응답시간을 조회합니다.

    Args:
        project_id: (선택) 특정 프로젝트만 조회. 없으면 전체.

    Returns:
        에이전트별 사용량 표 + 합계
    """
    teams = st.list_teams(project_id)
    if not teams:
        return "No agents found."

    total_in = total_out = 0
    total_cost = 0.0
    model_counts: dict[str, int] = {}

    lines = ["## 에이전트 사용량 리포트\n",
             "| 에이전트 | 팀 유형 | 모델 | 입력토큰 | 출력토큰 | 소요(초) | 비용($) | 상태 |",
             "|---------|--------|------|---------|---------|---------|--------|------|"]

    for t in sorted(teams, key=lambda x: x.created_at):
        in_tok = t.input_tokens or 0
        out_tok = t.output_tokens or 0
        cost = t.total_cost_usd or 0.0
        elapsed = t.elapsed_seconds or 0.0
        model = t.model_used or "-"

        total_in += in_tok
        total_out += out_tok
        total_cost += cost
        if model != "-":
            model_counts[model] = model_counts.get(model, 0) + 1

        # 모델명 축약
        short_model = model.replace("claude-", "").replace("-20251001", "")
        in_str = f"{in_tok:,}" if in_tok else "-"
        out_str = f"{out_tok:,}" if out_tok else "-"
        elapsed_str = f"{elapsed:.0f}" if elapsed else "-"
        cost_str = f"{cost:.4f}" if cost else "-"

        lines.append(
            f"| {t.name} | {t.team_type} | {short_model} "
            f"| {in_str} | {out_str} | {elapsed_str} | {cost_str} | {t.status} |"
        )

    lines.append("")
    lines.append(f"**합계** — 입력: {total_in:,} | 출력: {total_out:,} | 총 비용: ${total_cost:.4f}")
    if model_counts:
        model_summary = " / ".join(
            f"{m.replace('claude-','').replace('-20251001','')}: {c}회"
            for m, c in sorted(model_counts.items(), key=lambda x: -x[1])
        )
        lines.append(f"**모델 사용** — {model_summary}")

    return "\n".join(lines)


@mcp.tool()
def get_dashboard_url() -> str:
    """
    실시간 대시보드의 공개 Cloudflare Tunnel URL을 반환합니다.
    Notion embed 등에 사용할 수 있습니다.
    """
    if _DASHBOARD_PUBLIC_URL:
        return f"Dashboard URL: {_DASHBOARD_PUBLIC_URL}"
    return "Dashboard URL not available (cloudflared not installed or tunnel not started)"


if __name__ == "__main__":
    # Fix: Python 3.13 + Windows PIPE에서 stdout.flush() → OSError [Errno 22] 발생
    # anyio의 AsyncFile.flush에서 OSError를 suppress해서 연결 유지
    try:
        import anyio._core._fileio as _anyio_fio
        _orig_flush = _anyio_fio.AsyncFile.flush
        _orig_write = _anyio_fio.AsyncFile.write

        async def _safe_flush(self):
            try:
                await _orig_flush(self)
            except OSError:
                pass  # Windows PIPE flush 실패 무시

        async def _safe_write(self, b):
            try:
                return await _orig_write(self, b)
            except (OSError, BrokenPipeError) as _we:
                _logger.warning("anyio.AsyncFile.write 오류 (무시): %s", _we)
                return 0

        _anyio_fio.AsyncFile.flush = _safe_flush
        _anyio_fio.AsyncFile.write = _safe_write
        _logger.info("anyio.AsyncFile flush+write OSError-safe 패치 적용")
    except Exception as _e:
        _logger.warning("패치 실패 (무시): %s", _e)

    # atexit에서 서버 종료 시간 기록
    def _log_exit():
        _logger.info("=== server.py 종료 ===")
        logging.shutdown()
    atexit.register(_log_exit)

    # 백그라운드 스레드 시작 (mcp.run() 전에, __main__ 컨텍스트에서)
    threading.Thread(target=_start_dashboard_background, daemon=True).start()
    threading.Thread(target=_start_tg_listener, daemon=True).start()
    _logger.info("백그라운드 스레드 시작")

    try:
        _logger.info("mcp.run() 호출")
        mcp.run()
        _logger.info("mcp.run() 정상 종료 (stdin EOF - 클라이언트 연결 종료)")
    except SystemExit as _se:
        _logger.info("mcp.run() SystemExit: code=%s", _se.code)
    except Exception as _ex:
        _logger.critical("mcp.run() 예외:\n%s", traceback.format_exc())
