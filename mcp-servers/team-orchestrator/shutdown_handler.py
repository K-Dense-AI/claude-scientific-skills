"""종료 시 Asana + Notion 동기화"""
import json
import time
from pathlib import Path
from typing import Optional

# anthropic, requests는 shutdown 시에만 필요 → lazy import
# (top-level import 하면 MCP 서버 시작 시 패키지 없으면 전체 import 실패)


def _get_secrets() -> dict:
    secrets_path = Path.home() / ".claude" / "secrets.json"
    if secrets_path.exists():
        return json.loads(secrets_path.read_text(encoding="utf-8"))
    return {}


def _get_asana_pat() -> str:
    secrets = _get_secrets()
    pat = secrets.get("ASANA_PAT") or ""
    if not pat:
        import os
        pat = os.environ.get("ASANA_PAT", "")
    if not pat:
        raise RuntimeError("ASANA_PAT가 없습니다.")
    return pat


def _get_asana_api():
    """asana_api.py 헬퍼 모듈 로드 (단일 소스로 통합)."""
    import sys
    api_dir = str(Path.home() / ".claude" / "skills" / "asana-extended-api" / "scripts")
    if api_dir not in sys.path:
        sys.path.insert(0, api_dir)
    from asana_api import AsanaAPI
    return AsanaAPI(pat=_get_asana_pat())


def _get_anthropic_key() -> str:
    secrets = _get_secrets()
    key = secrets.get("ANTHROPIC_API_KEY") or ""
    if not key:
        import os
        key = os.environ.get("ANTHROPIC_API_KEY", "")
    if not key:
        raise RuntimeError("ANTHROPIC_API_KEY가 없습니다.")
    return key


# ── Asana 조회 ──────────────────────────────────────────────────────────────

def fetch_asana_tasks(workspace_gid: str = None) -> list[dict]:
    """asana_api.py 헬퍼를 통해 내 미완료 태스크 조회 (단일 소스 통합)."""
    api = _get_asana_api()
    return api.get_my_tasks(workspace_gid=workspace_gid)


def fetch_asana_workspaces() -> list[dict]:
    """asana_api.py 헬퍼를 통해 워크스페이스 목록 조회 (단일 소스 통합)."""
    api = _get_asana_api()
    return api.get_workspaces()


def create_asana_task(workspace_gid: str, name: str, notes: str, due_on: str = None) -> dict:
    """asana_api.py 헬퍼를 통해 태스크 생성 (단일 소스 통합)."""
    api = _get_asana_api()
    return api.create_task(workspace_gid, name, notes=notes, due_on=due_on)


def update_asana_task(task_gid: str, notes: str) -> dict:
    """asana_api.py 헬퍼를 통해 태스크 업데이트 (단일 소스 통합)."""
    api = _get_asana_api()
    return api.update_task(task_gid, notes=notes)


# ── 규칙 기반 동기화 플래너 (Claude API 불필요) ──────────────────────────────

def plan_sync(team_outputs: list[dict], asana_tasks: list[dict]) -> dict:
    """
    팀 산출물을 기반으로 Asana 동기화 계획 수립.
    Claude API 대신 규칙 기반으로 처리: 항상 새 태스크 생성.
    (기존 태스크 이름 중복 시 업데이트, 아니면 새로 생성)
    """
    existing_names = {t["name"].strip().lower(): t["gid"] for t in asana_tasks}
    actions = []

    for output in team_outputs:
        agent_name = output.get("name", "agent")
        agent_type = output.get("type", "")
        content = output.get("output", "").strip()
        if not content:
            continue

        task_name = f"[팀작업] {agent_name} 완료 보고"
        notes = f"에이전트: {agent_name} ({agent_type})\n\n{content[:2000]}"

        # 동일 이름 태스크가 이미 있으면 업데이트, 없으면 생성
        existing_gid = existing_names.get(task_name.strip().lower())
        if existing_gid:
            actions.append({
                "action": "update",
                "reason": "동일 이름 태스크 이미 존재",
                "task_name": task_name,
                "notes": notes,
                "asana_gid": existing_gid,
            })
        else:
            actions.append({
                "action": "create",
                "reason": "새 팀 작업 완료 보고",
                "task_name": task_name,
                "notes": notes,
            })

    return {
        "actions": actions,
        "summary": f"총 {len(actions)}개 태스크 동기화 계획 (규칙 기반)",
    }


def execute_sync(plan: dict, workspace_gid: str) -> list[str]:
    """동기화 계획 실행 → 결과 로그 반환"""
    logs = []
    for action in plan.get("actions", []):
        act = action["action"]
        if act == "create":
            result = create_asana_task(workspace_gid, action["task_name"], action["notes"])
            logs.append(f"[생성] {action['task_name']} → {result.get('permalink_url', result.get('gid', ''))}")
        elif act == "update":
            result = update_asana_task(action["asana_gid"], action["notes"])
            logs.append(f"[업데이트] {action['task_name']} (gid: {action['asana_gid']})")
        else:
            logs.append(f"[건너뜀] {action['task_name']}: {action['reason']}")
    return logs


def run_shutdown_sync(team_outputs: list[dict]) -> dict:
    """
    전체 종료 시퀀스:
    1. Asana 현재 태스크 조회
    2. Claude로 동기화 계획 수립
    3. 계획 실행
    4. 요약 반환
    """
    # 워크스페이스 첫 번째 사용
    workspaces = fetch_asana_workspaces()
    if not workspaces:
        return {"error": "Asana 워크스페이스 없음"}

    workspace_gid = workspaces[0]["gid"]
    current_tasks = fetch_asana_tasks(workspace_gid)

    # 산출물이 있는 팀만 필터
    meaningful_outputs = [t for t in team_outputs if t.get("output", "").strip()]
    if not meaningful_outputs:
        return {"summary": "산출물 없음, 동기화 건너뜀", "actions": []}

    plan = plan_sync(meaningful_outputs, current_tasks)
    logs = execute_sync(plan, workspace_gid)

    return {
        "summary": plan.get("summary", ""),
        "actions_count": len(plan.get("actions", [])),
        "logs": logs
    }
