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


# ── Claude 기반 동기화 플래너 ───────────────────────────────────────────────

def plan_sync(team_outputs: list[dict], asana_tasks: list[dict]) -> dict:
    """
    팀 산출물 + 현재 Asana 태스크를 Claude에게 주고
    무엇을 새로 만들지/업데이트할지/건너뛸지 결정
    """
    import anthropic
    client = anthropic.Anthropic(api_key=_get_anthropic_key())

    system = """당신은 프로젝트 관리 전문가입니다.
팀 작업 완료 후 Asana를 업데이트해야 합니다.
주어진 팀 산출물과 현재 Asana 태스크 목록을 비교하여
무엇을 해야 하는지 판단하세요.

반드시 아래 JSON 형식으로만 응답하세요:
{
  "actions": [
    {
      "action": "create|update|skip",
      "reason": "한 줄 이유",
      "task_name": "태스크 이름",
      "notes": "태스크 내용 (create/update 시)",
      "asana_gid": "업데이트 대상 gid (update 시만)"
    }
  ],
  "summary": "전체 동기화 요약"
}"""

    outputs_text = "\n\n".join([
        f"[{t['name']} ({t['type']})] {t['output'][:500]}..."
        for t in team_outputs
    ])

    asana_text = "\n".join([
        f"- [{t['gid']}] {t['name']}: {t.get('notes', '')[:100]}"
        for t in asana_tasks
    ])

    prompt = f"""## 팀 산출물
{outputs_text}

## 현재 Asana 태스크 (미완료)
{asana_text if asana_text else "(없음)"}

중복 생성 없이 필요한 업데이트만 계획해주세요."""

    response = client.messages.create(
        model="claude-sonnet-4-6",
        max_tokens=2000,
        system=system,
        messages=[{"role": "user", "content": prompt}]
    )

    raw = response.content[0].text.strip()
    # 코드블록 제거
    import re
    match = re.search(r"```(?:json)?\s*([\s\S]+?)\s*```", raw)
    if match:
        raw = match.group(1)
    else:
        # 중괄호 기준으로 JSON만 추출
        match = re.search(r"\{[\s\S]+\}", raw)
        if match:
            raw = match.group(0)
    return json.loads(raw)


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
