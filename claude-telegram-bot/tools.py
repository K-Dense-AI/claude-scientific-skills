"""
Tool implementations for Claude Telegram Bot.
Each tool is invoked via Claude API function calling.
"""

import json
import os
import subprocess
import sys
import time
import uuid
from pathlib import Path

SECRETS_FILE = Path.home() / ".claude" / "secrets.json"
ASANA_SCRIPT_DIR = (
    Path.home()
    / "claude-scientific-skills"
    / "scientific-skills"
    / "asana-extended-api"
    / "scripts"
)
ORCH_DIR = (
    Path.home()
    / "claude-scientific-skills"
    / "mcp-servers"
    / "team-orchestrator"
)


def _load_secrets() -> dict:
    return json.loads(SECRETS_FILE.read_text(encoding="utf-8"))


def _asana_request(method: str, path: str, data: dict = None) -> dict:
    import urllib.request, urllib.error

    secrets = _load_secrets()
    url = f"https://app.asana.com/api/1.0{path}"
    body = json.dumps({"data": data}).encode() if data else None
    req = urllib.request.Request(
        url,
        data=body,
        method=method,
        headers={
            "Authorization": f"Bearer {secrets['ASANA_PAT']}",
            "Content-Type": "application/json",
        },
    )
    try:
        resp = urllib.request.urlopen(req)
        content = resp.read().decode()
        return json.loads(content) if content.strip() else {}
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"Asana API {e.code}: {e.read().decode()[:300]}")


def _notion_request(method: str, path: str, body: dict = None) -> dict:
    import urllib.request, urllib.error

    secrets = _load_secrets()
    url = f"https://api.notion.com/v1{path}"
    data = json.dumps(body).encode() if body else None
    req = urllib.request.Request(
        url,
        data=data,
        method=method,
        headers={
            "Authorization": f"Bearer {secrets['NOTION_TOKEN']}",
            "Notion-Version": "2022-06-28",
            "Content-Type": "application/json",
        },
    )
    try:
        resp = urllib.request.urlopen(req)
        return json.loads(resp.read().decode())
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"Notion API {e.code}: {e.read().decode()[:300]}")


def _github_request(method: str, path: str, body: dict = None) -> dict:
    import urllib.request, urllib.error

    secrets = _load_secrets()
    url = f"https://api.github.com{path}"
    data = json.dumps(body).encode() if body else None
    req = urllib.request.Request(
        url,
        data=data,
        method=method,
        headers={
            "Authorization": f"Bearer {secrets['GITHUB_PAT']}",
            "Accept": "application/vnd.github+json",
            "Content-Type": "application/json",
        },
    )
    try:
        resp = urllib.request.urlopen(req)
        return json.loads(resp.read().decode())
    except urllib.error.HTTPError as e:
        raise RuntimeError(f"GitHub API {e.code}: {e.read().decode()[:300]}")


# ── Asana ──────────────────────────────────────────────────────────────────────

def asana_get_workspace() -> str:
    """사용자의 기본 워크스페이스 GID 반환"""
    result = _asana_request("GET", "/workspaces")
    workspaces = result.get("data", [])
    if not workspaces:
        raise RuntimeError("Asana 워크스페이스를 찾을 수 없습니다.")
    return workspaces[0]["gid"]


def asana_create_task(name: str, notes: str = "", project_gid: str = None, due_date: str = None) -> dict:
    """Asana 태스크 생성"""
    workspace = asana_get_workspace()
    data = {
        "name": name,
        "workspace": workspace,
        "assignee": "me",
    }
    if notes:
        data["html_notes"] = f"<body>{notes}</body>"
    if project_gid:
        data["projects"] = [project_gid]
    if due_date:
        data["due_on"] = due_date

    result = _asana_request("POST", "/tasks", data)
    task = result["data"]
    return {
        "gid": task["gid"],
        "name": task["name"],
        "url": f"https://app.asana.com/0/0/{task['gid']}/f",
    }


def asana_list_tasks(project_gid: str = None, limit: int = 10) -> list:
    """Asana 태스크 목록 조회"""
    fields = "name,due_on,completed,permalink_url,assignee_status"
    if project_gid:
        path = f"/projects/{project_gid}/tasks?opt_fields={fields}&limit={limit}"
    else:
        workspace = asana_get_workspace()
        path = (
            f"/tasks?assignee=me&workspace={workspace}"
            f"&opt_fields={fields}&limit={limit}&completed_since=now"
        )
    result = _asana_request("GET", path)
    return [
        {
            "gid": t["gid"],
            "name": t["name"],
            "due": t.get("due_on"),
            "done": t.get("completed", False),
            "url": t.get("permalink_url", f"https://app.asana.com/0/0/{t['gid']}/f"),
        }
        for t in result.get("data", [])
    ]


def asana_update_task(task_gid: str, name: str = None, notes: str = None,
                      due_date: str = None, completed: bool = None) -> dict:
    """Asana 태스크 업데이트"""
    data = {}
    if name is not None:
        data["name"] = name
    if notes is not None:
        data["html_notes"] = f"<body>{notes}</body>"
    if due_date is not None:
        data["due_on"] = due_date
    if completed is not None:
        data["completed"] = completed
    result = _asana_request("PUT", f"/tasks/{task_gid}", data)
    task = result["data"]
    return {
        "gid": task["gid"],
        "name": task["name"],
        "url": f"https://app.asana.com/0/0/{task['gid']}/f",
    }


# ── Notion ─────────────────────────────────────────────────────────────────────

def notion_search(query: str, limit: int = 5) -> list:
    """Notion 페이지/데이터베이스 검색"""
    result = _notion_request("POST", "/search", {"query": query, "page_size": limit})
    items = []
    for obj in result.get("results", []):
        title = ""
        if obj["object"] == "page":
            props = obj.get("properties", {})
            for v in props.values():
                if v.get("type") == "title":
                    title = "".join(t.get("plain_text", "") for t in v.get("title", []))
                    break
        elif obj["object"] == "database":
            title_obj = obj.get("title", [])
            title = "".join(t.get("plain_text", "") for t in title_obj)
        items.append({
            "id": obj["id"],
            "type": obj["object"],
            "title": title or "(제목 없음)",
            "url": obj.get("url", ""),
        })
    return items


def notion_create_page(parent_page_id: str, title: str, content: str = "") -> dict:
    """Notion 페이지 생성"""
    blocks = []
    if content:
        for line in content.split("\n")[:50]:
            if line.strip():
                blocks.append({
                    "object": "block",
                    "type": "paragraph",
                    "paragraph": {
                        "rich_text": [{"type": "text", "text": {"content": line}}]
                    },
                })
    body = {
        "parent": {"type": "page_id", "page_id": parent_page_id.replace("-", "")},
        "properties": {
            "title": {"title": [{"type": "text", "text": {"content": title}}]}
        },
        "children": blocks,
    }
    result = _notion_request("POST", "/pages", body)
    return {"id": result["id"], "url": result.get("url", "")}


def notion_get_page(page_id: str) -> dict:
    """Notion 페이지 정보 조회"""
    result = _notion_request("GET", f"/pages/{page_id.replace('-', '')}")
    props = result.get("properties", {})
    title = ""
    for v in props.values():
        if v.get("type") == "title":
            title = "".join(t.get("plain_text", "") for t in v.get("title", []))
            break
    return {"id": result["id"], "title": title, "url": result.get("url", "")}


# ── GitHub ─────────────────────────────────────────────────────────────────────

def github_create_issue(repo: str, title: str, body: str = "") -> dict:
    """GitHub 이슈 생성 (repo: owner/repo)"""
    result = _github_request("POST", f"/repos/{repo}/issues", {"title": title, "body": body})
    return {
        "number": result["number"],
        "title": result["title"],
        "url": result["html_url"],
    }


def github_list_issues(repo: str, state: str = "open", limit: int = 10) -> list:
    """GitHub 이슈 목록 조회"""
    result = _github_request("GET", f"/repos/{repo}/issues?state={state}&per_page={limit}")
    return [
        {
            "number": i["number"],
            "title": i["title"],
            "url": i["html_url"],
            "state": i["state"],
            "labels": [l["name"] for l in i.get("labels", [])],
        }
        for i in result
        if isinstance(i, dict)
    ]


def github_comment_issue(repo: str, issue_number: int, comment: str) -> dict:
    """GitHub 이슈에 댓글 추가"""
    result = _github_request(
        "POST",
        f"/repos/{repo}/issues/{issue_number}/comments",
        {"body": comment},
    )
    return {"id": result["id"], "url": result["html_url"]}


# ── File System ────────────────────────────────────────────────────────────────

def file_read(path: str) -> str:
    """로컬 파일 읽기"""
    p = Path(path).expanduser()
    if not p.exists():
        return f"파일 없음: {path}"
    if p.stat().st_size > 100_000:
        lines = p.read_text(encoding="utf-8", errors="replace").splitlines()
        return f"(파일이 큼 — 첫 100줄만)\n" + "\n".join(lines[:100])
    return p.read_text(encoding="utf-8", errors="replace")


def file_write(path: str, content: str) -> dict:
    """로컬 파일 쓰기"""
    p = Path(path).expanduser()
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(content, encoding="utf-8")
    return {"path": str(p), "bytes": len(content.encode("utf-8"))}


def file_list(directory: str = "~/lab-analyses") -> list:
    """디렉토리 목록 조회"""
    p = Path(directory).expanduser()
    if not p.exists():
        return [{"error": f"디렉토리 없음: {directory}"}]
    items = []
    for f in sorted(p.iterdir())[:40]:
        items.append({
            "name": f.name,
            "type": "dir" if f.is_dir() else "file",
            "size": f.stat().st_size if f.is_file() else None,
        })
    return items


def bash_run(command: str) -> str:
    """셸 명령 실행 (안전 필터 적용)"""
    dangerous_patterns = ["rm -rf", "del /f /s", "format ", "shutdown", "mkfs", "dd if="]
    for pattern in dangerous_patterns:
        if pattern.lower() in command.lower():
            return f"[차단] 위험 명령 포함: '{pattern}'"
    try:
        result = subprocess.run(
            command,
            shell=True,
            capture_output=True,
            text=True,
            timeout=30,
            encoding="utf-8",
            errors="replace",
        )
        output = (result.stdout + result.stderr).strip()
        return output[:3000] if output else "(출력 없음)"
    except subprocess.TimeoutExpired:
        return "[타임아웃] 30초 초과"
    except Exception as e:
        return f"[오류] {e}"


# ── Team Orchestrator ──────────────────────────────────────────────────────────

def agent_start(task: str, team_type: str = "general", agent_name: str = None) -> dict:
    """Claude Code 에이전트를 비동기 subprocess로 시작"""
    import tempfile

    DETACHED_PROCESS = 0x00000008
    CREATE_NEW_PROCESS_GROUP = 0x00000200

    run_agent_script = str(ORCH_DIR / "run_agent.py")
    if not Path(run_agent_script).exists():
        return {"error": "run_agent.py 없음. team-orchestrator가 설치되어 있는지 확인하세요."}

    proj_id = f"tg-{str(uuid.uuid4())[:6]}"
    agent_id = f"{proj_id}-1"
    name = agent_name or f"TG-{team_type}"

    tmp = tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", suffix=".txt",
        prefix=f"task_{agent_id}_", delete=False,
    )
    tmp.write(task)
    tmp.close()

    clean_env = {k: v for k, v in os.environ.items() if k != "CLAUDECODE"}

    try:
        subprocess.Popen(
            [sys.executable, run_agent_script,
             agent_id, team_type, proj_id, str(Path.home()), tmp.name, "solo"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            env=clean_env,
            creationflags=DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP,
            close_fds=True,
        )
        return {
            "agent_id": agent_id,
            "project_id": proj_id,
            "name": name,
            "status": "started",
            "message": f"에이전트 시작. agent_get_output('{agent_id}')로 결과 확인.",
        }
    except Exception as e:
        return {"error": f"에이전트 시작 실패: {e}"}


def agent_get_output(agent_id: str) -> dict:
    """에이전트 결과 조회"""
    import sqlite3

    orch_db = ORCH_DIR / "orchestrator.db"
    if not orch_db.exists():
        return {"error": "오케스트레이터 DB 없음. 에이전트를 먼저 시작하세요."}

    conn = sqlite3.connect(str(orch_db))
    conn.row_factory = sqlite3.Row
    row = conn.execute(
        "SELECT name, status, output, created_at FROM teams WHERE team_id=?",
        (agent_id,),
    ).fetchone()
    conn.close()

    if not row:
        return {"error": f"에이전트 '{agent_id}' 없음"}

    output = (row["output"] or "").strip()
    return {
        "agent_id": agent_id,
        "name": row["name"],
        "status": row["status"],
        "output": output[:2000] + ("..." if len(output) > 2000 else ""),
        "done": row["status"] in ("done", "failed"),
    }


def agent_list() -> list:
    """최근 에이전트 목록 조회"""
    import sqlite3

    orch_db = ORCH_DIR / "orchestrator.db"
    if not orch_db.exists():
        return [{"error": "오케스트레이터 DB 없음"}]

    conn = sqlite3.connect(str(orch_db))
    conn.row_factory = sqlite3.Row
    rows = conn.execute(
        "SELECT team_id, name, status, team_type, created_at FROM teams ORDER BY created_at DESC LIMIT 15"
    ).fetchall()
    conn.close()

    return [
        {
            "agent_id": r["team_id"],
            "name": r["name"],
            "status": r["status"],
            "type": r["team_type"],
        }
        for r in rows
    ]
