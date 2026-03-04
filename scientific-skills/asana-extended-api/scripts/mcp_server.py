"""
Asana Extended API — Local MCP Server
MCP 도구에 없는 Asana REST API(addProject, removeProject, batch 등)를 제공.

실행: python mcp_server.py  (Claude Code가 자동 실행)
환경변수: ASANA_PAT 필수
"""
import os
import json
import urllib.request
import urllib.error
import urllib.parse
from pathlib import Path
from mcp.server.fastmcp import FastMCP

ASANA_BASE = "https://app.asana.com/api/1.0"
_SECRETS_FILE = Path.home() / ".claude" / "secrets.json"

mcp = FastMCP("asana-extended")


def _get_pat():
    """PAT 로드 우선순위: 환경변수 → ~/.claude/secrets.json"""
    pat = os.environ.get("ASANA_PAT")
    if pat:
        return pat
    if _SECRETS_FILE.exists():
        try:
            secrets = json.loads(_SECRETS_FILE.read_text(encoding="utf-8"))
            pat = secrets.get("ASANA_PAT")
            if pat:
                return pat
        except Exception:
            pass
    raise RuntimeError(
        f"ASANA_PAT를 찾을 수 없습니다. "
        f"{_SECRETS_FILE} 에 {{\"ASANA_PAT\": \"2/xxxx...\"}} 저장하거나 "
        "환경변수를 설정하세요."
    )


def _request(method: str, path: str, data: dict | None = None) -> dict:
    url = f"{ASANA_BASE}{path}"
    body = json.dumps({"data": data}).encode() if data else None
    req = urllib.request.Request(url, data=body, method=method, headers={
        "Authorization": f"Bearer {_get_pat()}",
        "Content-Type": "application/json",
    })
    try:
        resp = urllib.request.urlopen(req)
        content = resp.read().decode()
        return json.loads(content) if content.strip() else {}
    except urllib.error.HTTPError as e:
        error_body = e.read().decode()
        return {"error": f"HTTP {e.code}", "detail": error_body[:500]}


@mcp.tool()
def asana_add_task_to_project(
    task_gid: str,
    project_gid: str,
    section_gid: str = "",
) -> str:
    """태스크를 프로젝트에 추가 (multi-homing). 원본 프로젝트에서 제거되지 않음.
    opt_silent=true 자동 적용 (알림 없음)."""
    data = {"project": project_gid}
    if section_gid:
        data["section"] = section_gid
    result = _request("POST", f"/tasks/{task_gid}/addProject?opt_silent=true", data)
    if "error" in result:
        return json.dumps(result, ensure_ascii=False)
    return f"태스크 {task_gid} → 프로젝트 {project_gid} 추가 완료"


@mcp.tool()
def asana_remove_task_from_project(
    task_gid: str,
    project_gid: str,
) -> str:
    """태스크를 프로젝트에서 제거. 태스크 자체는 삭제되지 않음."""
    result = _request("POST", f"/tasks/{task_gid}/removeProject", {"project": project_gid})
    if "error" in result:
        return json.dumps(result, ensure_ascii=False)
    return f"태스크 {task_gid} → 프로젝트 {project_gid}에서 제거 완료"


@mcp.tool()
def asana_batch_add_to_project(
    task_gids: list[str],
    project_gid: str,
    section_gid: str = "",
) -> str:
    """여러 태스크를 한 프로젝트에 일괄 추가. opt_silent=true 자동 적용."""
    success, fail = 0, 0
    errors = []
    for gid in task_gids:
        data = {"project": project_gid}
        if section_gid:
            data["section"] = section_gid
        result = _request("POST", f"/tasks/{gid}/addProject?opt_silent=true", data)
        if "error" in result:
            fail += 1
            errors.append(f"{gid}: {result['error']}")
        else:
            success += 1
    summary = f"완료: 성공 {success}건, 실패 {fail}건"
    if errors:
        summary += "\n실패 목록:\n" + "\n".join(errors)
    return summary


@mcp.tool()
def asana_create_section(
    project_gid: str,
    name: str,
) -> str:
    """프로젝트에 새 섹션 생성. 생성된 섹션 GID 반환."""
    result = _request("POST", f"/projects/{project_gid}/sections", {"name": name})
    if "error" in result:
        return json.dumps(result, ensure_ascii=False)
    section = result.get("data", {})
    return f"섹션 생성 완료: {section.get('name')} (GID: {section.get('gid')})"


@mcp.tool()
def asana_add_task_to_section(
    section_gid: str,
    task_gid: str,
) -> str:
    """태스크를 특정 섹션으로 이동."""
    result = _request("POST", f"/sections/{section_gid}/addTask", {"task": task_gid})
    if "error" in result:
        return json.dumps(result, ensure_ascii=False)
    return f"태스크 {task_gid} → 섹션 {section_gid} 이동 완료"


@mcp.tool()
def asana_add_tag_to_task(
    task_gid: str,
    tag_gid: str,
) -> str:
    """태스크에 태그 추가."""
    result = _request("POST", f"/tasks/{task_gid}/addTag", {"tag": tag_gid})
    if "error" in result:
        return json.dumps(result, ensure_ascii=False)
    return f"태스크 {task_gid}에 태그 {tag_gid} 추가 완료"


@mcp.tool()
def asana_remove_tag_from_task(
    task_gid: str,
    tag_gid: str,
) -> str:
    """태스크에서 태그 제거."""
    result = _request("POST", f"/tasks/{task_gid}/removeTag", {"tag": tag_gid})
    if "error" in result:
        return json.dumps(result, ensure_ascii=False)
    return f"태스크 {task_gid}에서 태그 {tag_gid} 제거 완료"


def _upload_multipart(path: str, fields: dict, files: dict) -> dict:
    """multipart/form-data로 파일 업로드."""
    boundary = b"----AsanaUploadBoundary7MA4YWxkTrZu0gW"
    body = b""
    for key, value in fields.items():
        body += b"--" + boundary + b"\r\n"
        body += f'Content-Disposition: form-data; name="{key}"\r\n\r\n'.encode()
        body += value.encode() + b"\r\n"
    for key, (filename, data, content_type) in files.items():
        body += b"--" + boundary + b"\r\n"
        body += f'Content-Disposition: form-data; name="{key}"; filename="{filename}"\r\n'.encode()
        body += f'Content-Type: {content_type}\r\n\r\n'.encode()
        body += data + b"\r\n"
    body += b"--" + boundary + b"--\r\n"
    url = f"{ASANA_BASE}{path}"
    req = urllib.request.Request(url, data=body, method="POST", headers={
        "Authorization": f"Bearer {_get_pat()}",
        "Content-Type": f"multipart/form-data; boundary={boundary.decode()}",
    })
    try:
        resp = urllib.request.urlopen(req)
        content = resp.read().decode()
        return json.loads(content) if content.strip() else {}
    except urllib.error.HTTPError as e:
        error_body = e.read().decode()
        return {"error": f"HTTP {e.code}", "detail": error_body[:500]}


@mcp.tool()
def asana_upload_attachment(task_gid: str, file_path: str) -> str:
    """로컬 파일을 Asana 태스크에 첨부 업로드. file_path는 절대경로."""
    import mimetypes
    p = Path(file_path)
    if not p.exists():
        return f"파일을 찾을 수 없습니다: {file_path}"
    content_type, _ = mimetypes.guess_type(str(p))
    if not content_type:
        content_type = "application/octet-stream"
    result = _upload_multipart(
        "/attachments",
        fields={"parent": task_gid},
        files={"file": (p.name, p.read_bytes(), content_type)},
    )
    if "error" in result:
        return json.dumps(result, ensure_ascii=False)
    att = result.get("data", {})
    return f"첨부 완료: {att.get('name')} (GID: {att.get('gid')})"


def _get_paginated(path: str, opt_fields: str, limit: int = 100) -> list:
    """페이지네이션 처리하여 전체 목록 반환."""
    results = []
    offset = None
    while True:
        url = f"{ASANA_BASE}{path}?opt_fields={opt_fields}&limit={limit}"
        if offset:
            url += f"&offset={offset}"
        result = _request("GET", url.replace(ASANA_BASE, ""))
        if "error" in result:
            return result
        results.extend(result.get("data", []))
        next_page = result.get("next_page")
        if not next_page:
            break
        offset = next_page.get("offset")
    return results


@mcp.tool()
def asana_get_tasks(
    project_id: str,
    opt_fields: str = "name,assignee.name,assignee.gid,gid,completed,completed_at,due_on,created_at,permalink_url",
    limit: int = 100,
) -> str:
    """프로젝트의 태스크 목록 조회 (페이지네이션 자동 처리).
    opt_fields 예시: name,assignee.name,gid,completed,due_on,completed_at,permalink_url"""
    result = _get_paginated(f"/projects/{project_id}/tasks", opt_fields, limit)
    if isinstance(result, dict) and "error" in result:
        return json.dumps(result, ensure_ascii=False)
    return json.dumps({"data": result}, ensure_ascii=False)


@mcp.tool()
def asana_get_task(
    task_gid: str,
    opt_fields: str = "name,assignee.name,gid,completed,completed_at,due_on,created_at,notes,permalink_url,projects.name",
) -> str:
    """태스크 상세 정보 조회. notes(본문), projects 등 포함 가능."""
    result = _request("GET", f"/tasks/{task_gid}?opt_fields={opt_fields}")
    if "error" in result:
        return json.dumps(result, ensure_ascii=False)
    return json.dumps(result, ensure_ascii=False)


@mcp.tool()
def asana_get_project_task_counts(
    project_id: str,
) -> str:
    """프로젝트 태스크 완료/미완료 수 조회."""
    result = _request("GET", f"/projects/{project_id}/task_counts?opt_fields=num_tasks,num_incomplete_tasks,num_completed_tasks")
    if "error" in result:
        return json.dumps(result, ensure_ascii=False)
    data = result.get("data", {})
    total = data.get("num_tasks", 0)
    completed = data.get("num_completed_tasks", 0)
    incomplete = data.get("num_incomplete_tasks", 0)
    rate = round(completed / total * 100, 1) if total > 0 else 0
    return json.dumps({"total": total, "completed": completed, "incomplete": incomplete, "rate": rate}, ensure_ascii=False)


@mcp.tool()
def asana_search_tasks_in_workspace(
    workspace_gid: str,
    text: str,
    opt_fields: str = "name,assignee.name,gid,completed,due_on,permalink_url,projects.name",
) -> str:
    """워크스페이스에서 키워드로 태스크 검색. Asana Free에서는 제한적."""
    path = f"/workspaces/{workspace_gid}/tasks/search?text={urllib.parse.quote(text)}&opt_fields={opt_fields}&limit=50"
    result = _request("GET", path)
    if "error" in result:
        return json.dumps(result, ensure_ascii=False)
    return json.dumps(result, ensure_ascii=False)


if __name__ == "__main__":
    mcp.run(transport="stdio")
