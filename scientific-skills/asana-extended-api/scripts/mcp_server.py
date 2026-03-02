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
from mcp.server.fastmcp import FastMCP

ASANA_BASE = "https://app.asana.com/api/1.0"

mcp = FastMCP("asana-extended")


def _get_pat():
    pat = os.environ.get("ASANA_PAT")
    if not pat:
        raise RuntimeError(
            "ASANA_PAT 환경변수가 설정되지 않았습니다. "
            "https://app.asana.com/0/my-apps 에서 PAT를 발급하세요."
        )
    return pat


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


if __name__ == "__main__":
    mcp.run(transport="stdio")
