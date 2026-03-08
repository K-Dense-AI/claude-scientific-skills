"""
Notion 작동 로그 연동 모듈 v2 — Project / Agent 분리 DB + 관계형 연결

구조:
  Orchestrator — Projects  (project_id, 설명, Repo, 총비용)
        ↑ relation
  Orchestrator — Agents    (agent_id, Role, Model, Task, Result, 토큰, 비용, ...)

전제 조건:
  ~/.claude/secrets.json 에 "NOTION_TOKEN" (Notion Internal Integration 토큰)
"""
import json
import re
import time
import threading
import requests
from pathlib import Path

# 완료 대기가 필요한 스레드 추적 (subprocess 종료 전 flush용)
_pending_threads: list[threading.Thread] = []
_pending_lock = threading.Lock()

_NOTION_PAGE_ID = "31cf91ac-a96f-8107-819b-ef1d4b05900a"
_API = "https://api.notion.com/v1"
_CACHE = Path(__file__).parent / ".notion_cache.json"

_PROJECT_DB_TITLE = "Orchestrator — Projects"
_AGENT_DB_TITLE   = "Orchestrator — Agents"

# 대시보드에 이미 존재하는 DB ID (하드코딩으로 검색/생성 우회)
_PROJECT_DB_ID = "31cf91ac-a96f-819e-813d-e1f3dce80bf7"
_AGENT_DB_ID   = "31cf91ac-a96f-81ae-bc8e-fa6bac46b3b7"

_KNOWN_REPOS = [
    "UDH_Clustering", "PeakPicker", "biosteam-tagatose",
    "Kinetic-modeling", "claude-scientific-skills",
]


# ── 기본 헬퍼 ────────────────────────────────────────────────────────────────

def _get_token() -> str:
    try:
        return json.loads(
            (Path.home() / ".claude" / "secrets.json").read_text(encoding="utf-8")
        ).get("NOTION_TOKEN", "")
    except Exception:
        return ""


def _headers(token: str) -> dict:
    return {
        "Authorization": f"Bearer {token}",
        "Notion-Version": "2022-06-28",
        "Content-Type": "application/json",
    }


def _load_cache() -> dict:
    try:
        return json.loads(_CACHE.read_text(encoding="utf-8")) if _CACHE.exists() else {}
    except Exception:
        return {}


def _save_cache(cache: dict):
    try:
        _CACHE.write_text(json.dumps(cache, indent=2, ensure_ascii=False), encoding="utf-8")
    except Exception:
        pass


def _detect_repos(text: str) -> list[str]:
    return [r for r in _KNOWN_REPOS if r.lower() in text.lower()]


def _extract_github_issue_urls(text: str) -> list[str]:
    """텍스트에서 모든 GitHub issues/pull URL 추출 (중복 제거)"""
    urls = re.findall(r'https://github\.com/[^\s\)>"\']+/(?:issues|pull)/\d+', text)
    return list(dict.fromkeys(url.rstrip('.,') for url in urls))


def _make_github_rt(urls: list[str]) -> list[dict]:
    """URL 목록 → Notion rich_text 클릭 가능 링크 배열"""
    rt = []
    for i, url in enumerate(urls):
        num = url.rstrip('/').split('/')[-1]
        if i > 0:
            rt.append({"type": "text", "text": {"content": ", "}})
        rt.append({"type": "text", "text": {"content": f"#{num}", "link": {"url": url}},
                   "annotations": {"color": "blue"}})
    return rt


def _short_model(model: str) -> str:
    return model.replace("claude-", "").replace("-20251001", "")


# ── DB 생성 / 검색 ───────────────────────────────────────────────────────────

def _search_db(token: str, title: str) -> str | None:
    resp = requests.post(
        f"{_API}/search", headers=_headers(token),
        json={"query": title, "filter": {"value": "database", "property": "object"}},
        timeout=10,
    )
    pid_no_dash = _NOTION_PAGE_ID.replace("-", "")
    for r in resp.json().get("results", []):
        parent = r.get("parent", {})
        # 제목 정확 일치 + 부모 페이지 일치
        db_title = "".join(
            t.get("plain_text", "") for t in r.get("title", [])
        )
        if (parent.get("type") == "page_id"
                and parent.get("page_id", "").replace("-", "") == pid_no_dash
                and db_title == title):
            return r["id"]
    return None


def _create_project_db(token: str) -> str:  # noqa: C901
    schema = {
        "parent": {"type": "page_id", "page_id": _NOTION_PAGE_ID},
        "title": [{"type": "text", "text": {"content": _PROJECT_DB_TITLE}}],
        "properties": {
            "Name":           {"title": {}},
            "Project ID":     {"rich_text": {}},
            "Status":         {"select": {"options": [
                {"name": "active",  "color": "yellow"},
                {"name": "done",    "color": "green"},
                {"name": "failed",  "color": "red"},
            ]}},
            "Description":    {"rich_text": {}},
            "Agent Count":    {"number": {"format": "number"}},
            "Total Cost ($)": {"number": {"format": "dollar"}},
            "Repo":           {"multi_select": {"options": [
                {"name": r, "color": c} for r, c in zip(
                    _KNOWN_REPOS, ["blue", "green", "orange", "purple", "gray"]
                )
            ]}},
            "Created":        {"date": {}},
        },
    }
    r = requests.post(f"{_API}/databases", headers=_headers(token), json=schema, timeout=10)
    r.raise_for_status()
    return r.json()["id"]


def _create_agent_db(token: str, project_db_id: str) -> str:
    schema = {
        "parent": {"type": "page_id", "page_id": _NOTION_PAGE_ID},
        "title": [{"type": "text", "text": {"content": _AGENT_DB_TITLE}}],
        "properties": {
            "Name":        {"title": {}},
            "Project":     {"relation": {"database_id": project_db_id, "single_property": {}}},
            "Agent ID":    {"rich_text": {}},
            "Role":        {"select": {"options": [
                {"name": "lead",   "color": "blue"},
                {"name": "worker", "color": "yellow"},
                {"name": "solo",   "color": "purple"},
            ]}},
            "Team Type":   {"select": {"options": []}},
            "Model":       {"select": {"options": [
                {"name": "sonnet-4-6", "color": "blue"},
                {"name": "haiku-4-5",  "color": "green"},
                {"name": "opus-4-6",   "color": "purple"},
            ]}},
            "Event":       {"select": {"options": [
                {"name": "agent_start",  "color": "yellow"},
                {"name": "agent_done",   "color": "green"},
                {"name": "agent_failed", "color": "red"},
            ]}},
            "Status":      {"select": {"options": [
                {"name": "running", "color": "yellow"},
                {"name": "done",    "color": "green"},
                {"name": "failed",  "color": "red"},
                {"name": "stuck",   "color": "gray"},
            ]}},
            "Elapsed (s)": {"number": {"format": "number"}},
            "Tokens In":   {"number": {"format": "number_with_commas"}},
            "Tokens Out":  {"number": {"format": "number_with_commas"}},
            "Cost ($)":    {"number": {"format": "dollar"}},
            "Task":        {"rich_text": {}},
            "Result":      {"rich_text": {}},
            "Timestamp":   {"date": {}},
        },
    }
    r = requests.post(f"{_API}/databases", headers=_headers(token), json=schema, timeout=10)
    r.raise_for_status()
    return r.json()["id"]


def _ensure_dbs(token: str) -> tuple[str, str]:
    # 하드코딩된 DB ID 우선 사용 (검색/생성 불필요)
    return _PROJECT_DB_ID, _AGENT_DB_ID


# ── Project DB 관리 ──────────────────────────────────────────────────────────

def _get_or_create_project_notion_id(token: str, project_id: str,
                                      project_db_id: str, description: str = "") -> str | None:
    """project_id → Notion 페이지 ID (없으면 Project DB에 생성)"""
    cache = _load_cache()
    key = f"proj:{project_id}"
    notion_id = cache.get(key)
    if notion_id:
        return notion_id

    now_iso = time.strftime("%Y-%m-%dT%H:%M:%S+09:00", time.localtime())
    repos = _detect_repos(description)
    props = {
        "Name":        {"title": [{"text": {"content": f"{project_id} — {description[:60]}"}}]},
        "Project ID":  {"rich_text": [{"text": {"content": project_id}}]},
        "Status":      {"select": {"name": "active"}},
        "Description": {"rich_text": [{"text": {"content": description[:2000]}}]},
        "Created":     {"date": {"start": now_iso}},
    }
    if repos:
        props["Repo"] = {"multi_select": [{"name": r} for r in repos]}

    r = requests.post(f"{_API}/pages", headers=_headers(token),
                      json={"parent": {"database_id": project_db_id}, "properties": props},
                      timeout=10)
    if r.status_code == 200:
        notion_id = r.json()["id"]
        cache[key] = notion_id
        _save_cache(cache)
        return notion_id
    return None


def _update_project_cost(token: str, project_notion_id: str, cost_delta: float):
    """Project 페이지의 Total Cost 증가 (read → add → write)"""
    try:
        r = requests.get(f"{_API}/pages/{project_notion_id}", headers=_headers(token), timeout=10)
        current = r.json().get("properties", {}).get("Total Cost ($)", {}).get("number") or 0.0
        requests.patch(
            f"{_API}/pages/{project_notion_id}", headers=_headers(token),
            json={"properties": {"Total Cost ($)": {"number": round(current + cost_delta, 6)}}},
            timeout=10,
        )
    except Exception:
        pass


# ── Agent DB 로그 ────────────────────────────────────────────────────────────

def _log_agent_sync(event: str, name: str, agent_id: str = "", project_id: str = "",
                    team_type: str = "", role: str = "solo", model: str = "",
                    status: str = "", elapsed: float = None,
                    tokens_in: int = None, tokens_out: int = None,
                    cost: float = None, task: str = "", result: str = ""):
    try:
        token = _get_token()
        if not token:
            return

        project_db_id, agent_db_id = _ensure_dbs(token)
        now_iso = time.strftime("%Y-%m-%dT%H:%M:%S+09:00", time.localtime())

        # Project Notion 페이지 확보
        project_notion_id = None
        if project_id:
            project_notion_id = _get_or_create_project_notion_id(
                token, project_id, project_db_id, task or name
            )

        props: dict = {
            "Name":      {"title": [{"text": {"content": name[:100]}}]},
            "Event":     {"select": {"name": event}},
            "Timestamp": {"date": {"start": now_iso}},
        }
        if project_notion_id:
            props["Project"] = {"relation": [{"id": project_notion_id}]}
        if agent_id:
            props["Agent ID"]  = {"rich_text": [{"text": {"content": agent_id}}]}
        if role:
            props["Role"]      = {"select": {"name": role}}
        if team_type:
            props["Team Type"] = {"select": {"name": team_type[:100]}}
        if model:
            props["Model"]     = {"select": {"name": _short_model(model)[:100]}}
        if status:
            props["Status"]    = {"select": {"name": status}}
        if elapsed is not None:
            props["Elap. (s)"] = {"number": round(elapsed, 1)}
        if tokens_in is not None:
            props["Tokens In"]  = {"number": tokens_in}
        if tokens_out is not None:
            props["Tokens Out"] = {"number": tokens_out}
        if cost is not None:
            props["Cost ($)"]   = {"number": round(cost, 6)}
        if task:
            props["Task"]   = {"rich_text": [{"text": {"content": task[:2000]}}]}
        if result:
            props["Result"] = {"rich_text": [{"text": {"content": result[:2000]}}]}
        # result + task에서 모든 GitHub 이슈 URL 수집 → rich_text 링크로 저장
        gh_urls = _extract_github_issue_urls((result or "") + " " + (task or ""))
        if gh_urls:
            props["GitHub Issue"] = {"rich_text": _make_github_rt(gh_urls)}

        requests.post(f"{_API}/pages", headers=_headers(token),
                      json={"parent": {"database_id": agent_db_id}, "properties": props},
                      timeout=10)

        # 완료 시 Project 비용 누적
        if project_notion_id and event == "agent_done" and cost:
            _update_project_cost(token, project_notion_id, cost)

    except Exception:
        pass


def _log_project_sync(project_id: str, description: str):
    try:
        token = _get_token()
        if not token:
            return
        project_db_id, _ = _ensure_dbs(token)
        _get_or_create_project_notion_id(token, project_id, project_db_id, description)
    except Exception:
        pass


# ── Agent 상태 실시간 동기화 ──────────────────────────────────────────────────

def _update_agent_status_sync(agent_id: str, status: str, elapsed: float = 0.0):
    """
    기존 Agent 행의 Status와 Elapsed 시간을 Notion에 업데이트 (실시간 폴링용).
    agent_id → Notion Agent DB에서 검색 → Status 및 Elapsed 필드 업데이트.
    """
    try:
        token = _get_token()
        if not token:
            return

        _, agent_db_id = _ensure_dbs(token)

        # Agent DB에서 해당 agent_id 행 검색
        resp = requests.post(
            f"{_API}/databases/{agent_db_id}/query",
            headers=_headers(token),
            json={
                "filter": {
                    "property": "Agent ID",
                    "rich_text": {"equals": agent_id}
                }
            },
            timeout=10,
        )
        resp.raise_for_status()
        results = resp.json().get("results", [])

        if not results:
            return  # 아직 생성되지 않은 에이전트 (agent_start 이벤트 대기 중)

        page_id = results[0]["id"]

        # Status와 Elapsed 필드 업데이트
        props = {
            "Status": {"select": {"name": status}},
            "Elapsed (s)": {"number": round(elapsed, 1)},
        }
        requests.patch(
            f"{_API}/pages/{page_id}",
            headers=_headers(token),
            json={"properties": props},
            timeout=10,
        )
    except Exception:
        pass


# ── 공개 API ─────────────────────────────────────────────────────────────────

def wait_for_pending(timeout: float = 10.0):
    """subprocess 종료 전 대기 중인 Notion 로그 스레드를 모두 flush한다."""
    with _pending_lock:
        threads = list(_pending_threads)
    for t in threads:
        t.join(timeout=timeout)


def _track(t: threading.Thread) -> threading.Thread:
    """스레드를 _pending_threads에 등록하고, 완료 시 자동 제거한다."""
    with _pending_lock:
        _pending_threads.append(t)

    def _remove():
        with _pending_lock:
            try:
                _pending_threads.remove(t)
            except ValueError:
                pass

    wrapper = threading.Thread(target=lambda: (t.join(), _remove()), daemon=True)
    wrapper.start()
    return t


def log_project(project_id: str, description: str):
    """프로젝트 생성 로그 (비동기, Project DB 행 생성)"""
    t = threading.Thread(target=_log_project_sync, args=(project_id, description), daemon=False)
    _track(t)
    t.start()


def update_agent_status(agent_id: str, status: str, elapsed: float = 0.0):
    """
    에이전트 상태를 Notion에 실시간으로 업데이트 (비동기, 폴링용).
    get_status() 호출 시 각 팀마다 호출되어 running 상태의 경과 시간을 동기화.
    """
    t = threading.Thread(
        target=_update_agent_status_sync,
        args=(agent_id, status, elapsed),
        daemon=False,
    )
    _track(t)
    t.start()


def log_event(event: str, name: str, agent_id: str = "", project_id: str = "",
              team_type: str = "", role: str = "solo", model: str = "",
              status: str = "", elapsed: float = None,
              tokens_in: int = None, tokens_out: int = None,
              cost: float = None, task: str = "", result: str = "",
              summary: str = ""):    # summary: 이전 버전 호환 파라미터
    """에이전트 이벤트 로그 (비동기, Agent DB 행 생성)"""
    # 이전 버전 summary → task / result 자동 매핑
    if summary:
        if event in ("agent_start", "project_start") and not task:
            task = summary
        elif event in ("agent_done", "agent_failed") and not result:
            result = summary

    # project_start는 Project DB만 처리
    if event == "project_start":
        log_project(project_id, task or name)
        return

    t = threading.Thread(
        target=_log_agent_sync,
        kwargs=dict(event=event, name=name, agent_id=agent_id, project_id=project_id,
                    team_type=team_type, role=role, model=model, status=status,
                    elapsed=elapsed, tokens_in=tokens_in, tokens_out=tokens_out,
                    cost=cost, task=task, result=result),
        daemon=False,
    )
    _track(t)
    t.start()
