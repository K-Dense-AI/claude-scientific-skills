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
import os
import re
import subprocess
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


# ── 주간 토큰 한도 조회 (ccusage) ────────────────────────────────────────────

_ccusage_cache: dict = {}
_ccusage_cache_time: float = 0.0


def _get_weekly_usage() -> dict | None:
    """ccusage CLI로 현재 5h 블록의 주간 한도 대비 사용량을 조회한다.
    결과를 60초간 캐시한다.
    Returns: {totalTokens, limit, percentUsed, costUSD, projectedUsage} or None
    """
    global _ccusage_cache, _ccusage_cache_time
    if time.time() - _ccusage_cache_time < 60 and _ccusage_cache:
        return _ccusage_cache

    try:
        result = subprocess.run(
            "ccusage blocks --recent --offline --token-limit max --json",
            capture_output=True, text=True, timeout=30, shell=True,
        )
        if result.returncode != 0:
            return None
        data = json.loads(result.stdout)
        blocks = data.get("blocks", [])
        active = next((b for b in blocks if b.get("isActive")), None)
        if not active and blocks:
            active = blocks[-1]  # 가장 최근 블록 사용
        if not active:
            return None

        tls = active.get("tokenLimitStatus", {})
        info = {
            "totalTokens": active.get("totalTokens", 0),
            "limit": tls.get("limit", 0),
            "percentUsed": round(tls.get("percentUsed", 0), 1),
            "costUSD": round(active.get("costUSD", 0), 2),
            "projectedUsage": tls.get("projectedUsage", 0),
        }
        _ccusage_cache = info
        _ccusage_cache_time = time.time()
        return info
    except Exception:
        return None


def _format_tokens(n: int) -> str:
    """173025844 → '173M'"""
    if n >= 1_000_000:
        return f"{n / 1_000_000:.0f}M"
    if n >= 1_000:
        return f"{n / 1_000:.0f}K"
    return str(n)


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

    # 페이지 본문에 프로젝트 정보 블록 추가
    children = [
        {
            "object": "block",
            "type": "heading_2",
            "heading_2": {"rich_text": [{"type": "text", "text": {"content": "Project Info"}}]},
        },
        {
            "object": "block",
            "type": "bulleted_list_item",
            "bulleted_list_item": {"rich_text": [
                {"type": "text", "text": {"content": "Project ID: "}, "annotations": {"bold": True}},
                {"type": "text", "text": {"content": project_id}},
            ]},
        },
        {
            "object": "block",
            "type": "bulleted_list_item",
            "bulleted_list_item": {"rich_text": [
                {"type": "text", "text": {"content": "Description: "}, "annotations": {"bold": True}},
                {"type": "text", "text": {"content": description[:500]}},
            ]},
        },
        {
            "object": "block",
            "type": "bulleted_list_item",
            "bulleted_list_item": {"rich_text": [
                {"type": "text", "text": {"content": "Created: "}, "annotations": {"bold": True}},
                {"type": "text", "text": {"content": now_iso}},
            ]},
        },
        {
            "object": "block",
            "type": "divider",
            "divider": {},
        },
        {
            "object": "block",
            "type": "heading_2",
            "heading_2": {"rich_text": [{"type": "text", "text": {"content": "Agent Results"}}]},
        },
        {
            "object": "block",
            "type": "paragraph",
            "paragraph": {"rich_text": [
                {"type": "text", "text": {"content": "(Results will be appended as agents complete)"}, "annotations": {"italic": True}},
            ]},
        },
    ]

    r = requests.post(f"{_API}/pages", headers=_headers(token),
                      json={"parent": {"database_id": project_db_id}, "properties": props,
                            "children": children},
                      timeout=10)
    if r.status_code == 200:
        notion_id = r.json()["id"]
        cache[key] = notion_id
        _save_cache(cache)
        return notion_id
    return None


def _append_agent_result_to_page(token: str, project_notion_id: str,
                                 agent_name: str, agent_id: str, status: str,
                                 elapsed: float, cost: float, result: str):
    """에이전트 완료 시 Project 페이지 본문에 결과 블록을 추가"""
    try:
        elapsed_str = f"{elapsed:.0f}s" if elapsed else "N/A"
        cost_str = f"${cost:.4f}" if cost else "N/A"
        status_icon = {"done": "v", "failed": "x"}.get(status, "?")

        blocks = [
            {
                "object": "block",
                "type": "heading_3",
                "heading_3": {"rich_text": [
                    {"type": "text", "text": {"content": f"[{status_icon}] {agent_name} ({agent_id})"}},
                ]},
            },
            {
                "object": "block",
                "type": "bulleted_list_item",
                "bulleted_list_item": {"rich_text": [
                    {"type": "text", "text": {"content": f"Status: {status} | Elapsed: {elapsed_str} | Cost: {cost_str}"}},
                ]},
            },
        ]

        # 결과 텍스트를 최대 1500자까지 paragraph 블록으로 추가
        if result:
            result_text = result[:1500]
            blocks.append({
                "object": "block",
                "type": "paragraph",
                "paragraph": {"rich_text": [
                    {"type": "text", "text": {"content": result_text}},
                ]},
            })

        requests.patch(
            f"{_API}/blocks/{project_notion_id}/children",
            headers=_headers(token),
            json={"children": blocks},
            timeout=15,
        )
    except Exception:
        pass


def _update_project_metrics(token: str, project_notion_id: str,
                            cost_delta: float = 0.0, elapsed_delta: float = 0.0,
                            tokens_delta: int = 0):
    """Project 페이지의 Total Cost / Total Elapsed / Total Tokens / Tok/s 증가 (read → add → write)"""
    try:
        r = requests.get(f"{_API}/pages/{project_notion_id}", headers=_headers(token), timeout=10)
        props = r.json().get("properties", {})
        current_cost = props.get("Total Cost ($)", {}).get("number") or 0.0
        current_elapsed = props.get("Total Elapsed (s)", {}).get("number") or 0.0
        current_tokens = props.get("Total Tokens", {}).get("number") or 0

        new_cost = round(current_cost + cost_delta, 6)
        new_elapsed = round(current_elapsed + elapsed_delta, 1)
        new_tokens = int(current_tokens + tokens_delta)
        new_toks = round(new_tokens / new_elapsed, 1) if new_elapsed > 0 else 0

        # 주간 한도 대비 사용률 계산 (Notion percent 포맷: 0.01 = 1%)
        weekly_pct = None
        usage = _get_weekly_usage()
        if usage and usage.get("limit") and new_tokens > 0:
            weekly_pct = round(new_tokens / usage["limit"], 6)

        update_props = {
            "Total Cost ($)": {"number": new_cost},
            "Total Elapsed (s)": {"number": new_elapsed},
            "Total Tokens": {"number": new_tokens},
            "Tok/s": {"number": new_toks},
        }
        if weekly_pct is not None:
            update_props["Weekly %"] = {"number": weekly_pct}
        requests.patch(
            f"{_API}/pages/{project_notion_id}", headers=_headers(token),
            json={"properties": update_props},
            timeout=10,
        )
    except Exception:
        pass


# ── Agent DB 로그 ────────────────────────────────────────────────────────────

def _log_agent_sync(event: str, name: str, agent_id: str = "", project_id: str = "",
                    team_type: str = "", role: str = "solo", model: str = "",
                    status: str = "", elapsed: float = None,
                    tokens_in: int = None, tokens_out: int = None,
                    cost: float = None, task: str = "", result: str = "",
                    mcp_tools_used: dict = None):
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

        # MCP 도구 사용량 추가
        final_result = result or ""
        if mcp_tools_used and isinstance(mcp_tools_used, dict):
            mcp_summary = "\n\n[MCP Tools Used]"
            for tool, tokens in sorted(mcp_tools_used.items()):
                mcp_summary += f"\n  {tool}: {tokens:,} tokens"
            final_result = (final_result + mcp_summary)[:2000]

        if final_result:
            props["Result"] = {"rich_text": [{"text": {"content": final_result}}]}
        # result + task에서 모든 GitHub 이슈 URL 수집 → rich_text 링크로 저장
        gh_urls = _extract_github_issue_urls((final_result or "") + " " + (task or ""))
        if gh_urls:
            props["GitHub Issue"] = {"rich_text": _make_github_rt(gh_urls)}

        requests.post(f"{_API}/pages", headers=_headers(token),
                      json={"parent": {"database_id": agent_db_id}, "properties": props},
                      timeout=10)

        # 완료 시 Project 메트릭 누적 + 페이지 본문에 결과 추가 + 프로젝트 상태 업데이트
        if project_notion_id and event in ("agent_done", "agent_failed"):
            _update_project_metrics(
                token, project_notion_id,
                cost_delta=cost or 0.0,
                elapsed_delta=elapsed or 0.0,
                tokens_delta=(tokens_in or 0) + (tokens_out or 0),
            )
            _append_agent_result_to_page(
                token, project_notion_id,
                agent_name=name, agent_id=agent_id,
                status="done" if event == "agent_done" else "failed",
                elapsed=elapsed or 0.0, cost=cost or 0.0,
                result=result,
            )
            _update_project_agent_count(token, project_notion_id, project_id)

        # 대시보드 헤더 자동 갱신
        if event in ("agent_done", "agent_failed", "agent_start"):
            _update_dashboard_header(token)

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

        # Status와 Elapsed 필드 업데이트 (Notion 실제 필드명: "Elap. (s)")
        props = {
            "Status": {"select": {"name": status}},
            "Elap. (s)": {"number": round(elapsed, 1)},
        }
        requests.patch(
            f"{_API}/pages/{page_id}",
            headers=_headers(token),
            json={"properties": props},
            timeout=10,
        )
    except Exception:
        pass


# ── 대시보드 헤더 자동 갱신 ──────────────────────────────────────────────────

def _count_agents_by_status(token: str) -> dict[str, int]:
    """Agent DB를 쿼리해서 Status별 카운트를 반환"""
    counts = {"running": 0, "done": 0, "failed": 0}
    try:
        # done 에이전트 카운트 (agent_done 이벤트만)
        for status_val in ("running", "done", "failed"):
            resp = requests.post(
                f"{_API}/databases/{_AGENT_DB_ID}/query",
                headers=_headers(token),
                json={
                    "filter": {"property": "Status", "select": {"equals": status_val}},
                    "page_size": 1,
                },
                timeout=10,
            )
            # Notion은 has_more + 총 개수를 직접 안 줌 → 간단히 쿼리
            # page_size=100으로 여러 번 페이징하는 대신, 마지막 업데이트 기준으로 추정
        # 더 효율적: filter 없이 전체를 한 번에 (최대 100)
        all_agents = []
        has_more = True
        start_cursor = None
        while has_more:
            body: dict = {"page_size": 100}
            if start_cursor:
                body["start_cursor"] = start_cursor
            resp = requests.post(
                f"{_API}/databases/{_AGENT_DB_ID}/query",
                headers=_headers(token), json=body, timeout=15,
            )
            data = resp.json()
            all_agents.extend(data.get("results", []))
            has_more = data.get("has_more", False)
            start_cursor = data.get("next_cursor")
            if len(all_agents) > 500:
                break  # 안전 제한

        for a in all_agents:
            s = a.get("properties", {}).get("Status", {}).get("select", {})
            sname = s.get("name", "") if s else ""
            if sname in counts:
                counts[sname] += 1
    except Exception:
        pass
    return counts


def _update_dashboard_header(token: str):
    """대시보드 페이지의 헤더 카운트를 갱신"""
    try:
        counts = _count_agents_by_status(token)
        now_str = time.strftime("%Y-%m-%d %H:%M", time.localtime())
        running = counts.get("running", 0)
        done = counts.get("done", 0)
        failed = counts.get("failed", 0)

        # Notion 페이지 내용에서 기존 헤더 찾아 교체
        new_header = (
            f"**🔄 마지막 업데이트**: {now_str} \\| "
            f"**실행 중**: {running}개 \\| "
            f"**완료**: {done}개 \\| "
            f"**실패**: {failed}개"
        )

        # Notion MCP update_content는 사용 불가 (동기 HTTP 직접 호출)
        # 대신 REST API로 페이지 블록을 업데이트
        # 페이지의 첫 번째 자식 블록(callout/quote)을 찾아서 교체
        resp = requests.get(
            f"{_API}/blocks/{_NOTION_PAGE_ID}/children?page_size=5",
            headers=_headers(token), timeout=10,
        )
        blocks = resp.json().get("results", [])

        for block in blocks:
            # quote 블록(>) 찾기 — 헤더가 여기 있음
            if block.get("type") == "quote":
                block_id = block["id"]
                # quote 블록은 rich_text 교체 (한 번에 덮어쓰기)
                rt = []
                parts = [
                    ("🔄 마지막 업데이트", False), (f": {now_str}", False),
                    (" | ", False),
                    ("실행 중", True), (f": {running}개", False),
                    (" | ", False),
                    ("완료", True), (f": {done}개", False),
                    (" | ", False),
                    ("실패", True), (f": {failed}개", False),
                ]
                # 주간 한도 사용량 추가
                usage = _get_weekly_usage()
                if usage and usage.get("limit"):
                    pct = usage["percentUsed"]
                    total = _format_tokens(usage["totalTokens"])
                    limit = _format_tokens(usage["limit"])
                    cost = usage["costUSD"]
                    parts.append((" | ", False))
                    parts.append(("📊 한도", True))
                    parts.append((f": {total}/{limit} ({pct}%) ${cost}", False))
                for text, bold in parts:
                    item = {"type": "text", "text": {"content": text}}
                    if bold:
                        item["annotations"] = {"bold": True}
                    rt.append(item)

                requests.patch(
                    f"{_API}/blocks/{block_id}",
                    headers=_headers(token),
                    json={"quote": {"rich_text": rt}},
                    timeout=10,
                )
                break
    except Exception:
        pass


def _update_project_agent_count(token: str, project_notion_id: str, project_id: str):
    """Project의 Agent Count와 Status를 업데이트 (전체 에이전트 완료 시 done으로)"""
    try:
        # 해당 프로젝트의 에이전트를 Project relation으로 필터 (agent_start + agent_done 모두 포함)
        resp = requests.post(
            f"{_API}/databases/{_AGENT_DB_ID}/query",
            headers=_headers(token),
            json={
                "filter": {
                    "property": "Project",
                    "relation": {"contains": project_notion_id}
                },
                "page_size": 100,
            },
            timeout=10,
        )
        agents = resp.json().get("results", [])

        # agent_id 기준으로 최신 이벤트만 추적 (agent_start → agent_done 순서)
        latest_by_agent: dict[str, str] = {}
        for a in agents:
            a_props = a.get("properties", {})
            aid_rt = a_props.get("Agent ID", {}).get("rich_text", [])
            aid = aid_rt[0].get("plain_text", "") if aid_rt else ""
            event_sel = a_props.get("Event", {}).get("select", {}) or {}
            event_name = event_sel.get("name", "")
            if aid:
                # agent_done/agent_failed가 agent_start보다 우선
                if event_name in ("agent_done", "agent_failed"):
                    latest_by_agent[aid] = event_name
                elif aid not in latest_by_agent:
                    latest_by_agent[aid] = event_name

        agent_count = len(latest_by_agent)
        props: dict = {"Agent Count": {"number": agent_count}}

        # 모든 고유 에이전트가 done/failed인지 확인
        all_terminal = (
            agent_count > 0
            and all(ev in ("agent_done", "agent_failed") for ev in latest_by_agent.values())
        )

        if all_terminal:
            props["Status"] = {"select": {"name": "done"}}

        requests.patch(
            f"{_API}/pages/{project_notion_id}",
            headers=_headers(token),
            json={"properties": props},
            timeout=10,
        )
    except Exception:
        pass


# ── 자동 폴링 (백그라운드 Notion Elapsed 동기화) ─────────────────────────────

_auto_sync_thread: threading.Thread | None = None
_auto_sync_stop = threading.Event()


def _auto_sync_loop(interval: float = 30.0):
    """
    백그라운드 스레드: running 에이전트의 Elapsed를 interval초마다 Notion에 동기화.
    get_status() 수동 호출 없이도 Notion DB가 실시간 갱신됨.
    """
    import importlib
    while not _auto_sync_stop.is_set():
        _auto_sync_stop.wait(interval)
        if _auto_sync_stop.is_set():
            break
        try:
            # state 모듈을 늦게 import (순환 import 방지)
            st = importlib.import_module("state")
            current_time = time.time()
            teams = st.list_teams()
            for t in teams:
                if t.status == "running" and t.started_at > 0:
                    elapsed = current_time - t.started_at
                    _update_agent_status_sync(t.team_id, "running", elapsed)
            # 대시보드 헤더도 갱신
            token = _get_token()
            if token:
                _update_dashboard_header(token)
        except Exception:
            pass


def start_auto_sync(interval: float = 30.0):
    """자동 동기화 시작 (MCP 서버 시작 시 호출)"""
    global _auto_sync_thread
    if _auto_sync_thread and _auto_sync_thread.is_alive():
        return  # 이미 실행 중
    _auto_sync_stop.clear()
    _auto_sync_thread = threading.Thread(
        target=_auto_sync_loop, args=(interval,), daemon=True
    )
    _auto_sync_thread.start()


def stop_auto_sync():
    """자동 동기화 중지"""
    _auto_sync_stop.set()
    if _auto_sync_thread:
        _auto_sync_thread.join(timeout=5)


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
    t.start()
    _track(t)


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
    t.start()
    _track(t)


def log_event(event: str, name: str, agent_id: str = "", project_id: str = "",
              team_type: str = "", role: str = "solo", model: str = "",
              status: str = "", elapsed: float = None,
              tokens_in: int = None, tokens_out: int = None,
              cost: float = None, task: str = "", result: str = "",
              summary: str = "", mcp_tools_used: dict = None):    # summary: 이전 버전 호환 파라미터
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
                    cost=cost, task=task, result=result, mcp_tools_used=mcp_tools_used),
        daemon=False,
    )
    t.start()
    _track(t)
