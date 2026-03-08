#!/usr/bin/env python3
"""
공유 메시지 보드 MCP 서버
- 병렬 실행 중인 모든 에이전트가 이 서버에 연결해 메시지를 주고받음
- 여러 인스턴스가 동시에 실행되어도 SQLite WAL 모드로 안전하게 공유
"""
import sys
import json
import time
import sqlite3
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from mcp.server.fastmcp import FastMCP

DB_PATH = Path(__file__).parent / "board.db"

mcp = FastMCP("shared-board")


def get_conn() -> sqlite3.Connection:
    conn = sqlite3.connect(str(DB_PATH), check_same_thread=False, timeout=10)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA busy_timeout=5000")
    conn.row_factory = sqlite3.Row
    return conn


def init_board():
    with get_conn() as conn:
        conn.executescript("""
        CREATE TABLE IF NOT EXISTS board_messages (
            msg_id     INTEGER PRIMARY KEY AUTOINCREMENT,
            from_agent TEXT NOT NULL,
            to_agent   TEXT,
            content    TEXT NOT NULL,
            timestamp  REAL NOT NULL,
            read_by    TEXT DEFAULT '[]'
        );

        CREATE TABLE IF NOT EXISTS active_agents (
            agent_id   TEXT PRIMARY KEY,
            team_type  TEXT,
            project_id TEXT,
            last_seen  REAL
        );
        """)


init_board()


@mcp.tool()
def register_agent(agent_id: str, team_type: str, project_id: str) -> str:
    """
    작업 시작 시 자신을 보드에 등록합니다. 반드시 가장 먼저 호출하세요.

    Args:
        agent_id: 이 에이전트의 고유 ID (팀 ID와 동일)
        team_type: 팀 유형 (research | code | writing | ops)
        project_id: 프로젝트 ID
    """
    with get_conn() as conn:
        conn.execute("""
        INSERT OR REPLACE INTO active_agents (agent_id, team_type, project_id, last_seen)
        VALUES (?, ?, ?, ?)
        """, (agent_id, team_type, project_id, time.time()))
    return f"[등록 완료] {agent_id} ({team_type}팀) 보드에 등록됨."


@mcp.tool()
def heartbeat(agent_id: str) -> str:
    """
    활성 상태를 갱신합니다. 긴 작업 중 주기적으로 호출하세요.

    Args:
        agent_id: 이 에이전트의 ID
    """
    with get_conn() as conn:
        conn.execute(
            "UPDATE active_agents SET last_seen=? WHERE agent_id=?",
            (time.time(), agent_id)
        )
    return "heartbeat ok"


@mcp.tool()
def post_message(from_agent: str, content: str, to_agent: str = None) -> str:
    """
    다른 에이전트 또는 전체에게 메시지를 전송합니다.
    질문, 중간 결과 공유, 답변 모두 이 도구로 전송하세요.

    Args:
        from_agent: 보내는 에이전트 ID (자신의 team_id)
        content: 전송할 내용
        to_agent: 수신 에이전트 ID. None이면 전체 브로드캐스트.
    """
    with get_conn() as conn:
        conn.execute("""
        INSERT INTO board_messages (from_agent, to_agent, content, timestamp)
        VALUES (?, ?, ?, ?)
        """, (from_agent, to_agent, content, time.time()))

    target = to_agent if to_agent else "전체"
    return f"[전송 완료] {from_agent} → {target}"


@mcp.tool()
def read_inbox(agent_id: str) -> str:
    """
    내 inbox의 읽지 않은 메시지를 모두 가져옵니다.
    나에게 직접 보낸 메시지 + 전체 브로드캐스트를 포함합니다.
    작업 중 주기적으로 호출해 다른 팀의 질문이나 공유 내용을 확인하세요.

    Args:
        agent_id: 이 에이전트의 ID
    """
    with get_conn() as conn:
        rows = conn.execute("""
        SELECT msg_id, from_agent, to_agent, content, read_by
        FROM board_messages
        WHERE (to_agent = ? OR to_agent IS NULL)
          AND from_agent != ?
        ORDER BY timestamp ASC
        """, (agent_id, agent_id)).fetchall()

        unread = [r for r in rows if agent_id not in json.loads(r["read_by"] or "[]")]

        if not unread:
            return "[inbox 비어있음]"

        for row in unread:
            read_by = json.loads(row["read_by"] or "[]")
            read_by.append(agent_id)
            conn.execute(
                "UPDATE board_messages SET read_by=? WHERE msg_id=?",
                (json.dumps(read_by), row["msg_id"])
            )

    lines = [f"[inbox: {len(unread)}개 새 메시지]\n"]
    for row in unread:
        target = f"→ {row['to_agent']}" if row["to_agent"] else "→ 전체"
        lines.append(f"[{row['from_agent']} {target}]\n{row['content']}\n")
    return "\n".join(lines)


@mcp.tool()
def list_active_agents(project_id: str = None) -> str:
    """
    현재 활성 에이전트 목록을 반환합니다.
    메시지를 보낼 대상 agent_id를 확인할 때 사용하세요.

    Args:
        project_id: 특정 프로젝트만 조회. None이면 전체.
    """
    cutoff = time.time() - 600
    with get_conn() as conn:
        if project_id:
            rows = conn.execute("""
            SELECT agent_id, team_type, last_seen FROM active_agents
            WHERE project_id=? AND last_seen > ?
            """, (project_id, cutoff)).fetchall()
        else:
            rows = conn.execute("""
            SELECT agent_id, team_type, last_seen FROM active_agents
            WHERE last_seen > ?
            """, (cutoff,)).fetchall()

    if not rows:
        return "활성 에이전트 없음."

    lines = [f"[활성 에이전트: {len(rows)}개]"]
    for r in rows:
        age = int(time.time() - r["last_seen"])
        lines.append(f"  {r['agent_id']} ({r['team_type']}팀) - {age}초 전 활성")
    return "\n".join(lines)


@mcp.tool()
def wait_for_workers(lead_id: str, worker_ids: list, timeout_per_check: int = 20, max_checks: int = 30) -> str:
    """
    리드 에이전트가 워커 완료를 대기합니다. 내부적으로 sleep을 수행합니다.
    모든 워커가 "[완료 보고]"를 보낼 때까지 주기적으로 inbox를 확인합니다.

    Args:
        lead_id: 리드 에이전트의 team_id
        worker_ids: 대기할 워커 ID 목록
        timeout_per_check: 각 체크 사이 대기 시간(초). 기본 20초.
        max_checks: 최대 체크 횟수. 기본 30회.

    Returns:
        수집된 워커 결과 목록 (완료된 것만)
    """
    import time as _time

    pending = set(worker_ids)
    collected = {}

    for check_num in range(max_checks):
        if not pending:
            break

        # inbox 확인
        with get_conn() as conn:
            rows = conn.execute("""
                SELECT msg_id, from_agent, content, read_by
                FROM board_messages
                WHERE (to_agent = ? OR to_agent IS NULL)
                  AND from_agent != ?
                ORDER BY timestamp ASC
            """, (lead_id, lead_id)).fetchall()

            unread = [r for r in rows if lead_id not in json.loads(r["read_by"] or "[]")]

            for row in unread:
                # read 표시
                read_by = json.loads(row["read_by"] or "[]")
                read_by.append(lead_id)
                conn.execute(
                    "UPDATE board_messages SET read_by=? WHERE msg_id=?",
                    (json.dumps(read_by), row["msg_id"])
                )
                # 완료 보고 감지
                if "[완료 보고]" in row["content"] and row["from_agent"] in pending:
                    collected[row["from_agent"]] = row["content"]
                    pending.discard(row["from_agent"])

        if not pending:
            break

        if check_num < max_checks - 1:
            _time.sleep(timeout_per_check)

    lines = [f"[대기 완료] {len(collected)}/{len(worker_ids)}명 워커 결과 수집\n"]
    for wid, content in collected.items():
        lines.append(f"--- {wid} ---\n{content[:1000]}\n")
    if pending:
        lines.append(f"[미완료 워커] {', '.join(pending)} — 시간 초과로 현재까지 결과로 진행")
    return "\n".join(lines)


@mcp.tool()
def spawn_worker(lead_id: str, sub_task: str, worker_type: str = None) -> str:
    """
    리드 에이전트가 워커 서브에이전트를 생성합니다.

    Args:
        lead_id: 리드 에이전트의 team_id
        sub_task: 워커에게 맡길 서브태스크 내용
        worker_type: 워커 유형 (None이면 리드와 같은 유형)

    Returns:
        생성된 워커 ID
    """
    import os
    import subprocess
    import tempfile
    import time
    import json

    DETACHED_PROCESS = 0x00000008
    CREATE_NEW_PROCESS_GROUP = 0x00000200
    CREATE_NO_WINDOW = 0x08000000
    PYTHON_BIN = sys.executable
    RUN_AGENT_SCRIPT = str(Path(__file__).parent / "run_agent.py")

    import state as st
    st.init_db()

    lead = st.get_team(lead_id)
    if not lead:
        return f"[오류] 리드 {lead_id}를 찾을 수 없습니다."

    all_teams = st.list_teams(lead.project_id)
    worker_count = sum(1 for t in all_teams if t.team_id.startswith(f"{lead_id}-w"))
    worker_num = worker_count + 1
    worker_id = f"{lead_id}-w{worker_num}"
    wtype = worker_type or lead.team_type
    worker_name = f"{lead.name}-W{worker_num}"

    worker = st.Team(
        team_id=worker_id, name=worker_name, team_type=wtype,
        status="pending", task=sub_task, project_id=lead.project_id,
        created_at=time.time(), depends_on=json.dumps([])
    )
    st.upsert_team(worker)

    tmp = tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", suffix=".txt",
        prefix=f"task_{worker_id}_", delete=False
    )
    tmp.write(sub_task)
    tmp.close()

    clean_env = {k: v for k, v in os.environ.items() if k != "CLAUDECODE"}
    workdir = str(Path.home())

    proc = subprocess.Popen(
        [PYTHON_BIN, RUN_AGENT_SCRIPT,
         worker_id, wtype, lead.project_id, workdir, tmp.name, "worker", lead_id],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        env=clean_env,
        creationflags=DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP | CREATE_NO_WINDOW,
        close_fds=True,
    )

    st.append_team_output(lead_id, f"[워커 생성: {worker_id}, PID={proc.pid}]\n")
    return f"워커 생성 완료: {worker_name} (ID: `{worker_id}`). read_inbox로 완료 보고를 기다리세요."


if __name__ == "__main__":
    mcp.run()
