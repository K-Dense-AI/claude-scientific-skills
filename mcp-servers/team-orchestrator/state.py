"""SQLite 기반 팀/태스크/메시지 상태 관리"""
import sqlite3
import json
import time
from pathlib import Path
from dataclasses import dataclass, asdict
from typing import Optional

DB_PATH = Path(__file__).parent / "orchestrator.db"


@dataclass
class Team:
    team_id: str
    name: str
    team_type: str          # research | code | writing | ops
    status: str             # pending | running | paused | done | failed
    task: str
    project_id: str
    created_at: float
    output: str = ""
    error: str = ""
    depends_on: str = "[]"  # JSON list of team names this team depends on
    input_tokens: int = 0
    output_tokens: int = 0
    model_used: str = ""
    started_at: float = 0.0
    elapsed_seconds: float = 0.0
    total_cost_usd: float = 0.0


@dataclass
class Message:
    msg_id: int
    team_id: str
    role: str               # ceo | team
    content: str
    timestamp: float
    consumed: bool = False  # 팀이 아직 읽지 않은 CEO 피드백


def get_conn() -> sqlite3.Connection:
    conn = sqlite3.connect(DB_PATH, check_same_thread=False)
    conn.row_factory = sqlite3.Row
    return conn


def init_db():
    with get_conn() as conn:
        # 기존 DB에 컬럼 추가 (없는 경우만)
        for col_sql in [
            "ALTER TABLE teams ADD COLUMN depends_on TEXT DEFAULT '[]'",
            "ALTER TABLE teams ADD COLUMN input_tokens INTEGER DEFAULT 0",
            "ALTER TABLE teams ADD COLUMN output_tokens INTEGER DEFAULT 0",
            "ALTER TABLE teams ADD COLUMN model_used TEXT DEFAULT ''",
            "ALTER TABLE teams ADD COLUMN started_at REAL DEFAULT 0",
            "ALTER TABLE teams ADD COLUMN elapsed_seconds REAL DEFAULT 0",
            "ALTER TABLE teams ADD COLUMN total_cost_usd REAL DEFAULT 0",
        ]:
            try:
                conn.execute(col_sql)
            except Exception:
                pass  # 이미 존재하면 무시
        conn.executescript("""
        CREATE TABLE IF NOT EXISTS teams (
            team_id         TEXT PRIMARY KEY,
            name            TEXT,
            team_type       TEXT,
            status          TEXT,
            task            TEXT,
            project_id      TEXT,
            created_at      REAL,
            output          TEXT DEFAULT '',
            error           TEXT DEFAULT '',
            depends_on      TEXT DEFAULT '[]',
            input_tokens    INTEGER DEFAULT 0,
            output_tokens   INTEGER DEFAULT 0,
            model_used      TEXT DEFAULT '',
            started_at      REAL DEFAULT 0,
            elapsed_seconds REAL DEFAULT 0,
            total_cost_usd  REAL DEFAULT 0
        );

        CREATE TABLE IF NOT EXISTS projects (
            project_id  TEXT PRIMARY KEY,
            description TEXT,
            created_at  REAL,
            status      TEXT DEFAULT 'active'
        );

        CREATE TABLE IF NOT EXISTS messages (
            msg_id    INTEGER PRIMARY KEY AUTOINCREMENT,
            team_id   TEXT,
            role      TEXT,
            content   TEXT,
            timestamp REAL,
            consumed  INTEGER DEFAULT 0
        );

        CREATE TABLE IF NOT EXISTS notifications (
            id        INTEGER PRIMARY KEY AUTOINCREMENT,
            team_id   TEXT,
            message   TEXT,
            timestamp REAL,
            read      INTEGER DEFAULT 0
        );
        """)


# ── Team CRUD ──────────────────────────────────────────────────────────────

def upsert_team(team: Team):
    with get_conn() as conn:
        conn.execute("""
        INSERT OR REPLACE INTO teams
            (team_id, name, team_type, status, task, project_id, created_at, output, error, depends_on,
             input_tokens, output_tokens, model_used, started_at, elapsed_seconds, total_cost_usd)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (team.team_id, team.name, team.team_type, team.status,
              team.task, team.project_id, team.created_at, team.output, team.error, team.depends_on,
              team.input_tokens, team.output_tokens, team.model_used,
              team.started_at, team.elapsed_seconds, team.total_cost_usd))


def get_team(team_id: str) -> Optional[Team]:
    with get_conn() as conn:
        row = conn.execute("SELECT * FROM teams WHERE team_id=?", (team_id,)).fetchone()
        if not row:
            return None
        return Team(**dict(row))


def list_teams(project_id: str = None) -> list[Team]:
    with get_conn() as conn:
        if project_id:
            rows = conn.execute("SELECT * FROM teams WHERE project_id=?", (project_id,)).fetchall()
        else:
            rows = conn.execute("SELECT * FROM teams ORDER BY created_at DESC").fetchall()
        return [Team(**dict(r)) for r in rows]


def update_team_status(team_id: str, status: str, output: str = None, error: str = None):
    with get_conn() as conn:
        if output is not None:
            conn.execute("UPDATE teams SET status=?, output=? WHERE team_id=?", (status, output, team_id))
        elif error is not None:
            conn.execute("UPDATE teams SET status=?, error=? WHERE team_id=?", (status, error, team_id))
        else:
            conn.execute("UPDATE teams SET status=? WHERE team_id=?", (status, team_id))


def append_team_output(team_id: str, chunk: str):
    with get_conn() as conn:
        conn.execute("UPDATE teams SET output = output || ? WHERE team_id=?", (chunk, team_id))


def get_team_output(team_id: str) -> str:
    with get_conn() as conn:
        row = conn.execute("SELECT output FROM teams WHERE team_id=?", (team_id,)).fetchone()
        return row["output"] if row else ""


def update_team_usage(team_id: str, input_tokens: int, output_tokens: int,
                      model_used: str, elapsed_seconds: float, total_cost_usd: float = 0.0):
    """에이전트 완료 시 토큰·모델·시간 정보를 저장"""
    with get_conn() as conn:
        conn.execute("""
        UPDATE teams SET input_tokens=?, output_tokens=?, model_used=?,
                         elapsed_seconds=?, total_cost_usd=?
        WHERE team_id=?
        """, (input_tokens, output_tokens, model_used, elapsed_seconds, total_cost_usd, team_id))


def set_team_started_at(team_id: str, started_at: float):
    """에이전트 실행 시작 시간 기록"""
    with get_conn() as conn:
        conn.execute("UPDATE teams SET started_at=? WHERE team_id=?", (started_at, team_id))


def cleanup_stale_agents(max_age_hours: float = 2.0) -> list[str]:
    """
    오래된 running/pending 에이전트를 'stuck'으로 업데이트하고 정리.
    max_age_hours: 이 시간 이상 running 상태로 있으면 stale로 처리 (기본 2시간)
    """
    cutoff = time.time() - max_age_hours * 3600
    with get_conn() as conn:
        rows = conn.execute("""
        SELECT team_id, name, status, created_at FROM teams
        WHERE status IN ('running', 'pending')
          AND created_at < ?
        """, (cutoff,)).fetchall()
        stale_ids = [r["team_id"] for r in rows]
        if stale_ids:
            placeholders = ",".join("?" * len(stale_ids))
            conn.execute(
                f"UPDATE teams SET status='stuck' WHERE team_id IN ({placeholders})",
                stale_ids
            )
    return stale_ids


# ── Project CRUD ───────────────────────────────────────────────────────────

def create_project(project_id: str, description: str):
    with get_conn() as conn:
        conn.execute("""
        INSERT OR REPLACE INTO projects (project_id, description, created_at, status)
        VALUES (?, ?, ?, 'active')
        """, (project_id, description, time.time()))


def get_project(project_id: str) -> Optional[dict]:
    with get_conn() as conn:
        row = conn.execute("SELECT * FROM projects WHERE project_id=?", (project_id,)).fetchone()
        return dict(row) if row else None


# ── Message Queue ──────────────────────────────────────────────────────────

def push_message(team_id: str, role: str, content: str):
    with get_conn() as conn:
        conn.execute("""
        INSERT INTO messages (team_id, role, content, timestamp, consumed)
        VALUES (?, ?, ?, ?, 0)
        """, (team_id, role, content, time.time()))


def pop_ceo_messages(team_id: str) -> list[str]:
    """팀이 아직 읽지 않은 CEO 피드백 가져오고 consumed 표시"""
    with get_conn() as conn:
        rows = conn.execute("""
        SELECT msg_id, content FROM messages
        WHERE team_id=? AND role='ceo' AND consumed=0
        ORDER BY timestamp ASC
        """, (team_id,)).fetchall()

        if rows:
            ids = [r["msg_id"] for r in rows]
            conn.execute(
                f"UPDATE messages SET consumed=1 WHERE msg_id IN ({','.join('?'*len(ids))})",
                ids
            )
        return [r["content"] for r in rows]


# ── Notifications ──────────────────────────────────────────────────────────

def add_notification(team_id: str, message: str):
    """에이전트 완료/실패 시 CEO 알림 추가"""
    with get_conn() as conn:
        conn.execute(
            "INSERT INTO notifications (team_id, message, timestamp) VALUES (?, ?, ?)",
            (team_id, message, time.time())
        )


def get_pending_notifications() -> list[dict]:
    """읽지 않은 알림 반환 후 자동 읽음 처리"""
    with get_conn() as conn:
        rows = conn.execute(
            "SELECT id, team_id, message, timestamp FROM notifications "
            "WHERE read=0 ORDER BY timestamp ASC"
        ).fetchall()
        if rows:
            ids = [r["id"] for r in rows]
            conn.execute(
                f"UPDATE notifications SET read=1 WHERE id IN ({','.join('?'*len(ids))})",
                ids
            )
        return [dict(r) for r in rows]


def get_conversation(team_id: str) -> list[dict]:
    """팀의 전체 대화 이력 (Claude API messages 형식)"""
    with get_conn() as conn:
        rows = conn.execute("""
        SELECT role, content FROM messages
        WHERE team_id=? ORDER BY timestamp ASC
        """, (team_id,)).fetchall()
        result = []
        for r in rows:
            api_role = "user" if r["role"] in ("ceo", "user") else "assistant"
            result.append({"role": api_role, "content": r["content"]})
        return result
