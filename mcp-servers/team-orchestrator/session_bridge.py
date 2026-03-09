"""
session_bridge.py
session_id <-> Telegram message_id 매핑 + inbox 관리용 공유 SQLite DB.

사용:
    from session_bridge import save_session, get_session_by_tg_msg_id, add_inbox, get_and_clear_inbox
"""
import sqlite3
from datetime import datetime
from pathlib import Path

DB = Path(__file__).parent / "session_bridge.db"


def _get_conn() -> sqlite3.Connection:
    conn = sqlite3.connect(str(DB))
    conn.row_factory = sqlite3.Row
    conn.executescript("""
        CREATE TABLE IF NOT EXISTS sessions (
            session_id  TEXT,
            cwd         TEXT,
            tg_msg_id   INTEGER,
            team_id     TEXT,
            created_at  TEXT
        );
        CREATE UNIQUE INDEX IF NOT EXISTS idx_sessions_tg_msg
            ON sessions(tg_msg_id);
        CREATE TABLE IF NOT EXISTS inbox (
            id          INTEGER PRIMARY KEY AUTOINCREMENT,
            session_id  TEXT,
            message     TEXT,
            delivered   INTEGER DEFAULT 0,
            created_at  TEXT
        );
    """)
    # team_id 컬럼 마이그레이션 (기존 DB 호환)
    try:
        conn.execute("SELECT team_id FROM sessions LIMIT 1")
    except sqlite3.OperationalError:
        conn.execute("ALTER TABLE sessions ADD COLUMN team_id TEXT")
        conn.commit()
    return conn


def save_session(session_id: str, cwd: str, tg_msg_id: int, team_id: str = None):
    """Telegram 알림 전송 후 msg_id -> session_id + team_id 매핑 저장."""
    conn = _get_conn()
    conn.execute(
        "INSERT OR REPLACE INTO sessions (session_id, cwd, tg_msg_id, team_id, created_at) VALUES (?,?,?,?,?)",
        (session_id, cwd, tg_msg_id, team_id, datetime.now().isoformat()),
    )
    conn.commit()
    conn.close()


def get_session_by_tg_msg_id(tg_msg_id: int) -> dict | None:
    """Telegram message_id로 세션 정보 조회."""
    conn = _get_conn()
    row = conn.execute(
        "SELECT * FROM sessions WHERE tg_msg_id=?", (tg_msg_id,)
    ).fetchone()
    conn.close()
    return dict(row) if row else None


def add_inbox(session_id: str, message: str):
    """특정 세션의 inbox에 메시지 저장."""
    conn = _get_conn()
    conn.execute(
        "INSERT INTO inbox (session_id, message, created_at) VALUES (?,?,?)",
        (session_id, message, datetime.now().isoformat()),
    )
    conn.commit()
    conn.close()


def get_and_clear_inbox(session_id: str) -> list[dict]:
    """세션 inbox 조회 후 delivered=1 표시 (한 번만 전달)."""
    conn = _get_conn()
    rows = conn.execute(
        "SELECT * FROM inbox WHERE session_id=? AND delivered=0 ORDER BY created_at",
        (session_id,),
    ).fetchall()
    if rows:
        ids = [r["id"] for r in rows]
        conn.execute(
            f"UPDATE inbox SET delivered=1 WHERE id IN ({','.join('?'*len(ids))})",
            ids,
        )
        conn.commit()
    conn.close()
    return [dict(r) for r in rows]
