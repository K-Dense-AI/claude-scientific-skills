"""
SQLite 기반 경량 이벤트 버스.
cross-process 이벤트 공유 (server.py <-> ws_server.py).
"""
import json
import time
import sqlite3
from pathlib import Path
from typing import Callable, Optional

DB_PATH = str(Path(__file__).parent / "orchestrator.db")

def _get_conn():
    conn = sqlite3.connect(DB_PATH, check_same_thread=False)
    conn.execute("""
        CREATE TABLE IF NOT EXISTS events (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            event_type TEXT NOT NULL,
            team_id TEXT,
            project_id TEXT,
            data TEXT,
            created_at REAL NOT NULL
        )
    """)
    conn.execute("CREATE INDEX IF NOT EXISTS idx_events_created_at ON events(created_at)")
    conn.commit()
    return conn

def publish(event_type: str, team_id: str = None, project_id: str = None, data: dict = None):
    """이벤트 발행. team_manager.py에서 호출."""
    conn = _get_conn()
    try:
        conn.execute(
            "INSERT INTO events (event_type, team_id, project_id, data, created_at) VALUES (?,?,?,?,?)",
            (event_type, team_id, project_id, json.dumps(data or {}), time.time())
        )
        conn.commit()
    finally:
        conn.close()

def get_events_since(since_id: int = 0, limit: int = 50) -> list[dict]:
    """since_id 이후의 새 이벤트 반환. ws_server.py 폴링용."""
    conn = _get_conn()
    try:
        rows = conn.execute(
            "SELECT id, event_type, team_id, project_id, data, created_at FROM events WHERE id > ? ORDER BY id ASC LIMIT ?",
            (since_id, limit)
        ).fetchall()
        return [
            {"id": r[0], "event_type": r[1], "team_id": r[2], "project_id": r[3],
             "data": json.loads(r[4] or "{}"), "created_at": r[5]}
            for r in rows
        ]
    finally:
        conn.close()

def get_latest_id() -> int:
    """현재 최신 이벤트 ID 반환."""
    conn = _get_conn()
    try:
        row = conn.execute("SELECT MAX(id) FROM events").fetchone()
        return row[0] or 0
    finally:
        conn.close()

def cleanup_old_events(max_age_seconds: int = 3600):
    """1시간 이상 된 이벤트 삭제."""
    conn = _get_conn()
    try:
        conn.execute("DELETE FROM events WHERE created_at < ?", (time.time() - max_age_seconds,))
        conn.commit()
    finally:
        conn.close()
