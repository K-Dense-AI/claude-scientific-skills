#!/usr/bin/env python3
"""
Team Orchestrator Monitor
Usage: python monitor.py [project_id] [interval_sec] [--no-board]
"""
import sys
import time
import sqlite3
from pathlib import Path

DIR = Path(__file__).parent
sys.path.insert(0, str(DIR))

import state as st

ORCH_DB = DIR / "orchestrator.db"
BOARD_DB = DIR / "board.db"

STATUS_ICON = {
    "running": ">>", "done": "OK", "failed": "!!", "pending": "..", "paused": "||",
}
TYPE_LABEL = {
    "research": "RESEARCH", "code": "CODE   ", "writing": "WRITING", "ops": "OPS    ",
}


def _safe(text: str, limit: int = 70) -> str:
    return text[:limit].encode("ascii", errors="replace").decode()


def get_board_messages(limit: int = 8) -> list:
    if not BOARD_DB.exists():
        return []
    try:
        conn = sqlite3.connect(str(BOARD_DB))
        conn.row_factory = sqlite3.Row
        rows = conn.execute("""
            SELECT from_agent, to_agent, content, timestamp
            FROM board_messages ORDER BY timestamp DESC LIMIT ?
        """, (limit,)).fetchall()
        conn.close()
        return [dict(r) for r in reversed(rows)]
    except Exception:
        return []


def get_active_agents() -> dict:
    if not BOARD_DB.exists():
        return {}
    try:
        conn = sqlite3.connect(str(BOARD_DB))
        conn.row_factory = sqlite3.Row
        cutoff = time.time() - 120
        rows = conn.execute("""
            SELECT agent_id, team_type, last_seen FROM active_agents
            WHERE last_seen > ? ORDER BY last_seen DESC
        """, (cutoff,)).fetchall()
        conn.close()
        return {r["agent_id"]: dict(r) for r in rows}
    except Exception:
        return {}


def show_status(project_id=None, show_board=True):
    teams = st.list_teams(project_id)
    now_str = time.strftime("%H:%M:%S")
    active_agents = get_active_agents()

    sep = "=" * 65
    dash = "-" * 65
    print(f"\n{sep}")
    proj_label = f"  project: {project_id}" if project_id else "  (all projects)"
    print(f"  Team Status  [{now_str}]{proj_label}")
    print(sep)

    if not teams:
        print("  (no active teams)")
    else:
        for t in teams:
            icon = STATUS_ICON.get(t.status, "??")
            type_label = TYPE_LABEL.get(t.team_type, t.team_type.upper().ljust(7))
            output_len = len(t.output) if t.output else 0
            alive = t.team_id in active_agents
            alive_mark = " [ACTIVE]" if alive else ""

            last_line = ""
            if t.output:
                lines = [l for l in t.output.strip().split("\n") if l.strip()]
                if lines:
                    last_line = _safe(lines[-1], 55)

            print(f"  [{icon}] [{type_label}] {_safe(t.name, 20)}  ({t.team_id}){alive_mark}")
            print(f"       status: {t.status:<10} | output: {output_len:>7,} chars")
            if last_line:
                print(f"       last:   {last_line}")
            print()

    # Board messages
    if show_board:
        msgs = get_board_messages(limit=8)
        if msgs:
            print(dash)
            print(f"  Board Messages (last {len(msgs)})")
            print(dash)
            for m in msgs:
                ts = time.strftime("%H:%M:%S", time.localtime(m["timestamp"]))
                to = f"-> {m['to_agent']}" if m["to_agent"] else "-> ALL"
                content = _safe(m["content"].replace("\n", " "), 65)
                print(f"  [{ts}] {_safe(m['from_agent'],12)} {to}")
                print(f"         {content}")
            print()

    # Summary
    if teams:
        counts = {}
        for t in teams:
            counts[t.status] = counts.get(t.status, 0) + 1
        summary = " | ".join(f"{s}:{n}" for s, n in counts.items())
        print(f"  Summary: {summary}")
    print(sep)


if __name__ == "__main__":
    project_id = None
    interval = 15
    show_board = True

    for arg in sys.argv[1:]:
        if arg == "--no-board":
            show_board = False
        elif arg.isdigit():
            interval = int(arg)
        else:
            project_id = arg

    print(f"Monitor started (interval={interval}s)")
    print("Stop: Ctrl+C  |  Hide board: --no-board")
    try:
        while True:
            show_status(project_id, show_board)
            time.sleep(interval)
    except KeyboardInterrupt:
        print("\nMonitor stopped.")
