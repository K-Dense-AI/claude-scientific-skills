#!/usr/bin/env python3
"""
에이전트 완료 폴링 스크립트.
CEO (Claude Code VSCode 세션)가 run_in_background로 실행 -> 완료 시 task-notification 발생.

사용법:
    python poll_until_done.py <agent_id> [interval=20]
"""
import sqlite3
import sys
import time
from pathlib import Path

DB = Path(__file__).parent / "orchestrator.db"
INTERVAL = int(sys.argv[2]) if len(sys.argv) > 2 else 20
agent_id = sys.argv[1]

print(f"[poll] {agent_id} 모니터링 시작 (간격 {INTERVAL}초)")

while True:
    try:
        if DB.exists() and DB.stat().st_size > 0:
            conn = sqlite3.connect(str(DB))
            row = conn.execute(
                "SELECT status, output FROM teams WHERE team_id=?", (agent_id,)
            ).fetchone()
            conn.close()
            if row and row[0] in ("done", "failed"):
                status = row[0]
                output = (row[1] or "").strip()[-300:]
                print(f"\nAgent done: {agent_id} [{status}]")
                print(f"--- output ---\n{output}")
                sys.exit(0)
    except sqlite3.OperationalError:
        pass  # DB 초기화 중 — 다음 폴링에서 재시도
    time.sleep(INTERVAL)
