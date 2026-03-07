#!/usr/bin/env python3
"""
Team Orchestrator Watchdog
- 5분마다 팀 상태 + claude.exe 수 체크
- running 팀이 20분 이상 claude.exe 없으면 STUCK 판정
- 로그: watchdog.log
- stuck 시 자동 재시작 시도
"""
import sys
import os
import sqlite3
import subprocess
import time
import json
from pathlib import Path
from datetime import datetime

DIR = Path(__file__).parent
DB_PATH = DIR / "orchestrator.db"
LOG_PATH = DIR / "watchdog.log"
CHECK_INTERVAL = 300        # 5분
STUCK_THRESHOLD = 1200      # 20분 (claude.exe 없는 상태 지속 시)

sys.path.insert(0, str(DIR))


def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with open(LOG_PATH, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def get_teams():
    if not DB_PATH.exists():
        return []
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    cur.execute("SELECT team_id, name, status, task FROM teams")
    rows = cur.fetchall()
    conn.close()
    return rows


def get_claude_count():
    try:
        result = subprocess.run(
            ["tasklist", "/FI", "IMAGENAME eq claude.exe", "/FO", "CSV"],
            capture_output=True, text=True, timeout=5
        )
        return sum(1 for l in result.stdout.split("\n") if "claude.exe" in l.lower())
    except Exception:
        return -1


def force_pending(team_id):
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    cur.execute("UPDATE teams SET status='pending' WHERE team_id=?", (team_id,))
    conn.commit()
    conn.close()


def restart_team(team_id, name, team_type, task, project_id):
    import state as st
    import team_manager as tm
    st.init_db()
    force_pending(team_id)
    team = st.get_team(team_id)
    if team:
        tm.start_team(team.team_id, team.name, team.team_type, team.task, team.project_id)
        log(f"  --> 재시작 완료: {name} ({team_id})")
    else:
        log(f"  --> 재시작 실패: 팀 없음 {team_id}")


def main():
    log("=" * 50)
    log("Watchdog 시작")
    log("=" * 50)

    # team_id -> 마지막으로 claude.exe 있었던 시각
    last_alive = {}

    while True:
        teams = get_teams()
        claude_count = get_claude_count()
        running = [(tid, name, status, task) for tid, name, status, task in teams if status == "running"]
        done = sum(1 for _, _, s, _ in teams if s == "done")
        failed = sum(1 for _, _, s, _ in teams if s == "failed")

        now = time.time()

        log(f"상태 | running:{len(running)} done:{done} failed:{failed} | claude.exe:{claude_count}개")

        for tid, name, status, task in running:
            task_short = (task[:40] + "...") if task and len(task) > 40 else (task or "")
            log(f"  [>] {name} ({tid}): {task_short}")

        if claude_count > 0:
            # claude.exe 살아있으면 alive 시각 갱신
            for tid, name, status, task in running:
                last_alive[tid] = now
        else:
            # claude.exe 없음 → stuck 여부 체크
            if running:
                log(f"  ! claude.exe 없음, running 팀 {len(running)}개 대기 중")
                for tid, name, status, task in running:
                    alive_at = last_alive.get(tid, now)
                    elapsed = now - alive_at
                    log(f"    {name}: {int(elapsed)}초 비활성")

                    if elapsed >= STUCK_THRESHOLD:
                        log(f"  [STUCK] {name} ({tid}) - {int(elapsed)}초 이상 비활성 -> 재시작 시도")
                        # DB에서 task, project_id 가져오기
                        conn = sqlite3.connect(DB_PATH)
                        cur = conn.cursor()
                        cur.execute("SELECT team_type, project_id FROM teams WHERE team_id=?", (tid,))
                        row = cur.fetchone()
                        conn.close()
                        if row:
                            team_type, project_id = row
                            try:
                                restart_team(tid, name, team_type, task, project_id)
                            except Exception as e:
                                log(f"  [ERROR] 재시작 오류: {e}")
            else:
                log("  모든 팀 완료 또는 없음. Watchdog 종료.")
                break

        time.sleep(CHECK_INTERVAL)


if __name__ == "__main__":
    main()
