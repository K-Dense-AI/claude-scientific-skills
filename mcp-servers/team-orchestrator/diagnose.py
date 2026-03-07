#!/usr/bin/env python3
"""
Team Orchestrator MCP 진단 스크립트
실행: python diagnose.py
"""
import sys
import os
import subprocess
import sqlite3
import json
from pathlib import Path

DIR = Path(__file__).parent
DB_PATH = DIR / "orchestrator.db"

SEP = "-" * 60

def section(title):
    print(f"\n{SEP}")
    print(f"  {title}")
    print(SEP)

def check_imports():
    section("1. Import 체크")
    modules = ["mcp.server.fastmcp", "state", "task_router", "team_manager", "shutdown_handler"]
    sys.path.insert(0, str(DIR))
    for mod in modules:
        try:
            __import__(mod)
            print(f"  [OK] {mod}")
        except Exception as e:
            print(f"  [FAIL] {mod}: {e}")

def check_db():
    section("2. DB 상태")
    if not DB_PATH.exists():
        print(f"  [FAIL] DB 없음: {DB_PATH}")
        return
    print(f"  [OK] DB: {DB_PATH}")
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    # 팀 목록
    cur.execute("SELECT team_id, name, status, task FROM teams ORDER BY created_at DESC LIMIT 20")
    rows = cur.fetchall()
    print(f"\n  팀 현황 ({len(rows)}개):")
    for team_id, name, status, task in rows:
        task_short = (task[:60] + "...") if task and len(task) > 60 else (task or "")
        status_icon = {"running": "[>]", "pending": "[ ]", "paused": "[~]", "done": "[v]", "failed": "[x]"}.get(status, "[?]")
        print(f"    {status_icon} {team_id} | {name} | {status}")
        print(f"       작업: {task_short}")
    conn.close()

def check_claude_exe():
    section("3. Claude.exe 프로세스")
    try:
        result = subprocess.run(
            ["tasklist", "/FI", "IMAGENAME eq claude.exe", "/FO", "CSV"],
            capture_output=True, text=True, timeout=5
        )
        lines = [l for l in result.stdout.strip().split("\n") if "claude.exe" in l.lower()]
        if lines:
            print(f"  [OK] claude.exe 실행 중: {len(lines)}개")
            for l in lines:
                print(f"    {l.strip()}")
        else:
            print("  [INFO] claude.exe 프로세스 없음 (팀 작업 중 아님)")
    except Exception as e:
        print(f"  [FAIL] tasklist 실행 오류: {e}")

def check_mcp_config():
    section("4. MCP 설정 (.claude.json)")
    config_path = Path.home() / ".claude.json"
    if not config_path.exists():
        print(f"  [FAIL] 파일 없음: {config_path}")
        return
    try:
        data = json.loads(config_path.read_text(encoding="utf-8"))
        servers = data.get("mcpServers", {})
        if "team-orchestrator" in servers:
            cfg = servers["team-orchestrator"]
            print(f"  [OK] team-orchestrator 등록됨")
            print(f"    command: {cfg.get('command')}")
            args = cfg.get("args", [])
            print(f"    args: {args}")
            # server.py 존재 여부
            if args:
                sp = Path(args[-1])
                exists = sp.exists()
                print(f"    server.py: {'[OK] 존재' if exists else '[FAIL] 없음'} -> {sp}")
        else:
            print("  [FAIL] team-orchestrator 설정 없음")
        print(f"\n  등록된 MCP 서버 목록:")
        for name in servers:
            print(f"    - {name}")
    except Exception as e:
        print(f"  [FAIL] 설정 파싱 오류: {e}")

def check_server_syntax():
    section("5. server.py 문법 체크")
    server_path = DIR / "server.py"
    try:
        result = subprocess.run(
            ["python", "-m", "py_compile", str(server_path)],
            capture_output=True, text=True, timeout=10
        )
        if result.returncode == 0:
            print(f"  [OK] 문법 정상")
        else:
            print(f"  [FAIL] 문법 오류:\n{result.stderr}")
    except Exception as e:
        print(f"  [FAIL] 실행 오류: {e}")

def check_pending_teams():
    section("6. pending 팀 start 가능 여부")
    if not DB_PATH.exists():
        print("  [SKIP] DB 없음")
        return
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    cur.execute("SELECT team_id, name, task FROM teams WHERE status='pending'")
    rows = cur.fetchall()
    conn.close()
    if not rows:
        print("  [INFO] pending 팀 없음")
    else:
        print(f"  [INFO] pending 팀 {len(rows)}개 (MCP 재연결 후 start_team으로 시작 가능):")
        for tid, name, task in rows:
            task_short = (task[:55] + "...") if task and len(task) > 55 else (task or "")
            print(f"    - {tid} | {name} | {task_short}")

def resume_pending_with_python():
    """MCP 없이 직접 팀 프로세스 실행 (옵션)"""
    section("7. MCP 없이 pending 팀 직접 실행 (선택)")
    print("  MCP 연결 없이 팀을 직접 시작하려면 아래 명령 실행:")
    if not DB_PATH.exists():
        print("  [SKIP] DB 없음")
        return
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    cur.execute("SELECT team_id, name, task, project_id FROM teams WHERE status='pending'")
    rows = cur.fetchall()
    conn.close()
    if not rows:
        print("  [INFO] pending 팀 없음")
        return
    for tid, name, task, proj_id in rows:
        print(f"\n  팀: {name} ({tid})")
        print(f"  python diagnose.py --start {tid}")

def main():
    print("\n" + "=" * 60)
    print("  Team Orchestrator MCP 진단")
    print("=" * 60)

    # --start 옵션: 직접 팀 실행
    if len(sys.argv) >= 3 and sys.argv[1] == "--start":
        team_id = sys.argv[2]
        sys.path.insert(0, str(DIR))
        import state as st
        import team_manager as tm
        st.init_db()
        team = st.get_team(team_id)
        if not team:
            print(f"팀 {team_id} 없음")
            return
        if team.status != "pending":
            # force to pending
            st.update_team_status(team_id, "pending")
            team = st.get_team(team_id)
        print(f"팀 '{team.name}' 시작 중...")
        tm.start_team(team.team_id, team.name, team.team_type, team.task, team.project_id)
        print("  프로세스 실행됨. 잠시 후 상태 확인하세요.")
        return

    # --start-all 옵션: pending 팀 전부 실행
    if len(sys.argv) >= 2 and sys.argv[1] == "--start-all":
        sys.path.insert(0, str(DIR))
        import state as st
        import team_manager as tm
        st.init_db()
        teams = [t for t in st.list_teams() if t.status == "pending"]
        if not teams:
            print("pending 팀 없음")
            return
        for team in teams:
            print(f"시작: {team.name} ({team.team_id})")
            tm.start_team(team.team_id, team.name, team.team_type, team.task, team.project_id)
        print(f"\n총 {len(teams)}개 팀 시작됨")
        return

    check_imports()
    check_db()
    check_claude_exe()
    check_mcp_config()
    check_server_syntax()
    check_pending_teams()
    resume_pending_with_python()

    print(f"\n{'=' * 60}")
    print("  진단 완료")
    print(f"{'=' * 60}\n")

if __name__ == "__main__":
    main()
