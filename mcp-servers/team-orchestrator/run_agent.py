#!/usr/bin/env python3
"""
Agent SDK 실행 진입점 — detached subprocess로 실행됨
사용법: python run_agent.py <team_id> <team_type> <project_id> <workdir> <task_file>
"""
import sys
import asyncio
from pathlib import Path

DIR = Path(__file__).parent
sys.path.insert(0, str(DIR))

import state as st


def main():
    if len(sys.argv) < 6:
        print("Usage: run_agent.py <team_id> <team_type> <project_id> <workdir> <task_file> [mode] [lead_id]")
        sys.exit(1)

    team_id = sys.argv[1]
    team_type = sys.argv[2]
    project_id = sys.argv[3]
    workdir = sys.argv[4]
    task_file = sys.argv[5]
    mode = sys.argv[6] if len(sys.argv) > 6 else "lead"  # lead | worker | solo
    lead_id = sys.argv[7] if len(sys.argv) > 7 else None

    st.init_db()

    try:
        task = Path(task_file).read_text(encoding="utf-8")
    except Exception as e:
        st.update_team_status(team_id, "failed", error=f"task 파일 읽기 실패: {e}")
        sys.exit(1)

    import agent_runner
    if mode == "solo":
        asyncio.run(agent_runner.run_solo_agent(team_id, team_type, task, project_id, workdir))
    else:
        is_lead = mode != "worker"
        asyncio.run(agent_runner.run_agent(team_id, team_type, task, project_id, workdir,
                                           is_lead=is_lead, lead_id=lead_id))


if __name__ == "__main__":
    main()
