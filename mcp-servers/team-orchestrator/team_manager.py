"""팀 생명주기 관리 - Agent SDK + Detached subprocess 방식"""
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Optional

import state as st

PYTHON_BIN = sys.executable
RUN_AGENT_SCRIPT = str(Path(__file__).parent / "run_agent.py")
DEFAULT_WORKDIR = str(Path.home())

DETACHED_PROCESS = 0x00000008
CREATE_NEW_PROCESS_GROUP = 0x00000200
CREATE_NO_WINDOW = 0x08000000

_pids: dict[str, int] = {}
_task_files: dict[str, str] = {}
_project_feedback_timeouts: dict[str, int] = {}


def _resolve_dependents(completed_team_id: str, feedback_timeout: int = 10):
    """완료된 팀에 의존하던 pending 팀을 자동 시작"""
    if feedback_timeout > 0:
        time.sleep(feedback_timeout)

    completed = st.get_team(completed_team_id)
    if not completed:
        return

    all_teams = st.list_teams(completed.project_id)
    name_to_team = {t.name: t for t in all_teams}

    for team in all_teams:
        if team.status != "pending" or not team.depends_on:
            continue
        try:
            deps = json.loads(team.depends_on)
        except Exception:
            deps = []
        if not deps:
            continue

        all_deps_done = all(
            name_to_team.get(dep) and name_to_team[dep].status == "done"
            for dep in deps
        )
        if all_deps_done:
            st.append_team_output(
                team.team_id,
                f"\n[의존성 해결 - 자동 시작: {', '.join(deps)} 완료]\n"
            )
            ft = _project_feedback_timeouts.get(team.project_id, 10)
            start_team(team.team_id, team.name, team.team_type, team.task,
                       team.project_id, feedback_timeout=ft)


def _launch_detached_agent(team_id: str, team_type: str, task: str,
                           project_id: str, workdir: str) -> int:
    """task를 임시 파일에 저장 후 run_agent.py를 detached process로 실행. PID 반환."""
    # task를 파일로 저장 (긴 내용도 안전하게 전달)
    tmp = tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", suffix=".txt",
        prefix=f"task_{team_id}_", delete=False
    )
    tmp.write(task)
    tmp.close()
    _task_files[team_id] = tmp.name

    clean_env = {k: v for k, v in os.environ.items() if k != "CLAUDECODE"}

    proc = subprocess.Popen(
        [PYTHON_BIN, RUN_AGENT_SCRIPT,
         team_id, team_type, project_id, workdir, tmp.name],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        env=clean_env,
        creationflags=DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP | CREATE_NO_WINDOW,
        close_fds=True,
    )
    return proc.pid


def start_team(team_id: str, name: str, team_type: str, task: str,
               project_id: str, feedback_timeout: int = 10,
               depends_on: list = None,
               workdir: str = DEFAULT_WORKDIR) -> str:
    _project_feedback_timeouts[project_id] = feedback_timeout

    team = st.Team(
        team_id=team_id,
        name=name,
        team_type=team_type,
        status="pending",
        task=task,
        project_id=project_id,
        created_at=time.time(),
        depends_on=json.dumps(depends_on or []),
    )
    st.upsert_team(team)

    pid = _launch_detached_agent(team_id, team_type, task, project_id, workdir)
    _pids[team_id] = pid

    st.append_team_output(team_id, f"[{team_type.upper()} 에이전트 시작 PID={pid}]\n")
    return team_id


def send_feedback(team_id: str, message: str):
    team = st.get_team(team_id)
    if not team:
        raise ValueError(f"팀 {team_id}를 찾을 수 없습니다.")

    # 메시지 큐에 CEO 피드백 저장 (에이전트가 read_inbox로 읽을 수 있도록)
    st.push_message(team_id, "ceo", message)

    # 완료된 팀이면 새 프로세스로 재실행
    if team.status == "done":
        ft = _project_feedback_timeouts.get(team.project_id, 10)
        prev = st.get_team_output(team_id) or ""
        context = prev[-1000:] if len(prev) > 1000 else prev
        continuation_task = (
            f"[이전 작업 결과 요약]\n{context}\n\n"
            f"[CEO 추가 지시]\n{message}"
        )
        st.update_team_status(team_id, "pending")
        start_team(team_id, team.name, team.team_type, continuation_task,
                   team.project_id, feedback_timeout=ft)


def set_project_feedback_timeout(project_id: str, timeout: int):
    _project_feedback_timeouts[project_id] = timeout


def pause_team(team_id: str):
    """팀 상태를 paused로 표시 (detached 프로세스는 다음 체크포인트에서 감지)"""
    st.update_team_status(team_id, "paused")


def get_status_all(project_id: Optional[str] = None) -> list[dict]:
    teams = st.list_teams(project_id)
    return [
        {
            "team_id": t.team_id,
            "name": t.name,
            "type": t.team_type,
            "status": t.status,
            "task": t.task[:100] + "..." if len(t.task) > 100 else t.task,
            "has_output": bool(t.output),
            "pid": _pids.get(t.team_id),
        }
        for t in teams
    ]
