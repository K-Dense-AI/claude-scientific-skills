"""
team-orchestrator 직접 테스트
MCP 서버 없이 core 로직만 검증:
  1. DB 초기화
  2. 프로젝트 생성 + 팀 배정 (task_router)
  3. 팀 실행 (team_manager - claude -p subprocess)
  4. 상태 / 출력 조회
  5. CEO 피드백 전송
"""
import asyncio
import sys
import io
import time
from pathlib import Path

# cp949 콘솔 한국어 출력 문제 방지
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

sys.path.insert(0, str(Path(__file__).parent))

import state as st
import task_router
import team_manager

# ── 1. DB 초기화 ────────────────────────────────────────────────────────────
print("=== 1. DB 초기화 ===")
st.init_db()
print("  orchestrator.db OK")

# ── 2. 프로젝트 생성 + 팀 배정 ─────────────────────────────────────────────
PROJECT_DESC = (
    "biosteam-tagatose 레포(C:/Users/Jahyun/biosteam-tagatose)의 "
    "README.md가 있는지 확인하고, 없으면 프로젝트 개요를 한 문단으로 요약한 "
    "README.md를 작성하라."
)

print("\n=== 2. 프로젝트 분석 + 팀 배정 ===")
plan = task_router.analyze_and_split(PROJECT_DESC)
print(f"  요약: {plan['project_summary']}")
print(f"  팀 수: {len(plan['teams'])}")
for t in plan["teams"]:
    print(f"    [{t['type']}] {t['name']}: {t['task'][:60]}...")

# ── 3. 프로젝트 + 팀 등록 및 실행 ──────────────────────────────────────────
import uuid

project_id = str(uuid.uuid4())[:8]
st.create_project(project_id, PROJECT_DESC)
print(f"\n=== 3. 프로젝트 시작: {project_id} ===")

teams_created = []
for i, team_def in enumerate(plan["teams"]):
    team_id = f"{project_id}-{i+1}"
    name = team_def["name"]
    team_type = team_def["type"]
    task = team_def["task"]

    team_manager.start_team(team_id, name, team_type, task, project_id)
    teams_created.append(team_id)
    print(f"  [{team_type}] {name} → {team_id} 시작")

# ── 4. 실행 대기 + 상태 폴링 ────────────────────────────────────────────────
print("\n=== 4. 실행 중 (최대 600초 대기) ===")

loop = asyncio.get_event_loop()

async def wait_and_poll():
    for _ in range(120):  # 5초 × 120 = 600초
        await asyncio.sleep(5)
        statuses = team_manager.get_status_all(project_id)
        all_done = all(s["status"] in ("done", "failed") for s in statuses)
        for s in statuses:
            print(f"  [{s['team_id']}] {s['status']} | out={s['has_output']}", flush=True)
        if all_done:
            print("  → 모든 팀 완료!", flush=True)
            return
    print("  → 타임아웃")

loop.run_until_complete(wait_and_poll())

# ── 5. 출력 확인 ────────────────────────────────────────────────────────────
print("\n=== 5. 팀 산출물 ===")
for team_id in teams_created:
    output = st.get_team_output(team_id)
    print(f"\n[{team_id}] 출력 ({len(output)}자):")
    print(output[:500] + ("..." if len(output) > 500 else ""))

print("\n=== 테스트 완료 ===")
