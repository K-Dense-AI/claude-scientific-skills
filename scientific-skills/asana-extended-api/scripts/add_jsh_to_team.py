"""
장서현 미완료 태스크 → 팀원관리 프로젝트 multi-homing
실행: python add_jsh_to_team.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from asana_api import AsanaAPI

PROJECT_TEAM = "1213473127047670"  # 팀원관리

# 장서현 미완료 태스크 15건
TASKS = [
    ("1209564722654322", "Xylulose kinase (XylB and XKS1)"),
    ("1209564722654320", "Xylose isomerase (XylA)"),
    ("1211779567725587", "ManUA Substrate condition"),
    ("1205475236174590", "GalUA Substrate condition"),
    ("1212538124551483", "비효소적 배경 반응 검증"),
    ("1213009032968949", "UDH activity on L-GulUA"),
    ("1213260420488362", "Mannaric acid standard"),
    ("1212636872667469", "NADH:product 불일치 확인"),
    ("1213236563380174", "Inhibition of HMP/Furfural"),
    ("1210268641271868", "Structure differences between uronates"),
    ("1210268641271866", "Substrate specificity of UDHs"),
    ("1211960508955745", "실험조건 설정"),
    ("1211145723763555", "7. Cell stock 제작"),
    ("1209564722654324", "Ribulose-phosphate 3-epimerase (Rpe)"),
    ("1210270078280992", "Uronate dehydrogenase (UDH)"),
]

if __name__ == "__main__":
    api = AsanaAPI()
    task_gids = [gid for gid, _ in TASKS]

    print(f"장서현 태스크 {len(task_gids)}건 → 팀원관리 프로젝트 추가 시작\n")
    result = api.batch_add_to_project(task_gids, PROJECT_TEAM)

    print(f"\n완료: 성공 {result['success']}건, 실패 {result['fail']}건")
    if result['fail'] > 0:
        for d in result['details']:
            if d.get('error'):
                print(f"  실패: {d}")
