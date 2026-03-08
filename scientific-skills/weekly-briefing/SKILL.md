# 주간 브리핑 스킬

매주 초(월요일) 실행하여 이번 주 신경 써야 할 사항들을 한눈에 정리한다.

## 실행 순서

```
[0] 공휴일 확인 (CEO 직접)
[1] Notion 이력 참조 (CEO MCP — notion-fetch)
[2] 병렬 실행:
    ├── bash: python collect_all.py  → Asana + GitHub + OneDrive JSON 저장
    └── CEO MCP: gcal_list_events    → Calendar 조회
[3] 결과 파일 읽기 + Calendar 결과 종합
[4] 브리핑 작성 + 출력
[5] Notion 업데이트 (CEO MCP — notion-update-page)
```

## [2] 스크립트 실행

```bash
C:/Users/Jahyun/anaconda3/python.exe \
  "C:/Users/Jahyun/claude-scientific-skills/scientific-skills/weekly-briefing/scripts/collect_all.py"
```

결과 파일 (C:\Users\Jahyun\lab-analyses\temp\):
- `briefing_asana.json`   — overdue / due_this_week / new_this_week
- `briefing_github.json`  — 레포별 최근 커밋
- `briefing_onedrive.json` — 최근 7일 변경 파일

### scripts/ 디렉토리 구성

```
weekly-briefing/scripts/
├── collect_all.py        — 통합 실행 진입점 (위 3개 + 아래 2개 호출)
├── collect_asana.py      — Asana 태스크 수집 (asana-extended-api 연동)
├── collect_github.py     — GitHub PR/커밋/이슈 수집
├── collect_onedrive.py   — OneDrive 최근 변경 파일 수집
├── collect_orders.py     — 시약 주문 현황 수집
├── collect_reagents.py   — 시약 재고 현황 수집
└── google_auth.py        — Google API 인증 헬퍼
```

### Asana 연동 (asana-extended-api 스킬)

`collect_asana.py`는 `asana-extended-api/scripts/asana_api.py`의 `AsanaAPI` 클래스를 직접 import하여 사용한다:

```python
import sys
sys.path.insert(0, r"C:\Users\Jahyun\claude-scientific-skills\scientific-skills\asana-extended-api\scripts")
from asana_api import AsanaAPI
from pathlib import Path
import json

secrets = json.loads(Path("~/.claude/secrets.json").expanduser().read_text())
api = AsanaAPI(token=secrets["ASANA_PAT"])
```

수집 항목: `overdue`, `due_this_week`, `new_this_week` (created 7일 이내)

### GitHub 연동

`collect_github.py`는 `GITHUB_PAT`으로 REST API를 호출하여 PR/이슈/커밋을 수집한다:

```python
secrets = json.loads(Path("~/.claude/secrets.json").expanduser().read_text())
headers = {"Authorization": f"Bearer {secrets['GITHUB_PAT']}"}
# GET https://api.github.com/repos/{owner}/{repo}/commits?since={7일전}&per_page=10
# GET https://api.github.com/repos/{owner}/{repo}/pulls?state=open
# GET https://api.github.com/repos/{owner}/{repo}/issues?state=open&since={7일전}
```

### Google Drive 연동 (gdrive MCP)

브리핑 실행 시 CEO가 직접 MCP로 조회 (스크립트 미포함):

```
mcp__gdrive__search(query="modifiedTime > '{7일전}' order by modifiedTime desc", pageSize=10)
```

결과에서 드라이브 공유 파일 중 최근 변경 항목을 브리핑에 포함한다.

## [2] Calendar MCP (CEO 직접)

```
gcal_list_events(
  calendarId="primary",
  timeMin=<이번주 월요일 00:00 KST>,
  timeMax=<이번주 일요일 23:59 KST>
)
```

## [5] Notion 업데이트

- 페이지: Ribose Team 연구 Overview (315f91ac-a96f-8136-9f05-cf398f64b801)
- 형식: 주간 브리핑 날짜 제목으로 하위 페이지 생성

---

## 브리핑 출력 형식

```markdown
# 주간 브리핑 — YYYY-MM-DD

## Asana
- 마감 초과 (즉시): N건 — [태스크명](permalink_url) 형식으로 링크 포함
- 이번 주 마감: N건
- 신규: N건

## 이번 주 일정
- [요일] [시간] [일정명]

## OneDrive 최근 변경 (7일)
- 호서대: [파일명] (날짜)
- 고려대: [파일명] (날짜)

## Google Drive 최근 변경 (7일)
- [파일명] (날짜, 공유자)

## GitHub 활동
- [레포명]: 커밋 N건 (마지막: YYYY-MM-DD), 오픈 PR N건, 오픈 이슈 N건
```

## 주의사항

- Calendar: `gcal_list_events` calendarId는 `primary` (jahyunlee082@gmail.com)
- Asana `asana_search_tasks`는 Premium 전용 — `get_tasks` 사용
- Notion 하위 페이지 생성 시 parent page_id: `315f91aca96f81369f05cf398f64b801`
