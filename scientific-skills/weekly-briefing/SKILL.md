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
- 마감 초과 (즉시): N건
- 이번 주 마감: N건
- 신규: N건

## 이번 주 일정
- [요일] [시간] [일정명]

## OneDrive 최근 변경 (7일)
- 호서대: [파일명] (날짜)
- 고려대: [파일명] (날짜)

## GitHub 활동
- [레포명]: 커밋 N건 (마지막: YYYY-MM-DD)
```

## 주의사항

- Calendar: `gcal_list_events` calendarId는 `primary` (jahyunlee082@gmail.com)
- Asana `asana_search_tasks`는 Premium 전용 — `get_tasks` 사용
- Notion 하위 페이지 생성 시 parent page_id: `315f91aca96f81369f05cf398f64b801`
