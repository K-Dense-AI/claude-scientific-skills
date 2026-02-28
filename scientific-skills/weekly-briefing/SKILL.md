# 주간 브리핑 스킬

매주 초(월요일) 실행하여 이번 주 신경 써야 할 사항들을 한눈에 정리한다.

## 실행 순서

memory 파일 `~/.claude/projects/C--Users-Jahyun-Desktop/memory/weekly-briefing.md`를 먼저 읽어 소스 URL 및 경로를 확인한다.

아래 7개 소스를 **병렬로** 최대한 동시에 조회하고, 결과를 통합해 한 번에 브리핑을 출력한다.

## 모델 분배 전략 (토큰 효율화)

| 소스 | 모델 | 방식 | 이유 |
|------|------|------|------|
| [1] Asana 태스크 | **Haiku** | Task subagent | MCP 호출 + 결과 반환만 |
| [2] Google Calendar | **Haiku** | Task subagent | MCP 호출 + 결과 반환만 |
| [3] FBE 주문목록 | **Sonnet** | Task subagent | Chrome JS + 데이터 파싱 |
| [4] FBE 시약목록 | **Sonnet** | Task subagent | Chrome JS + 데이터 파싱 |
| [5] OneDrive 호서대 | **Haiku** | Task subagent | `ls -lt` 실행 + 결과 반환 |
| [6] OneDrive 고려대 | **Haiku** | Task subagent | `ls -lt` 실행 + 결과 반환 |
| [7] GitHub | **Haiku** | Task subagent | 웹페이지 읽기만 |
| **최종 브리핑 작성** | **Opus** (메인) | 직접 | 종합 분석 + 정리 |

- Haiku: 단순 데이터 수집 (MCP 호출, bash 명령, 웹 읽기)
- Sonnet: Chrome 브라우저 조작 + 데이터 파싱이 필요한 작업
- Opus: 최종 종합 분석 및 브리핑 출력 (메인 에이전트)

---

## 소스별 조회 방법

### [1] Asana (Asana MCP)
- `asana_get_tasks`로 project=`1213473127047670` (팀원관리) 조회 (search_tasks는 Premium 전용이라 사용 불가)
- opt_fields로 `name,assignee,due_on,completed,created_at` 포함
- 마감 임박(3일 이내), 마감 초과(overdue), 이번 주 신규 태스크로 분류

### [2] Google Calendar (Calendar MCP)
- `gcal_list_events`로 이번 주 월~일 일정 조회 (calendarId="primary")
- 공유 캘린더가 있다면 함께 조회
- 실험 일정, 미팅, 발표, 마감일 표시

### [3] FBE 주문목록 (Chrome JS)
- `https://docs.google.com/spreadsheets/d/1MD9pugKxYEjk8NX0Dr-YSHAigSFfSmOuQbr02qNQeGc/edit` 로 이동
- 현재 연도 탭(`시약 및 소모품_2026` 등) 접근
- JavaScript로 시트 데이터 읽기:
  ```javascript
  // Google Sheets API v4 사용 (로그인 상태에서 가능)
  // 또는 화면에서 직접 읽기
  ```
- F열(주문요청일) 최근 2주 / K열 주문 상태 미처리 건 확인

### [4] FBE 시약목록 (Chrome JS)
- `https://docs.google.com/spreadsheets/d/1_Vn0HKqmCpqU6YAH_Odzcv74H7Y_t4sG6CPVS7wDwho/edit` 로 이동
- 현재 연도(2026) 검수 컬럼이 비어있는 시약 확인
- 탭별(`209 시약`, `서관 226` 등) 스캔

### [5] OneDrive - 호서대학교 (로컬 파일)
- **주의: Git Bash의 `find -newermt`는 Windows에서 작동하지 않음**
- PowerShell 스크립트 사용: `powershell -File "C:/Users/Jahyun/Desktop/scan_onedrive.ps1"`
- 또는 연구 관련 폴더만 `ls -lt`로 개별 조회:
  ```bash
  ls -lt "/c/Users/Jahyun/OneDrive - 호서대학교/바탕 화면 [Labtop]/" | head -15
  ls -lt "/c/Users/Jahyun/OneDrive - 호서대학교/바탕 화면 [Labtop]/Team meeting/Journal Club/" | head -10
  ls -lt "/c/Users/Jahyun/OneDrive - 호서대학교/바탕 화면 [Labtop]/Team meeting/research presentation/" | head -10
  ls -lt "/c/Users/Jahyun/OneDrive - 호서대학교/바탕 화면 [Labtop]/resutls/" | head -10
  ```
- 연구 관련 폴더: `바탕 화면 [Labtop]/` (실험정리, graph, 논문), `FBE 공유폴더/`

### [6] OneDrive - 고려대학교 (로컬 파일)
- 연구 관련 폴더만 `ls -lt`로 개별 조회:
  ```bash
  ls -lt "/c/Users/Jahyun/OneDrive - 고려대학교/저장소/" | head -15
  ls -lt "/c/Users/Jahyun/OneDrive - 고려대학교/저장소/8. D-Gal to MA and D-tagatose/" | head -10
  ls -lt "/c/Users/Jahyun/OneDrive - 고려대학교/저장소/9. 한미해조류_해수부/" | head -10
  ```
- 연구 관련 폴더: `저장소/` (프로젝트별 method 파일, 결과 데이터)

### [7] GitHub (Chrome → GitHub 웹)
- `https://github.com/jahyunlee00299?tab=repositories` 접속
- 주요 활성 레포 커밋 날짜 확인:
  - UDH_Clustering, Kinetic-modeling-and-optimization, biosteam-tagatose
- 지난 주 이후 업데이트된 레포 표시

---

## 브리핑 출력 형식

```markdown
# 주간 브리핑 — [월요일 날짜]

## Asana
- **마감 초과 (즉시 처리 필요)**: N건
  - [태스크명] (마감: X일 초과)
- **이번 주 마감**: N건
  - [태스크명] (D-N)
- **이번 주 신규**: N건

## 이번 주 일정
- [요일] [시간] [일정명]
- ...

## 주문 현황
- 처리 대기 중: N건
  - [품명] (요청일: YYMMDD, 요청자: ...)
- 이번 주 주문 예정: N건

## 시약 재고 확인 필요
- 2026 검수 미완료: N건 (우선순위 높은 것 위주)
  - [시약명] (위치: ...)

## OneDrive 최근 변경 (지난 7일)
- 호서대: [파일명] (수정일)
- 고려대: [파일명] (수정일)

## GitHub 활동
- 커밋한 레포: [레포명] (X일 전)
- 오래된 레포(주의): [레포명] (마지막 커밋: N주 전)
```

---

## 주의사항
- Google Sheets 접근 시 Chrome 브라우저가 Google 계정에 로그인되어 있어야 함
- **OneDrive 파일 스캔**: Git Bash `find -newermt`는 Windows에서 작동하지 않음. `ls -lt`로 연구 폴더별 개별 조회하거나 PowerShell 스크립트 사용
- `asana_search_tasks`는 Premium 전용 — `asana_get_tasks` (project 기반) 사용
- Asana multi-homing 작업 시 항상 `opt_silent=true` 파라미터 사용
