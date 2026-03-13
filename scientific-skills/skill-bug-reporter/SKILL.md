---
name: skill-bug-reporter
description: >
  스킬 실행 중 발생하는 오류(cp949 인코딩, 경로 오류, 형식 오류 등)를
  Haiku 에이전트가 자동으로 감지하고, 수정 브랜치+PR 자동 생성 → 깃 커밋을 통해 추적한다.
  같은 (스킬, 오류유형) 조합이 2회 이상이면 자동으로 수정 에이전트를 호출한다.

  다음 상황에서 반드시 이 스킬을 사용할 것:
  - 스킬 스크립트 실행 시 UnicodeDecodeError / UnicodeEncodeError (cp949, utf-8 등)
  - 스킬 스크립트에 잘못된 경로(하드코딩, 상대경로 오류 등)
  - SKILL.md 코드 블록 형식 오류
  - print() 구문의 이모지로 인한 Windows cp49 출력 깨짐
  - 스킬에서 생성한 스크립트가 FileNotFoundError / ModuleNotFoundError 발생
  - Claude가 코드 분석 중 버그/문제를 언어적으로 식별한 경우:
    * "버그가 있습니다", "두 가지 버그", "N가지 버그" 등 버그 개수 명시
    * "문제 파악했습니다", "문제가 있습니다", "문제를 발견" 등 문제 발견 선언
    * "코드를 확인합니다", "코드를 살펴보겠습니다" 등 의심 기반 코드 점검 선언
    * "오류 발견", "잘못된 로직", "그룹키 버그", "표시 오류" 등 구체적 버그 지칭
license: MIT license
metadata:
    skill-author: lab-internal
---

# Skill Bug Reporter (스킬 버그 자동 기록 + PR/브랜치 자동 생성)

스킬 관련 오류를 Haiku 에이전트가 MD 파일에 기록하고,
**같은 오류가 2회 이상 반복되면 수정 에이전트를 자동 호출**한다.

## 기록 파일

```
C:\Users\Jahyun\.claude\logs\skill_bugs.md
```

없으면 자동 생성. 있으면 맨 아래에 append.

---

## 오류 분류표

| 오류 유형 | 예시 증상 |
|----------|----------|
| **encoding** | UnicodeEncodeError cp949, 이모지 출력 깨짐 |
| **path** | FileNotFoundError, 하드코딩 경로 오류 |
| **format** | SKILL.md 형식 오류, 코드블록 누락 |
| **import** | ModuleNotFoundError, ImportError |
| **runtime** | 스킬 스크립트 실행 중 예외 |
| **logic** | Claude가 코드 리뷰/분석 중 언어적으로 버그·문제를 식별한 경우. 실행 없이도 기록. 예: "두 가지 버그가 있습니다", "그룹키 버그", "문제 파악했습니다", "코드를 확인합니다" |

---

## 실행 흐름

```
오류 감지
    |
    v
[Phase 0] Branch/PR 자동 생성
    |
    |-- feature 브랜치 생성: fix/{skill_name}-{error_type}-{timestamp}
    |-- 브랜치 체크아웃
    |
    v
[Phase 1] Haiku 에이전트 — 기록 + 중복 카운트
    |
    |-- 해당 (스킬, 유형) 첫 번째 오류?
    |   --> skill_bugs.md에 행 append + PR 링크 포함, 종료
    |
    \-- 2회 이상 반복?
        --> [Phase 2] Sonnet 수정 에이전트 자동 호출
            |
            v
        실제 파일 수정 후
        PR draft -> ready로 전환
        skill_bugs.md에 [FIXED] 표시 + PR 링크 업데이트
```

---

## Phase 0: Branch/PR 자동 생성

```bash
# 현재 브랜치: feat/team-orchestrator-v2 또는 main
git checkout -b fix/{skill_name}-{error_type}-{timestamp}
  # 예: fix/weekly-briefing-encoding-20260313-143052

# PR 자동 생성 (draft 상태, 구체적 설명 포함)
gh pr create \
  --title "[skill-bug-reporter] {skill_name}: {error_type} 자동 수정" \
  --body "## 오류 정보

- **스킬:** {skill_name}
- **오류 유형:** {error_type}
- **오류 메시지:** {error_message}
- **발생 파일:** {file_path}

## 예상 수정사항
[Phase 2의 수정 에이전트가 구체적 수정 내용 추가]

자동 생성된 PR입니다. 준비 완료 후 main으로 merge합니다." \
  --draft

# 출력에서 PR URL 추출 (예: https://github.com/jahyunlee00299/claude-scientific-skills/pull/42)
```

**중요사항:**
- PR 생성 후 즉시 skill_bugs.md 기록에 PR URL 포함
- Phase 2에서 수정 후 `gh pr ready {PR_NUMBER}` 로 draft → ready 전환

---

## Phase 1: 기록 에이전트 프롬프트

```
Agent(subagent_type="general-purpose", model="haiku", prompt="""
## 작업
C:\\Users\\Jahyun\\.claude\\logs\\skill_bugs.md 를 읽고 버그를 기록하라.
파일 없으면 헤더부터 생성.

## 입력
- 오류_유형: {error_type}  # encoding/path/format/import/runtime/logic 중 하나
- 스킬_이름: {skill_name}
- 오류_메시지: {error_message}  # logic 유형이면 Claude가 발한 진단 문장 그대로
- 발생_파일: {file_path}       # logic 유형이면 분석 중이던 파일 경로 (없으면 "unknown")
- PR_URL: {pr_url}              # Phase 0에서 생성된 PR URL
- 날짜: 오늘 날짜 (YYYY-MM-DD)

## logic 유형 특이사항
- 실제 Python 예외가 없어도 기록 가능
- error_message에는 Claude가 발화한 진단 문장 전체를 그대로 기입
  예) "두 가지 버그가 있습니다: 그룹키 버그(t0 보정값 사용), 제목 표시 오류(보정값 표시)"
- file_path는 분석 대상 파일; 알 수 없으면 "unknown" 입력
- 나머지 절차(중복 카운트, append)는 동일

## 절차
1. skill_bugs.md 읽기
2. 같은 (스킬_이름, 오류_유형) 조합이 기존에 몇 번 기록됐는지 카운트
3. 파일에 새 행 append:
   | {날짜} | {skill_name} | {error_type} | {한줄_증상} | `{file_path}` | {카운트}회째 | [PR]({pr_url}) |

## 파일 없을 때 헤더
# Skill Bug Log

| 날짜 | 스킬 | 유형 | 증상 요약 | 발생 파일 | 횟수 | PR 링크 | 상태 |
|------|------|------|----------|----------|------|--------|------|

## 반환값
- count: 이번 포함 누적 횟수 (int)
- pr_url: PR URL
- file_path: skill_bugs.md 경로
""")
```

**카운트 >= 2이면 즉시 Phase 2 호출.**

---

## Phase 2: 수정 에이전트 프롬프트

```
Agent(subagent_type="general-purpose", model="sonnet", prompt="""
## 역할
스킬 파일의 반복 오류를 직접 수정하고 PR을 업데이트하는 에이전트.
아래 스킬들의 규칙을 적용하여 수정한다:
- code-excellence: 수정 후 코드 품질·중복·간결성 재검토
- code-validator: 수정된 파일 실행 검증 (가능한 경우 pytest 또는 python -c import)
- conda-env-manager: import 오류 시 환경 진단 및 패키지 설치 확인
- solid-principles: 구조적 문제(경로 하드코딩 등)는 SOLID 원칙 기준으로 리팩토링
- git-workflow-manager: 수정 완료 후 변경 파일 git add + commit (메시지: "fix({skill_name}): {error_type} 반복 오류 자동 수정")

## 대상
- 스킬: {skill_name}
- 오류 유형: {error_type}
- 발생 파일: {file_path}
- 오류 메시지: {error_message}
- PR URL: {pr_url}
- skill_bugs.md 누적 횟수: {count}회
- fix 브랜치: fix/{skill_name}-{error_type}-{timestamp}

## 오류 유형별 수정 방법

### encoding (cp949 / 이모지)
1. {file_path} 열기
2. `import sys` 줄 다음에 삽입:
   import io
   sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
3. print() 내 이모지 전수 교체:
   ✅->[OK]  ❌->[ERROR]  ⚠️->[WARN]  🔄->[...]  📊->[CHART]  🔬->[SCI]
4. 같은 스킬의 다른 .py 파일도 동일 패턴 검색 후 일괄 수정 (Grep 활용)
5. code-validator로 수정 파일 import 검증

### path (경로 오류)
1. {file_path} 에서 하드코딩 경로 전수 탐색 (Grep: r'[A-Z]:\\\\' 또는 '../')
2. 교체:
   SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
   output_path = os.path.join(SCRIPT_DIR, ...)
3. solid-principles의 단일책임 원칙: 경로 상수는 파일 상단에 모아서 정의
4. code-validator로 경로 존재 여부 검증

### import (ModuleNotFoundError)
1. conda-env-manager로 현재 환경에 패키지 설치 여부 진단
2. 미설치 시 설치 명령 실행 또는 {skill_name}/SKILL.md 설치 단계에 추가
3. 조건부 import 패턴 적용 (ImportError try/except + 안내 메시지)

### format (SKILL.md 형식 오류)
1. {skill_name}/SKILL.md 열어 YAML frontmatter `---` 닫힘 확인
2. 코드블록 ``` 짝 불일치 탐색 후 수정
3. code-excellence로 전체 문서 구조 재검토

### runtime (일반 실행 오류)
1. 오류 traceback 분석 → 근본 원인 파악
2. code-excellence 규칙으로 방어 코드 최소한으로 추가
3. code-validator로 수정 후 재실행 검증

### logic (Claude 진단 버그)
1. error_message(진단 문장)에서 버그 항목 목록 추출
2. {file_path} 열어 해당 로직 확인
3. 버그별 수정:
   - 그룹키 오류: 보정값 대신 공칭값(nominal) 기준으로 그룹핑하도록 수정
   - 표시값 오류: label/title에 원본 파라미터(공칭값) 변수로 교체
   - 기타: 진단 문장 기반으로 직접 분석 후 최소 변경으로 수정
4. code-validator로 수정 후 검증 (가능한 경우 python -c import 또는 pytest)
5. 수정 불가능한 경우(file_path="unknown") → skill_bugs.md에 [PENDING] 표시 후 종료

## 완료 후
1. git add + commit (메시지: "fix({skill_name}): {error_type} 반복 오류 자동 수정")
2. git push origin fix/{skill_name}-{error_type}-{timestamp}
3. PR body 업데이트 + gh pr ready {PR_NUMBER} (draft → ready 전환)
4. skill_bugs.md 해당 행 상태 [FIXED] 로 업데이트
5. CEO에게: [PR]({pr_url}) + 1줄 요약 + commit hash 보고
""")
```

---

## 오류 유형별 빠른 수정 참조

### encoding
```python
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
```
이모지 대체표: `✅→[OK]` `❌→[ERROR]` `⚠️→[WARN]` `🔄→[...]` `📊→[CHART]` `🔬→[SCI]`

### path
```python
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
output_path = os.path.join(SCRIPT_DIR, "output", "result.csv")
```

### import
SKILL.md 설치 단계에 `pip install {package}` 추가.

---

## 사용 예

```
# encoding 유형
1회: lab-viz 이모지
  -> Phase 0: fix/lab-viz-encoding-20260313-143052 브랜치 + PR #42 자동 생성
  -> Phase 1: skill_bugs.md 기록 (PR 링크 포함), 종료

2회: lab-viz 이모지 재발
  -> Phase 0: 기존 fix 브랜치에 새 커밋 추가
  -> Phase 1: skill_bugs.md 업데이트 (카운트 2회째)
  -> Phase 2: Sonnet 수정 에이전트 자동 호출
             -> lab-viz/scripts/plot.py 실제 수정
             -> PR #42 body 업데이트 + draft → ready 전환
             -> git push + commit
             -> skill_bugs.md [FIXED] 표시

# logic 유형 (신규)
Claude: "두 가지 버그가 있습니다: 그룹키 버그(t0 보정값 사용), 제목 표시 오류(보정값 표시)"
1회:
  -> Phase 0: fix/kinetic-bo-pipeline-logic-20260313-150000 브랜치 + PR #43 생성
  -> Phase 1: skill_bugs.md 기록, 종료

2회: 같은 버그 재발견
  -> Phase 2: 수정 에이전트 자동 호출 + PR ready + [FIXED]
```

**CEO 보고:**
```
[1회] [PR #42](https://github.com/jahyunlee00299/claude-scientific-skills/pull/42) — lab-viz encoding 버그 기록

[2회] [PR #42](https://github.com/jahyunlee00299/claude-scientific-skills/pull/42) 자동 수정 완료
      - lab-viz/scripts/plot.py (UTF-8 래퍼 + 이모지 교체)
      - commit: a1b2c3d4
```
