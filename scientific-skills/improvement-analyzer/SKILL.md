---
name: improvement-analyzer
description: 사이클 종료 시 개선점 도출 도구. (1) P1~P4 에러·회귀 패턴 분석, (2) 병목 구간 식별, (3) upstream 변경사항 확인, (4) 사이클 로그(docs/logs/cycle-N.md) 생성, (5) 다음 사이클 최적화 제안. P5 전담.
---

# Improvement Analyzer (개선점 분석) — P5 전담

사이클 종료 시 에러 패턴을 분석하고 다음 사이클을 개선하는 스킬.

## 실행 방식

**팀 생성 금지. 직접 순차 실행한다.**

Step 1–3은 파일 읽기+분석이므로 에이전트 없이 직접 처리한다.
upstream 확인만 Bash로 `git fetch upstream && git log upstream/main..HEAD --oneline` 실행.

## 워크플로우

```
P1~P4 전 과정 데이터 입력
    ↓
Step 1 + Step 2 + Step 3: 직접 순차 실행
    ↓
Step 4: 사이클 로그 생성 (docs/logs/cycle-N.md)
Step 5: 다음 사이클 개선 제안
    ↓
다음 사이클 P1 입력에 반영
```

---

## Step 1: 에러·회귀 패턴 분석

수집 대상: P3 에러 로그, P4 판정 결과, 회귀 이력, 디버그 루프 횟수.
반복 에러 유형, 회귀 원인 비율(코드 vs 요건), 에러 집중 모듈을 분석한다.

## Step 2: 병목 구간 식별

각 Phase의 반복 횟수를 기대치와 비교. P3 디버그 3회, P4 회귀 1회 이상이면 병목.

## Step 3: upstream 변경사항 확인

```bash
git fetch upstream
git log upstream/main..HEAD --oneline
```

## Step 4: 사이클 로그 생성

`docs/logs/cycle-{N}.md`에 P1~P4 전 과정 요약, 미해결 에러, 개선 메모를 기록.

## Step 5: 다음 사이클 개선 제안

문제 → 원인 → 제안 → 적용 Phase 형식으로 구체적 개선안 제시.

---

## 사전 확인 (모르면 물어보기)

작업 시작 전 다음 항목이 불명확하면 **추측하지 말고 사용자에게 질문한다**:
- 분석 대상 사이클 범위 (현재 사이클만? 전체?)
- 에러 로그 위치
- upstream 레포 URL/브랜치

## 완료 체크리스트 (끝나면 대조 보고)

```
✅ 에러·회귀 패턴 분석 완료
✅ 병목 구간 식별 완료
✅ upstream 변경사항 확인 완료
✅ 사이클 로그 docs/logs/cycle-N.md 생성
✅ 개선 제안 N건 도출
⬜ [누락 항목 있으면 여기에 명시]
```

## 종료 조건

사이클 로그 저장, 미해결 목록 작성, requirements.txt 갱신, 커밋 `[P5] log: 사이클 N 완료`.

## 스킬 간 연동

- session-continuity (컨텍스트 보존)
- research-discussion Record (디스커션 기록)
- git-workflow-manager (upstream 추적·파일 정리)
