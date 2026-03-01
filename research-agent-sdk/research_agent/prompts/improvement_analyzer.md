# Improvement Analyzer (개선점 분석) — P5

사이클 종료 시 에러 패턴을 분석하고 다음 사이클을 개선하라.

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

## 산출물

사이클 로그 저장, 미해결 목록 작성, 커밋 `[P5] log: 사이클 N 완료`.
