# Code Implementer (코드 구현) — P2

P1.5 매핑 문서 기반으로 Fork 레포에서 코드를 구현하라.

## 전제 조건

- P1.5 매핑 문서(docs/reference_mapping.md)가 존재
- Fork 레포가 clone되어 있음
- upstream remote가 설정되어 있음

## Step 1: 환경 세팅

Python 프로젝트 구조 표준 적용. pyproject.toml에 ruff/mypy/pytest 설정.

## Step 2: 브랜치 생성 & 작업 계획

main → dev → exp/{name} → fix/{issue}
매핑 문서의 "수정 필요" 항목을 의존성 순으로 분해.

## Step 3: 코드 수정 원칙

1. 레퍼런스 코드의 작동하는 구조를 유지 → 불필요한 리팩토링 금지
2. 변경 최소화 원칙 → 한 번에 하나의 수정 항목만 작업
3. SOLID 원칙 적용 → 새로 작성하는 코드에만 적용
4. 타입 힌트 추가 → 새로 작성/수정하는 함수에 필수

수정 흐름: 코드 수정 → ruff check → mypy → pytest → 커밋

## Step 4: 커밋 규칙

[P2] feat/fix/refactor/test/docs: 설명

## Step 5: 품질 게이트

ruff, mypy, pytest를 병렬 실행하여 모두 통과 확인.

## 산출물

수정된 코드 + 통과한 품질 게이트 + 커밋 로그 → P3로 전달
