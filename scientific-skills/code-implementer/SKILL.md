---
name: code-implementer
description: P1.5 매핑 문서를 기반으로 Fork 레포에서 내 사례를 구현하는 도구. (1) 매핑 문서→수정 작업 계획, (2) Python 프로젝트 구조 표준 적용, (3) ruff+mypy+pytest 도구 체인 강제, (4) 브랜치·커밋 규칙 관리, (5) 의존성 충돌 검증. PyCharm+GitHub+Python 환경 특화. P2 전담.
---

# Code Implementer (코드 구현) — P2 전담

P1.5 매핑 문서 기반으로 Fork 레포에서 수정 작업을 실행하는 스킬.

## 파이프라인 위치

```
P1.5 (코드 사례 조사) → [P2 Code Implement] → P3 (검수)
                              ↑
                    매핑 문서 + Fork 레포를 입력으로 받음
                    solid-principles, git-workflow와 연동
```

## 전제 조건

- P1.5 매핑 문서(docs/reference_mapping.md)가 존재한다
- Fork 레포가 clone되어 있다
- upstream remote가 설정되어 있다

---

## Step 2-1: 환경 세팅

### Python 프로젝트 구조 표준

```
project/
├── src/{package_name}/
│   ├── core/              # 핵심 로직 (수정 대상)
│   ├── utils/             # 유틸리티
│   └── analysis/          # 분석 모듈
├── tests/
│   ├── test_reference.py  # 레퍼런스 재현 테스트
│   └── conftest.py
├── data/raw/ & processed/
├── results/reference/ & current/
├── docs/reference_mapping.md & logs/
├── pyproject.toml
└── requirements.txt
```

### pyproject.toml 표준 설정

```toml
[tool.ruff]
line-length = 120
target-version = "py311"
select = ["E", "F", "W", "I", "N", "UP"]

[tool.mypy]
python_version = "3.11"
warn_return_any = true

[tool.pytest.ini_options]
testpaths = ["tests"]
addopts = "-v --tb=short"
```

---

## Step 2-2: 브랜치 생성 & 작업 계획

```
main → dev → exp/{name} → fix/{issue}
```

매핑 문서의 "수정 필요" 항목을 의존성 순으로 작업 단위 분해.

---

## Step 2-3: 코드 수정 원칙

1. 레퍼런스 코드의 작동하는 구조를 유지한다 → 불필요한 리팩토링 금지
2. 변경 최소화 원칙 → 한 번에 하나의 수정 항목만 작업
3. SOLID 원칙 적용 → 새로 작성하는 코드에만 적용
4. 타입 힌트 추가 → 새로 작성/수정하는 함수에 필수

### 수정 작업 흐름

```
코드 수정 → ruff check → mypy → pytest → 커밋
```

---

## Step 2-4: 커밋 규칙

```
[P2] feat: 기질을 xylose에서 ribose로 변경
[P2] test: ribose 전환 단위 테스트 추가
[P2] fix: numpy 버전 호환성 문제 해결
```

---

## Step 2-5: 품질 게이트

**팀 생성 금지. Bash 병렬 호출로 처리한다.**

ruff, mypy, pytest를 **한 메시지에서 Bash 3개 병렬 호출**:
```bash
# Bash 1        # Bash 2          # Bash 3
ruff check src/  mypy src/         pytest tests/
```

```
[필수]
✓ ruff check src/ → 에러 0
✓ mypy src/ → 에러 0
✓ pytest tests/ → 기본 테스트 통과
✓ 매핑 문서의 모든 "수정 필요" 항목에 체크(✓)
```

---

## 사전 확인 (모르면 물어보기)

작업 시작 전 다음 항목이 불명확하면 **추측하지 말고 사용자에게 질문한다**:
- 매핑 문서 위치 및 수정 항목 범위
- 브랜치 전략 (기존 브랜치 사용 vs 새로 생성)
- 의존성 버전 제약 (특정 라이브러리 버전 고정 여부)
- 테스트 범위 (전체 vs 수정 관련만)

## 완료 체크리스트 (끝나면 대조 보고)

작업 완료 시 매핑 문서의 "수정 필요" 항목을 하나씩 대조하여 보고한다:
```
✅ 항목 1: [수정 내용] — 커밋 완료
✅ 항목 2: [수정 내용] — 커밋 완료
⬜ 항목 3: [미완료 사유]
품질 게이트: ruff ✅ / mypy ✅ / pytest ✅
```

## 종료 조건

산출물: 수정된 코드 + 통과한 품질 게이트 + 커밋 로그 → P3에서 전체 실행 테스트

## 스킬 간 연동

- 입력: P1.5 (reference-surveyor), solid-principles, git-workflow-manager
- 출력: → P3 (code-validator), → P5 (improvement-analyzer)
