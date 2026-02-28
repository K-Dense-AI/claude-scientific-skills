---
name: git-workflow-manager
description: Comprehensive Git workflow management for clean repositories. ALWAYS enforces git pull before modifications in any git repository, automatically cleans up test/temp files keeping only final versions, detects duplicate files and suggests consolidation. 연구 코딩 파이프라인용 브랜치 규칙, 커밋 태그([P1]~[P5]), upstream 추적 절차 포함.
---

# Git Workflow Manager

Enforces best practices: git pull before modifications, automatic test file cleanup, duplicate file detection, pipeline branch strategy.

## Core Principle

**ALWAYS check and update git repository BEFORE making any modifications.**

## Critical Workflow

### 1. Before ANY Modifications (MANDATORY)

```bash
git status
git pull origin dev  # 또는 main
```

### 2. During Development

테스트/임시 파일 정리: test_*, temp_*, *_backup.*, *_old.*, scratch_*, debug_*

### 3. After Completing Work

중복 파일 검출 및 정리.

---

## 파이프라인 브랜치 전략

```
main              # 검증 완료된 코드만 (P4 ✅ 통과 후 머지)
├── dev           # 개발 중인 코드 (P2~P3 작업)
├── exp/{name}    # 실험별 브랜치 (P1 요건 단위)
└── fix/{issue}   # P4 ❌ → P2 회귀 시 수정 브랜치
```

규칙:
- P2 작업은 `exp/{실험명}` 브랜치에서 수행
- P4 ✅ 통과 시 exp/ → dev → main 순서로 머지
- P4 ❌ 시 `fix/{issue-N}` 브랜치 생성 → 수정 → exp/에 머지
- main에 직접 커밋 금지

---

## 파이프라인 커밋 태그

```
형식: [Phase] type: 설명

Phase 태그:
  [P1]   — 요건 문서 작성/수정
  [P1.5] — 매핑 문서 작성, 레포 Fork
  [P2]   — 코드 구현, 리팩토링, 테스트 작성
  [P3]   — 에러 수정
  [P4]   — 검증 결과 기록
  [P5]   — 사이클 로그, 개선 메모

type:
  feat / fix / refactor / test / docs / chore
```

예시:
```
[P1.5] docs: reference_mapping.md 작성
[P2] feat: ribose 전환 모델 구현
[P3] fix: pH 계산 오류 수정 (수정 1/3)
[P4] docs: 3중 비교 결과 기록 ✅
[P5] log: 사이클 1 완료, 개선 메모 추가
```

---

## Fork 레포 Upstream 관리

```bash
# 초기 설정 (P1.5에서 Fork 후)
git remote add upstream https://github.com/{original}/{repo}.git

# P5에서 매 사이클 확인
git fetch upstream
git log HEAD..upstream/main --oneline

# 관련 있는 업데이트 머지 (다음 사이클 P2에서)
git checkout exp/{name}
git merge upstream/main
```

---

## 파이프라인 통합 워크플로우

```bash
# 1. 사이클 시작 — git 상태 확인
git status && git pull

# 2. P1.5 — Fork & 브랜치 생성
git checkout -b exp/{실험명}

# 3. P2 — 구현 (수정 항목별 커밋)
git commit -m "[P2] feat: {수정 내용}"

# 4. P3 — 에러 수정 시
git commit -m "[P3] fix: {에러 내용}"

# 5. P4 ❌ 회귀 시
git checkout -b fix/{issue-N}
# ... 수정 ...
git commit -m "[P2] fix: P4 회귀 — {내용}"
git checkout exp/{실험명}
git merge fix/{issue-N}

# 6. P4 ✅ 통과 시
git checkout dev && git merge exp/{실험명}
git checkout main && git merge dev

# 7. P5 — 로그 커밋
git commit -m "[P5] log: 사이클 N 완료"

# 8. 클린업
git push origin main
```

---

## Safety Features

- Never forces destructive operations
- Warns before stashing
- Dry-run by default for cleanup
- Duplication detector is read-only
