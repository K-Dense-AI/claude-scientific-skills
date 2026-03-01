# Git Workflow Manager

Git 워크플로우 관리: pull, 브랜치 전략, 커밋 태그, 클린업을 수행하라.

## 핵심 원칙

수정 전에 반드시 git status && git pull 실행.

## 브랜치 전략

```
main              # 검증 완료된 코드만 (P4 ✅ 통과 후 머지)
├── dev           # 개발 중인 코드 (P2~P3 작업)
├── exp/{name}    # 실험별 브랜치 (P1 요건 단위)
└── fix/{issue}   # P4 ❌ → P2 회귀 시 수정 브랜치
```

## 커밋 태그

```
[P1] — 요건 문서
[P1.5] — 매핑 문서, Fork
[P2] — 코드 구현
[P3] — 에러 수정
[P4] — 검증 결과
[P5] — 사이클 로그
type: feat / fix / refactor / test / docs / chore
```

## Upstream 관리

```bash
git remote add upstream https://github.com/{original}/{repo}.git
git fetch upstream
git log HEAD..upstream/main --oneline
```

## 안전 규칙

- 파괴적 작업 금지 (force push 등)
- stash 전 경고
- 클린업은 dry-run 기본
