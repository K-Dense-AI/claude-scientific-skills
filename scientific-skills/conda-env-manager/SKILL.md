---
name: conda-env-manager
description: Conda environment management for scientific Python projects. Diagnose environments (list, verify existence, detect conda/pip conflicts), scan project imports via AST, identify missing packages, resolve version conflicts (numba-numpy etc.), create/repair environments with conda-forge priority, and configure VS Code interpreter settings. Use when setting up new projects, migrating from PyCharm, troubleshooting import errors, or auditing environment health.
license: MIT license
metadata:
    skill-author: Jahyun Lee
---

# Conda Environment Manager

## Overview

프로젝트별 conda 환경의 생성, 진단, 복구를 자동화하는 스킬.
PyCharm → VS Code 이전, 환경 재구축, 의존성 충돌 해결 등에 사용.

## When to Use

- 프로젝트에서 `ImportError` / `ModuleNotFoundError` 발생 시
- 새 프로젝트 환경 세팅 시
- PyCharm에서 VS Code로 이전 시
- conda/pip 혼용으로 패키지 충돌 발생 시
- 환경이 삭제/손상되었는지 확인할 때

## Workflow

```
Step 1: 환경 진단 (check_env.py --list)
    → 모든 conda 환경 + 디스크 존재 여부 + Python 버전
    ↓
Step 2: 프로젝트 스캔 (scan_imports.py <dir> [env])
    → AST로 import 추출 → 설치된/누락된 패키지 분류
    ↓
Step 3: 충돌 체크 (check_env.py --check <env>)
    → conda/pip 중복 감지 + 알려진 버전 충돌 체크
    ↓
Step 4: 설치/수정
    → 누락 패키지 설치 (conda-forge 우선 → pip 폴백)
    → 버전 충돌 수정
    ↓
Step 5: VS Code 연결 (선택)
    → .vscode/settings.json에 인터프리터 경로 설정
```

## Scripts

### `scripts/scan_imports.py`

프로젝트 디렉토리의 `.py` 파일에서 third-party import를 AST로 추출.

```bash
# 기본 스캔
conda run -n base python scripts/scan_imports.py /path/to/project

# conda 환경과 비교하여 누락 패키지 확인
conda run -n base python scripts/scan_imports.py /path/to/project env_name
```

**기능:**
- stdlib vs third-party 자동 구분
- try/except 내 optional import 식별
- import 이름 → 패키지명 매핑 (PIL→Pillow, sklearn→scikit-learn 등)
- conda-forge 전용 패키지 표시 (rdkit, openmm 등)

### `scripts/check_env.py`

conda 환경 상태 진단.

```bash
# 전체 환경 목록 + 존재 여부
conda run -n base python scripts/check_env.py --list

# 특정 환경 진단 (중복/충돌 체크)
conda run -n base python scripts/check_env.py --check biosteam
```

**기능:**
- 등록된 환경 vs 실제 디스크 존재 확인
- conda/pip 동일 패키지 중복 감지
- 알려진 버전 충돌 자동 체크 (numba↔numpy 등)

### `references/known_conflicts.yml`

알려진 버전 충돌 패턴 DB. check_env.py가 참조.

## Installation Rules

패키지 설치 시 아래 우선순위를 따른다:

### 1. conda-forge 우선

```bash
# 좋은 예: conda-forge에서 pre-built 바이너리 설치
conda install -n env_name -c conda-forge rdkit numpoly

# 나쁜 예: pip로 C 확장 패키지 빌드 시도 (Windows에서 실패 가능)
pip install rdkit  # ❌ 빌드 실패
```

### 2. pip는 --no-deps로 최소한만

conda 패키지와 pip 패키지가 같은 의존성을 공유하면 충돌 발생.
pip는 필요한 패키지만 `--no-deps`로 설치.

```bash
# conda로 기반 설치 후
conda install -n env_name -c conda-forge numpoly

# pip로 상위 패키지만 (의존성은 conda가 관리)
conda run -n env_name pip install biorefineries --no-deps
```

### 3. conda/pip 중복 방지

```bash
# 확인: 같은 패키지가 conda와 pip 양쪽에 있으면
conda run -n env_name pip show numpoly  # pip 버전 확인
conda list -n env_name numpoly          # conda 버전 확인

# 해결: pip 버전 제거 → conda 버전만 유지
conda run -n env_name pip uninstall numpoly -y
conda install -n env_name -c conda-forge numpoly --force-reinstall
```

## Known Issues

### Windows에서 C 확장 빌드 실패

**증상:** `Microsoft Visual C++ 14.0 or greater is required`
**원인:** Build Tools 미설치
**해결:** conda-forge에서 pre-built 바이너리 사용

### numba ↔ numpy 버전 불일치

**증상:** `Numba needs NumPy 2.3 or less`
**원인:** pip가 numpy 최신(2.4+) 설치, numba는 2.3까지만 지원
**해결:** `pip install "numpy<=2.3"`

### conda/pip 유령 패키지

**증상:** conda install 후에도 이전 pip 버전이 import됨
**원인:** pip의 dist-info가 conda 설치를 가림
**해결:** `pip uninstall <pkg>` → `conda install --force-reinstall`

### Claude Code에서 conda activate 불가

**증상:** `conda activate` 명령이 동작하지 않음
**원인:** Claude Code bash가 `~/.bashrc`를 로드하지 않음
**해결:** `conda run -n env_name <command>` 사용

## VS Code Integration

프로젝트별 인터프리터를 `.vscode/settings.json`에 설정:

```json
{
    "python.defaultInterpreterPath": "C:\\Users\\Jahyun\\anaconda3\\envs\\biosteam\\python.exe"
}
```

또는 VS Code 하단 상태바 클릭 → `Python: Select Interpreter` → 환경 선택.
