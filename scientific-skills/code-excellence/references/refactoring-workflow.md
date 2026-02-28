# Refactoring Workflow - 5단계 리팩토링 가이드

기존 프로젝트를 체계적으로 분석하고 개선하기 위한 상세 워크플로우입니다.

## Overview

이 워크플로우는 **5단계**로 구성됩니다:
1. **분석** - 현재 상태 파악
2. **구조 리팩토링** - 파일/디렉토리 정리
3. **SOLID 검증** - solid-principles 스킬 호출
4. **중복 제거** - 공통 패턴 추출
5. **통합 & 테스트** - 최종 검증

**중요**: 각 단계는 순차적으로 진행하며, 사용자 승인 없이 다음 단계로 넘어가지 않습니다.

---

## Phase 1: 분석 (Analysis)

### 목표
현재 프로젝트의 상태를 완전히 이해하고 문제점을 식별합니다.

### Step 1.1: 프로젝트 구조 분석

**사용 도구**: Glob, Bash (ls, find 등)

```bash
# 모든 Python 파일 찾기
find . -name "*.py" -type f

# 파일 크기 확인
find . -name "*.py" -exec wc -l {} \; | sort -rn

# 디렉토리 구조 확인
tree -L 3
```

**체크리스트**:
- [ ] 모든 Python 파일의 위치 파악
- [ ] 각 파일의 줄 수 확인
- [ ] 중복된 파일명 찾기
- [ ] `__init__.py` 파일 존재 여부
- [ ] `requirements.txt`, `setup.py` 존재 여부
- [ ] 테스트 파일 위치 파악

**출력 예시**:
```
프로젝트 구조 분석 결과:

파일 개수: 45개
총 줄 수: 12,345줄

디렉토리 구조:
├── analysis.py (1,200줄) 🔴
├── utils.py (800줄) 🔴
├── visualizer.py (650줄) 🟡
├── pipeline/
│   ├── loader.py (400줄)
│   └── processor.py (500줄)
└── tests/ (흩어져 있음)

문제점:
🔴 2개 파일이 800줄 초과
🟡 1개 파일이 500줄 초과
❌ __init__.py 파일 없음
❌ requirements.txt 없음
```

### Step 1.2: 코드 품질 분석

**사용 도구**: Read, Grep

#### 파일 크기 분석

```python
# 분석 스크립트 예시
def analyze_file_sizes(files):
    """파일 크기 분석"""
    critical = []  # > 800줄
    moderate = []  # > 500줄
    good = []      # < 500줄

    for file in files:
        lines = count_lines(file)
        if lines > 800:
            critical.append((file, lines))
        elif lines > 500:
            moderate.append((file, lines))
        else:
            good.append((file, lines))

    return critical, moderate, good
```

**출력 예시**:
```
파일 크기 분석:

🔴 Critical (800줄 초과):
  - analysis.py: 1,200줄
  - utils.py: 800줄

🟡 Moderate (500줄 초과):
  - visualizer.py: 650줄
  - pipeline/processor.py: 550줄

🟢 Good (500줄 이하):
  - 41개 파일
```

#### 클래스 크기 분석

```python
def analyze_class_sizes(files):
    """클래스 크기 분석"""
    for file in files:
        classes = extract_classes(file)
        for cls in classes:
            lines = count_class_lines(cls)
            if lines > 450:
                print(f"🔴 {file}:{cls.name} = {lines}줄 (분할 권장)")
            elif lines > 300:
                print(f"🟡 {file}:{cls.name} = {lines}줄 (분할 검토)")
```

**출력 예시**:
```
클래스 크기 분석:

🔴 Critical (450줄 초과):
  - analysis.py:DataAnalyzer = 600줄
  - utils.py:DataProcessor = 500줄

🟡 Moderate (300줄 초과):
  - visualizer.py:Plotter = 350줄
```

#### SOLID 위반 사전 검사

**주의**: 이것은 개요만 확인합니다. 상세 검증은 Phase 3에서 solid-principles 스킬이 수행합니다.

```python
# 간단한 패턴 검사
patterns = {
    'god_class': r'class \w+Manager:',        # Manager 클래스
    'direct_import': r'from \w+ import \w+Database',  # 구체 클래스 import
    'isinstance': r'isinstance\(',            # 타입 체크
}
```

**출력 예시**:
```
SOLID 위반 가능성 (개요):

⚠️ God 클래스 패턴:
  - utils.py: class DataManager
  - api.py: class ApplicationManager

⚠️ 구체 클래스 직접 사용:
  - service.py: from db import MySQLDatabase

⚠️ isinstance 타입 체크:
  - processor.py: 15개 발견

Note: 상세 검증은 Phase 3에서 solid-principles 스킬이 수행합니다.
```

#### 코드 중복 분석

```python
def find_duplicate_code(files):
    """중복 코드 패턴 찾기"""
    # 유사한 함수명 찾기
    # 반복되는 로직 패턴 찾기
    pass
```

**출력 예시**:
```
코드 중복 분석:

🔴 중복 패턴 발견:
  - load_data() 함수가 3개 파일에 반복
  - validate_input() 로직이 5개 파일에 유사하게 반복
  - JSON 파싱 코드가 여러 곳에 반복
```

### Step 1.3: 개선 계획 생성

모든 분석 결과를 종합하여 **마크다운 리포트** 생성

```markdown
# 프로젝트 리팩토링 계획

## 현재 상태 요약

- 총 파일: 45개
- 총 줄 수: 12,345줄
- Critical 이슈: 5개
- Moderate 이슈: 8개

## 발견된 문제점

### 🔴 Critical Issues

1. **파일 크기 초과** (2개)
   - analysis.py: 1,200줄 → 3개 파일로 분할 필요
   - utils.py: 800줄 → 2개 파일로 분할 필요

2. **God 클래스** (2개)
   - DataAnalyzer: 600줄 → 여러 클래스로 분할 필요
   - DataProcessor: 500줄 → 여러 클래스로 분할 필요

3. **프로젝트 구조 부재**
   - __init__.py 파일 없음
   - src/ 디렉토리 없음
   - tests/ 구조화 안됨

### 🟡 Moderate Issues

1. **파일 크기 주의** (2개)
2. **클래스 크기 주의** (3개)
3. **코드 중복** (다수)

## 제안 디렉토리 구조

```
src/
├── __init__.py
├── analysis/
│   ├── __init__.py
│   ├── structure_analyzer.py  (from analysis.py)
│   ├── noise_analyzer.py      (from analysis.py)
│   └── comparison_analyzer.py (from analysis.py)
├── processing/
│   ├── __init__.py
│   ├── loader.py              (from utils.py)
│   └── transformer.py         (from utils.py)
└── visualization/
    ├── __init__.py
    └── plotter.py             (from visualizer.py)
```

## 리팩토링 단계

### Phase 2: 구조 리팩토링
- 디렉토리 구조 생성
- 파일 분할 및 이동
- import 경로 업데이트

### Phase 3: SOLID 검증
- solid-principles 스킬 호출
- SOLID 원칙 상세 검증
- 리팩토링 패턴 적용

### Phase 4: 중복 제거
- 공통 유틸리티 추출
- 기본 클래스 생성

### Phase 5: 통합 & 테스트
- CLI 진입점 생성
- import 검증
- 문서 업데이트

## 추정 작업량

- Phase 2: 중간 (파일 분할 및 이동)
- Phase 3: 중간 (SOLID 리팩토링)
- Phase 4: 낮음 (중복 제거)
- Phase 5: 낮음 (통합)

## 승인 요청

위 계획을 진행해도 될까요?
```

**사용자 승인 대기** → 승인 후 Phase 2로 진행

---

## Phase 2: 구조 리팩토링 (Structure Refactoring)

### 목표
파일과 디렉토리를 논리적으로 재구성하고 크기 제한을 준수합니다.

### Step 2.1: 디렉토리 구조 생성

```bash
# 새 디렉토리 구조 생성
mkdir -p src/{analysis,processing,visualization,utils}
mkdir -p src/{analysis,processing,visualization,utils}/__init__.py
mkdir -p tests/{unit,integration}
```

### Step 2.2: 파일 분할

#### 예시: analysis.py (1,200줄) → 3개 파일로 분할

**Before**:
```python
# analysis.py (1,200줄)
class StructureAnalyzer:  # 400줄
    ...

class NoiseAnalyzer:      # 400줄
    ...

class ComparisonAnalyzer: # 400줄
    ...
```

**After**:
```python
# src/analysis/structure_analyzer.py (400줄)
class StructureAnalyzer:
    ...

# src/analysis/noise_analyzer.py (400줄)
class NoiseAnalyzer:
    ...

# src/analysis/comparison_analyzer.py (400줄)
class ComparisonAnalyzer:
    ...

# src/analysis/__init__.py
from .structure_analyzer import StructureAnalyzer
from .noise_analyzer import NoiseAnalyzer
from .comparison_analyzer import ComparisonAnalyzer

__all__ = ['StructureAnalyzer', 'NoiseAnalyzer', 'ComparisonAnalyzer']
```

### Step 2.3: Import 경로 업데이트

```python
# Before
from analysis import StructureAnalyzer

# After
from src.analysis import StructureAnalyzer
```

**체크리스트**:
- [ ] 모든 import 문 업데이트
- [ ] 순환 의존성 확인
- [ ] 각 파일이 독립적으로 import 가능한지 확인

### Step 2.4: 파일 크기 최적화

**800줄 초과 파일 처리**:

```python
# 분할 전: data_processor.py (1,000줄)

# 분할 후:
# src/processing/
# ├── loader.py       (250줄) - 데이터 로딩
# ├── validator.py    (250줄) - 데이터 검증
# ├── transformer.py  (250줄) - 데이터 변환
# └── exporter.py     (250줄) - 데이터 내보내기
```

**450줄 초과 클래스 처리**:

```python
# 분할 전: class DataProcessor (600줄)

# 분할 후:
class DataLoader:      # 200줄
    """데이터 로딩 담당"""
    ...

class DataValidator:   # 200줄
    """데이터 검증 담당"""
    ...

class DataTransformer: # 200줄
    """데이터 변환 담당"""
    ...
```

### Step 2.5: 패키지 파일 생성

#### requirements.txt 생성

```python
def generate_requirements(src_dir):
    """import 문에서 requirements.txt 생성"""
    imports = set()
    for file in find_python_files(src_dir):
        imports.update(extract_imports(file))

    # 표준 라이브러리 제외
    external = filter_external_packages(imports)

    with open('requirements.txt', 'w') as f:
        for package in sorted(external):
            f.write(f"{package}\n")
```

#### setup.py 생성

```python
from setuptools import setup, find_packages

setup(
    name='your-project',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        # requirements.txt와 동기화
    ],
)
```

#### README.md 업데이트

```markdown
# Project Name

## Project Structure

```
src/
├── analysis/      # 데이터 분석 모듈
├── processing/    # 데이터 처리 모듈
└── visualization/ # 시각화 모듈
```

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```python
from src.analysis import StructureAnalyzer
```
```

---

## Phase 3: SOLID 원칙 검증 및 적용

### 목표
solid-principles 스킬을 호출하여 코드 품질을 검증하고 개선합니다.

### 사용자 알림

```
이제 **solid-principles 스킬**을 사용하여 코드 품질을 검증하고 개선하겠습니다.

solid-principles 스킬이 다음을 수행합니다:
✅ SOLID 위반 자동 감지
✅ 리팩토링 패턴 제안
✅ Protocol/Interface 정의
✅ 의존성 주입 적용
✅ 계층화된 구조 제안
```

### solid-principles 스킬 호출

**Note**: 이 단계는 solid-principles 스킬이 담당합니다.
code-excellence 스킬은 결과를 받아서 다음 단계로 진행합니다.

**검증 항목**:
1. 단일 책임 원칙 (SRP)
2. 개방-폐쇄 원칙 (OCP)
3. 리스코프 치환 원칙 (LSP)
4. 인터페이스 분리 원칙 (ISP)
5. 의존성 역전 원칙 (DIP)

**예상 출력**:
```
SOLID 검증 결과:

✅ 단일 책임 원칙: 대부분 준수
⚠️ 개방-폐쇄 원칙: 3개 위반 발견
⚠️ 의존성 역전 원칙: 5개 위반 발견

권장 리팩토링:
1. UserService에 Protocol 정의 필요
2. DatabaseConnection을 주입하도록 변경
3. Factory 패턴 도입 권장
```

---

## Phase 4: 코드 중복 제거 (Deduplication)

### 목표
반복되는 코드를 찾아서 공통 패턴으로 추출합니다.

### Step 4.1: 중복 코드 식별

```python
# 중복 패턴 예시

# File 1: loader_a.py
def load_data(path):
    with open(path, 'r') as f:
        data = json.load(f)
    return data

# File 2: loader_b.py
def load_json(filepath):
    with open(filepath, 'r') as f:
        content = json.load(f)
    return content

# File 3: parser.py
def read_json(file_path):
    with open(file_path, 'r') as f:
        result = json.load(f)
    return result
```

### Step 4.2: 공통 유틸리티 추출

```python
# src/utils/file_io.py
def load_json(path: str) -> dict:
    """JSON 파일 로드 (공통 유틸리티)"""
    with open(path, 'r') as f:
        return json.load(f)

# 다른 파일들에서 사용
from src.utils.file_io import load_json

data = load_json('data.json')
```

### Step 4.3: 기본 클래스 생성

```python
# 중복: 여러 Analyzer 클래스가 비슷한 구조

# Before
class StructureAnalyzer:
    def __init__(self, data):
        self.data = data
        self.results = {}

    def analyze(self):
        # 분석 로직
        pass

    def get_results(self):
        return self.results

class NoiseAnalyzer:
    def __init__(self, data):
        self.data = data
        self.results = {}

    def analyze(self):
        # 분석 로직
        pass

    def get_results(self):
        return self.results

# After: 기본 클래스 추출
class BaseAnalyzer:
    """모든 Analyzer의 기본 클래스"""

    def __init__(self, data):
        self.data = data
        self.results = {}

    def analyze(self):
        """서브클래스에서 구현"""
        raise NotImplementedError

    def get_results(self):
        return self.results

class StructureAnalyzer(BaseAnalyzer):
    def analyze(self):
        # 구조 분석 로직만
        pass

class NoiseAnalyzer(BaseAnalyzer):
    def analyze(self):
        # 노이즈 분석 로직만
        pass
```

### Step 4.4: 템플릿 메서드 패턴

```python
class DataProcessor:
    """템플릿 메서드 패턴"""

    def process(self, data):
        """공통 처리 흐름"""
        loaded = self.load(data)
        validated = self.validate(loaded)
        transformed = self.transform(validated)
        return self.export(transformed)

    def load(self, data):
        """서브클래스에서 구현"""
        raise NotImplementedError

    def validate(self, data):
        """서브클래스에서 구현"""
        raise NotImplementedError

    def transform(self, data):
        """서브클래스에서 구현"""
        raise NotImplementedError

    def export(self, data):
        """서브클래스에서 구현"""
        raise NotImplementedError
```

---

## Phase 5: 통합 & 테스트 (Integration & Testing)

### 목표
모든 변경사항을 통합하고 정상 작동을 검증합니다.

### Step 5.1: CLI 진입점 생성

```python
# cli.py
import argparse
from src.analysis import StructureAnalyzer
from src.processing import DataLoader
from src.visualization import Plotter

def main():
    parser = argparse.ArgumentParser(description='Data Analysis Tool')
    parser.add_argument('command', choices=['analyze', 'process', 'visualize'])
    parser.add_argument('--input', required=True)
    parser.add_argument('--output')

    args = parser.parse_args()

    if args.command == 'analyze':
        analyzer = StructureAnalyzer(args.input)
        results = analyzer.analyze()
        print(results)

    elif args.command == 'process':
        loader = DataLoader()
        data = loader.load(args.input)
        # ...

    elif args.command == 'visualize':
        plotter = Plotter()
        plotter.plot(args.input)

if __name__ == '__main__':
    main()
```

### Step 5.2: Import 검증

```python
# test_imports.py
def test_all_imports():
    """모든 import가 작동하는지 테스트"""

    # 주요 모듈 import
    from src.analysis import StructureAnalyzer
    from src.processing import DataLoader
    from src.visualization import Plotter

    # 순환 의존성 체크
    import sys
    import importlib

    modules = [
        'src.analysis.structure_analyzer',
        'src.analysis.noise_analyzer',
        'src.processing.loader',
        'src.processing.transformer',
    ]

    for module in modules:
        try:
            importlib.import_module(module)
            print(f"✅ {module}")
        except ImportError as e:
            print(f"❌ {module}: {e}")
```

### Step 5.3: 기능 테스트

```bash
# 기본 기능이 작동하는지 테스트
python cli.py analyze --input data.csv
python -m pytest tests/
```

### Step 5.4: 문서 업데이트

#### 아키텍처 문서

```markdown
# Architecture.md

## Project Structure

```
src/
├── analysis/      # 데이터 분석 모듈
│   ├── structure_analyzer.py  # 구조 분석
│   ├── noise_analyzer.py      # 노이즈 분석
│   └── comparison_analyzer.py # 비교 분석
├── processing/    # 데이터 처리 모듈
│   ├── loader.py      # 데이터 로딩
│   └── transformer.py # 데이터 변환
└── visualization/ # 시각화 모듈
    └── plotter.py     # 플로팅
```

## Design Patterns Used

- **Repository Pattern**: data access layer
- **Factory Pattern**: object creation
- **Strategy Pattern**: algorithm selection

## SOLID Principles

All modules follow SOLID principles:
- Single Responsibility: Each class has one clear purpose
- Dependency Inversion: Depends on abstractions (Protocols)
```

#### Docstring 업데이트

```python
class StructureAnalyzer:
    """
    구조 분석을 수행하는 클래스.

    이 클래스는 데이터의 구조적 패턴을 분석하여
    통계 정보와 인사이트를 제공합니다.

    Attributes:
        data: 분석할 데이터
        results: 분석 결과

    Example:
        >>> analyzer = StructureAnalyzer(data)
        >>> results = analyzer.analyze()
    """
    ...
```

---

## Checklist: 전체 리팩토링 완료 확인

### Phase 1: 분석
- [ ] 프로젝트 구조 파악 완료
- [ ] 파일/클래스 크기 분석 완료
- [ ] SOLID 위반 사전 검사 완료
- [ ] 코드 중복 분석 완료
- [ ] 개선 계획 생성 및 사용자 승인

### Phase 2: 구조 리팩토링
- [ ] 디렉토리 구조 생성 완료
- [ ] 800줄 초과 파일 분할 완료
- [ ] 450줄 초과 클래스 분할 완료
- [ ] 모든 import 경로 업데이트 완료
- [ ] `__init__.py` 파일 생성 완료
- [ ] requirements.txt 생성 완료
- [ ] setup.py 생성 완료

### Phase 3: SOLID 검증
- [ ] solid-principles 스킬 호출 완료
- [ ] SOLID 위반 수정 완료
- [ ] Protocol/Interface 정의 완료
- [ ] 의존성 주입 적용 완료

### Phase 4: 중복 제거
- [ ] 중복 코드 식별 완료
- [ ] 공통 유틸리티 추출 완료
- [ ] 기본 클래스 생성 완료

### Phase 5: 통합 & 테스트
- [ ] CLI 진입점 생성 완료
- [ ] 모든 import 테스트 통과
- [ ] 순환 의존성 없음
- [ ] 기능 테스트 통과
- [ ] 문서 업데이트 완료

### 최종 검증
- [ ] 모든 파일 < 800줄
- [ ] 모든 클래스 < 450줄
- [ ] SOLID 원칙 준수
- [ ] 중복 코드 제거
- [ ] 테스트 통과
- [ ] 문서 완비

---

## Common Pitfalls (주의사항)

### 1. Import 경로 놓치기
```python
# ❌ 업데이트 안 된 import
from analysis import StructureAnalyzer

# ✅ 업데이트된 import
from src.analysis import StructureAnalyzer
```

### 2. 순환 의존성
```python
# ❌ 순환 의존성
# a.py
from b import ClassB

# b.py
from a import ClassA

# ✅ 해결: 공통 인터페이스 사용
# protocols.py
class InterfaceA(Protocol): ...

# a.py
from protocols import InterfaceA

# b.py
from protocols import InterfaceA
```

### 3. __init__.py 누락
```python
# ✅ 모든 패키지에 __init__.py 필요
src/
├── __init__.py          # 필수!
├── analysis/
│   └── __init__.py      # 필수!
└── processing/
    └── __init__.py      # 필수!
```

### 4. 기존 기능 손상
```python
# ⚠️ 리팩토링 전 테스트 작성
def test_existing_functionality():
    """기존 기능이 유지되는지 확인"""
    # 리팩토링 전후 동일한 결과가 나와야 함
    pass
```

---

## Success Criteria

리팩토링이 성공적으로 완료되었는지 확인:

✅ **구조**
- 논리적 디렉토리 구조
- 모든 패키지에 `__init__.py`
- 명확한 모듈 분리

✅ **크기**
- 모든 파일 < 800줄
- 모든 클래스 < 450줄
- 읽을 때 25000 토큰 미만

✅ **품질**
- SOLID 원칙 준수
- 중복 코드 제거
- Protocol/Interface 정의됨

✅ **기능**
- 모든 import 작동
- 순환 의존성 없음
- 기존 기능 유지

✅ **문서**
- README 업데이트
- Architecture 문서
- Docstring 완비

리팩토링 완료! 🎉
