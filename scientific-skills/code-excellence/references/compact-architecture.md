# Compact Architecture - 실용적 아키텍처 가이드

새 코드 작성 시 실용적이고 유지보수 가능한 아키텍처를 위한 상세 가이드입니다.

## Core Philosophy

**실용성 우선**: 완벽한 설계보다 간단하고 동작하는 솔루션
**점진적 개선**: 필요할 때 복잡도를 추가
**명확성 우선**: 영리한 코드보다 이해하기 쉬운 코드

## File Size Guidelines

### 권장 크기

```python
# Python files: 목표 500줄, 최대 800줄
# 예시: 적절한 크기의 파일

# data_processor.py (약 450줄)
"""
Data processing module with multiple related functions.
Each function is focused and clear.
"""

class DataLoader:           # 150줄
    """데이터 로딩 담당"""
    pass

class DataValidator:        # 150줄
    """데이터 검증 담당"""
    pass

class DataTransformer:      # 150줄
    """데이터 변환 담당"""
    pass

# 총 450줄: 3개의 관련된 클래스, 하나의 관심사 (데이터 처리)
```

### 분할 신호

**다음의 경우 파일 분할 고려:**
- 파일이 800줄 초과
- 5개 이상의 클래스
- 20개 이상의 함수
- 3개 이상의 서로 다른 관심사

```python
# ❌ BAD: 너무 많은 관심사가 섞임
# application.py (1200줄)

class UserAuthentication: ...  # 인증
class DatabaseConnection: ...  # 데이터베이스
class EmailService: ...        # 이메일
class ReportGenerator: ...     # 리포트
class PaymentProcessor: ...    # 결제
# 너무 많은 관련 없는 것들!

# ✅ GOOD: 관심사별로 분할
# auth/authentication.py (300줄)
class UserAuthentication: ...

# db/connection.py (250줄)
class DatabaseConnection: ...

# notifications/email.py (200줄)
class EmailService: ...

# reports/generator.py (400줄)
class ReportGenerator: ...

# payments/processor.py (350줄)
class PaymentProcessor: ...
```

## Class Size Guidelines

### 권장 크기

**일반 클래스**: 목표 225줄, 최대 450줄
**데이터 클래스**: 최대 150줄
**복잡한 로직 클래스**: 목표 300줄, 최대 600줄

### 좋은 예시

```python
# ✅ GOOD: 적절한 크기의 클래스 (약 200줄)
class ImageProcessor:
    """이미지 처리 클래스 - 단일 책임"""

    def __init__(self, config: Config):
        self.config = config
        self.filters = {}

    def load_image(self, path: str) -> Image:
        """이미지 로드 (10줄)"""
        ...

    def resize(self, image: Image, width: int, height: int) -> Image:
        """크기 조정 (15줄)"""
        ...

    def crop(self, image: Image, box: tuple) -> Image:
        """이미지 자르기 (12줄)"""
        ...

    def rotate(self, image: Image, angle: float) -> Image:
        """이미지 회전 (10줄)"""
        ...

    def adjust_brightness(self, image: Image, factor: float) -> Image:
        """밝기 조정 (15줄)"""
        ...

    def adjust_contrast(self, image: Image, factor: float) -> Image:
        """대비 조정 (15줄)"""
        ...

    def apply_filter(self, image: Image, filter_name: str) -> Image:
        """필터 적용 (20줄)"""
        ...

    def save_image(self, image: Image, path: str) -> None:
        """이미지 저장 (12줄)"""
        ...

    # 총 약 200줄
    # 단일 책임: 이미지 처리
    # 응집도 높음: 모든 메서드가 이미지와 관련
    # 이해하기 쉬움: 각 메서드가 명확
```

### 나쁜 예시

```python
# ❌ BAD: God 클래스 (600줄)
class Application:
    """너무 많은 책임"""

    def __init__(self):
        self.db = Database()
        self.auth = Auth()
        self.email = Email()
        self.cache = Cache()
        # ... 더 많은 의존성

    def handle_http_request(self): ...     # 100줄
    def manage_database(self): ...          # 80줄
    def send_notifications(self): ...       # 70줄
    def generate_reports(self): ...         # 90줄
    def process_payments(self): ...         # 85줄
    def manage_users(self): ...             # 75줄
    def handle_webhooks(self): ...          # 100줄

    # 총 600줄
    # 문제: 여러 책임, 낮은 응집도, 이해하기 어려움
```

## Method Size Guidelines

### 기본 원칙

**메서드가 100줄을 넘으면**, 더 작게 나눌 수 있는지 검토하세요.

### 좋은 예시

```python
# ✅ GOOD: 복잡한 로직을 작은 메서드로 분할
class DataAnalyzer:

    def analyze_dataset(self, data: pd.DataFrame) -> dict:
        """전체 분석 프로세스 (15줄)"""
        cleaned_data = self._clean_data(data)
        validated_data = self._validate_data(cleaned_data)
        results = self._compute_statistics(validated_data)
        insights = self._generate_insights(results)

        return {
            'results': results,
            'insights': insights,
            'metadata': self._get_metadata(data)
        }

    def _clean_data(self, data: pd.DataFrame) -> pd.DataFrame:
        """데이터 정제 (25줄)"""
        ...

    def _validate_data(self, data: pd.DataFrame) -> pd.DataFrame:
        """데이터 검증 (30줄)"""
        ...

    def _compute_statistics(self, data: pd.DataFrame) -> dict:
        """통계 계산 (40줄)"""
        ...

    def _generate_insights(self, results: dict) -> list:
        """인사이트 생성 (35줄)"""
        ...

    def _get_metadata(self, data: pd.DataFrame) -> dict:
        """메타데이터 추출 (20줄)"""
        ...

    # 각 메서드가 50줄 이하로 명확하고 이해하기 쉬움
```

### 나쁜 예시

```python
# ❌ BAD: 하나의 거대한 메서드 (250줄)
class DataProcessor:

    def process(self, data):
        """모든 것을 한 번에 (250줄)"""
        # 데이터 정제 (50줄)
        for row in data:
            if row['value'] is None:
                row['value'] = 0
            # ... 더 많은 정제 로직

        # 데이터 검증 (50줄)
        errors = []
        for row in data:
            if row['value'] < 0:
                errors.append(row)
            # ... 더 많은 검증 로직

        # 통계 계산 (70줄)
        total = 0
        count = 0
        for row in data:
            total += row['value']
            count += 1
            # ... 더 많은 계산 로직

        # 인사이트 생성 (50줄)
        insights = []
        if total > 1000:
            insights.append("High total")
            # ... 더 많은 인사이트 로직

        # 결과 포맷팅 (30줄)
        result = {
            'data': data,
            'stats': {'total': total, 'count': count},
            'insights': insights
        }
        # ... 더 많은 포맷팅 로직

        return result

        # 문제: 너무 길고, 여러 일을 하고, 이해하기 어려움
```

## Class Design Principles

### 1. 단일 책임 원칙 (Single Responsibility)

**하나의 클래스는 하나의 이유로만 변경되어야 합니다.**

```python
# ✅ GOOD: 단일 책임
class UserRepository:
    """사용자 데이터 영속성만 담당"""

    def save(self, user: User) -> None:
        """사용자 저장"""
        ...

    def find_by_id(self, user_id: int) -> Optional[User]:
        """ID로 사용자 찾기"""
        ...

    def delete(self, user_id: int) -> None:
        """사용자 삭제"""
        ...

# ❌ BAD: 여러 책임
class UserManager:
    """너무 많은 책임"""

    def save_to_database(self, user): ...      # 데이터베이스
    def send_welcome_email(self, user): ...    # 이메일
    def validate_user_data(self, user): ...    # 검증
    def generate_user_report(self, user): ...  # 리포트
    # 4가지 서로 다른 이유로 변경될 수 있음!
```

### 2. 의존성 주입 (Dependency Injection)

**구체 클래스가 아닌 추상화에 의존하세요.**

```python
# ✅ GOOD: Protocol과 의존성 주입
from typing import Protocol

class Database(Protocol):
    def save(self, data: dict) -> None: ...
    def load(self, id: int) -> dict: ...

class UserService:
    def __init__(self, db: Database):
        self.db = db  # 추상화에 의존

    def create_user(self, user_data: dict) -> None:
        self.db.save(user_data)

# 테스트나 다른 구현으로 쉽게 교체 가능
class MySQLDatabase:
    def save(self, data: dict) -> None: ...
    def load(self, id: int) -> dict: ...

class PostgreSQLDatabase:
    def save(self, data: dict) -> None: ...
    def load(self, id: int) -> dict: ...

# ❌ BAD: 구체 클래스에 직접 의존
class UserService:
    def __init__(self):
        self.db = MySQLDatabase()  # 구체 클래스에 직접 의존
        # 테스트나 변경이 어려움
```

### 3. 상속보다 컴포지션 (Composition over Inheritance)

```python
# ✅ GOOD: 컴포지션
class Engine:
    def start(self): ...
    def stop(self): ...

class Wheels:
    def rotate(self): ...

class Car:
    def __init__(self):
        self.engine = Engine()  # 컴포지션
        self.wheels = Wheels()  # 컴포지션

    def drive(self):
        self.engine.start()
        self.wheels.rotate()

# ❌ AVOID: 깊은 상속 계층
class Vehicle: ...
class LandVehicle(Vehicle): ...
class MotorizedVehicle(LandVehicle): ...
class Car(MotorizedVehicle): ...
# 계층이 너무 깊고 복잡함
```

## Project Structure Example

### 좋은 프로젝트 구조

```
project/
├── src/
│   ├── __init__.py
│   ├── domain/              # 비즈니스 로직 (의존성 없음)
│   │   ├── __init__.py
│   │   ├── entities/        # 엔티티 클래스
│   │   │   ├── user.py      (150줄)
│   │   │   └── order.py     (200줄)
│   │   └── protocols/       # 인터페이스 정의
│   │       ├── repository.py (100줄)
│   │       └── service.py    (80줄)
│   │
│   ├── application/         # 유스케이스 (domain에 의존)
│   │   ├── __init__.py
│   │   ├── user_service.py  (300줄)
│   │   └── order_service.py (350줄)
│   │
│   ├── infrastructure/      # 구현 세부사항 (domain에 의존)
│   │   ├── __init__.py
│   │   ├── database/
│   │   │   ├── __init__.py
│   │   │   ├── connection.py (200줄)
│   │   │   └── repositories.py (400줄)
│   │   └── external/
│   │       ├── email_client.py (150줄)
│   │       └── payment_client.py (250줄)
│   │
│   └── utils/              # 유틸리티 함수
│       ├── __init__.py
│       ├── validators.py   (200줄)
│       └── formatters.py   (150줄)
│
├── tests/
│   ├── unit/               # 단위 테스트
│   │   ├── test_user.py    (400줄)
│   │   └── test_order.py   (450줄)
│   └── integration/        # 통합 테스트
│       └── test_api.py     (600줄)
│
├── cli.py                  # CLI 진입점 (200줄)
├── requirements.txt
├── setup.py
└── README.md
```

## When to Split Files/Classes

### 파일 분할 시점

```python
# 시나리오 1: 파일이 800줄 초과
# analysis.py (1200줄) → 분할

# 분할 후:
# analysis/
# ├── __init__.py
# ├── structure_analyzer.py  (400줄)
# ├── noise_analyzer.py      (400줄)
# └── comparison_analyzer.py (400줄)

# 시나리오 2: 관련 없는 클래스들
# utilities.py (600줄)
# class EmailSender: ...
# class DataValidator: ...
# class FileParser: ...

# 분할 후:
# utils/
# ├── email.py      (200줄)
# ├── validators.py (200줄)
# └── parsers.py    (200줄)
```

### 클래스 분할 시점

```python
# 시나리오 1: 클래스가 450줄 초과
class DataProcessor:  # 600줄
    def load(self): ...
    def validate(self): ...
    def transform(self): ...
    def analyze(self): ...
    def export(self): ...

# 분할 후:
class DataLoader:      # 150줄
    def load(self): ...

class DataValidator:   # 150줄
    def validate(self): ...

class DataTransformer: # 150줄
    def transform(self): ...

class DataAnalyzer:    # 150줄
    def analyze(self): ...

# 시나리오 2: 여러 책임을 가진 클래스
class UserManager:  # 400줄
    def authenticate(self): ...  # 인증
    def save_to_db(self): ...    # 저장
    def send_email(self): ...    # 알림

# 분할 후:
class AuthService:     # 150줄
    def authenticate(self): ...

class UserRepository:  # 150줄
    def save_to_db(self): ...

class NotificationService: # 100줄
    def send_email(self): ...
```

## Common Patterns

### 1. Factory Pattern

```python
class UserFactory:
    """사용자 객체 생성 담당"""

    @staticmethod
    def create_regular_user(email: str, name: str) -> User:
        return User(email=email, name=name, role='user')

    @staticmethod
    def create_admin_user(email: str, name: str) -> User:
        return User(email=email, name=name, role='admin')
```

### 2. Strategy Pattern

```python
class PaymentProcessor(Protocol):
    def process(self, amount: float) -> bool: ...

class CreditCardProcessor:
    def process(self, amount: float) -> bool: ...

class PayPalProcessor:
    def process(self, amount: float) -> bool: ...

class PaymentService:
    def __init__(self, processor: PaymentProcessor):
        self.processor = processor

    def make_payment(self, amount: float) -> bool:
        return self.processor.process(amount)
```

### 3. Repository Pattern

```python
class UserRepository(Protocol):
    def find_by_id(self, id: int) -> Optional[User]: ...
    def save(self, user: User) -> None: ...
    def delete(self, id: int) -> None: ...

class SQLUserRepository:
    def __init__(self, db: Database):
        self.db = db

    def find_by_id(self, id: int) -> Optional[User]:
        # SQL 구현
        ...

    def save(self, user: User) -> None:
        # SQL 구현
        ...
```

## Best Practices Summary

### DO ✅

- Protocol/Interface 먼저 정의
- 의존성 주입 사용
- 작고 집중된 클래스 작성 (225줄 목표)
- 작고 명확한 메서드 작성 (100줄 이하 권장)
- 단일 책임 원칙 준수
- 상속보다 컴포지션 사용
- 코드 중복 제거

### DON'T ❌

- God 클래스 만들기
- 구체 클래스에 직접 의존
- 깊은 상속 계층
- 거대한 메서드 (100줄 초과)
- 여러 책임을 가진 클래스
- 조기 최적화
- 불필요한 추상화

## Checklist

새 코드 작성 시 확인사항:

- [ ] 파일이 800줄 이하인가?
- [ ] 각 클래스가 450줄 이하인가?
- [ ] 각 메서드가 100줄 이하인가? (초과 시 분할 검토)
- [ ] 각 클래스가 단일 책임을 가지는가?
- [ ] Protocol/Interface를 정의했는가?
- [ ] 의존성 주입을 사용했는가?
- [ ] 코드가 이해하기 쉬운가?
- [ ] 테스트하기 쉬운가?
