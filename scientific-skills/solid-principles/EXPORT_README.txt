===============================================================================
                    SOLID Principles Skill - Export Package
===============================================================================

버전: 1.0
생성일: 2025-11-28
작성자: Claude (Anthropic)
라이선스: MIT

===============================================================================
패키지 내용
===============================================================================

solid-principles/
├── SKILL.md                      # 메인 스킬 정의 (필수)
├── reference.md                  # 상세 예제 및 패턴 가이드
├── EXPORT_README.txt            # 이 파일
└── scripts/
    ├── solid_analyzer.py        # 코드 분석 도구
    └── solid_templates.py       # 템플릿 생성기


===============================================================================
빠른 설치 (Quick Install)
===============================================================================

Windows:
────────
1. 이 폴더를 복사:
   xcopy solid-principles C:\Users\<사용자명>\.claude\skills\solid-principles\ /E /I

2. 설치 확인:
   dir C:\Users\<사용자명>\.claude\skills\solid-principles


Mac/Linux:
──────────
1. 이 폴더를 복사:
   cp -r solid-principles ~/.claude/skills/

2. 실행 권한 부여:
   chmod +x ~/.claude/skills/solid-principles/scripts/*.py

3. 설치 확인:
   ls -la ~/.claude/skills/solid-principles


===============================================================================
스킬 설명
===============================================================================

이 스킬은 코드 작성 및 리뷰 시 SOLID 설계 원칙을 자동으로 적용합니다.

SOLID 5대 원칙:
───────────────
S - Single Responsibility Principle    (단일 책임 원칙)
O - Open/Closed Principle              (개방-폐쇄 원칙)
L - Liskov Substitution Principle      (리스코프 치환 원칙)
I - Interface Segregation Principle    (인터페이스 분리 원칙)
D - Dependency Inversion Principle     (의존성 역전 원칙)


===============================================================================
자동 감지 기능
===============================================================================

이 스킬은 다음 위반 사항을 자동으로 감지합니다:

[S] Single Responsibility
  ✗ 100줄 이상의 클래스
  ✗ 여러 책임을 가진 클래스
  ✗ 혼재된 관심사 (DB + UI + 검증)
  ✓ Extract Class 패턴 제안

[O] Open/Closed
  ✗ 긴 if/elif 체인 (3개 이상)
  ✗ 타입별 switch 문
  ✓ Strategy Pattern, Factory Pattern 제안

[L] Liskov Substitution
  ✗ 빈 메서드 오버라이드 (pass)
  ✗ NotImplementedError 발생
  ✓ 상속 계층 재설계 제안

[I] Interface Segregation
  ✗ 15개 이상의 메서드를 가진 인터페이스
  ✗ 불필요한 메서드 강제 구현
  ✓ 인터페이스 분리 제안

[D] Dependency Inversion
  ✗ __init__에서 직접 인스턴스 생성
  ✗ 구체 클래스에 대한 의존성
  ✓ Dependency Injection 패턴 제안


===============================================================================
사용 방법
===============================================================================

방법 1: 자동 활성화 (권장)
──────────────────────────
단순히 Claude Code와 대화하면 자동으로 활성화됩니다:

예시:
  "결제 시스템 설계해줘"
  "이 클래스 리팩토링해줘"
  "UserManager 클래스의 설계 문제 분석해줘"
  "의존성 주입 패턴으로 변경해줘"


방법 2: 코드 분석 도구 직접 실행
─────────────────────────────────
터미널에서:

# 코드 분석
python ~/.claude/skills/solid-principles/scripts/solid_analyzer.py your_file.py

출력 예시:
  ================================================================================
  SOLID Analysis: your_file.py
  ================================================================================

  [S] Single Responsibility Violations:
  --------------------------------------------------------------------------------

  🔴 Line 10 in UserManager
     Problem: Class has 150 lines (>100). Likely has multiple responsibilities.
     Fix: Apply Extract Class pattern to separate concerns.

  [D] Dependency Inversion Violations:
  --------------------------------------------------------------------------------

  🔴 Line 25 in UserService
     Problem: Direct instantiation of MySQLDatabase in __init__. Depends on concrete class.
     Fix: Inject dependency via __init__ parameter with Protocol type.

  ================================================================================
  Total violations: 2
  ================================================================================


방법 3: 템플릿 생성기 실행
──────────────────────────
터미널에서:

python ~/.claude/skills/solid-principles/scripts/solid_templates.py

대화형 메뉴:
  === SOLID Template Generator ===

  Choose a template:
    1. Protocol
    2. Dependency Injection
    3. Strategy Pattern
    4. Repository Pattern
    5. Factory Pattern
    6. Layered Architecture

  Enter choice (1-6):


===============================================================================
활용 예시
===============================================================================

예시 1: 설계 리뷰
─────────────────
사용자: "이 UserManager 클래스 리뷰해줘"

Claude: [solid-principles 스킬 자동 활성화]
        "이 클래스는 여러 SOLID 원칙을 위반하고 있습니다:

        1. [S] 단일 책임 원칙 위반: 데이터 저장, 이메일 발송, 검증을 모두 담당
        2. [D] 의존성 역전 원칙 위반: MySQLDatabase를 직접 생성

        리팩토링 제안:
        - User (엔티티)
        - UserRepository (데이터 저장)
        - EmailService (이메일)
        - UserValidator (검증)
        - UserService (조율)"


예시 2: 새 기능 설계
────────────────────
사용자: "다양한 결제 방법을 지원하는 시스템 설계해줘"

Claude: [solid-principles 스킬 자동 활성화]
        "Open/Closed 원칙을 따르는 설계를 제안합니다:

        ```python
        from typing import Protocol

        class PaymentMethod(Protocol):
            def process(self, amount: float) -> None: ...

        class CreditCardPayment: ...
        class PayPalPayment: ...

        class PaymentProcessor:
            def __init__(self, payment_method: PaymentMethod):
                self._payment_method = payment_method
        ```

        새 결제 방법 추가 시 기존 코드 수정 불필요합니다."


예시 3: 리팩토링
────────────────
사용자: "이 코드 SOLID 원칙으로 리팩토링해줘"

Claude: [solid-principles 스킬 자동 활성화]
        "다음 패턴을 적용하여 리팩토링하겠습니다:

        1. Extract Class: OrderProcessor → Order, OrderValidator, OrderRepository
        2. Dependency Injection: 의존성을 생성자로 주입
        3. Protocol 사용: 구체 클래스 대신 추상화에 의존

        리팩토링된 코드..."


===============================================================================
권장 프로젝트 구조
===============================================================================

이 스킬은 다음 계층 구조를 권장합니다:

project/
├── domain/                    # 비즈니스 로직 (의존성 없음)
│   ├── entities/             # 도메인 엔티티
│   │   └── user.py
│   └── protocols/            # 인터페이스 정의
│       └── repository.py
│
├── application/               # 유즈케이스 (domain에만 의존)
│   └── services/
│       └── user_service.py
│
└── infrastructure/            # 구현체 (domain 구현)
    ├── repositories/         # 데이터 접근
    │   └── user_repository.py
    └── external/             # 외부 서비스
        └── email_service.py

의존성 흐름:
  infrastructure → domain ← application
       ↓              ↑
       └──────────────┘
    (구현)         (사용)


===============================================================================
지원 언어
===============================================================================

주요 언어: Python
지원 버전: Python 3.8+

다른 언어:
- 원칙은 언어 독립적이므로 Java, C#, TypeScript 등에도 적용 가능
- 코드 분석 도구는 Python 전용


===============================================================================
문제 해결
===============================================================================

문제: 스킬이 활성화되지 않음
해결:
  ✓ SKILL.md 파일 위치 확인: ~/.claude/skills/solid-principles/SKILL.md
  ✓ YAML frontmatter 형식 확인
  ✓ Claude Code 재시작

문제: 스크립트 실행 오류
해결:
  # Python 버전 확인
  python --version  # 3.8 이상 필요

  # 권한 설정 (Mac/Linux)
  chmod +x ~/.claude/skills/solid-principles/scripts/*.py

문제: Protocol 타입 오류
해결:
  Python 3.8 미만인 경우:
  - typing_extensions 설치: pip install typing_extensions
  - 또는 Python 업그레이드


===============================================================================
커스터마이징
===============================================================================

프로젝트별 맞춤 설정:

1. SKILL.md의 description 수정
   - 사용하는 프레임워크 명시 (Django, FastAPI 등)
   - 프로젝트 특정 용어 추가

2. reference.md에 프로젝트 패턴 추가
   - 팀의 코딩 스타일
   - 자주 사용하는 패턴

3. scripts/에 커스텀 도구 추가
   - 프레임워크별 분석기
   - 프로젝트 템플릿 생성기


===============================================================================
팀과 공유
===============================================================================

방법 1: Git 저장소
─────────────────
# 프로젝트에 추가
cp -r ~/.claude/skills/solid-principles .claude/skills/

# 커밋 및 푸시
git add .claude/skills/solid-principles
git commit -m "Add SOLID principles skill"
git push

팀원들은 git pull 후 자동으로 사용 가능


방법 2: 압축 파일
─────────────────
이 폴더 전체를 압축하여 공유

받는 사람: 압축 해제 후 ~/.claude/skills/에 복사


===============================================================================
업데이트 내역
===============================================================================

Version 1.0 (2025-11-28)
  - 초기 릴리스
  - SOLID 5대 원칙 자동 감지
  - 코드 분석 도구 (solid_analyzer.py)
  - 템플릿 생성기 (solid_templates.py)
  - 상세 참조 문서 (reference.md)


===============================================================================
라이선스 및 저작권
===============================================================================

라이선스: MIT License

Copyright (c) 2025 Anthropic

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


===============================================================================
연락처
===============================================================================

이슈 및 피드백: https://github.com/anthropics/claude-code/issues
문서: https://code.claude.com/docs/en/skills.md
커뮤니티: https://github.com/anthropics/claude-code/discussions

생성: Claude Code
작성자: Claude (Anthropic)
===============================================================================
