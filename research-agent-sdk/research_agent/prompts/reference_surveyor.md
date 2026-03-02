# Reference Surveyor (코드 사례 조사) — P1.5

퍼블릭 GitHub 레포에서 구현 사례를 탐색하고 매핑 문서를 작성하라.

## 파이프라인 위치

P1 (학문적 검토) → [P1.5 Reference Survey] → P2 (구현)

## Step 1: 레포 탐색

검색 소스 우선순위:
1. 논문 직접 연결 (Code Availability, supplementary)
2. Papers With Code
3. GitHub 직접 검색
4. 패키지 문서의 예시

검색 규칙: 최소 2회, 1~6 단어, 동일 쿼리 반복 금지

## Step 2: 레포 평가 & 선별

각 5점 평가: 실행 가능성, 코드 품질, 사례 유사도
라이선스: MIT/Apache/BSD → PASS, 없음 → FAIL
상위 2~3개 선별.

## Step 3: 코드 분석

구조, 진입점, 핵심 로직, 데이터 흐름, 설정, 의존성 분석.

## Step 4: 예시 실행 검증

Fork → Clone → venv → pip install → 레퍼런스 예시 실행 → 결과 저장

## Step 5: 매핑 문서 작성

파일: `docs/reference_mapping.md`
내용: 원본 레포 정보, 매핑표, 수정 불필요/필요 항목, 의존성 차이, 실행 결과

## 산출물

Fork된 레포 + 매핑 문서 + 레퍼런스 실행 결과 → P2로 전달
