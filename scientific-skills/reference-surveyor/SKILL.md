---
name: reference-surveyor
description: 퍼블릭 GitHub 레포에서 구현 사례를 조사하고 내 사례에 매핑하는 도구. (1) 논문/방법론의 공개 구현체 탐색, (2) 레포 평가(실행가능성·품질·유사도·라이선스), (3) 핵심 코드 패턴 분석, (4) 예시 실행 검증, (5) 매핑 문서 작성. P1.5 전담.
---

# Reference Surveyor (코드 사례 조사) — P1.5 전담

퍼블릭 레포 탐색, 코드 분석, 매핑 문서 작성을 위한 스킬.

## 파이프라인 위치

```
P1 (학문적 검토) → [P1.5 Reference Survey] → P2 (구현)
```

## 실행 방식

**팀 생성 금지. WebSearch + Bash 직접 실행으로 처리한다.**

- 레포 탐색: WebSearch 병렬 호출 (Papers With Code + GitHub 동시 검색)
- 레포 평가: Bash로 `gh repo view` 등 직접 확인
- 예시 실행: Bash로 clone → venv → pip install → 실행

## 전제 조건

- P1 요건 문서가 완성되어 있다
- 검색 키워드가 준비되어 있다

---

## Step 1.5-1: 레포 탐색 (팀 병렬)

### 검색 소스 우선순위

1. 논문 직접 연결 (Code Availability, supplementary)
2. Papers With Code
3. GitHub 직접 검색
4. 패키지 문서의 예시
5. 커뮤니티 (Stack Overflow, Reddit)

검색 규칙: 최소 2회, 1~6 단어, 동일 쿼리 반복 금지

---

## Step 1.5-2: 레포 평가 & 선별

각 5점 평가: 실행 가능성, 코드 품질, 사례 유사도
라이선스: MIT/Apache/BSD → PASS, 없음 → FAIL
Python: 3.10+ → PASS

상위 2~3개 선별.

---

## Step 1.5-3: 코드 분석

```
[구조]     프로젝트 디렉토리 구조는?
[진입점]   main 실행 흐름은?
[핵심 로직] 내가 필요한 알고리즘이 어느 파일, 어느 함수에?
[데이터 흐름] 입력 → 변환 → 출력의 데이터 타입과 형식은?
[설정]     파라미터를 어떻게 관리하는가?
[의존성]   어떤 라이브러리를 쓰는가?
```

---

## Step 1.5-4: 예시 실행 검증

Fork → Clone → venv → pip install → 레퍼런스 예시 실행 → 결과 저장

```
results/reference/
├─ reference_output.json
├─ reference_run_log.txt
└─ reference_params.yaml
```

---

## Step 1.5-5: 매핑 문서 작성

파일: `docs/reference_mapping.md`

핵심 내용:
- 원본 레포 정보 (URL, Fork, 커밋 해시, 라이선스)
- 매핑표 (레퍼런스 위치 ↔ 내 사례 위치 ↔ 수정 내용)
- 수정 불필요 / 수정 필요 항목 분류
- 의존성 차이
- 레퍼런스 실행 결과

---

## 사전 확인 (모르면 물어보기)

작업 시작 전 다음 항목이 불명확하면 **추측하지 말고 사용자에게 질문한다**:
- 검색 대상 방법론/알고리즘 범위
- 언어/프레임워크 선호 (Python만? R도 OK?)
- 라이선스 제약 (GPL 허용 여부)
- Fork 필요 여부 (읽기만 vs 실제 수정 예정)

## 완료 체크리스트 (끝나면 대조 보고)

```
✅ 레퍼런스 레포 N개 탐색, M개 선별
✅ 각 레포 평가 점수 (실행가능성/품질/유사도)
✅ 예시 실행 검증: 통과/실패
✅ 매핑 문서 docs/reference_mapping.md 작성
✅ Fork + upstream remote 설정
✅ 라이선스 확인
⬜ [누락 항목 있으면 여기에 명시]
```

## 종료 조건

산출물: Fork된 레포 + 매핑 문서 + 레퍼런스 실행 결과 → P2

## 스킬 간 연동

- 입력: P1 요건 문서, S3 검색 결과
- 출력: → P2 (code-implementer), → P3 (code-validator), → P4 (3중 비교)
