먼저 ~/.claude/commands/research-commons.md 파일을 읽고 공통 규칙(인용 스타일, Storage 키 체계, 스킬 간 연동, 트리거 분기)을 숙지한 후 아래 스킬을 실행하세요.

---
name: research-assistant
description: 연구 논문 검색, 문헌 리뷰 자동화, 사실관계 확인, 메소드 참조 검색 도구. (1) 키워드 기반 논문 검색 및 요약, (2) 체계적 문헌 리뷰 수행, (3) 실험 데이터/주장의 사실관계 검증, (4) 논문/제조사 프로토콜/데이터시트 기반 메소드 참조 검색, (5) 연구 동향 분석에 사용. 효소공학, 바이오전환, 희귀당 생산 등 생명공학 연구에 최적화.
---

# Research Assistant (연구 보조)

논문 검색, 문헌 리뷰, 사실관계 확인, 메소드 참조 검색을 위한 연구 보조 스킬.

## 핵심 기능

1. **논문 검색 & 요약** — 키워드 기반 논문 탐색, 핵심 내용 요약
2. **체계적 문헌 리뷰** — 주제별 논문 분류, 연구 동향/갭 분석
3. **사실관계 확인** — 주장/수치/조건의 문헌 기반 검증
4. **메소드 참조 검색** — 논문 실험법, 분석 조건, 제조사 프로토콜/데이터시트 탐색

## 모드 1: 논문 검색 (Search)

### 다중 소스 검색 전략

검색은 반드시 **2회 이상** 실행하여 커버리지를 확보한다.

**필수 검색 (2회):**
```
쿼리 1: {핵심 키워드} site:pubmed.ncbi.nlm.nih.gov
쿼리 2: {핵심 키워드} {추가 키워드/동의어} review OR recent
```

**선택 검색 (필요 시):**
```
쿼리 3: {키워드} site:biorxiv.org
쿼리 4: {키워드} site:nature.com OR site:sciencedirect.com
쿼리 5: {EC번호 또는 기질명} enzyme characterization
```

### DOI 자동 보완: 핵심 논문 5편 이내에서 web_fetch로 PubMed DOI 추출

### 인용 규칙: APA 7th Edition, DOI 필수

---

## 모드 2: 체계적 문헌 리뷰 (Review)

최소 검색 3회, 목표 논문 10편 이상. 분류 → 동향 분석 → 갭 분석.

---

## 모드 3: 사실관계 확인 (Fact-Check)

### 판정 기준
- ✅ 확인됨: 2개 이상 독립 출처에서 확인
- ⚠️ 부분 확인: 수치는 확인되나 조건이 다름
- ❌ 불일치: 문헌 데이터와 명확히 다름
- ❓ 확인 불가: 관련 문헌을 찾을 수 없음

---

## 모드 4: 메소드 참조 검색 (Method Search)

### 3-카테고리 병렬 검색

**[A] 논문 메소드**: Materials & Methods에서 반응/분석 조건 추출
**[B] 분석/키트 프로토콜**: HPLC/GC 컬럼 제조사, 분석 키트 프로토콜
**[C] 제조사 데이터시트**: 효소/시약 최적 조건, 보관법, 안정성

### experiment-hub 연동 (자동 호출)

프로토콜 작성(experiment-hub M1) 시 자동으로 메소드 검색을 수행하여 조건의 근거를 확보한다.

---

## 스킬 간 연동

> 상세 연동 흐름도, 키 매핑, 트리거 분기 규칙은 **COMMONS.md** 참조.

### 본 스킬의 연동 포인트
- research-discussion 해석/비교 중 수치·주장 검증 → 팩트체크(M3) 호출
- research-discussion 체인 모드(M5) Step 2 → 팩트체크 자동 삽입
- experiment-hub 프로토콜 작성 시 → 메소드 검색(M4) 호출
- experiment-hub 실험 제안 시 → 논문 검색(M1) 호출
- 수집된 참고문헌(research:references) → manuscript-writer References 섹션에 저널 포맷으로 변환 제공
- manuscript-writer References 부족 감지 시 → 추가 논문 검색(M1) 트리거

---

## Storage 키
- `research:search:{timestamp}` — 검색 결과
- `research:review:{topic}` — 리뷰 보고서
- `research:factcheck:{timestamp}` — 팩트체크 결과
- `research:references:{id}` — 수집된 참고문헌 DB
- `research:method:{timestamp}` — 메소드 검색 결과
- `research:method:protocol:{protocol_id}` — 프로토콜 연결 메소드 참조

### manuscript-writer 연동

**참고문헌 DB → References 변환:**
```
research:references (APA 7th)
    ↓
manuscript-writer가 ms:refs로 조회
    ↓
저널 스타일에 맞게 변환 (Nature/ACS/Elsevier/Angewandte)
    ↓
본문 인용 마커도 동시 변환: (Author, Year) → [1] 등
```

**추가 검색 트리거:**
```
manuscript-writer 작성 중 References 부족 감지
  → "이 주장을 뒷받침하는 문헌이 부족합니다"
  → research-assistant M1(검색) 자동 호출
  → 결과를 research:references에 추가
  → ms:refs 업데이트
```

## 주의사항

- 검색 결과는 웹 검색 기반이므로 전체 논문 본문 접근은 제한적
- 팩트체크 시 수치뿐 아니라 실험 조건 차이도 반드시 분석
- 검색은 최소 2회 이상 실행하여 단일 소스 의존 방지
