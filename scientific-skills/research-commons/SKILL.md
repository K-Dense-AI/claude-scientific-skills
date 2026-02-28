---
name: research-commons
description: 4개 연구 스킬(experiment-hub, research-discussion, research-assistant, manuscript-writer)의 공통 규칙. 인용 스타일, Storage 키 체계, 스킬 간 연동, 트리거 분기, 공통 주의사항을 정의.
---

# Research Commons (공통 규칙)

experiment-hub, research-discussion, research-assistant, manuscript-writer 4개 스킬이 공유하는 규칙.

---

## 인용 규칙 (APA 7th Edition)

### 본문 인용
- 저자 1명: (Kim, 2024)
- 저자 2명: (Kim & Park, 2024)
- 저자 3명 이상: (Kim et al., 2024)
- 직접 언급: Kim et al. (2024)에 따르면...

### 필수 요건
- 모든 참고문헌에 DOI 링크(https://doi.org/...) 필수 포함
- DOI 미확인 시 PubMed 링크 또는 URL 대체 제공
- DOI 누락 논문은 web_fetch로 PubMed 페이지에서 보완 시도 (핵심 5편 이내)
- 문헌 내용을 언급할 때 반드시 (저자, 연도) 표기
- 답변 끝에 참고문헌 목록(References) 포함, 각 항목에 DOI 명시

### 참고문헌 형식
```
References
──────────────────────────────────────────────────
Author, A. B., & Author, C. D. (Year). Article title. Journal Name, Volume(Issue), Pages. https://doi.org/xx.xxxx/xxxxxxx
```

### 저작권 준수
- 논문 내용 직접 인용 최소화, 요약/해석 중심
- 검색 결과는 웹 기반이므로 전체 본문 접근 제한적

---

## 스킬 간 데이터 연동

Claude Code에서는 파일 기반으로 데이터를 공유한다. 키 접두어로 데이터 출처를 식별:
- `exp:` — experiment-hub (실험 데이터)
- `disc:` — research-discussion (디스커션)
- `research:` — research-assistant (문헌/팩트체크)
- `ms:` — manuscript-writer (원고)

주요 연동 흐름: 실험결과 → 디스커션 해석 → 원고 작성, 리비전 → 추가실험 프로토콜

---

## 4-스킬 연동 요약

| 방향 | 트리거 | 흐름 |
|------|--------|------|
| 실험→디스커션 | 결과 기록 완료 | exp:log → 해석(M1) |
| 디스커션→실험 | 액션 아이템 도출 | disc:action → 프로토콜(M1) |
| 디스커션→연구보조 | 수치 검증/문헌 비교 필요 | 팩트체크/검색 → 판정 반환 |
| 실험→연구보조 | 프로토콜 작성(자동) | 메소드 검색(M4) → 참조 반환 |
| 디스커션→원고 | Discussion 구조화 완료 | disc:conclusion → ms:draft |
| 실험→원고 | Results/Figure 작성 | exp:result+viz → ms:draft |
| 원고→실험 | 리뷰어 추가실험 요청 | ms:revision → exp:protocol |

---

## 트리거 분기 규칙 (겹침 해소)

"비교" / "해석" 트리거가 SKILL-1과 SKILL-2에서 겹칠 수 있으므로 다음 규칙으로 분기한다.

### "해석" 관련 트리거 분기

| 사용자 표현 | 분기 대상 | 근거 |
|------------|----------|------|
| "왜 이런 결과?", "원인이 뭐야?", "변인별 해석" | SKILL-1 M6 (양상분석) | 정량적 분석 — 추세선, 통계, 패턴 분류 |
| "이 결과 어떻게 해석?", "의미가 뭐야?", "메커니즘이 뭐야?" | SKILL-2 M1 (결과해석) | 정성적 분석 — 메커니즘, 맥락, Interaction Matrix |
| "결과 해석부터 Discussion까지" | SKILL-2 M5 (체인) | 해석→팩트체크→비교→구조화 파이프라인 |

**원칙**: 데이터/통계 중심이면 SKILL-1, 의미/맥락/논문 중심이면 SKILL-2.

### "비교" 관련 트리거 분기

| 사용자 표현 | 분기 대상 | 근거 |
|------------|----------|------|
| "이전 실험이랑 비교", "실험 간 차이", "어떤 조건이 효과?" | SKILL-1 M7 (실험비교) | 내 실험 데이터 간 비교 |
| "문헌이랑 비교", "다른 연구랑 어떻게 달라?", "수율이 높은 편?" | SKILL-2 M2 (비교분석) | 내 결과 vs 문헌 데이터 비교 |

**원칙**: "내 실험끼리" = SKILL-1, "문헌과" = SKILL-2.

### "논문 작성" 관련 트리거 분기

| 사용자 표현 | 분기 대상 | 근거 |
|------------|----------|------|
| "Discussion 어떻게 쓸까?", "논의 구조 잡아줘" | SKILL-2 M3 (구조화) | 논리 흐름 설계 단계 |
| "Discussion 써줘", "원고 작성", "초안 작성" | manuscript-writer | 실제 원고 텍스트 생성 |
| "문법 체크", "에러 확인", "교정해줘" | manuscript-writer S2 | 프루프리딩 |
| "레퍼런스 정리", "참고문헌 포맷 변환" | manuscript-writer (refs) | 저널 스타일 변환 |
| "리비전 대응", "리뷰어 답변" | manuscript-writer (revision) | 포인트별 대응 |
| "추가실험 필요" (리비전 중) | manuscript-writer → exp:protocol | 연동 트리거 |

**원칙**: 구조/논리 설계 = SKILL-2, 실제 텍스트 생성/교정 = manuscript-writer.

---

## 학술 Figure 표기법 & Nomenclature

Figure/그래프 생성, 원고 교정, 시각화 시 반드시 적용. experiment-hub M5, manuscript-writer S2 등 모든 시각화/교정 작업에 해당.

### 1. 분자생물학 명명법 (Gene vs Protein)

| 구분 | 표기 | 폰트 | 예시 |
|------|------|------|------|
| **유전자(Gene)** | 이탤릭 | *italic* | *PsXYL1*, *lacZ*, *NOX* |
| **단백질(Protein)** | 로만(정체) | regular | PsXR, LacZ, NOX |
| **플라스미드/벡터** | 소문자 p + 로만 | regular | pETDuet-1, pUC19, pACYC |
| **종명(Species)** | 이탤릭 + 축약 | *italic* | *E. coli*, *S. stipitis*, *L. pentosus* |
| **종명(첫 등장)** | 이탤릭 + 전체 | *italic* | *Escherichia coli*, *Scheffersomyces stipitis* |

**핵심 원칙**: DNA/RNA 수준 = 이탤릭, 단백질 수준 = 로만. Construct diagram에서 유전자명 = 이탤릭.

### 2. 단백질 Engineering 변이 표기법

| 표기법 | 형식 | 예시 | 사용 맥락 |
|--------|------|------|----------|
| **단일문자(Single-letter)** | [WT잔기][위치][변이잔기] | E223A | Figure, 본문 (권장) |
| **삼문자(Three-letter)** | [WT잔기][위치][변이잔기] | Glu223Ala | 본문, 임상유전학 |
| **변이체 명명** | 단백질명-변이 | PsXR-E223A | Figure 제목 |
| **다중 변이** | 슬래시 구분 | PsXR-E223A/K270M | 본문 |
| **Codon 변이** | WT코돈→변이코돈 | GAA→GCA | Sanger 검증 figure |

**주의**: 단일문자 표기에서 아미노산 코드는 항상 대문자 (E, A, K, M 등). 위치 번호는 mature protein 기준 (signal peptide 제외 후 번호).

### 3. Construct 표기법

| 요소 | 형식 | 예시 |
|------|------|------|
| **벡터::유전자** | 벡터명::*Gene1*–*Gene2* | pETDuet-1::*PsXYL1*–*PsFDH* |
| **내성 마커** | 약어 + 위첨자 R | Amp^R, Cm^R, Kan^R |
| **복제원점** | 정식 명칭 | ColE1, p15A, pSC101 |
| **발현 요소** | P_{프로모터}, T_{터미네이터} | P_T7, T_T7 |
| **태그** | 약어 + 위치 | N-His₆, C-S-tag |
| **MCS** | 번호 표기 | MCS1, MCS2 (Duet 벡터) |

### 4. Figure 기술 표준

| 항목 | 기준 | 비고 |
|------|------|------|
| **DPI** | 최소 300, 권장 600 | 인쇄용 |
| **폰트** | Arial 또는 Helvetica | sans-serif 필수 |
| **최소 폰트 크기** | 6 pt (축소 후 기준) | 저널 가이드라인 확인 |
| **단일 컬럼 폭** | 90 mm (3.54 in) | Nature/Elsevier/ACS 공통 |
| **2단 컬럼 폭** | 180 mm (7.09 in) | 저널별 차이 있음 |
| **패널 라벨** | 대문자 볼드 A, B, C | 좌상단 배치 |
| **색상 모드** | CMYK (인쇄) / RGB (온라인) | 색맹 고려 팔레트 권장 |
| **PDF fonttype** | Type 42 (pdf.fonttype: 42) | 텍스트 편집 가능 유지 |
| **파일 형식** | PDF (벡터), PNG 600 DPI (래스터) | TIFF는 호환성 이슈 가능 |

### 5. Genetic Construct Diagram (SBOL Visual 표준)

Figure에 유전자 발현 construct를 그릴 때 SBOL Visual 심볼 사용:

| 요소 | 심볼 | 설명 |
|------|------|------|
| **CDS (유전자)** | 화살표 오각형 (Pentagon) | 전사 방향 표시, 유전자명 이탤릭 |
| **Promoter** | 꺾인 화살표 (Bent arrow) | 전사 시작 방향 표시 |
| **Terminator** | T자 모양 | 전사 종결 위치 |
| **RBS** | 반원 (Semicircle) | 리보솜 결합 부위 |
| **Tag** | 작은 사각형 | His₆, S-tag, FLAG 등 |
| **Backbone** | 수평선 | 벡터 골격 |
| **변이 부위** | 빨간 별표 (★) | SDM/engineering 부위 |

**색상 규칙**: 유전자별 고유 색상 일관 유지. 동일 유전자는 모든 construct에서 동일 색상.

### 6. Chromatogram / Sequencing Figure 규칙

| 항목 | 표준 |
|------|------|
| **Trace 색상** | G=검정, A=초록, T=빨강, C=파랑 (ABI 표준) |
| **변이 코돈** | 점선 박스로 강조 |
| **WT vs Mutant** | 상하 비교 배치 (WT 위, Mutant 아래) |
| **코돈 표기** | 변이 코돈 위에 3문자 표시 (GAA, GCA 등) |
| **아미노산 표기** | 코돈 위에 이탤릭 표시: Glu (E), Ala (A) |
| **읽기 방향** | 5'→3' 화살표 하단 표시 |

### 7. Sanger Sequencing 자동 검증 워크플로우

SDM (Site-Directed Mutagenesis) 등 변이 확인을 위한 Sanger 시퀀싱 결과를 자동으로 검증하는 프로토콜.

#### 기본 원칙
- **전수 검사**: 디렉토리 내 모든 `.ab1` 파일을 자동 스캔하여 빠짐없이 검증
- **패턴 기반 탐지**: 고정 위치 대신 **고유 flanking sequence anchor** 사용 (위치 독립적)
- **양방향 확인**: Forward/Reverse 프라이머 모두 검증하여 cross-validation
- **품질 우선 선택**: Figure 표시용 샘플은 변이 부위 Phred quality 기준 최고 품질 선택

#### 자동 탐지 패턴 설계
```python
import re
from Bio.Seq import Seq

# 변이 부위 양쪽의 고유한 DNA 서열을 anchor로 사용
# 예: PsXYL1 E223A → ...TCTTTCGTT[GAA/GCA]TTGAACC...
E223_RE = re.compile(r'TCTTTCGTT(...)TTGAACC')

def detect_mutation(seq_str):
    """Forward + RevComp 양방향 검색."""
    m = E223_RE.search(seq_str)
    if m:
        return m.group(1), m.start() + 10, 'fwd'  # codon, center_pos, strand
    rc = str(Seq(seq_str).reverse_complement())
    m = E223_RE.search(rc)
    if m:
        return m.group(1), None, 'rev'
    return None, None, None
```

#### Anchor 설계 규칙
| 항목 | 기준 |
|------|------|
| **Left anchor 길이** | 8-12 bp (codon 직전 invariant 서열) |
| **Right anchor 길이** | 6-10 bp (codon 직후 invariant 서열) |
| **고유성** | 해당 유전자 내에서 유일하게 매칭 (false positive 방지) |
| **변이 포함 금지** | Anchor 자체에 변이 부위를 포함하지 않음 |

#### 파일명 파싱 규칙
```
{Construct}_{Type}_{Clone#}_{Primer}.ab1
예: XRF_E_1_pET-upstream.ab1
    ├── Construct: XRF
    ├── Type: E (E223A mutant), S (Standard/WT)
    ├── Clone#: 1
    └── Primer: pET-upstream
```

#### Figure 표시 샘플 선택
- Forward-strand 읽기 방향 reads만 사용 (trace가 자연 방향으로 표시)
- Reverse primer reads는 검증용으로만 사용 (trace가 complement로 표시되어 혼동)
- 동일 조건 다수 클론 중 변이 부위 Phred quality가 가장 높은 것 선택

#### 검증 보고서 필수 항목
1. **전체 파일 목록**: 스캔된 모든 .ab1 파일 나열
2. **개별 결과**: 각 파일별 코돈, 상태(WT/MUT/N/A), 방향, 품질
3. **클론 요약**: 클론별 종합 판정 (CONFIRMED / NOT DETECTED)
4. **최종 집계**: `E223A verified: n/N clones`
5. **Figure 반영**: 검증 결과를 Figure에 stamp로 표기

### 8. Matplotlib 공통 설정 (Python)

```python
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 7,
    'axes.linewidth': 0.5,
    'figure.dpi': 600,
    'savefig.dpi': 600,
    'pdf.fonttype': 42,      # Type 42 (TrueType)
    'ps.fonttype': 42,
    'mathtext.default': 'regular',
})
```

**LaTeX 표기 팁**:
- 이탤릭 유전자명: `r'$\it{PsXYL1}$'`
- 위첨자: `r'Amp$^R$'`
- 아래첨자: `r'His$_6$'`
- Unicode ᴿ, ₆ 등은 폰트 호환 문제 → LaTeX mathtext 권장

---

## 공통 주의사항

- 모든 기록/분석은 공식 연구 노트를 대체하지 않음 (보조 도구)
- 해석/제안은 참고용이며, 최종 판단은 연구자가 수행
- 문헌 비교 시 실험 조건 차이를 반드시 고려
- 사실관계 확인 시 단일 출처에만 의존하지 않음 (최소 2개 독립 출처)
- persistent storage 데이터는 정기적으로 외부 백업 권장
- 민감한 미발표 데이터는 shared 모드 사용 금지
- 추세선 외삽(데이터 범위 밖 예측)은 신뢰도 낮음을 항상 명시
- 시각화 데이터 추출(그래프→수치)은 근사값이므로 원본 확인 권장
