# Inverse PCR Primer Design (iPCR 프라이머 설계)

원형 플라스미드에 대해 Inverse PCR을 이용한 삽입(insertion), 치환(substitution), 결실(deletion) 조작을 위한 프라이머를 설계하는 스킬.

---

## 개요

Inverse PCR은 원형 플라스미드 전체를 증폭하면서, 프라이머 5' tail에 삽입/치환/결실 서열을 포함시켜 선형 PCR 산물을 만든 후, KLD(Kinase-Ligase-DpnI) 또는 Gibson Assembly로 재원형화하는 방법이다.

**핵심 스크립트**: `C:/Users/Jahyun/PycharmProjects/pythonProject1/inverse_pcr_atg_insertion.py`
**보조 스크립트**: `C:/Users/Jahyun/PycharmProjects/pythonProject1/primer_quality_checker.py` (프라이머 품질 검증)

---

## 방법론 선택 가이드

작업 시작 전 아래 표로 적절한 방법을 먼저 결정한다.

| 상황 | 권장 방법 | 이유 |
|------|----------|------|
| 삽입 ≤ 6 nt | **NEB Q5 SDM (KLD)** | back-to-back, 5분 KLD, 효율 최고 |
| 삽입 6~100 nt | **NEB Q5 SDM (KLD)** | 삽입서열을 양쪽에 반씩 배치 |
| 삽입 > 100 nt | **Gibson Assembly** | 긴 서열은 별도 합성 + Gibson |
| 치환 (substitution) | **NEB Q5 SDM (KLD)** | mismatch를 F primer 중앙에 배치 |
| 결실 (deletion) | **NEB Q5 SDM (KLD)** | back-to-back non-mutagenic primer |
| Overlap 기반 cloning과 통합 | **Overlapping primer** | In-Fusion/Gibson과 동일 워크플로우 |

### NEB Q5 SDM (back-to-back) 방식과 현재 스크립트(overlapping)의 차이

| 항목 | NEB Q5 SDM (back-to-back) | 현재 스크립트 (overlapping) |
|------|--------------------------|--------------------------|
| Primer overlap | 없음 (back-to-back) | 18 bp |
| Amplification | Exponential | Exponential |
| 후처리 | KLD 5분 (PNK+Ligase+DpnI) | T4 PNK + T4 Ligase 또는 Gibson |
| Tm 계산 도구 | NEBaseChanger 권장 | BioPython Tm_NN (saltcorr=7) |
| 삽입 설계 | ≤6 nt: F primer 5'에 전부 / >6 nt: 양쪽 반씩 | overlap 영역 내 양쪽에 배분 |
| 효율 | 높음 | 중간 |

**→ 현재 스크립트는 overlapping 방식 구현. NEB Q5 SDM 모드는 향후 추가 예정.**

---

## 1단계: 입력 정보 수집

사용자에게 다음 정보를 순서대로 확인한다. 괄호 안은 기본값.

| 항목 | 설명 | 기본값 |
|------|------|--------|
| **DNA 파일 경로** | SnapGene .dna 파일의 절대 경로 | 없음 (필수) |
| **대상 유전자 이름** | CDS feature 이름 (예: PsFDHV9) | 없음 (필수) |
| **조작 유형** | insertion / substitution / deletion | insertion |
| **삽입/치환 서열** | 삽입할 DNA 서열 (예: AT, ATG, 전체 코돈) | 없음 (필수) |
| **삽입 위치 조정 (N)** | `ins_pos = gene_start - N` (아래 로직 참조) | 0 |
| **목표 Tm** | Annealing portion 기준 Tm (degC) | 61.0 |
| **중합효소** | PCR 중합효소 (Q5, Phusion 등) | Q5 |
| **Overlap 길이** | F/R 프라이머 5' tail의 overlap 길이 (bp) | 18 |
| **최소/최대 annealing 길이** | Annealing portion 탐색 범위 (bp) | 18 / 35 |

### 삽입 위치 조정 로직 (ins_pos)

```
ins_pos = gene_start - N
```

- **N = 0** (기본): 유전자 시작 위치 바로 앞에 삽입. 삽입 서열이 독립적으로 완전한 경우.
- **N = 1**: 유전자 시작 위치에서 1bp 업스트림. 업스트림 마지막 염기를 삽입 서열과 결합하여 새 코돈을 완성하는 경우.
  - 예: 업스트림 끝이 `...G`이고 `AT`를 삽입하면, `AT + G(업스트림)` = `ATG` (개시코돈)
  - 이 경우 FDH reading frame이 pos(ins_pos+1)부터 시작하여 원래 frame 유지

**사용자에게 확인할 핵심 질문**: "삽입 서열이 업스트림의 마지막 염기와 결합하여 코돈을 완성하나요? (예: AT + 업스트림 G = ATG)"
- "예" -> N = 1
- "아니오" -> N = 0

---

## 2단계: 스크립트 설정 확인 및 수정

`inverse_pcr_atg_insertion.py` 파일의 설정 변수를 사용자 입력에 맞게 수정한다.

### 수정 대상 변수 (파일 상단 설정 블록)

```python
# 파일 경로
DNA_FILE = r"<사용자 입력 .dna 파일 경로>"

# 유전자 및 서열
TARGET_GENE  = "<사용자 입력 유전자 이름>"
INSERT_SEQ   = "<사용자 입력 삽입 서열>"

# Tm 및 프라이머 파라미터
TARGET_TM    = <목표 Tm>       # degC
OVERLAP_LEN  = <overlap 길이>  # bp
MIN_LEN      = <최소 길이>     # bp
MAX_LEN      = <최대 길이>     # bp
```

### ins_pos 로직 수정 (main 함수 내)

N 값에 따라 `ins_pos` 라인을 수정한다:

```python
# N = 0 (기본, 유전자 시작 바로 앞 삽입)
ins_pos = fdh_start

# N = 1 (업스트림 1bp 포함, 예: AT + G = ATG)
ins_pos = fdh_start - 1
```

### 중합효소 조건 수정 (Q5가 아닌 경우)

Q5 외 중합효소 사용 시 버퍼 조건을 변경한다:

```python
# Phusion 예시
Q5_Na    = 50    # mM
Q5_Mg    = 1.5   # mM
Q5_DNAC1 = 500   # nM
Q5_DNTPS = 0.2   # mM
```

---

## 3단계: 스크립트 실행

```bash
cd C:/Users/Jahyun/PycharmProjects/pythonProject1
python inverse_pcr_atg_insertion.py
```

### 정상 출력 예시

스크립트는 다음을 출력한다:
1. Template 정보 (파일명, 길이, circular/linear)
2. 대상 유전자 위치 및 현재 시작 코돈 상태
3. 삽입 위치 컨텍스트 (전후 서열)
4. Forward primer (iPCR_ATGins_F): 전체 서열, tail/insertion/annealing 분해, Tm, GC%, homopolymer, 3' clamp
5. Reverse primer (iPCR_ATGins_R): 전체 서열, tail/insertion/annealing 분해, Tm, GC%, homopolymer, 3' clamp
6. Overlap 상보성 검증
7. PCR 조건 (annealing 온도, extension 시간)

### 결과 파일

- `C:/Users/Jahyun/PycharmProjects/pythonProject1/inverse_pcr_ATG_primers.txt`
- 파일명은 `OUTPUT_TXT` 변수로 변경 가능

---

## 4단계: 결과 검증 체크리스트

스크립트 실행 후 다음 항목을 반드시 확인한다:

### 필수 확인 (FAIL 시 재설계 필요)

- [ ] **Overlap 상보성**: `Complementary: YES` 확인. NO이면 프라이머 설계 오류.
- [ ] **Reading frame 보존**: 삽입 후 단백질 서열이 의도한 대로인지 확인.
  - ATG 삽입의 경우: `M` + 원래 단백질 서열 (첫 아미노산부터) 이 유지되는지 확인
  - 삽입 길이가 3의 배수가 아닌 경우 frameshift 주의
- [ ] **Tm 적절성**: 목표 Tm +/-2 degC 이내. 범위 밖이면 annealing 길이 조정.

### 권장 확인 (WARNING 수준)

- [ ] **Homopolymer 없음**: 4bp 이상 연속 동일 염기 없음 확인. 있으면 시퀀싱 오류 위험. SDM 코돈 변경으로 해소 가능 (발현 효율 영향 무시 가능).
- [ ] **3' GC clamp**: 마지막 염기가 G 또는 C. 약한 3' end는 priming 효율 저하.
- [ ] **GC 함량**: Annealing portion의 GC%가 40-60% 범위 내.
- [ ] **프라이머 총 길이**: 일반적으로 45bp 이하 권장. 합성 오류율은 길이에 비례.

### Reading frame 확인 방법

스크립트 출력에서 다음을 확인:
```
ATG codon: [삽입서열](inserted) + [X](posNN, upstream) = ATG
FDH reading frame: posNN = XYZ (AA) ...
```

삽입 후 예상 단백질 서열:
- N=1 케이스: `M` (ATG 완성) + 원래 단백질 첫 코돈부터 reading frame 유지
- N=0 케이스: `M` (ATG 완전 삽입) + 원래 단백질 서열 그대로 추가

---

## 5단계: 추가 품질 검증 (선택)

설계된 프라이머의 심화 품질 분석이 필요하면 `primer_quality_checker.py`의 로직을 참고하여 다음을 수행:

- **Hairpin 분석**: 6bp 이상의 self-complementary stem 탐지 (dG < -4 kcal/mol 이면 FAIL)
- **3' Self-dimer**: 3' 말단 6bp가 프라이머 자체 서열과 상보적인지 확인
- **Heterodimer**: F/R 프라이머 간 3' 상보성 확인

---

## 주의사항

### 삽입 위치 로직 상세

1. **업스트림 G 활용 ATG 완성** (`ins_pos = fdh_start - 1`):
   - 벡터 서열에서 유전자 바로 앞의 마지막 염기가 G인 경우
   - `AT`만 삽입하면 `AT + G(기존)` = `ATG` 개시코돈 완성
   - 이 방식은 삽입 길이를 최소화하고, reading frame을 자연스럽게 보존
   - 예시: `...GGATCC-G-GCTAAA...` -> `...GGATCC-[AT]-G-GCTAAA...` -> ATG + GCT(Ala) + AAA(Lys)...

2. **완전 ATG 삽입** (`ins_pos = fdh_start`):
   - `ATG` 전체를 삽입
   - 업스트림 서열에 의존하지 않음
   - 3bp 삽입이므로 reading frame 보존 (3의 배수)

### SDM 코돈 선택 — E. coli 발현 최적화 vs Homopolymer 회피

SDM(치환/삽입)으로 아미노산을 도입할 때, **동의 코돈(synonymous codon)** 선택이 primer 품질에 영향을 줄 수 있다.

#### 주요 원칙

1. **E. coli 최적 코돈 ≠ 항상 최선**
   - E. coli K-12 codon usage 기준으로 선호 코돈이 primer 내 homopolymer run을 유발할 수 있음
   - 예: Lys(K) — AAA(선호, ~76%) vs AAG(차선, ~24%)
   - AAA를 쓰면 인접 서열에 따라 Ax5~6 run 발생 → primer 합성 불량, PCR 실패 위험

2. **단일 코돈 SDM에서의 발현 차이는 무시 가능**
   - 코돈 1개 변경의 발현량 차이: 실험적으로 구분 불가 수준 (<5%)
   - 문제가 되는 경우: 여러 희귀 코돈이 3개 이상 cluster될 때
   - 전체 유전자 최적화가 목표라면 gene synthesis (IDT/GenScript Codon Optimization Tool) 사용 권장

3. **Homopolymer 발생 구조 파악**
   - prefix에 mutation 코돈 포함 → 코돈 마지막 염기 + gene_seq 첫 염기 = junction
   - junction에서 동일 염기 연속이 발생하는지 반드시 확인
   - 예: prefix 끝 `...AAA` + gene 시작 `AAC...` = `AAAAA` (Ax5) → FAIL

#### 코돈 선택 판단 흐름

```
Homopolymer 발생? (≥5bp 연속 동일 염기)
├── YES → 동의 코돈으로 변경 (마지막 or 첫 번째 염기 변경)
│         발현 효율 차이 무시 가능 → 바로 변경
└── NO  → E. coli 선호 코돈 우선 사용
```

#### E. coli 자주 참조하는 코돈 쌍

| AA | 선호 코돈 | 차선 코돈 | Homopolymer 유발 시 |
|----|-----------|-----------|---------------------|
| Lys (K) | AAA (76%) | **AAG** | AAA → AAG |
| Glu (E) | GAA (68%) | **GAG** | GAA → GAG |
| Gln (Q) | CAG (65%) | **CAA** | — |
| Arg (R) | CGT (38%) | CGC (40%) | context 의존 |
| Ile (I) | ATT (51%) | ATC (39%) | ATT → ATC |

---

### Post-PCR 처리

- **KLD 반응**: T4 PNK (5' 인산화) + T4 DNA Ligase (blunt-end ligation) + DpnI (template 제거)
- **대안**: Gibson Assembly (overlap이 충분한 경우)
- KLD 후 형질전환 -> Colony PCR -> 시퀀싱 확인

### 일반적 문제 해결

| 문제 | 원인 | 해결 |
|------|------|------|
| CDS feature 미발견 | 유전자 이름 불일치 | SnapGene에서 feature 이름 정확히 확인 |
| Tm 목표 미달 | Annealing 구간의 AT-rich 서열 | MAX_LEN 증가 또는 TARGET_TM 하향 |
| Overlap 상보성 NO | ins_pos 또는 overlap 분할 로직 오류 | k1/k2 값 및 서열 수동 확인 |
| 3' clamp 실패 | 마지막 염기가 A 또는 T | 1bp 연장하여 G/C로 끝나게 조정 |

---

## 스크립트 구조 요약

`inverse_pcr_atg_insertion.py`의 주요 함수:

| 함수 | 역할 |
|------|------|
| `parse_snapgene()` | SnapGene .dna 바이너리 파일에서 서열, circular 여부, feature 목록 추출 |
| `calc_tm()` | Q5 버퍼 조건에서 nearest-neighbor Tm 계산 (Biopython `Tm_NN`) |
| `gc_pct()` | GC 함량 (%) 계산 |
| `homopolymer_run()` | 4bp 이상 homopolymer run 탐지 |
| `gc_clamp_ok()` | 3' 말단이 G/C인지 확인 |
| `design_primer()` | 주어진 anchor 위치에서 Tm 목표에 도달할 때까지 연장하여 최적 annealing 구간 결정 |
| `main()` | 전체 워크플로우: 파일 로드 -> 유전자 탐색 -> 삽입 위치 결정 -> overlap 분할 -> F/R 프라이머 설계 -> 검증 -> 출력 |

### Overlap 설계 원리

```
F primer: 5'-[upstream_tail(k1 bp)]-[INSERT_SEQ]-[annealing(FDH 방향, Tm 기준)]---3'
R primer: 5'-[FDH_RC(k2 bp)]-[INSERT_SEQ_RC]-[annealing(upstream 방향, Tm 기준)]---3'

Overlap region = upstream_tail + INSERT_SEQ + FDH_short
                 (k1 bp)        (삽입 길이)   (k2 bp)
                 총 OVERLAP_LEN bp

k1 = (OVERLAP_LEN - len(INSERT_SEQ)) // 2
k2 = OVERLAP_LEN - len(INSERT_SEQ) - k1
```

F primer의 5' overlap과 R primer의 5' overlap은 상보적(complementary)이어야 하며, 이를 통해 PCR 산물의 blunt end가 KLD 반응으로 재원형화된다.
