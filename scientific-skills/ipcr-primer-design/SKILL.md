Base directory for this skill: C:\Users\Jahyun\.claude\skills\ipcr-primer-design

# iPCR & 클로닝 프라이머 설계 (UDH Primer Design Suite)

원형 플라스미드 mutagenesis(치환·결실), RE 클로닝, Colony PCR 스크리닝, E. coli 발현 최적화 분석, Macrogen 주문서 생성까지 지원하는 통합 프라이머 설계 스킬.

---

## 개요

**패키지 위치**: `C:/Users/Jahyun/PycharmProjects/UDH Clustering/src/primer_design/`
**MCP 서버**: `mcp_server.py` (FastMCP, stdio transport)

### MCP 서버 활성화 확인

`.claude/settings.json`에 등록 필요:

```json
{
  "mcpServers": {
    "primer-design": {
      "command": "python",
      "args": ["-m", "src.primer_design.mcp_server"],
      "cwd": "C:/Users/Jahyun/PycharmProjects/UDH Clustering"
    }
  }
}
```

---

## 기능 요약

| 기능 | 조작 유형 | MCP Tool / 클래스 |
|------|----------|------------------|
| 아미노산 치환 (SDM) | Mutagenesis | `iPCRSubstDesigner` |
| 코돈 결실 | Mutagenesis | `iPCRDelDesigner` |
| RE 클로닝 프라이머 | Cloning | `design_re_cloning_primers` |
| RE 쌍 추천 | Cloning | `recommend_re_pair` |
| Colony PCR 프라이머 | Screening | `suggest_colony_pcr` |
| E. coli 발현 최적화 분석 | Analysis | `analyze_expression` |
| Reading frame 검증 | Analysis | `check_reading_frame_tool` |
| 벡터 목록 조회 | Info | `list_vectors` |
| 제한효소 목록 조회 | Info | `list_restriction_enzymes` |
| Macrogen 주문서 (XLSX) | Output | `generate_macrogen_order` |

---

## 방법론 선택 가이드

| 상황 | 권장 방법 |
|------|----------|
| 아미노산 치환 (SDM, 1–수십 nt) | **iPCRSubstDesigner** |
| 코돈 결실 (in-frame 또는 frameshift) | **iPCRDelDesigner** |
| 유전자를 발현 벡터에 삽입 | **design_re_cloning_primers** |
| 어떤 RE 조합이 최적인지 모를 때 | **recommend_re_pair** 먼저 실행 |
| 클로닝 후 콜로니 확인 | **suggest_colony_pcr** |
| 발현 전 CDS 품질 사전 확인 | **analyze_expression** |
| 삽입 > 100 nt (ATG 삽입 포함) | Gibson Assembly 고려 |

---

## MCP Tool 상세

### Tool 1: `design_re_cloning_primers`

유전자를 발현 벡터에 클로닝하기 위한 RE 클로닝 프라이머를 설계한다. 내부적으로 reading frame 검증, 발현 viability 체크, SnapGene 파일 생성까지 수행한다.

**파라미터**

| 파라미터 | 타입 | 기본값 | 설명 |
|---------|------|--------|------|
| `insert_seq` | str | 필수 | Insert CDS 서열 (ATG ~ stop) |
| `re_5prime` | str | 필수 | 5' 제한효소 이름 (예: `"BamHI"`, `"NdeI"`, `"NheI"`) |
| `re_3prime` | str | 필수 | 3' 제한효소 이름 (예: `"XhoI"`, `"NotI"`, `"HindIII"`) |
| `vector_name` | str\|None | None | 벡터 이름 — 지정 시 reading frame 자동 검증 (예: `"pET-28a(+)"`) |
| `include_start_codon` | bool | True | F 프라이머에 ATG 포함 여부 (NdeI/NcoI는 RE site 자체에 ATG 포함됨) |
| `include_stop_codon` | bool | False | R 프라이머에 stop codon 추가 여부 (True → C-terminal 태그 차단됨) |
| `target_tm` | float | 62.0 | annealing 영역 목표 Tm (°C) |
| `gene_name` | str | `"Insert"` | 리포트 및 SnapGene feature 레이블 |
| `output_dir` | str\|None | None | PNG 리포트 + .dna 파일 저장 디렉토리 |

**프라이머 구조**

```
F: 5'-[protection(4-6 bp)]-[RE5 site]-[annealing(~20 bp)]---3'
R: 5'-[protection(4-6 bp)]-[RE3 site]-[RC(stop)?]-[RC annealing(~20 bp)]---3'

Tm 계산: annealing 영역만 (RE tail은 template에 결합하지 않으므로 Tm에 미기여)
Protection bases: 6-cutter → 4 bp, 8-cutter(NotI 등) → 6 bp (GCGC... 패턴, NEB 가이드라인)
```

**자동 경고 항목**

- Insert 내부에 사용한 RE site 존재 시: `"CRITICAL: RE cuts within insert"`
- NdeI/NcoI 사용 + `include_start_codon=True`: `"RE site itself provides start codon"`
- Compatible overhang 쌍 (BamHI+BglII, EcoRI+MfeI, SalI+XhoI 등): 양방향 삽입 경고
- Blunt-end enzyme (EcoRV): 방향성 없음 경고
- 프라이머 길이 > 60 nt, homopolymer, 3' GC clamp 부재

**반환 주요 키**

| 키 | 설명 |
|----|------|
| `f_full`, `r_full` | 최종 프라이머 서열 (tail + annealing) |
| `f_tail`, `r_tail` | protection + RE site + spacer |
| `f_ann`, `r_ann` | annealing 영역만 |
| `f_tm`, `r_tm` | Annealing 영역 Tm (°C) |
| `anneal_temp` | 권장 PCR annealing 온도 = min(Tm) + 1 (°C) |
| `f_qc`, `r_qc` | Hairpin/homodimer QC (PASS/WARNING/FAIL) |
| `frame_check` | Reading frame 검증 결과 (vector_name 지정 시) |
| `expression_check` | 발현 viability 분석 결과 — `verdict`: PASS/WARNING/FAIL |
| `report_image_path` | PNG 리포트 경로 (output_dir 지정 시) |
| `snapgene_path` | SnapGene .dna 파일 경로 (output_dir 지정 시) |
| `internal_re_sites_5/3` | Insert 내부 RE site 위치 목록 |
| `warnings` | 경고 메시지 목록 |

**`expression_check` 상세**

발현 viability 종합 판정으로, 다음 6가지를 분석한다:

1. **벡터 발현 context**: promoter, inducer, RBS, 발현 타입 (예: T7/lac, IPTG, cytoplasmic high-level)
2. **Fusion protein 토폴로지**: N/C-terminal 태그 정보, 태그-insert 연결 서열 (linker AA)
3. **발현 단백질 전체 서열**: MW 계산, N/C linker 포함한 fusion protein 길이
4. **Insert 내부 RE site**: 사용한 제한효소가 insert를 자르는지 확인 (`CRITICAL` 판정)
5. **Premature stop codon**: ORF 내 조기 종결 코돈 위치
6. **CDS 품질**: CAI, 희귀 코돈 빈도, 균주 추천, signal peptide

```
verdict = FAIL    → 치명적 문제 (내부 RE site, premature stop, frameshift, 길이 오류)
verdict = WARNING → 경미한 문제 (낮은 CAI, QC WARNING)
verdict = PASS    → 문제 없음
```

**사용 예시**

```python
result = design_re_cloning_primers(
    insert_seq="ATGCGTAACCTGGCG...TAA",
    re_5prime="NdeI",
    re_3prime="XhoI",
    vector_name="pET-28a(+)",
    include_start_codon=True,
    include_stop_codon=False,   # C-His6 태그 유지
    gene_name="AtUdh",
    output_dir="C:/Users/Jahyun/results/primers",
)
print(result["f_full"])           # 최종 F 프라이머
print(result["expression_check"]["verdict"])   # PASS/WARNING/FAIL
print(result["frame_check"]["frame_report"])   # reading frame 요약
```

---

### Tool 2: `recommend_re_pair`

Insert와 벡터를 입력하면 최적 RE 조합을 점수 순으로 추천한다. 벡터 MCS에 있는 모든 RE 쌍을 시험하여 아래 기준으로 점수를 매긴다.

**점수 체계**

| 기준 | 점수 |
|------|------|
| Double digest 가능 (동일 버퍼 또는 CutSmart ≥ 75%) | +3 |
| 양쪽 모두 CutSmart 100% | +2 |
| Insert-vector reading frame 호환 | +2 |
| HF variant 존재 | +1 |
| Insert 내부 RE site 존재 시 해당 쌍 제외 | — |

**반환값**: 점수 내림차순 목록, 각 항목에 `re_5prime`, `re_3prime`, `buffer`, `double_digest_ok`, `frame_ok`, `score`, `reason` 포함

```python
recs = recommend_re_pair(
    insert_seq="ATGCGT...",
    vector_name="pET-28a(+)",
    prefer_hf=True,
)
# [{"re_5prime": "BamHI", "re_3prime": "XhoI", "score": 8, ...}, ...]
```

---

### Tool 3: `suggest_colony_pcr`

클로닝 후 콜로니 스크리닝을 위한 PCR 조건을 제안한다. Macrogen 표준 프라이머 데이터베이스를 기반으로 벡터 타입에 맞는 유니버설 프라이머를 자동 매칭한다.

**파라미터**

| 파라미터 | 설명 |
|---------|------|
| `vector_name` | 벡터 이름 (fuzzy matching 지원 — `"pET28a"`, `"pETDuet MCS1"` 등) |
| `insert_length_bp` | Insert 길이 (bp) |

**반환값**

- F/R 프라이머 이름 및 서열
- 각 프라이머 Tm
- 예상 band size: insert 포함 vs 빈 벡터
- 권장 annealing 온도

---

### Tool 4: `analyze_expression`

CDS 서열 입력 시 E. coli 발현에 영향을 줄 수 있는 요소를 종합 분석한다.

**반환값 상세**

| 항목 | 설명 |
|------|------|
| `basic_info` | 단백질 길이 (aa), MW (kDa), GC% |
| `cai` | Codon Adaptation Index (0–1, E. coli K-12 기준) |
| `rare_codons` | 희귀 코돈 빈도 (%), cluster 위치 (연속 3개 이상 시 경고) |
| `signal_peptide` | Signal peptide 예측 (rule-based, 분비 발현 확인용) |
| `map_removal` | Met aminopeptidase에 의한 N-terminal Met 제거 예측 |
| `strain_recommendation` | 추천 발현 균주: BL21, Rosetta (tRNA 보충), Rosetta 2 |
| `warnings` | CAI < 0.2, 희귀 코돈 cluster, 비표준 서열 등 |

```python
result = analyze_expression(cds_seq="ATGCGT...TAA")
print(result["cai"])                    # 예: 0.72
print(result["rare_codons"]["rare_codon_pct"])   # 예: 8.3%
print(result["strain_recommendation"]["primary_strain"])  # 예: "BL21(DE3)"
```

---

### Tool 5: `check_reading_frame_tool`

RE 클로닝 전략의 reading frame 호환성을 클로닝 전에 미리 검증한다. `design_re_cloning_primers`에서 자동 호출되지만, 독립적으로도 사용 가능.

**파라미터**

| 파라미터 | 기본값 | 설명 |
|---------|--------|------|
| `vector_name` | 필수 | 벡터 이름 |
| `re_5prime` | 필수 | 5' 제한효소 |
| `re_3prime` | 필수 | 3' 제한효소 |
| `insert_has_atg` | True | Insert 자체 ATG 포함 여부 |
| `insert_has_stop` | False | Insert 자체 stop codon 포함 여부 |
| `insert_cds_bp` | None | Insert 길이 (선택, 3의 배수 검증용) |

**반환값**

- `in_frame_5prime` / `in_frame_3prime`: True/False
- `topology`: N-terminal 및 C-terminal 태그 구조 다이어그램
- `linker_5prime_aa` / `linker_3prime_aa`: 태그-insert 연결 아미노산 서열
- `frame_report`: 사람이 읽기 쉬운 reading frame 요약 텍스트
- `warnings`: 프레임 불일치, C-terminal 태그 차단 등

---

### Tool 6: `list_vectors`

지원 벡터 목록과 MCS의 RE site, 태그 정보를 조회한다.

**지원 벡터**

| 벡터 | 태그 | 특징 |
|------|------|------|
| `pET-21a(+)` | C-His6 | T7/lac, IPTG |
| `pET-28a(+)` | N-His6+T7tag, C-His6 | T7/lac, IPTG |
| `pMAL-c6T` | N-MBP-TEV | Ptac, amylose 정제 |
| `pETDuet-1:MCS1` | N-His6 | T7/lac, co-expression |
| `pETDuet-1:MCS2` | S-tag (optional) | co-expression with MCS1 |
| `pACYCDuet-1:MCS1` | N-His6 | CmR, ColA ori, pET 호환 |
| `pACYCDuet-1:MCS2` | S-tag (optional) | CmR, co-expression |

---

### Tool 7: `list_restriction_enzymes`

지원 제한효소 목록, recognition site, overhang 타입, cut 위치 조회.

**지원 효소 예시**: BamHI(-HF), XhoI, NdeI, NheI(-HF), NotI(-HF), NcoI(-HF), EcoRI(-HF), HindIII(-HF), SalI(-HF), KpnI(-HF), BglII, EcoRV(-HF), MfeI(-HF), SacI(-HF), AvrII, FseI, AscI, PacI 등

---

### Tool 8: `generate_macrogen_order`

설계된 프라이머 목록으로 Macrogen Oligo Order 형식 XLSX 파일을 생성한다.

**파라미터**

| 파라미터 | 기본값 | 설명 |
|---------|--------|------|
| `primers` | 필수 | `[{"name": "...", "sequence": "..."}, ...]` |
| `project_name` | `"primer_order"` | 파일명 prefix |
| `output_dir` | None (cwd) | 저장 디렉토리 |

**반환값**: `file_path`, `total_primers`, `total_length_nt`, `estimated_cost_krw`

```python
result = generate_macrogen_order(
    primers=[
        {"name": "iPCR_AtUdh_T166D_F", "sequence": "ctcctgtGATccggaaccc..."},
        {"name": "iPCR_AtUdh_T166D_R", "sequence": "ggttccggATCacaggag..."},
    ],
    project_name="UDH_Phase2",
    output_dir="C:/Users/Jahyun/results/primer_order",
)
print(result["file_path"])       # XLSX 파일 경로
print(result["estimated_cost_krw"])  # 예상 합성 비용
```

---

## Python 직접 사용 — Mutagenesis 모드

MCP 서버 없이 Python 클래스를 직접 import하여 사용.

### 치환 (Substitution) — `iPCRSubstDesigner`

```python
from src.primer_design import iPCRSubstDesigner

designer = iPCRSubstDesigner()   # Q5 buffer 기본 (Na=50mM, Mg=2mM, dnac1=500nM)
result = designer.design(
    seq=template_seq,            # 전체 플라스미드 서열
    subst_pos=100,               # 치환 시작 위치 (0-indexed)
    old_seq="ACG",               # 원래 서열 (검증용, 불일치 시 ValueError)
    new_seq="GAT",               # 치환 후 서열
    target_tm=61.0,
    overlap_len=18,
    min_len=18,
    max_len=35,
)
```

**설계 원리**

```
overlap = up_tail(k1 bp) + new_seq + dn_tail(k2 bp)

F: 5'-[up_tail]-[new_seq]-[dn_tail]-[annealing_downstream]-3'
R: 5'-RC(dn_tail)-RC(new_seq)-RC(up_tail)-[annealing_upstream]-3'

Tm 계산: tail이 template과 일치하므로 tail + annealing 전체 effective binding 사용
k1 = (overlap_len - len(new_seq)) // 2
k2 = overlap_len - len(new_seq) - k1
```

**Hairpin 회피 2-pass 로직**

1. Pass 1: 최단 + GC clamp 기준으로 annealing 영역 결정
2. QC verdict = FAIL이면 Pass 2 자동 재시도: 모든 후보를 PASS > WARNING > FAIL 순 랭킹 → 최적 annealing 선택
3. Pass 2에서도 FAIL이면 `"F hairpin: stem in dn_tail+template boundary, annealing length adjustment cannot resolve"` 경고 발생 → `overlap_len` 조정 필요

---

### 결실 (Deletion) — `iPCRDelDesigner`

```python
from src.primer_design import iPCRDelDesigner

designer = iPCRDelDesigner()
result = designer.design(
    seq=template_seq,
    del_start=30,    # 결실 시작 (0-indexed, inclusive)
    del_end=33,      # 결실 끝 (0-indexed, exclusive) → 3 bp 결실
    target_tm=61.0,
    overlap_len=18,
    min_len=18,
    max_len=35,
    cds_start=0,     # reading frame 판단용 CDS 시작 위치 (선택)
)
```

**설계 원리**

```
원본: ...AAATTT [DEL_REGION] GGGCCC...
F: 5'-[left_before_del(k1 bp)]-[annealing_downstream]-3'
R: 5'-RC[right_after_del(k2 bp)]-[annealing_upstream]-3'
overlap = left_before_del + right_after_del  (결실 영역 건너뜀)

Tm 계산: annealing 영역만 (tail은 결실 영역 건너편이라 초기 cycle에서 template과 연속 결합 불가)
```

**Frameshift 판정 로직**

| 조건 | 결과 |
|------|------|
| `del_len % 3 == 0` | In-frame deletion — 경고만 (`"X bp = Y codons removed"`) |
| `del_len % 3 != 0` | `frameshift_warning: True` |
| 1 bp 결실 + `cds_start` 제공 | 코돈 내 위치(1st/2nd/3rd)까지 판단 |
| 1 bp 결실 + `cds_start` 미제공 | `"provide cds_start to check codon position"` 경고 |

---

### SnapGene 파일 로드

```python
from src.primer_design import parse_snapgene

seq, is_circular, features = parse_snapgene(
    r"C:\Users\Jahyun\OneDrive - 고려대학교\저장소\...\vector.dna"
)
# seq: str          전체 서열
# is_circular: bool 원형 여부
# features: list    [{name, start, end, strand, type}, ...]
```

---

## 결과 검증 체크리스트

### 필수 (FAIL 시 재설계)

- [ ] **Overlap 상보성**: `overlap_verified: True` (mutagenesis 모드)
- [ ] **Reading frame**: 치환/결실 후 ORF가 의도대로인지 확인
  - 치환: `len(new_seq) == len(old_seq)` 이면 frame 유지
  - 결실: `frameshift_warning: False` 이거나, frameshift가 의도적인지 확인
- [ ] **RE cloning**: `expression_check.verdict != "FAIL"` (내부 RE site, premature stop 없음)
- [ ] **Overlap 검증**: `overlap_verified: True`

### 권장 (WARNING 수준)

- [ ] **QC 판정**: PASS 또는 WARNING 허용 — FAIL이면 재설계

| 판정 | 기준 |
|------|------|
| FAIL | Hairpin/homodimer Tm > anneal_temp - 10°C |
| WARNING | Hairpin/homodimer Tm > anneal_temp - 20°C |
| PASS | Hairpin/homodimer Tm ≤ anneal_temp - 20°C 또는 구조체 없음 |

- [ ] **Homopolymer**: 4 bp 이상 연속 동일 염기 없음 → 동의 코돈 변경으로 해소
- [ ] **3' GC clamp**: 마지막 염기가 G 또는 C
- [ ] **GC%**: Annealing 영역 40–60%
- [ ] **프라이머 길이**: 60 nt 이하 권장

---

## 주의사항

### SDM 코돈 선택 — E. coli 발현 최적화 vs Homopolymer 회피

```
Homopolymer 발생 (≥4 bp 연속)?
├── YES → 동의 코돈 변경 (단일 코돈 변경은 발현 영향 < 5% — 바로 변경)
└── NO  → E. coli 선호 코돈 우선
```

| AA | 선호 코돈 | 차선 코돈 | Homopolymer 발생 시 |
|----|-----------|-----------|---------------------|
| Lys (K) | AAA (76%) | AAG | AAA → AAG |
| Glu (E) | GAA (68%) | GAG | GAA → GAG |
| Ile (I) | ATT (51%) | ATC | ATT → ATC |
| Arg (R) | CGT/CGC | context 의존 | — |

### Post-PCR 처리

**KLD 반응** (NEB #M0554):

```
1 µL PCR product  +  1 µL 10x KLD Enzyme Mix
5 µL 2x KLD Reaction Buffer  +  3 µL nuclease-free water
→ 25°C, 5 min → Transform 5 µL
```

KLD = T4 PNK (5' 인산화) + T4 DNA Ligase (blunt-end ligation) + DpnI (template 제거)

### 일반적 문제 해결

| 문제 | 원인 | 해결 |
|------|------|------|
| `ValueError: Mismatch` | `old_seq`가 `seq[subst_pos]`와 불일치 | SnapGene에서 위치 재확인 |
| Tm 목표 미달 | AT-rich annealing 서열 | `max_len` 증가 또는 `target_tm` 하향 |
| `overlap_verified: False` | overlap 분할 오류 | `k1/k2` 수동 확인 |
| F/R QC FAIL 지속 | Hairpin이 tail-template 경계에 위치 | `overlap_len` 조정 |
| `primer3 not available` | primer3 미설치 | `pip install primer3-py` |
| `Unknown RE` | 제한효소 이름 오류 | `list_restriction_enzymes()`로 정확한 이름 확인 |
| `expression_check FAIL` | Insert 내 RE site 또는 premature stop | Insert 서열 재확인 또는 SDM으로 site 제거 |

---

## 스크립트 구조 요약

```
src/primer_design/
├── subst_primer_mode.py         # iPCRDesignerBase + iPCRSubstDesigner (치환)
├── del_primer_mode.py           # iPCRDelDesigner (결실)
├── restriction_cloning_mode.py  # RestrictionCloningDesigner (RE 클로닝)
├── colony_pcr_mode.py           # ColonyPCRDesigner (colony PCR)
├── expression_analyzer.py       # ExpressionAnalyzer (CAI, 희귀 코돈, 균주 추천)
├── order_sheet.py               # PrimerOrderSheet → Macrogen XLSX
├── snapgene_parser.py           # parse_snapgene() — .dna 바이너리 파서
├── snapgene_writer.py           # write_cloning_construct() — .dna 생성
├── vector_registry.py           # EXPRESSION_VECTORS, RESTRICTION_ENZYMES DB
├── vector_dna_config.py         # .dna 파일 경로 매핑
├── cloning_report.py            # PNG 리포트 생성
├── mcp_server.py                # FastMCP server (8 tools)
└── clients/
    ├── udh_variant_mutagenesis.py   # UDH variant mutagenesis 실사용 예시
    └── psxr_cofactor_mutagenesis.py # PsXR cofactor mutagenesis 예시
```

### `iPCRDesignerBase` 공통 유틸리티

| 메서드 | 역할 |
|--------|------|
| `calc_tm(seq)` | Nearest-neighbor Tm (Owczarzy 2008, saltcorr=7, Mg2+ direct) |
| `gc_pct(seq)` | GC% 계산 |
| `gc_clamp_ok(seq)` | 3' 말단 G/C 확인 |
| `homopolymer_run(seq, n=4)` | 4 bp 이상 homopolymer run 탐지 |
| `check_primer(seq, anneal_temp)` | Hairpin + homodimer Tm-based QC (primer3 사용) |
| `check_heterodimer(seq1, seq2)` | F-R heterodimer dG, Tm |
| `_design_annealing(...)` | Tm 목표 달성하는 최적 annealing 영역 탐색 (2-pass hairpin 회피) |
