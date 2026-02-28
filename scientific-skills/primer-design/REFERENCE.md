# Primer Design Suite — 상세 레퍼런스

## 신규 환경 설치

```bash
git clone https://github.com/jahyunlee00299/claude-scientific-skills ~/.claude/claude-scientific-skills
cp -r ~/.claude/claude-scientific-skills/scientific-skills/primer-design ~/.claude/skills/primer-design
pip install -r ~/.claude/skills/primer-design/requirements.txt
# ~/.claude/settings.json에 MCP 등록 (SKILL.md 참조)
```

> `vector_dna_config.py`의 OneDrive 경로는 각 컴퓨터에 맞게 수동 수정 필요.

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
```

---

### Tool 3: `suggest_colony_pcr`

| 파라미터 | 설명 |
|---------|------|
| `vector_name` | 벡터 이름 (fuzzy matching 지원) |
| `insert_length_bp` | Insert 길이 (bp) |

반환: F/R 프라이머 이름·서열, Tm, 예상 band size, 권장 annealing 온도

---

### Tool 4: `analyze_expression`

| 항목 | 설명 |
|------|------|
| `basic_info` | 단백질 길이 (aa), MW (kDa), GC% |
| `cai` | Codon Adaptation Index (0–1, E. coli K-12 기준) |
| `rare_codons` | 희귀 코돈 빈도 (%), cluster 위치 |
| `signal_peptide` | Signal peptide 예측 |
| `map_removal` | N-terminal Met 제거 예측 |
| `strain_recommendation` | 추천 발현 균주 |

---

### Tool 5: `check_reading_frame_tool`

| 파라미터 | 기본값 | 설명 |
|---------|--------|------|
| `vector_name` | 필수 | 벡터 이름 |
| `re_5prime` | 필수 | 5' 제한효소 |
| `re_3prime` | 필수 | 3' 제한효소 |
| `insert_has_atg` | True | Insert 자체 ATG 포함 여부 |
| `insert_has_stop` | False | Insert 자체 stop codon 포함 여부 |
| `insert_cds_bp` | None | Insert 길이 (선택) |

반환: `in_frame_5prime/3prime`, `topology`, `linker_aa`, `frame_report`, `warnings`

---

### Tool 6: `list_vectors`

| 벡터 | 태그 | 특징 |
|------|------|------|
| `pET-21a(+)` | C-His6 | T7/lac, IPTG |
| `pET-28a(+)` | N-His6+T7tag, C-His6 | T7/lac, IPTG |
| `pMAL-c6T` | N-MBP-TEV | Ptac, amylose 정제 |
| `pETDuet-1:MCS1` | N-His6 | T7/lac, co-expression |
| `pETDuet-1:MCS2` | S-tag (optional) | co-expression with MCS1 |
| `pACYCDuet-1:MCS1` | N-His6 | CmR, ColA ori, pET 호환 |
| `pACYCDuet-1:MCS2` | S-tag (optional) | CmR, co-expression |

### Tool 7: `list_restriction_enzymes`

지원: BamHI(-HF), XhoI, NdeI, NheI(-HF), NotI(-HF), NcoI(-HF), EcoRI(-HF), HindIII(-HF), SalI(-HF), KpnI(-HF), BglII, EcoRV(-HF), MfeI(-HF), SacI(-HF), AvrII, FseI, AscI, PacI 등

### Tool 8: `generate_macrogen_order`

| 파라미터 | 기본값 | 설명 |
|---------|--------|------|
| `primers` | 필수 | `[{"name": "...", "sequence": "..."}, ...]` |
| `project_name` | `"primer_order"` | 파일명 prefix |
| `output_dir` | None | 저장 디렉토리 |

반환: `file_path`, `total_primers`, `total_length_nt`, `estimated_cost_krw`

---

## Python 직접 사용 — Mutagenesis

### 치환 — `iPCRSubstDesigner`

```python
from src.primer_design import iPCRSubstDesigner

designer = iPCRSubstDesigner()
result = designer.design(
    seq=template_seq,
    subst_pos=100,
    old_seq="ACG",
    new_seq="GAT",
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
Tm: tail + annealing 전체 effective binding
k1 = (overlap_len - len(new_seq)) // 2
k2 = overlap_len - len(new_seq) - k1
```

Hairpin 회피: 2-pass 로직 (Pass 1 실패 시 모든 후보 PASS>WARNING>FAIL 랭킹)

### 결실 — `iPCRDelDesigner`

```python
from src.primer_design import iPCRDelDesigner

designer = iPCRDelDesigner()
result = designer.design(
    seq=template_seq,
    del_start=30,
    del_end=33,
    target_tm=61.0,
    overlap_len=18,
    min_len=18,
    max_len=35,
    cds_start=0,
)
```

Frameshift: `del_len % 3 == 0` → in-frame, 아니면 `frameshift_warning: True`

### SnapGene 파일 로드

```python
from src.primer_design import parse_snapgene
seq, is_circular, features = parse_snapgene(r"path\to\vector.dna")
```

---

## 주의사항

### SDM 코돈 선택

```
Homopolymer ≥4 bp? → YES: 동의 코돈 변경 (발현 영향 <5%) / NO: E. coli 선호 코돈
```

| AA | 선호 | 차선 | Homopolymer 시 |
|----|------|------|---------------|
| Lys | AAA (76%) | AAG | AAA→AAG |
| Glu | GAA (68%) | GAG | GAA→GAG |
| Ile | ATT (51%) | ATC | ATT→ATC |

### Post-PCR: KLD 반응 (NEB #M0554)

```
1µL PCR + 1µL KLD Enzyme Mix + 5µL KLD Buffer + 3µL water → 25°C 5min → Transform
```

### 문제 해결

| 문제 | 해결 |
|------|------|
| `ValueError: Mismatch` | SnapGene에서 위치 재확인 |
| Tm 미달 | `max_len` 증가 또는 `target_tm` 하향 |
| `overlap_verified: False` | `k1/k2` 수동 확인 |
| QC FAIL 지속 | `overlap_len` 조정 |
| `primer3 not available` | `pip install primer3-py` |
| `expression_check FAIL` | Insert 서열 재확인, SDM으로 site 제거 |

---

## 스크립트 구조

```
src/primer_design/
├── subst_primer_mode.py         # iPCRDesignerBase + iPCRSubstDesigner
├── del_primer_mode.py           # iPCRDelDesigner
├── restriction_cloning_mode.py  # RestrictionCloningDesigner
├── colony_pcr_mode.py           # ColonyPCRDesigner
├── expression_analyzer.py       # ExpressionAnalyzer
├── order_sheet.py               # Macrogen XLSX
├── snapgene_parser.py           # .dna 파서
├── snapgene_writer.py           # .dna 생성
├── vector_registry.py           # 벡터/RE DB
├── vector_dna_config.py         # .dna 경로 매핑
├── cloning_report.py            # PNG 리포트
├── mcp_server.py                # FastMCP server (8 tools)
└── clients/
    ├── udh_variant_mutagenesis.py
    └── psxr_cofactor_mutagenesis.py
```
