# Task: Primer Design

> 이 파일을 데이터 파일 없이 claude.ai에 첨부하고 전송하면 Claude가 Python으로 프라이머를 설계합니다.

---

## 설계 모드 선택 / Select Mode

- [ ] **SDM** — 아미노산 치환 (iPCR back-to-back mutagenesis)
- [ ] **RE Cloning** — 제한효소 클로닝 프라이머

둘 다 필요하면 두 섹션 모두 작성하세요.

---

## 공통 입력 / Common Input

### Template (주형) 서열

```
ATGAAAGCAATTTTCGTACTGAAAGGTTTTGTTGGTTTTCTTGCATTTATAATGTATCGT
TTATTTAATCTGTTTAAAATGGTATCAAATCGTAAAGGTATTGAAAGTTCTTTAGGTGGT
ACAGTTATGGCTTCAGCAATCGGTCGTGGTGTTGGTGCATTTGGTTTTTTAGCAGAACGT
```

> 코딩 서열(CDS) 전체를 붙여넣어 주세요. ATG 시작 권장.

---

## SDM 입력 (아미노산 치환)

| 항목 | 값 |
|------|----|
| 치환할 아미노산 위치 (1-based) | `42` |
| 현재 아미노산 | `Cys (C)` |
| 바꿀 아미노산 | `Ala (A)` |
| 선호 코돈 | (비워두면 E. coli 최적 코돈 자동 선택) |
| Tm 목표 | `60°C` |
| Overlap 길이 | `15 nt` (기본값) |

---

## RE Cloning 입력

| 항목 | 값 |
|------|----|
| 삽입 유전자명 | `GDH` |
| 목적 벡터 | `pET-28a` |
| Forward RE | `NcoI` (비워두면 자동 추천) |
| Reverse RE | `XhoI` (비워두면 자동 추천) |
| His-tag 포함 (N-terminal) | `Yes` |
| Stop codon 포함 | `Yes` |
| ATG 유지 여부 | `Yes` (NcoI 사용 시 ATG 포함) |

### 벡터 MCS 서열 (선택)

```
...CCATGGATCCGAATTCGAGCTCGGTACCCGGGGATCCTCTAGA...
   NcoI                          XhoI
```

> 비워두면 Claude가 표준 pET-28a MCS로 처리합니다.

---

## 출력 요청 / Requested Output

- [ ] 프라이머 서열 5'→3' (Forward / Reverse)
- [ ] 각 프라이머 Tm (°C)
- [ ] Amplicon 크기 (bp)
- [ ] Overlap 상보성 검증 (SDM)
- [ ] Reading frame 확인 (RE Cloning)
- [ ] QC 요약 (PASS / WARNING / FAIL)
- [ ] Macrogen 주문표 형식으로 출력

---

## 알고리즘 가이드 (Claude에게)

아래 규칙을 따라 Python으로 프라이머를 설계하고 결과를 출력하세요.

### Tm 계산
```
단순 공식 (기본): Tm = 64.9 + 41 * (GC_count - 16.4) / length
단, length < 14 nt: Tm = 2*(A+T) + 4*(G+C)
```

### SDM — iPCR back-to-back 설계 원칙
```
1. 주형에서 target_position 번째 아미노산의 코돈을 찾는다 (0-based index: (pos-1)*3)
2. Fw primer: 변이 코돈부터 시작 → 3' 방향으로 연장 (총 18-25 nt, Tm≥58°C)
3. Rv primer: 변이 코돈 직전 끝 → 5' 방향(주형의 반대방향)으로 연장 (reverse complement)
4. 두 프라이머 5' 말단이 서로 등을 맞대는 구조 (back-to-back)
5. 5' overhang: 각 프라이머 5' 끝에 15 nt 상보적 서열 추가 (blunt-end ligation용)
6. overlap region의 reverse complement 확인: 두 5' overhang이 정확히 상보적이어야 함
```

### RE Cloning — 프라이머 설계 원칙
```
Forward primer: 5'-[6nt padding]-[RE1 site]-[추가 nt for frame]-[20nt 유전자 시작]-3'
Reverse primer: 5'-[6nt padding]-[RE2 site]-[RC of 20nt 유전자 끝 (stop 전/후)]-3'

주요 RE site 서열:
  NcoI:  CCATGG  (ATG 내장)
  NdeI:  CATATG  (ATG 내장)
  XhoI:  CTCGAG
  BamHI: GGATCC
  EcoRI: GAATTC
  HindIII: AAGCTT
  NheI:  GCTAGC
  XbaI:  TCTAGA

Reading frame 검증: RE site 이후 첫 번째 코돈이 벡터 ATG와 frame이 맞는지 확인
```

### E. coli 최적 코돈표 (주요)
```python
ECOLI_PREFERRED_CODONS = {
    'A': 'GCG', 'R': 'CGT', 'N': 'AAC', 'D': 'GAT', 'C': 'TGC',
    'Q': 'CAG', 'E': 'GAA', 'G': 'GGT', 'H': 'CAT', 'I': 'ATT',
    'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTT', 'P': 'CCG',
    'S': 'AGC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTT',
    '*': 'TAA',
}
```

### QC 기준
```
PASS:    Tm 55-68°C, 3' 말단 GC clamp (G 또는 C로 끝), hairpin 없음
WARNING: Tm 50-54°C 또는 68-72°C, 또는 3' 말단 A/T
FAIL:    Tm < 50°C 또는 > 72°C, 또는 3' 말단 AAAA/TTTT
```

### Macrogen 주문표 출력 형식
```
| # | Name | Sequence (5'→3') | Scale | Purification |
|---|------|-----------------|-------|--------------|
| 1 | GDH_F | ATGCATGCA... | 0.2 µmol | Desalt |
| 2 | GDH_R | TTACATGCA... | 0.2 µmol | Desalt |
```
