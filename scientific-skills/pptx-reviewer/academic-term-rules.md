# Academic Term Rules — 생명공학/생화학 학술 표기 규칙

> pptx-reviewer, peer-review 등 학술 검토 스킬에서 공통 참조.
> 각 규칙은 [자동 탐지 가능] 또는 [육안 검토 필요]로 분류.

---

## 1. 생물 학명 표기 (Organism/Strain names)

**규칙**
- 속·종명(*Genus species*): **이탤릭**, 속명 대문자 시작, 종명 소문자
  - 첫 언급: *Escherichia coli*, 이후 축약: *E. coli*
  - PPT에서 이탤릭 미적용 → 댓글로 지적
- 균주 코드(strain code): 이탤릭 **불필요** — BL21(DE3), K-12, ATCC 13032
- 종명 단독 사용 금지: "cerinus" → "*Gluconobacter cerinus*"

**흔한 오류**
| 잘못된 표기 | 올바른 표기 |
|---|---|
| E.coli (공백 없음) | *E. coli* |
| Bacillus subtilis (비이탤릭) | *Bacillus subtilis* |
| G. cerinus (첫 언급) | *Gluconobacter cerinus* |
| gluconobacter cerinus | *Gluconobacter cerinus* |

---

## 2. 유전자 및 단백질/효소 표기

**유전자명**
- 원핵생물: **이탤릭, 소문자** — *xylA*, *gldA*, *ptsG*
- 돌연변이/allele: *xylA*⁺, *xylA*⁻ (superscript)
- Locus tag: 비이탤릭 — BSU12345

**단백질/효소명**
- 이탤릭 **불필요**, 대문자 시작 또는 약어 대문자
  - XylA (단백질), GDH (glucose dehydrogenase), FDH (formate dehydrogenase)
- 기원 species prefix (약어):

| Prefix | Organism |
|---|---|
| Bs | *Bacillus subtilis* |
| Ec | *Escherichia coli* |
| Ps / Pf | *Pseudomonas* sp. |
| Gc | *Gluconobacter cerinus* |
| Ao | *Aspergillus oryzae* |

- 예: BsGDH, PsFDH, EcXylA — 첫 등장 시 전체명 병기 필수
  - "BsGDH (glucose dehydrogenase from *B. subtilis*)"

**효소 이름 표기**
- 일반명: 소문자 — glucose dehydrogenase, formate dehydrogenase
- 약어: 대문자 — GDH, FDH, XylA, PNP
- EC number: EC 1.1.1.47 형식

---

## 3. 조효소·화학물질 표기

**조효소 (Cofactors) — 항상 대문자**
- NAD⁺, NADH, NADP⁺, NADPH (위첨자 + 필수)
- ATP, ADP, AMP
- FAD, FADH₂
- PQQ (pyrroloquinoline quinone)
- CoA, Acetyl-CoA

**흔한 오류**
| 잘못된 표기 | 올바른 표기 |
|---|---|
| NAD+, NADP+ (일반 +) | NAD⁺, NADP⁺ (superscript) |
| nadph | NADPH |
| AcP (informal) | acetyl phosphate (AcP) — 첫 등장 시 전체명 |

**화학식 아래첨자**
- CO₂, H₂O, H₂O₂, O₂ — 아래첨자 필수 (PPT에서 미적용 → 댓글)
- 이온: Mg²⁺, Ca²⁺, Fe³⁺, NH₄⁺, Pi (inorganic phosphate)

**기질/시약 약어** — 첫 등장 시 전체명 + 괄호 안 약어
- D-xylose (D-Xyl), triazole carboxamide (TCA), acetyl phosphate (AcP)
- 5-fluorouracil (5-FU), dihydroxyacetone (DHA)

---

## 4. 단위 표기 [자동 탐지 가능]

**부피**
| 잘못된 표기 | 올바른 표기 |
|---|---|
| ul, uL | µL |
| ml | mL |
| l | L |

**농도**
| 잘못된 표기 | 올바른 표기 |
|---|---|
| uM, UM | µM |
| nM (단독, 단위 앞 공백 없음) | nM (숫자와 공백: "50 nM") |

**온도**
| 잘못된 표기 | 올바른 표기 |
|---|---|
| ℃ (Unicode 기호) | °C (degree + C 조합) |
| 30°C (공백 없음) | 30 °C |
| C (단독) | °C |

**시간·속도**
| 잘못된 표기 | 올바른 표기 |
|---|---|
| hr, hrs | h |
| min. | min |
| RPM | rpm |

**숫자-단위 공백 규칙**: 숫자와 단위 사이 반드시 공백
- "50mM" → "50 mM", "30℃" → "30 °C", "6h" → "6 h"

---

## 5. 효소 Kinetics 표기

**파라미터 표기** (논문에서 이탤릭, PPT는 이탤릭 권장)
| 기호 | 올바른 표기 | 주의 |
|---|---|---|
| Michaelis constant | *K*m 또는 Km | KM, km, kM 모두 비표준 |
| turnover number | *k*cat 또는 kcat | Kcat 오류 |
| maximum velocity | *V*max 또는 Vmax | vmax, VMAX 오류 |
| catalytic efficiency | *k*cat/*K*m | kcat/Km |
| Hill coefficient | *n*H | nH |

**[S] 범위 체크리스트**
- Km 결정 실험: [S] 범위가 반드시 0.2×Km ~ 10×Km 포함 여부 확인
- [S] 범위 < Km → Km 과대 추정 위험, 댓글 필수

---

## 6. 통계 표기 [자동 탐지 가능]

**기술통계**
| 잘못된 표기 | 올바른 표기 |
|---|---|
| mean ± standard deviation | mean ± SD |
| n=3 (공백 없음) | n = 3 |
| n=3 (반복 종류 미기재) | n = 3 independent experiments |
| ± std | ± SD |

**추론통계**
- p-value: p < 0.05 (이탤릭 p 권장, 부등호 앞뒤 공백)
  - "p<0.05" → "p < 0.05"
- 검정 방법 명시: "Student's t-test", "one-way ANOVA + Tukey's HSD"
- Error bar 정체 명시 필수: SD, SEM, 95% CI 중 무엇인지

**Figure에서 체크**
- [ ] Error bar 있는가?
- [ ] n수 있는가?
- [ ] 통계 유의성 기호 (*p < 0.05, **p < 0.01, ***p < 0.001) 있는가?
- [ ] 통계 검정 방법 caption 또는 Methods에 있는가?

---

## 7. Figure Caption 규칙

**필수 구성 요소**
1. **번호**: "Figure 1." 또는 "Fig. 1." — 번호 없는 "Figure." 금지
2. **제목**: Bold, 첫 단어만 대문자 (sentence case)
3. **패널 설명**: (a), (b) 또는 (A), (B) — 대소문자 일관성
4. **실험 조건**: 기질 농도, 온도, pH, 시간, rpm 포함 (self-contained)
5. **정규화 기준**: "Relative amount" 사용 시 기준값 명시
   - ✗ "Relative amount of product"
   - ✓ "Relative yield normalized to maximum (pH 7.0 = 1.0)"
6. **n수 및 오차**: "Data are mean ± SD (n = 3 independent experiments)"
7. **약어**: 캡션 내 첫 등장 시 전체명 병기

**패턴 오류**
| 잘못된 표기 | 올바른 표기 |
|---|---|
| Figure. (번호 없음) | Figure 1. |
| pH 3-5 (hyphen) | pH 3–5 (en dash) |
| 25-45 °C | 25–45 °C |
| mean ± standard deviation | mean ± SD |

---

## 8. 문장 부호 — Dash 구분 [자동 탐지 가능]

| 기호 | 용도 | 예시 |
|---|---|---|
| `-` hyphen | 복합어 연결 | cell-free, sucrose-based, dose-dependent |
| `–` en dash | 범위 표시 | pH 3–5, 25–45 °C, 2020–2025 |
| `—` em dash | 삽입구, 강조 | The result—surprisingly—showed... |

**흔한 오류**: 범위에 hyphen(-) 사용 → en dash(–)로 수정 필요

---

## 9. 플라스미드·분자생물학 표기

**플라스미드명**: 이탤릭 권장
- pET-28a, pUC19, pBAD/His

**제한효소명**: 첫 글자 대문자, 나머지 소문자 (속명 이탤릭 불필요)
- EcoRI, BamHI, NdeI, XhoI — *Eco*RI 형식도 허용

**PCR 관련**
- primer (소문자), PCR (대문자)
- Tm (melting temperature) — 이탤릭 T, subscript m

---

## 10. 반응 조건 표기 표준 (Methods/Caption)

**표준 순서**: 기질 → 버퍼 → 온도 → 시간 → 교반속도
```
Reactions were performed with [substrate] ([conc.] mM) and [co-substrate] ([conc.] mM)
in [buffer name] ([conc.] mM, pH [X.X]) at [temp.] °C for [time] h with shaking at [rpm] rpm.
```

**버퍼 명**: 소문자 (고유명사 제외)
- sodium phosphate buffer, citrate buffer, Tris-HCl buffer
- "Tris": 대문자 (trade name)

---

## 11. pptx-reviewer 자동 탐지 대상 패턴

아래 패턴을 python-pptx 추출 텍스트에서 regex로 탐지:

```python
import re

TYPO_PATTERNS = [
    # 단위 오류
    (r'\bul\b',           'µL'),
    (r'\buL\b',           'µL'),
    (r'\buM\b',           'µM'),
    (r'\bml\b',           'mL'),
    (r'\bhr\b',           'h'),
    (r'\bhrs\b',          'h'),
    (r'\bRPM\b',          'rpm'),
    # 온도
    (r'℃',               '°C'),
    (r'(\d)°C',          r'\1 °C'),   # 공백 없음
    (r'(\d)mM',          r'\1 mM'),   # 공백 없음
    # 통계
    (r'n=(\d)',          r'n = \1'),  # 공백 없음
    (r'mean ± standard deviation', 'mean ± SD'),
    # dash
    (r'pH (\d)-(\d)',    r'pH \1–\2'),  # hyphen → en dash (범위)
    (r'(\d+)-(\d+) °C', r'\1–\2 °C'),
    # 생화학
    (r'\bNAD\+\b',       'NAD⁺'),
    (r'\bNADP\+\b',      'NADP⁺'),
    (r'supertanant',     'supernatant'),
    (r'seperati',        'separati'),
    (r'recombinant\b',   'recombinant'),  # 확인용
]

MANUAL_REVIEW = [
    "학명 이탤릭 여부 확인",
    "유전자명 이탤릭·소문자 확인 (xylA, gldA 등)",
    "단백질/효소 약어 첫 등장 시 전체명 병기 여부",
    "Figure 번호 존재 여부",
    "캡션 내 실험 조건 완비 여부 (온도, pH, 시간, 기질 농도)",
    "Relative amount 정규화 기준 명시 여부",
    "Error bar 정의 명시 여부 (SD/SEM/CI)",
    "n = 3 independent experiments 표기 여부",
    "Km 결정 실험 [S] 범위 적절성",
]
```
