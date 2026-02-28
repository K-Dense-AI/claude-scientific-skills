먼저 ~/.claude/commands/research-commons.md 파일을 읽고 공통 규칙(인용 스타일, Storage 키 체계, 스킬 간 연동, 트리거 분기)을 숙지한 후 아래 스킬을 실행하세요.

---
name: manuscript-writer
description: Academic manuscript writing, proofreading, and research support tool. Use when (1) analyzing manuscript structure and logic flow, (2) checking grammar, terminology, abbreviations, figure/table references, formatting errors, (3) searching literature and managing references, (4) writing cover letters, highlights, or reviewer responses, (5) checking equations and cross-section connectivity. Supports journal-specific styles (Nature, Angewandte, Elsevier, ACS). Integrates with experiment-hub, research-discussion, research-assistant via commons.
---

# Manuscript Writer

Academic manuscript research, writing, and proofreading skill.

> 공통 인용 규칙, Storage 키 체계, 연동 흐름은 **COMMONS.md** 참조.

## Workflow

```
Step 1: Initial Analysis & Structure
    ↓
Step 2: Section-by-Section Error Check
    ↓
Step 3: Cross-Section Connectivity Check
```

## Step 1: Initial Setup & Analysis

### 1.1 Ask Target Journal
Before analysis, ask user for target journal to load appropriate style guide.

### 1.2 Analyze Manuscript Structure
Read provided file/text and generate structure report:

```
═══════════════════════════════════════════════════
           MANUSCRIPT STRUCTURE ANALYSIS
═══════════════════════════════════════════════════

1. BASIC INFO
   ├─ Title: [title]
   ├─ Target Journal: [journal]
   └─ Manuscript Type: [Research Article / Review / Letter]

2. ABSTRACT SUMMARY
   ├─ Background: [1-2 sentences]
   ├─ Objective: [research aim]
   ├─ Methods: [key methods]
   ├─ Results: [main findings]
   └─ Conclusion: [conclusion]

3. LOGIC FLOW
   Problem/Gap → Hypothesis → Approach → Findings → Impact

4. SECTION BREAKDOWN
   ├─ Introduction: [gap, objective]
   ├─ Methods: [key techniques]
   ├─ Results: [findings with figures]
   └─ Discussion: [interpretation, comparison]

5. KEY CLAIMS & EVIDENCE
   [Claim] ↔ [Supporting Figure/Table]

6. FIGURES & TABLES OVERVIEW
   [List all with brief description]

7. NOVELTY STATEMENT
   [One sentence summary of contribution]
═══════════════════════════════════════════════════
```

## Step 2: Section-by-Section Error Check

Analyze each section separately. Report errors by category:

### Error Categories

**[GRAMMAR]**
- Subject-verb agreement
- Tense consistency (Methods: past, Results: past, Discussion: present/past)
- Articles (a/an/the)
- Prepositions
- Parallel structure
- **Avoid nominalizations** (use strong verbs)
  - ❌ Estimate → Estimation | ✅ "We estimated" not "Estimation was performed"
  - ❌ Decide → Decision | ✅ "We decided" not "Decision was made"
  - ❌ Confirm → Confirmation | ✅ "We confirmed" not "Confirmation was carried out"
  - ❌ Correlate → Correlation | ✅ "The data correlated" not "Correlation was observed"

**[TERMINOLOGY & ABBREVIATION]**
- First use: full name (abbreviation) format
  - Example: D-glyceraldehyde 3-phosphate (D-GAP)
- Define separately in Abstract and main text
- Consistent abbreviation (D-GAP, not GAP or d-GAP)
- Enzyme: regular for protein, italic for gene
- Species: E. coli after first full mention

**[REFERENCE]**
- Space before reference: "reported [1]." or "reported [1],"
- No duplicate references
- All claims need supporting references
- Check journal format: [1] or (1) or superscript

**[FIGURE/TABLE/SCHEME]**
- Sequential numbering (Figure 1 → 2 → 3)
- Mention order = number order
- All figures mentioned in text
- Lowercase subpanel: Figure 1a, Figure 1b
- SI format: Figure S1

**[FORMATTING]**
- Number + unit spacing: 10 mM, 50 °C, 100 rpm
- Percentage: 50% (no space)
- Range: en-dash 10–20 (not hyphen 10-20)
- Temperature: 25 °C
- pH: pH 7.0
- Time: 24 h, 30 min, 60 s
- Decimal: 0.1 (not .1)

**[EQUATION]**
- Sequential numbering: Eq. 1, Eq. 2
- All equations referenced in text
- All variables defined at first use
- Unit consistency

**[COMPLEXITY]**
- Sentence length (break if >40 words)
- Remove unnecessary modifiers
- One paragraph = one topic

### Section Error Report Format

```
═══════════════════════════════════════════════════
         SECTION ERROR REPORT: [SECTION NAME]
═══════════════════════════════════════════════════

[GRAMMAR]
├─ Line X: [error] → [correction]
└─ Line Y: [error] → [correction]

[TERMINOLOGY & ABBREVIATION]
├─ Line X: [error] → [correction]
└─ Line Y: [error] → [correction]

[REFERENCE]
...

[FIGURE/TABLE/SCHEME]
...

[FORMATTING]
...

[EQUATION]
...

[COMPLEXITY]
...

───────────────────────────────────────────────────
Summary: X errors found
  Grammar: X | Terminology: X | Reference: X
  Figure: X | Formatting: X | Equation: X | Complexity: X
═══════════════════════════════════════════════════
```

Analyze in order:
1. Abstract
2. Introduction
3. Methods
4. Results
5. Discussion
6. Conclusion
7. Supplementary Information

## Step 3: Cross-Section Connectivity Check

After all sections analyzed, check connectivity:

```
═══════════════════════════════════════════════════
           CROSS-SECTION CONNECTIVITY CHECK
═══════════════════════════════════════════════════

[LOGIC FLOW]
├─ Introduction gap → Methods approach connection
├─ Methods technique → Results data presence
└─ Results findings → Discussion interpretation

[TERMINOLOGY CONSISTENCY]
├─ Same abbreviations across all sections
└─ Enzyme/compound names consistent

[FIGURE/TABLE CROSS-CHECK]
├─ All figures referenced in appropriate sections
├─ SI figures referenced in main text
└─ No orphan figures/tables

[CLAIM-EVIDENCE ALIGNMENT]
├─ All claims have supporting data
└─ No unsupported conclusions

[ABSTRACT ↔ MAIN TEXT]
├─ Numbers match (yield, conversion, etc.)
└─ Conclusions align

[SI ↔ MAIN TEXT]
├─ Data consistency
└─ All SI items referenced

───────────────────────────────────────────────────
Connectivity Score: X/Y passed
Critical Issues: X
═══════════════════════════════════════════════════
```

## Additional Features

### Cover Letter
Journal-specific cover letter generation.

### Highlights (Elsevier)
- 3-5 bullet points
- Max 85 characters each

### Revision Support
- Point-by-point response format
- Track changes guidance

### Literature Search
- Keyword-based search (5-100 papers)
- Filter: recent + high-citation + preprint
- Gap analysis support

## Journal Style References

Supported journal styles:
- Nature
- Angewandte Chemie
- Elsevier journals
- ACS journals

---

## 스킬 간 연동

> 상세 연동 흐름도, 키 매핑, 트리거 분기 규칙은 **COMMONS.md** 참조.

### 입력 연동 (다른 스킬 → manuscript-writer)

| 소스 | 대상 | 용도 |
|------|------|------|
| `disc:conclusion` | Discussion 섹션 초안 | 구조화된 논의를 원고 Discussion으로 변환 |
| `research:references` | References 섹션 | 수집된 참고문헌을 저널 포맷으로 변환 |
| `exp:result` | Results 섹션 서술 | 실험 데이터/조건을 정확히 서술 |
| `exp:viz` | Figure/Table 목록 | 시각화 데이터를 Figure 범례에 연결 |
| `disc:chain` | 종합 분석 활용 | 체인 모드 결과(해석→팩트체크→비교→구조화)를 원고에 반영 |

### 출력 연동 (manuscript-writer → 다른 스킬)

| 소스 | 대상 | 트리거 |
|------|------|--------|
| `ms:revision` | `exp:protocol` | 리뷰어 추가실험 요청 감지 시 프로토콜 자동 생성 |
| `ms:errlog` | `disc:action` | 교차검증(S3)에서 논리 에러 발견 시 액션 아이템 생성 |
| `ms:refs` (부족) | `research:search` | References 부족 시 research-assistant 추가 검색 트리거 |

### 섹션별 자동 연동

- **Discussion 작성**: disc:conclusion(7단계 구조) + 팩트체크 결과 + 문헌 비교 통합 → 초안 생성. 데이터 없으면 "research-discussion M3 먼저 실행" 안내.
- **Results 작성**: exp:result(조건/결과) + exp:viz(그래프) → Figure 매핑 → 저널 스타일 서술.
- **References 변환**: APA 7th 기본 → 저널별 자동 변환 (Nature/ACS/Elsevier/Angewandte).

### 리비전 자동 감지

| 리뷰어 패턴 | 연동 액션 |
|------------|----------|
| "추가 실험 필요" | → experiment-hub M1 프로토콜 생성 |
| "통계 분석 부족" | → experiment-hub M6 재분석 |
| "문헌 비교 부족" | → research-assistant M1 추가 검색 |
| "논리 비약" | → research-discussion M4 액션 아이템 |

---

## 주의사항

- 원고 초안은 보조 도구, 최종 텍스트는 연구자 검토/수정
- 저널 스타일 변환 후 수동 검증 필요

사용자 요청: $ARGUMENTS
