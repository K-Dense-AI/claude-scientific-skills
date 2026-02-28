---
name: journal-presentation-maker
description: Automated academic presentation creation and review. Use when (1) creating journal club presentations, literature review slides, or research overview presentations, or (2) reviewing lab meeting/journal club PPTs for graph/table accuracy, scientific content validity, and speaker notes consistency. Searches academic databases (PubMed, Crossref, bioRxiv, arXiv), generates professional PowerPoint slides with Arial font and Korean speaker notes. Also reviews existing PPTs via python-pptx analysis and 호서대 OneDrive PowerPoint Online visual inspection.
---

# Journal Presentation Maker

## Overview

Create professional academic presentations from research literature. Search multiple databases (peer-reviewed journals and preprints), intelligently filter papers by relevance and impact, extract key content, and generate structured slides with proper citation tracking.

## Workflow

### Step 1: Paper Search and Selection

Use web_search and web_fetch tools to search academic literature:

1. **Search strategy**:
   - Use web_search with targeted keywords for academic databases
   - Query multiple sources: PubMed, Google Scholar, bioRxiv/medRxiv, arXiv
   - Example: `web_search("enzyme cascade rare sugar biosynthesis PubMed")`

2. **Smart filtering** (apply automatically):
   - **Recent papers (2023-2025)**: Include regardless of citations
   - **Medium recent (2020-2022)**: Include if well-cited (50+ citations)
   - **Classic papers (pre-2020)**: Include only high-impact (200+ citations)
   - **Preprints**: Always include regardless of date/citations
   
3. **Retrieve paper details**:
   - Use web_fetch to access DOI links and get full metadata
   - Extract: title, authors, year, journal, abstract, figures
   - Target: 5-100 papers (user specifies range)

4. **Present to user**:
   ```
   1. [Title] (2024, Nature, cited: 245)
      Authors et al.
      DOI: 10.1038/xxxxx
   
   2. [Title] (2025, bioRxiv - PREPRINT)
      Authors et al.
      DOI: 10.1101/xxxxx
   
   3. [Title] (2018, Cell, cited: 523) ★ High-impact classic
      Authors et al.
      DOI: 10.1016/xxxxx
   ```

5. **User selects papers** to include in presentation

### Step 2: Content Extraction and Organization

1. **Fetch full paper content**:
   - Use web_fetch on DOI links
   - Extract: abstract, methods, results, discussion, tables
   - **CRITICAL**: Extract ALL figure captions completely

2. **Figure extraction strategy** (4-tier approach):

   **우선순위**: Tier 0 (로컬 PDF) → Tier 1 (PMC) → Tier 2 (Chrome 스크린샷) → Tier 3 (Placeholder)

   **Tier 0: 로컬 PDF 직접 추출 (최우선 시도 — 유료 저널에 가장 효과적)**

   사용자가 로컬 PDF 파일을 제공한 경우 PyMuPDF(fitz)로 직접 이미지를 추출합니다.

   ```python
   import fitz  # PyMuPDF (pip install pymupdf)
   import os

   doc = fitz.open("paper.pdf")

   # 1. 각 페이지별 이미지 목록 확인
   for page_num in range(len(doc)):
       page = doc[page_num]
       images = page.get_images(full=True)
       for img in images:
           xref = img[0]
           info = doc.extract_image(xref)
           w, h = info["width"], info["height"]
           print(f"Page {page_num+1}: xref={xref}, {w}x{h}, {info['ext']}")
           # 작은 이미지(로고, 아이콘 등)는 건너뜀
           if w < 500 or h < 300:
               continue

   # 2. Figure 이미지 저장
   doc.extract_image(xref)["image"]  # bytes로 반환
   with open(f"fig{n}.{ext}", "wb") as f:
       f.write(info["image"])
   ```

   **Figure caption 추출 (텍스트 기반):**
   ```python
   for page_num in range(len(doc)):
       page = doc[page_num]
       text = page.get_text("text")
       # "FIGURE N" 또는 "Fig. N" 패턴으로 caption 위치 찾기
       idx = text.find("FIGURE 1")
       if idx >= 0:
           caption = text[idx:idx+600]  # 약 500자
   ```

   **python-pptx로 PPT에 삽입:**
   ```python
   from pptx.util import Inches
   from PIL import Image

   # 이미지 원본 크기 확인 (비율 유지)
   with Image.open(fig_path) as img:
       orig_w, orig_h = img.size

   # 슬라이드 가용 영역에 맞게 비율 유지 스케일링
   max_w = Inches(12.5)
   max_h = Inches(3.5)
   scale = min(max_w / orig_w, max_h / orig_h)
   fig_w = orig_w * scale
   fig_h = orig_h * scale

   # 슬라이드 가로 중앙 정렬
   fig_left = (slide_width - fig_w) / 2
   fig_top  = slide_height - fig_h - Inches(0.5)  # 하단에서 여백

   slide.shapes.add_picture(fig_path, fig_left, fig_top, fig_w, fig_h)
   ```

   **장점:**
   - 인터넷 접속 불필요
   - 유료 저널(Wiley, Elsevier, ACS 등) 논문에서도 완벽 작동
   - 원본 고해상도 이미지 그대로 추출 (2000px+ 해상도 일반적)
   - PDF에 embedded된 모든 figure를 정확히 추출
   - 캡션 텍스트도 함께 추출 가능

   **주의사항:**
   - 작은 이미지(너비 < 500px 또는 높이 < 300px)는 로고/아이콘이므로 건너뜀
   - 한 페이지에 여러 이미지가 있으면 크기로 필터링
   - PDF 페이지와 논문 Figure 번호를 수동으로 매핑 필요
     (예: Page 2 = Figure 1, Page 3 = Figure 2 등)
   - PyMuPDF 필요: `pip install pymupdf`
   - Pillow 필요 (크기 계산): `pip install pillow`

   **Tier 1: PMC Open Access (로컬 PDF 없는 경우 첫 번째 시도)**
   - DOI → PMCID 변환: `web_fetch("https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids={DOI}&format=json")`
   - PMC 논문 페이지 접근: `web_fetch("https://pmc.ncbi.nlm.nih.gov/articles/{PMCID}/")`
   - Figure 이미지 URL 추출: PMC에서 figure는 보통 `/pmc/articles/PMC.../figure/fig1/` 형태
   - 캡션 텍스트 함께 추출
   - Open Access 논문은 figure 직접 접근 가능

   **Tier 2: Chrome MCP 브라우저 캡처 (Tier 1 실패 시)**
   - 브라우저 도구 사용하여 논문 페이지 내비게이션:
     ```
     1. tabs_context_mcp → 탭 ID 확인
     2. navigate(url=논문_URL, tabId=탭ID)
     3. find(query="Figure 1" 또는 "Fig. 1", tabId=탭ID) → figure 요소 위치 확인
     4. scroll_to(ref=figure_ref, tabId=탭ID) → figure로 스크롤
     5. computer(action="screenshot", tabId=탭ID) → 전체 화면 캡처
     6. computer(action="zoom", region=[x0,y0,x1,y1], tabId=탭ID) → figure 영역만 확대 캡처
     ```
   - 캡처된 스크린샷을 슬라이드에 삽입할 이미지로 활용
   - 반드시 figure caption과 출처 표기 포함

   **Tier 3: Placeholder + 사용자 안내 (Tier 1, 2 모두 실패 시)**
   - 슬라이드에 figure placeholder 박스 생성:
     ```html
     <div class="figure-placeholder">
       <p>[Figure X를 여기에 삽입하세요]</p>
       <p class="citation">Source: Author et al. (Year), Figure X</p>
       <p class="instruction">다운로드: DOI_URL → Figure X 클릭 → 이미지 저장</p>
     </div>
     ```
   - Speaker notes에 figure 다운로드 방법 상세 안내

3. **Figure Caption PPT 삽입 모범 사례:**

   1. **Caption 위치**: Figure 바로 아래 (Inches(0.06) 여백), 캡션 박스 높이 Inches(0.75)
   2. **Caption 폰트**: **9.5pt**, 이탤릭, MID_GRAY(`#888888`) 색상
   3. **Caption 내용**: 논문 원본 캡션 전체 (축약 금지) + 출처 표기
   4. **출처 형식**: 캡션 마지막 줄 `(Author et al., Journal Year)` — 이 줄은 bold + DARK_GRAY
   5. **줄 나누기**: 캡션 텍스트를 `\n`으로 분리 후 line별 paragraph로 추가
   6. **정렬**: LEFT 정렬 (CENTER보다 가독성 우수)

   **예시 캡션 구현:**
   ```python
   CAPTION_H = Inches(0.75)
   cap_top = fig_top + fig_h + Inches(0.06)
   cap_box = slide.shapes.add_textbox(
       Inches(0.3), cap_top,
       slide_width - Inches(0.6), CAPTION_H)
   tf = cap_box.text_frame
   tf.word_wrap = True
   lines = caption_text.split("\n")
   first = True
   for line in lines:
       p = tf.paragraphs[0] if first else tf.add_paragraph()
       first = False
       p.alignment = PP_ALIGN.LEFT
       run = p.add_run()
       run.text = line
       run.font.size = Pt(9.5)
       run.font.italic = True
       run.font.color.rgb = RGBColor(0x88, 0x88, 0x88)
       if "et al." in line:          # 출처 줄은 bold + 진한 색
           run.font.bold = True
           run.font.color.rgb = RGBColor(0x33, 0x33, 0x33)
       run.font.name = "Arial"
   ```

4. **Figure attribution (필수)**:
   - 모든 figure에 출처 표기: "(Figure X from [ref_num])"
   - 원본 캡션 완전 포함
   - Open Access가 아닌 경우 "Adapted from" 또는 "Reprinted with permission" 표기
   - Speaker notes에 저작권 정보 기록

5. **Organize content by section**:
   - **Introduction**: Background, motivation, research gap with in-text citations [1,2]
   - **Methods**: Experimental design, techniques, materials
   - **Results**: Key findings with figures and original captions
   - **Discussion**: Interpretation, implications, limitations with citations
   - **Conclusion**: Main takeaways, future directions

6. **Reference numbering system**:
   - Assign reference numbers in order of first appearance: [1], [2], [3]...
   - Track first mention of each paper in content
   - Maintain consistent numbering throughout presentation

7. **Track sources meticulously**:
   ```json
   {
     "slide_content": "...",
     "ref_numbers": [1, 3, 5],
     "sources": [
       {
         "ref_num": 1,
         "paper_idx": 0,
         "author": "Kim et al.",
         "year": 2024,
         "doi": "10.1038/xxxxx",
         "figures": ["Figure 2B", "Figure 3A"],
         "figure_captions": {
           "Figure 2B": "Original caption text...",
           "Figure 3A": "Original caption text..."
         }
       }
     ]
   }
   ```

8. **Figure recommendation**:
   - Identify most impactful figures from each paper
   - Prioritize: key results, mechanisms, summary figures
   - Include complete original captions
   - Suggest 1-2 figures per major finding

### Step 2.5: Academic Notation Enforcement

콘텐츠를 슬라이드에 배치하기 전에, 아래 학술 표기 규칙을 반드시 적용:

#### 이탤릭체 규칙
- **유전자명**: 종에 따라 구분
  - 인간 유전자: 이탤릭 + 전체 대문자 → *BRCA1*, *TP53*
  - 마우스 유전자: 이탤릭 + 첫글자만 대문자 → *Brca1*, *Tp53*
  - 단백질: 이탤릭 아님 → BRCA1 (protein), Brca1 (mouse protein)
- **학명**: 항상 이탤릭 → *Escherichia coli*, 이후 축약 *E. coli*
- **라틴어**: *in vitro*, *in vivo*, *in situ*, *de novo*, *et al.* — 이탤릭 처리
- HTML 구현: `<em>` 또는 `<i>` 태그 사용

#### 화학식 표기
- **아래첨자** (원자 수): H₂O, CO₂, CH₃OH, C₆H₁₂O₆
- **위첨자** (이온 전하): Na⁺, Ca²⁺, Fe³⁺, SO₄²⁻, PO₄³⁻
- HTML 구현: `<sub>`, `<sup>` 태그
  - 예: `H<sub>2</sub>O`, `Fe<sup>3+</sup>`
- **반응식**: 기질 → 생성물 (→ 사용, 화살표 유니코드 U+2192)
- **효소 위치 표기**: EC 번호 포함 (예: EC 5.3.1.5, D-xylose isomerase)

#### 단위 표기
- SI 단위 준수, 숫자와 단위 사이 공백:
  - 농도: 5 μM, 10 mM, 0.1 M
  - 분자량: 45 kDa, 2.5 MDa
  - 온도: 37 °C (° 앞에 공백)
  - 부피: 100 μL, 5 mL
  - 시간: 2 h, 30 min, 10 s
  - pH: pH 7.4 (pH 이탤릭 아님, 공백 후 숫자)
- 잘못된 예: `5μM` → 올바른 예: `5 μM`

#### 통계/수치 표기
- p-value: *p* < 0.05 (이탤릭 소문자 p)
- 샘플 수: *n* = 10 (이탤릭 소문자 n)
- 표준편차: mean ± SD
- 수율/전환율: 95% yield, 85% conversion
- Michaelis-Menten 파라미터: K<sub>m</sub>, V<sub>max</sub>, *k*<sub>cat</sub>

#### 약어 규칙
- 첫 등장 시 반드시 풀어쓰기: "adenosine triphosphate (ATP)"
- 이후 약어만 사용: "ATP"
- 슬라이드별로 약어 목록 관리하여 누락 방지

### Step 3: Presentation Generation

#### python-pptx 폰트 크기 기준표 (13.33" × 7.5" 와이드스크린)

| 요소 | 폰트 크기 | 스타일 | 색상 |
|------|-----------|--------|------|
| 콘텐츠 슬라이드 헤더 | **Pt(26)** | bold | WHITE (DARK_BLUE 배경) |
| 피겨 슬라이드 헤더 | **Pt(26)** | bold | WHITE (DARK_BLUE 배경) |
| 섹션 소제목 (box 내) | **Pt(18)** | bold | ORANGE / MID_BLUE |
| 피겨 슬라이드 소제목 | **Pt(14)** | bold + italic | MID_BLUE |
| 본문 불릿 (주요 슬라이드) | **Pt(16)** | regular | DARK_GRAY |
| 키포인트 불릿 (피겨 슬라이드) | **Pt(14)** | regular | DARK_GRAY |
| 테이블 헤더 | **Pt(12)** | bold | WHITE (DARK_BLUE 배경) |
| 테이블 셀 | **Pt(12)** | regular | DARK_GRAY |
| 피겨 캡션 | **Pt(9.5)** | italic | MID_GRAY(`#888888`) |
| 푸터 텍스트 | **Pt(9–10)** | regular | LIGHT_GRAY |
| 슬라이드 조건/주석 | **Pt(12)** | regular | MID_BLUE |

> **헤더 바 높이**: Inches(1.1), margin_left Inches(0.4), margin_top Inches(0.12)
> **피겨 슬라이드 소제목 위치**: top Inches(1.15), height Inches(0.5)
> **키포인트 시작 위치**: top Inches(1.7), height Inches(1.3) (3개 불릿 기준)

#### 피겨 전용 슬라이드 레이아웃 템플릿

논문 Figure를 별도 슬라이드에 배치할 때 아래 레이아웃을 사용:

```python
def make_figure_slide(prs, title, subtitle, key_points, fig_path, caption_text):
    """
    피겨 전용 슬라이드 레이아웃 (13.33" × 7.5"):
      0.00" ─ DARK_BLUE 헤더 바 (1.1") ── 제목 Pt(26) WHITE bold
      1.15" ─ 소제목 (0.5") ──────────── Pt(14) MID_BLUE bold italic
      1.70" ─ 키포인트 (1.3") ─────────── ▸ 불릿 3개 Pt(14) DARK_GRAY
      3.10" ─ 피겨 이미지 (가변) ────────── 중앙 정렬, 비율 유지
      fig_top+fig_h+0.06" ─ 캡션 (0.75") ─ Pt(9.5) italic MID_GRAY
      7.12" ─ 푸터 (0.38") ──────────────── Pt(9) LIGHT_GRAY
    """
    DARK_BLUE = RGBColor(0x1A, 0x35, 0x5E)
    MID_BLUE  = RGBColor(0x2E, 0x5E, 0x9B)
    DARK_GRAY = RGBColor(0x33, 0x33, 0x33)
    MID_GRAY  = RGBColor(0x88, 0x88, 0x88)
    WHITE     = RGBColor(0xFF, 0xFF, 0xFF)
    W = prs.slide_width
    H = prs.slide_height

    sl = prs.slides.add_slide(prs.slide_layouts[6])  # blank
    sl.background.fill.solid(); sl.background.fill.fore_color.rgb = WHITE

    # 1. 헤더 바
    bar = sl.shapes.add_shape(1, Inches(0), Inches(0), W, Inches(1.1))
    bar.fill.solid(); bar.fill.fore_color.rgb = DARK_BLUE
    bar.line.fill.background()
    tf = bar.text_frame; tf.word_wrap = True
    tf.margin_left = Inches(0.4); tf.margin_top = Inches(0.12)
    p = tf.paragraphs[0]; p.alignment = PP_ALIGN.LEFT
    run = p.add_run(); run.text = title
    run.font.size = Pt(26); run.font.bold = True
    run.font.color.rgb = WHITE; run.font.name = "Arial"

    # 2. 소제목
    txb = sl.shapes.add_textbox(Inches(0.4), Inches(1.15), W - Inches(0.8), Inches(0.5))
    txb.text_frame.word_wrap = True
    p2 = txb.text_frame.paragraphs[0]; p2.alignment = PP_ALIGN.LEFT
    run2 = p2.add_run(); run2.text = subtitle
    run2.font.size = Pt(14); run2.font.bold = True; run2.font.italic = True
    run2.font.color.rgb = MID_BLUE; run2.font.name = "Arial"

    # 3. 키포인트
    kp_box = sl.shapes.add_textbox(Inches(0.4), Inches(1.7), W - Inches(0.8), Inches(1.3))
    kp_box.text_frame.word_wrap = True
    first = True
    for kp in key_points:
        p3 = kp_box.text_frame.paragraphs[0] if first else kp_box.text_frame.add_paragraph()
        first = False
        p3.space_before = Pt(4)
        run3 = p3.add_run(); run3.text = "  ▸  " + kp
        run3.font.size = Pt(14); run3.font.color.rgb = DARK_GRAY; run3.font.name = "Arial"

    # 4. 피겨 삽입 (비율 유지 스케일링)
    from PIL import Image as PILImage
    with PILImage.open(fig_path) as img:
        orig_w, orig_h = img.size
    CAPTION_H = Inches(0.75); FOOTER_H = Inches(0.38)
    avail_h = H - FOOTER_H - CAPTION_H - Inches(3.1) - Inches(0.1)
    avail_w = W - Inches(0.6)
    scale = min(avail_w / orig_w, avail_h / orig_h)
    fig_w = orig_w * scale; fig_h = orig_h * scale
    fig_left = (W - fig_w) / 2
    fig_top  = Inches(3.1) + (avail_h - fig_h) / 2
    sl.shapes.add_picture(fig_path, fig_left, fig_top, fig_w, fig_h)

    # 5. 캡션
    cap_box = sl.shapes.add_textbox(
        Inches(0.3), fig_top + fig_h + Inches(0.06),
        W - Inches(0.6), CAPTION_H)
    cap_box.text_frame.word_wrap = True
    first_cap = True
    for line in caption_text.split("\n"):
        p4 = cap_box.text_frame.paragraphs[0] if first_cap else cap_box.text_frame.add_paragraph()
        first_cap = False
        run4 = p4.add_run(); run4.text = line
        run4.font.size = Pt(9.5); run4.font.italic = True
        run4.font.color.rgb = MID_GRAY; run4.font.name = "Arial"
        if "et al." in line:
            run4.font.bold = True
            run4.font.color.rgb = DARK_GRAY

    return sl


def insert_slide_at(prs, new_slide, position):
    """prs의 position 인덱스 위치에 new_slide를 삽입"""
    slides = prs.slides._sldIdLst
    new_elem = slides[-1]
    slides.remove(new_elem)
    slides.insert(position, new_elem)
```

#### 후처리 스크립트 패턴 (Post-processing Workflow)

콘텐츠 생성과 피겨 삽입, 노트 작성을 분리하면 유지보수가 쉽다:

```
1. make_ppt.py       → 11개 콘텐츠 슬라이드 생성 (베이스 PPTX)
2. fix_figures.py    → 피겨 전용 슬라이드 삽입 (베이스 → v2, 15개 슬라이드)
3. update_notes.py   → 슬라이드 노트 전면 업데이트 (v2 인플레이스 수정)
```

각 스크립트는 독립적으로 실행 가능하며, 노트만 재작성하거나 피겨만 교체하는 작업이 용이함.

> **Windows 인코딩 주의**: 콘솔 출력에 한글/특수문자 사용 시 반드시 첫 줄에 추가:
> ```python
> import sys; sys.stdout.reconfigure(encoding='utf-8')
> ```

1. **Plan presentation structure** (maximum 30 slides):
   - Title slide
   - Introduction (2-4 slides)
   - Methods overview (1-3 slides)  
   - Results (8-12 slides, figure-heavy with original captions)
   - Discussion (2-4 slides)
   - Conclusion (1 slide)
   - References (1-3 slides)

2. **Create CSS styling** (see [references/styling_guide.md](references/styling_guide.md)):
   ```css
   /* shared-styles.css */
   :root {
     --primary-dark: #1a365d;
     --primary-medium: #2c5282;
     --text-base: 20px;
     --text-xl: 32px;
     --font-family: 'Arial', sans-serif;  /* REQUIRED: Use Arial */
   }
   
   body, h1, h2, h3, h4, h5, h6, p, li {
     font-family: var(--font-family);
   }
   ```

3. **Generate HTML for each slide**:
   - Use Arial font for all text
   - Include in-slide citation: [1,2,3]
   - Place detailed citation in speaker notes
   
   **For Introduction/Discussion slides**:
   - Include claims with reference numbers: "Recent studies show... [1,2]"
   - Track which papers support each claim

   **For Results slides**:
   - Include figures with complete original captions
   - Add reference at end of caption: "(Figure 2 from [3])"

   **학술 표기 적용 예시 (Step 2.5 규칙 반영)**:
   ```html
   <!-- 화학식 + 학명 -->
   <p>D-glucose를 <em>E. coli</em> 유래 xylose isomerase (EC 5.3.1.5)로
   D-fructose로 전환 (K<sub>m</sub> = 15 mM, <em>k</em><sub>cat</sub> = 120 s<sup>−1</sup>)</p>

   <!-- 유전자/단백질 구분 -->
   <p><em>xylA</em> 유전자가 코딩하는 XylA 단백질의 활성을 측정</p>

   <!-- 통계 표기 -->
   <p>수율 85 ± 3% (<em>n</em> = 3, <em>p</em> &lt; 0.01)</p>

   <!-- 단위 (숫자와 단위 사이 공백) -->
   <p>반응 조건: 50 mM Tris-HCl (pH 7.5), 37 °C, 2 h</p>
   ```

4. **Create schematic diagrams for each section**:
   ```html
   <!-- Use SVG or simple HTML/CSS diagrams -->
   <div class="schematic">
     <!-- Summarize key concepts from multiple papers -->
     <!-- Show workflow, mechanism, or relationships -->
     <!-- Add caption explaining the synthesis -->
   </div>
   <p class="caption">
     Schematic summary of [topic] integrating findings from [1,2,5].
   </p>
   ```

5. **MANDATORY: Add comprehensive speaker notes to EVERY slide**:

   Speaker notes are REQUIRED for all slides. Each slide must include:

   **For Title Slide**:
   ```
   [인사 및 소개]
   - 논문 제목, 저자, 저널 정보
   - 논문의 중요성 및 선정 이유
   - 발표 순서 개요
   ```

   **For Content Slides**:
   ```
   [슬라이드 주제]
   - 핵심 포인트 상세 설명
   - 배경 지식 및 맥락
   - 발표 시 강조할 내용

   [출처]
   - [1] Kim et al. (2024), DOI: 10.1038/xxxxx
     → Figure 2B (caption: "...")
   - [3] Lee et al. (2023), DOI: 10.1016/xxxxx
     → Figure 1A (caption: "...")
   ```

   **For Discussion/Conclusion Slides**:
   ```
   [토론 포인트]
   - 각 질문에 대한 답변 가이드
   - 예상 질문 및 대응
   - 개인적 의견/비평 (강점, 약점, 한계)
   ```

   **Speaker Notes Content Guidelines**:
   - 한국어로 작성 (발표자 편의)
   - 슬라이드당 **400–800자** 분량 (너무 짧으면 발표 중 막힘)
   - 첫 줄에 `【슬라이드 제목 — 섹션 설명】` 형식의 헤더 추가
   - 발표 시 말할 내용을 대본처럼 작성
   - 전문 용어 설명 포함
   - 청중 예상 질문과 답변 준비

   **권장 노트 구조 (피겨 슬라이드):**
   ```
   【Figure N — 그림 제목】

   (a) 패널 설명:
     • 핵심 내용 1
     • 핵심 내용 2

   (b) 패널 설명:
     • 결과 수치 (예: 79±5% yield, 95% ee)
     • 비교 포인트

   질문 포인트: '예상 질문?' → 답변 가이드
   ```

   **권장 노트 구조 (콘텐츠 슬라이드):**
   ```
   【슬라이드 N — 섹션 이름】

   [핵심 내용]
     • 첫 번째 포인트: 상세 설명
     • 두 번째 포인트: 배경 지식

   [실험/방법 조건]
     • 조건 1, 조건 2

   Note: 강조 포인트 또는 청중이 놓치기 쉬운 내용
   ```

   **update_notes.py 패턴** — 노트를 별도 스크립트로 일괄 업데이트:
   ```python
   import sys; sys.stdout.reconfigure(encoding='utf-8')
   from pptx import Presentation

   NOTES = {
       1: "【타이틀 슬라이드 — 발표 오프닝】\n\n...",
       2: "【Figure 1 — 연구 배경】\n\n(a) ...\n(b) ...",
       # ...슬라이드 번호(1-based): 노트 텍스트
   }

   prs = Presentation("presentation.pptx")
   for slide_idx, note_text in NOTES.items():
       notes_tf = prs.slides[slide_idx - 1].notes_slide.notes_text_frame
       notes_tf.clear()
       notes_tf.text = note_text
   prs.save("presentation.pptx")
   ```

6. **Convert HTML to PowerPoint** using html2pptx:
   ```javascript
   const pptxgen = require("pptxgenjs");
   const { html2pptx } = require("@ant/html2pptx");
   
   const pptx = new pptxgen();
   pptx.layout = "LAYOUT_16x9";
   
   // Process each HTML slide (max 30 slides)
   await html2pptx(pptx, "slide_01.html");
   await html2pptx(pptx, "slide_02.html");
   
   // Add notes to slides with reference details
   pptx.slides[0].addNotes("Speaker notes with [ref] sources...");
   
   await pptx.writeFile("presentation.pptx");
   ```

7. **Create References slide**:
   - List papers in order of first citation: [1], [2], [3]...
   - Use APA or Nature format
   - Include DOI links
   - Example:
     ```
     1. Kim J, et al. (2024). Title. Nature 123:456-789. DOI: 10.1038/xxxxx
     2. Lee S, et al. (2023). Title. Cell 456:123-456. DOI: 10.1016/xxxxx
     ```

### Step 4: Quality Checks

1. **Citation completeness**:
   - Every factual claim has reference number [1,2,3]
   - References numbered in order of first appearance
   - All figures attributed with original captions
   - References slide includes all cited papers in order

2. **Content quality**:
   - Maximum 30 slides total
   - Max 6-8 bullet points per slide
   - Max 40-50 words per slide (excluding title and captions)
   - Figures with complete original captions included
   - Section summary schematics present

3. **Technical validation**:
   - All fonts are Arial
   - All DOI links functional
   - Author names and years accurate
   - Reference numbers consistent throughout
   - Figure captions complete and attributed

4. **Speaker Notes validation (MANDATORY)**:
   - EVERY slide must have speaker notes
   - Notes written in Korean for presenter convenience
   - Each note includes:
     - Detailed explanation of slide content
     - Background context and terminology
     - Source citations with DOI
     - Talking points for presentation
   - Title slide: introduction and overview
   - Content slides: explanations and sources
   - Discussion slides: Q&A preparation
   - Final slide: summary and closing remarks

5. **Academic notation validation**:
   - 유전자명/학명에 이탤릭 적용 확인 (`<em>` 또는 `<i>` 태그)
   - 화학식의 subscript/superscript 정확성 (`<sub>`, `<sup>`)
   - 단위와 숫자 사이 공백 확인 (5 μM, not 5μM)
   - 약어 첫 등장 시 풀어쓰기 확인
   - 효소명에 EC 번호 포함 확인
   - K<sub>m</sub>, V<sub>max</sub> 등 kinetic 파라미터 표기 정확성
   - *p*-value, *n* 등 통계 기호 이탤릭 확인
   - Figure 출처 표기 완전성 (모든 외부 figure에 attribution 존재)

### Step 5: Presentation Review (PPT 검토)

생성된 PPT 파일을 검토하여 학술 표기, 레이아웃, 인용 정확성을 점검한다.
사용자가 기존 PPT를 검토 요청하는 경우에도 이 단계를 독립적으로 실행할 수 있다.

#### 5-1. PPT 열기 (2가지 방법)

**방법 A: python-pptx 프로그래밍 분석 (로컬 파일 — 권장)**

로컬 PPT 파일을 python-pptx로 직접 분석. 브라우저 불필요, 텍스트 기반 정밀 검사에 최적.

```python
from pptx import Presentation
from pptx.util import Pt

prs = Presentation("presentation.pptx")
for i, slide in enumerate(prs.slides):
    for shape in slide.shapes:
        if shape.has_text_frame:
            for para in shape.text_frame.paragraphs:
                for run in para.runs:
                    text = run.text
                    font = run.font
                    # 학술 표기 검증:
                    # - 이탤릭 여부: font.italic
                    # - 아래첨자: font.subscript (화학식 H2O 등)
                    # - 위첨자: font.superscript (이온 전하 등)
                    # - 폰트: font.name (Arial 확인)
                    # - 크기: font.size
```

검사 항목:
- 유전자명/학명에 `font.italic = True` 적용 여부
- 화학식에 `font.subscript`/`font.superscript` 적용 여부
- 숫자와 단위 사이 공백 (정규식 `r'\d[μmMkK]'` 매칭)
- 폰트가 Arial인지 확인
- 슬라이드당 텍스트 양 (40-50 단어 초과 여부)

**방법 B: 호서대 OneDrive + PowerPoint Online (시각적 검토)**

OneDrive에 동기화된 PPT 파일을 PowerPoint Online에서 열어 시각적으로 검토.
호서대 OneDrive는 브라우저에 이미 로그인되어 있어 별도 인증 불필요.

```
1. tabs_context_mcp → 탭 ID 확인 (없으면 createIfEmpty=true)
2. navigate(url="https://visionhoseo-my.sharepoint.com/personal/20141204_365_hoseo_edu/Documents/", tabId=탭ID)
3. OneDrive 폴더 탐색 → PPT 파일 클릭 → PowerPoint Online 자동 열림
4. 슬라이드 패널에서 클릭으로 슬라이드 이동
5. computer(action="screenshot", tabId=탭ID) → 슬라이드 캡처
6. computer(action="zoom", region=[x0,y0,x1,y1], tabId=탭ID) → 세부 영역 확대
```

- OneDrive 로컬 경로: `C:\Users\Jahyun\OneDrive - 호서대학교\`
- 이 폴더에 PPT를 넣으면 자동으로 OneDrive 웹에 동기화됨
- PPT 클릭 시 PowerPoint Online이 자동으로 열림 (Google Slides 업로드 불필요)

**방법 선택 기준**:
- 텍스트/표기 오류 정밀 분석 → **방법 A** (python-pptx)
- 레이아웃/Figure 렌더링 시각 확인 → **방법 B** (OneDrive + PowerPoint Online)
- 최적: **A + B 병행** — python-pptx로 표기 오류 검출 후, PowerPoint Online으로 시각 확인

#### 5-2. 슬라이드별 검토 체크리스트

각 슬라이드 스크린샷을 분석하여 다음 항목을 점검:

**A. 레이아웃 및 가독성**:
- [ ] 텍스트가 슬라이드 영역을 벗어나지 않는지 (overflow 확인)
- [ ] 폰트 크기가 충분히 큰지 (헤더 Pt(26), 본문 불릿 Pt(16), 테이블 Pt(12), 캡션 Pt(9.5))
- [ ] 여백이 적절한지 (텍스트가 가장자리에 붙어있지 않은지)
- [ ] 그림과 텍스트 간 겹침이 없는지
- [ ] 전체적인 시각적 균형과 정렬

**B. 학술 표기 정확성** (Step 2.5 기준):
- [ ] 이탤릭체: 유전자명(*xylA*), 학명(*E. coli*), 라틴어(*in vitro*)
- [ ] 화학식: 아래첨자/위첨자 (H₂O가 H2O로 표시되지 않는지)
- [ ] 단위: 숫자와 단위 사이 공백 (5 μM, not 5μM)
- [ ] 통계 기호: *p*, *n* 이탤릭 여부
- [ ] 효소 kinetic 파라미터: K_m이 K<sub>m</sub>으로 정확히 렌더링되는지

**C. Figure 품질**:
- [ ] Figure가 실제로 삽입되어 있는지 (placeholder만 남아있지 않은지)
- [ ] Figure 해상도가 충분한지 (깨지거나 흐릿하지 않은지)
- [ ] Figure caption이 완전히 표시되는지 (잘림 없는지)
- [ ] 출처 표기: "(Figure X from [ref_num])" 존재 여부

**D. 인용 일관성**:
- [ ] 본문 내 인용 번호 [1,2,3]가 표시되는지
- [ ] References 슬라이드의 번호와 본문 인용이 일치하는지
- [ ] DOI가 포함되어 있는지

**E. 전체 흐름**:
- [ ] 슬라이드 순서가 논리적인지 (Introduction → Methods → Results → Discussion → Conclusion)
- [ ] 각 섹션 전환이 자연스러운지
- [ ] 전체 슬라이드 수가 적절한지 (15-30장)

#### 5-3. 검토 결과 리포트

검토 완료 후 사용자에게 다음 형식으로 결과 보고:

```markdown
## PPT 검토 결과

### 전체 요약
- 총 슬라이드: X장
- 검토 통과: Y항목 / 전체 Z항목
- 수정 필요: N건

### 슬라이드별 이슈

**슬라이드 3 (Methods)**:
- ❌ "H2O" → "H₂O" (아래첨자 누락)
- ❌ "E. coli" → "*E. coli*" (이탤릭 누락)
- ⚠️ Figure 1 해상도 낮음 — 원본 재삽입 권장

**슬라이드 7 (Results)**:
- ❌ "5μM" → "5 μM" (단위 공백 누락)
- ❌ Figure 출처 표기 없음 → "(Figure 2A from [3])" 추가 필요

**슬라이드 12 (Discussion)**:
- ✅ 모든 항목 통과

### 수정 우선순위
1. 🔴 높음: 학술 표기 오류 (X건) — 발표 신뢰도에 직접 영향
2. 🟡 중간: Figure 품질/출처 (X건) — 저작권 및 가독성
3. 🟢 낮음: 레이아웃 미세 조정 (X건) — 선택적 개선
```

#### 5-4. 자동 수정 제안

검토에서 발견된 이슈에 대해 수정 방안 제시:

- **학술 표기 오류**: 수정된 텍스트를 구체적으로 제시 (before → after)
- **Figure 이슈**: 대체 이미지 소스 또는 재캡처 방법 안내
- **레이아웃 문제**: CSS/HTML 수정 코드 제시 (PPT 재생성 시 반영)
- **인용 불일치**: 누락된 인용 번호 및 위치 명시

사용자 승인 시, 이슈를 반영하여 PPT를 재생성하거나 부분 수정 가능.

### Step 6: Lab Meeting Presentation Review (랩미팅 발표 검토)

랩미팅 PPT를 종합적으로 검토한다. 그래프/테이블의 정확성, 슬라이드 내용의 과학적 타당성,
발표 노트의 적합성을 점검하고, 제안사항 및 오류를 보고한다.
Journal Club PPT 검토(Step 5)와 독립적으로 실행 가능.

#### 6-1. PPT 파일 로드

Step 5-1과 동일한 방법 사용:
- **방법 A (python-pptx)**: 로컬 파일 프로그래밍 분석 — 텍스트/수치/노트 정밀 검사
- **방법 B (호서대 OneDrive)**: PowerPoint Online 시각적 검토 — 그래프/레이아웃 확인
- **최적**: A + B 병행

랩미팅 PPT 일반 경로:
```
C:\Users\Jahyun\OneDrive - 호서대학교\바탕 화면 [Labtop]\Team meeting\
├── research presentation\    ← 연구 발표 PPT
├── Journal Club\             ← Journal Club PPT (Step 5 대상)
└── 20XX team meeting.pptx    ← 연간 팀미팅 PPT
```

#### 6-2. 그래프 및 Figure 검토

**A. 그래프 품질 점검** (PowerPoint Online 시각적 확인):
```
각 그래프/차트에 대해:
1. computer(action="screenshot", tabId=탭ID) → 슬라이드 캡처
2. computer(action="zoom", region=[x0,y0,x1,y1], tabId=탭ID) → 그래프 영역 확대
3. 아래 체크리스트 항목 점검
```

- [ ] **축 레이블**: X축, Y축 레이블이 존재하고, 단위가 표기되어 있는지
  - 올바른 예: "Conversion (%)", "Time (h)", "Concentration (mM)"
  - 잘못된 예: 축 레이블 없음, 단위 누락
- [ ] **축 범위**: 데이터에 적합한 범위인지 (불필요한 여백, 잘린 데이터 없는지)
- [ ] **범례 (Legend)**: 여러 데이터셋이 있을 때 범례가 존재하고 구분 가능한지
- [ ] **데이터 포인트**: 에러바(error bar)가 있는지, 에러바 유형 명시 (SD, SEM, 95% CI)
- [ ] **그래프 유형 적합성**: 데이터 특성에 맞는 그래프 유형인지
  - 시간 경과 → line chart, 비교 → bar chart, 분포 → box plot, 상관관계 → scatter plot
- [ ] **폰트 크기**: 그래프 내 텍스트가 발표 시 가독 가능한 크기인지 (≥ 14pt 권장)
- [ ] **색상 구분**: 흑백 출력 시에도 구분 가능한지 (색약 대응)
- [ ] **Figure 번호/캡션**: 논문 figure 인용 시 출처가 명확한지

**B. 테이블 점검**:
- [ ] **헤더**: 열/행 헤더가 명확하고 단위가 표기되어 있는지
- [ ] **정렬**: 숫자는 우측 정렬, 텍스트는 좌측 정렬
- [ ] **유효 숫자**: 일관된 소수점 자릿수 사용 (예: 모두 소수점 2자리)
- [ ] **통계 표기**: p-value, n 수 등이 테이블 내 또는 각주에 표기되어 있는지
- [ ] **강조**: 핵심 결과가 굵은 글씨나 색상으로 강조되어 있는지
- [ ] **가독성**: 너무 많은 행/열로 인해 읽기 어렵지 않은지 (≤ 8열 권장)

**C. 그래프 수치 교차 확인** (python-pptx):
```python
# 슬라이드 내 텍스트에서 수치 추출
import re
for slide in prs.slides:
    for shape in slide.shapes:
        if shape.has_text_frame:
            text = shape.text_frame.text
            # 수율, 전환율 등 핵심 수치 추출
            numbers = re.findall(r'(\d+\.?\d*)\s*(%|mM|μM|kDa|°C|h|min)', text)
            # 발표 노트의 수치와 비교
```
- 슬라이드 본문의 수치와 발표 노트의 수치가 일치하는지 확인
- 그래프 캡션에 언급된 수치와 본문 설명이 모순되지 않는지 확인

#### 6-3. 내용 과학적 타당성 점검

**A. 논리적 흐름**:
- [ ] 연구 목적/가설이 명확하게 제시되어 있는지
- [ ] 실험 설계가 가설을 검증하기에 적합한지
- [ ] 결과 해석이 데이터에 기반하는지 (과대 해석 없는지)
- [ ] 대조군(control) 실험이 포함되어 있는지
- [ ] 결론이 결과에서 논리적으로 도출되는지

**B. 실험 조건 명시**:
- [ ] 반응 조건 (온도, pH, 시간, 농도) 표기 여부
- [ ] 사용 균주/세포주 명시 여부
- [ ] 분석 방법 (HPLC, GC-MS, SDS-PAGE 등) 명시 여부
- [ ] 반복 실험 횟수 (*n*) 명시 여부

**C. 기존 문헌과의 비교**:
- [ ] 자신의 결과를 기존 연구와 비교했는지
- [ ] 차이점/유사점에 대한 논의가 있는지
- [ ] 참고 문헌 인용이 적절한지

**D. 흔한 과학적 오류 체크**:
- [ ] 상관관계를 인과관계로 해석하지 않았는지
- [ ] 통계적 유의성 (*p* < 0.05)과 생물학적 의의를 구분했는지
- [ ] 단일 실험 결과를 일반화하지 않았는지
- [ ] 음성 결과(negative results)를 적절히 다루었는지

#### 6-4. 발표 노트 (Speaker Notes) 점검

**A. 노트 존재 여부** (python-pptx):
```python
for i, slide in enumerate(prs.slides):
    notes_slide = slide.notes_slide if slide.has_notes_slide else None
    if notes_slide:
        notes_text = notes_slide.notes_text_frame.text
        # 노트 내용 분석
    else:
        print(f"⚠️ 슬라이드 {i+1}: 발표 노트 없음")
```

**B. 노트 품질 체크리스트**:
- [ ] **모든 슬라이드에 노트 존재**: 빈 노트 슬라이드 검출
- [ ] **노트-슬라이드 일치**: 노트 내용이 해당 슬라이드 내용과 관련 있는지
  - 슬라이드 제목/키워드가 노트에도 등장하는지 교차 확인
  - 노트가 다른 슬라이드의 내용을 잘못 참조하지 않는지
- [ ] **수치 일관성**: 노트에서 언급한 수치가 슬라이드 그래프/텍스트와 일치하는지
  ```python
  # 슬라이드 본문 수치 vs 노트 수치 비교
  slide_numbers = extract_numbers(slide_text)
  notes_numbers = extract_numbers(notes_text)
  mismatches = find_mismatches(slide_numbers, notes_numbers)
  ```
- [ ] **충분한 분량**: 각 노트가 최소 100자 이상인지 (너무 짧으면 발표 시 어려움)
- [ ] **발표 어조**: 구어체로 발표하기 적합한지 (논문 문장 그대로 복붙 아닌지)
- [ ] **예상 질문**: Discussion/Results 슬라이드 노트에 예상 질문과 답변이 준비되어 있는지

**C. 노트 내용 정확성 검증**:
- 노트에서 주장하는 내용이 슬라이드의 데이터와 모순되지 않는지
- 노트에서 언급한 문헌 정보 (저자, 연도)가 References와 일치하는지
- 노트에서 설명하는 메커니즘/방법이 실제 실험과 일치하는지

#### 6-5. 랩미팅 검토 결과 리포트

```markdown
## 랩미팅 PPT 검토 결과

### 전체 요약
- 총 슬라이드: X장
- 그래프/Figure: Y개 (이슈 N건)
- 테이블: Z개 (이슈 N건)
- 발표 노트: 존재 X장 / 누락 Y장
- 수정 필요: 총 N건

### 1. 그래프/테이블 이슈

**슬라이드 5 (Results - Conversion Graph)**:
- ❌ Y축 레이블 누락 → "Conversion (%)" 추가 필요
- ❌ 에러바 없음 → SD 또는 SEM 에러바 추가 권장
- ⚠️ X축 범위 0-100이나 데이터는 0-50 구간에 집중 → 범위 조정 고려

**슬라이드 8 (Results - Kinetics Table)**:
- ❌ K_m 값 소수점 불일치: "15.2 mM" vs "15.23 mM" (같은 테이블 내)
- ⚠️ 유효숫자 통일 필요

### 2. 내용 과학적 이슈

**슬라이드 10 (Discussion)**:
- ⚠️ "전환율이 크게 향상되었다" → 구체적 수치 없이 정성적 표현만 사용
- ❌ 대조군 결과 미제시 → 비교 근거 불충분
- 💡 제안: "WT 대비 2.3배 향상 (85% vs 37%)" 형태로 정량적 비교 추가

**슬라이드 12 (Conclusion)**:
- ⚠️ 결론이 Results에서 보여준 데이터 범위를 넘어서는 일반화 포함
- 💡 제안: "본 조건에서" 등 범위 한정 표현 추가

### 3. 발표 노트 이슈

**슬라이드 3 (Methods)**:
- ❌ 발표 노트 없음 → 실험 조건 설명 노트 추가 필요

**슬라이드 7 (Results)**:
- ❌ 노트의 수치 "90% 전환율"이 슬라이드 그래프의 "85%"와 불일치
- ⚠️ 노트가 논문 문장 그대로 복붙 → 구어체로 재작성 권장

**슬라이드 11 (Discussion)**:
- ✅ 예상 질문 3개와 답변 준비됨 — 양호
- 💡 제안: "왜 이 효소를 선택했는지"에 대한 질문도 추가 권장

### 4. 종합 제안사항

#### 🔴 반드시 수정 (발표 전 필수)
1. 그래프 축 레이블 누락 (슬라이드 5, 9)
2. 슬라이드-노트 수치 불일치 (슬라이드 7)
3. 대조군 데이터 미제시 (슬라이드 10)

#### 🟡 수정 권장 (품질 향상)
4. 에러바 추가 (슬라이드 5, 6, 9)
5. 정량적 비교 표현 보강 (슬라이드 10, 12)
6. 발표 노트 없는 슬라이드 보완 (슬라이드 3, 4)

#### 🟢 선택적 개선 (시간 여유 시)
7. 테이블 유효숫자 통일 (슬라이드 8)
8. 그래프 축 범위 최적화 (슬라이드 5)
9. 예상 질문 추가 (슬라이드 11)

### 5. 발표 준비도 평가
- 슬라이드 완성도: ★★★☆☆ (3/5)
- 데이터 신뢰도: ★★★★☆ (4/5)
- 발표 노트 준비: ★★☆☆☆ (2/5)
- 예상 질문 대비: ★★★☆☆ (3/5)
- **종합**: 수정사항 반영 후 발표 가능
```

#### 6-6. 자동 수정 지원

검토 결과에 기반하여 수정 제안:

- **그래프 이슈**: 축 레이블, 범례 추가 방법 안내 (PowerPoint Online 편집 모드 활용)
- **수치 불일치**: 정확한 수치를 양쪽 (슬라이드 + 노트)에 통일하여 제시
- **발표 노트 보완**: 누락된 슬라이드에 대해 노트 초안 자동 생성
  ```python
  # python-pptx로 노트 추가/수정
  from pptx import Presentation
  prs = Presentation("labmeeting.pptx")
  slide = prs.slides[2]  # 노트 없는 슬라이드
  notes_slide = slide.notes_slide
  notes_slide.notes_text_frame.text = "생성된 발표 노트..."
  prs.save("labmeeting_reviewed.pptx")
  ```
- **과학적 표현 개선**: before → after 형태로 구체적 수정안 제시
- **예상 질문 생성**: 연구 내용 기반으로 교수/동료가 물어볼 만한 질문 3-5개 자동 생성

## Example Usage

**User**: "Make a presentation about enzyme cascades for rare sugar biosynthesis"

**Claude**:
1. Searches PubMed, bioRxiv for "enzyme cascade rare sugar biosynthesis"
2. Finds 45 papers, filters to 15 most relevant
3. Presents list to user with citations and impact metrics
4. User selects 8 papers
5. Extracts content from all 8 papers including:
   - All figure captions
   - Key claims with context
   - Methods and results details
6. Assigns reference numbers in order of first mention: [1]-[8]
7. Generates 25-slide presentation with:
   - Title slide
   - Background (3 slides) - with claims cited as [1,2]
   - Enzyme systems (5 slides) - figures with original captions
   - Schematic summary (1 slide) - integrating findings from [1,3,5]
   - Cascade optimization (4 slides) - figures and data
   - Kinetic analysis (3 slides) - figures with captions from [6,7]
   - Applications (3 slides) - with citations [2,4,8]
   - Conclusion (1 slide)
   - References (2 slides) - papers listed [1]-[8] in order
8. All text in Arial font
9. EVERY slide has comprehensive Korean speaker notes:

   **Example - Title Slide Notes**:
   ```
   [인사 및 소개]
   안녕하세요. 오늘 발표할 논문은 "Enzyme Cascades for Rare Sugar Biosynthesis"입니다.

   [논문의 중요성]
   - 희귀당 생산의 새로운 패러다임 제시
   - 산업적 응용 가능성 높음

   [발표 순서]
   Background → Methods → Results → Discussion → Conclusion
   ```

   **Example - Content Slide Notes**:
   ```
   [효소 캐스케이드 메커니즘]
   이 슬라이드에서는 3단계 효소 반응의 메커니즘을 설명합니다.

   - 첫 번째 효소 (Isomerase): D-glucose를 D-fructose로 전환
   - 두 번째 효소 (Epimerase): C3 위치의 입체화학 반전
   - 세 번째 효소 (Reductase): 케톤기를 알코올로 환원

   [출처]
   - [1] Kim et al. (2024), DOI: 10.1038/xxxxx
     → Figure 2B (caption: "Enzymatic cascade for L-ribose production...")
   - [3] Park et al. (2023), DOI: 10.1016/xxxxx
     → Figure 1A (caption: "Reaction scheme showing...")

   [발표 포인트]
   - 각 효소의 기질 특이성 강조
   - 전체 수율 85%의 의미 설명
   ```

## Resources

### references/

- **[api_reference.md](references/api_reference.md)**: Database APIs, search strategies, citation metrics, data formats for PubMed/Crossref/preprint servers
- **[styling_guide.md](references/styling_guide.md)**: Academic presentation design patterns, color schemes, typography, layout templates, citation display formats

**When to read**: 
- Read `api_reference.md` when implementing paper search
- Read `styling_guide.md` before creating slide HTML

### scripts/

- **search_papers.py**: Template for paper search implementation (placeholder for actual API integration)
- **generate_presentation.py**: Template for presentation generation workflow (demonstrates content organization and note tracking)

**Note**: Scripts are templates showing structure. Actual implementation uses web_search and web_fetch tools combined with pptx skill workflow.
