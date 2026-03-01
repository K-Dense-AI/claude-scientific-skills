# 시각화 스킬 현황 분석

> 분석 일시: 2026-03-02
> 검색 경로: `c:\Users\Jahyun\claude-scientific-skills\scientific-skills\`

---

## 1. scientific-visualization

### SKILL.md 내용 요약
- **역할**: 학술 논문 출판용 figure 제작을 위한 메타 스킬 (matplotlib/seaborn/plotly를 조율)
- **핵심 기능**:
  - 저널별(Nature, Science, Cell, PLOS 등) 규격에 맞는 figure 생성
  - Multi-panel 레이아웃, 유의성 마커, 에러바, 색맹 안전 팔레트 지원
  - PDF/EPS/TIFF 등 출판 포맷 내보내기
  - 저널별 스타일 프리셋 (`configure_for_journal()`)
  - seaborn 통합 섹션이 매우 상세 (FacetGrid, clustermap, 통계적 시각화 등)
  - Plotly는 간략히 언급 (정적 이미지 내보내기 정도)
- **워크플로우**: Plan → Configure (저널 설정) → Create → Verify (사이즈 체크) → Export → Review
- **통계적 엄격성 강조**: 에러바(SD/SEM/CI), 샘플 수(n), 유의성 마커(*/***), 개별 데이터 포인트 표시
- **피해야 할 사항 10가지**: 작은 폰트, JPEG 사용, 빨강-녹색 조합, 저해상도, 단위 누락, 3D 효과, chart junk, 절단된 축, 불일치 스타일, 에러바 미표시

### references/ 파일 목록 및 핵심
| 파일 | 핵심 내용 |
|------|----------|
| `publication_guidelines.md` | 해상도(300-1200 DPI), 파일 포맷(Vector 우선), 타이포그래피(Arial/Helvetica, 6-12pt), 색상 사용 원칙, 레이아웃 규칙 |
| `color_palettes.md` | Okabe-Ito, Wong 팔레트 RGB 명세, Sequential/Diverging 컬러맵 추천, 접근성 테스트 절차 |
| `journal_requirements.md` | Nature/Science/Cell/PLOS 등 저널별 해상도/사이즈/포맷/폰트 세부 규격 |
| `matplotlib_examples.md` | 10개 완전한 코드 예제 (line plot, bar, heatmap, multi-panel, 저널별 figure 등) |

### assets/ 파일 목록 및 핵심
| 파일 | 핵심 내용 |
|------|----------|
| `color_palettes.py` | Okabe-Ito, Wong, Paul Tol (Bright/Muted/Light/High Contrast) 팔레트 Python 상수, Sequential/Diverging 컬러맵 리스트, `apply_palette()` 헬퍼 함수 |
| `publication.mplstyle` | 범용 출판 품질 matplotlib 스타일 파일 |
| `nature.mplstyle` | Nature 저널 특화 스타일 파일 |
| `presentation.mplstyle` | 발표/포스터용 큰 폰트 스타일 파일 |

### scripts/ 파일 목록
| 파일 | 핵심 내용 |
|------|----------|
| `figure_export.py` | `save_publication_figure()`, `save_for_journal()`, `check_figure_size()` 유틸리티 |
| `style_presets.py` | `apply_publication_style()`, `configure_for_journal()`, `set_color_palette()` 프리셋 |

### 현재 커버리지 분석
- **강점**: 출판용 figure 워크플로우가 매우 체계적. 저널 규격 + 접근성 + 통계적 엄격성을 모두 커버
- **부족한 점**:
  - Table(표) 작성 가이드라인 부재 (학술 논문에서 Figure만큼 중요)
  - Scheme(반응 메커니즘 등) 관련 가이드 부재
  - Supplementary figure 가이드라인 없음
  - 학술 figure numbering/captioning 모범 사례 부족
  - Graphical Abstract 제작 가이드 없음
  - TOC(Table of Contents) graphics 가이드 없음

---

## 2. matplotlib

### SKILL.md 내용 요약
- **역할**: Python 기초 시각화 라이브러리, 세밀한 커스터마이징이 필요할 때 사용
- **핵심 기능**:
  - pyplot(MATLAB 스타일) vs Object-Oriented(fig, ax) 두 가지 인터페이스 설명
  - Figure → Axes → Artist → Axis 계층 구조 설명
  - 주요 플롯 유형: line, scatter, bar, histogram, heatmap, contour, box, violin
  - 3D 시각화 (mpl_toolkits.mplot3d)
  - subplot_mosaic, GridSpec 등 고급 레이아웃
  - rcParams 커스터마이징, 스타일시트 사용법
  - 저장 포맷(PNG/PDF/SVG) 및 DPI 설정
- **Best Practice**: OO 인터페이스 권장, constrained_layout 사용, 컬러맵 선택 가이드, 성능 최적화

### references/ 파일 목록 및 핵심
| 파일 | 핵심 내용 |
|------|----------|
| `plot_types.md` | 모든 플롯 유형의 완전한 카탈로그 + 코드 예제 + 사용 사례 |
| `styling_guide.md` | 상세한 스타일링 옵션, 컬러맵, 커스터마이징 가이드 |
| `api_reference.md` | 핵심 클래스/메서드 레퍼런스 |
| `common_issues.md` | 흔한 문제 해결 가이드 |

### scripts/ 파일 목록
| 파일 | 핵심 내용 |
|------|----------|
| `plot_template.py` | 다양한 플롯 유형 시작 템플릿 |
| `style_configurator.py` | 인터랙티브 스타일 설정 유틸리티 |

### assets/ 파일
- **없음** (assets 디렉토리 자체가 없음)

### 현재 커버리지 분석
- **강점**: matplotlib 핵심 기능을 포괄적으로 다룸. 초보자~고급 사용자 모두 활용 가능
- **부족한 점**:
  - 출판 품질 관련 내용은 scientific-visualization으로 위임되어 있음
  - 학술 특화 figure 유형(Kaplan-Meier, volcano plot, Manhattan plot 등) 부재
  - Animation 관련 내용이 매우 간략

---

## 3. plotly

### SKILL.md 내용 요약
- **역할**: 인터랙티브 시각화 라이브러리. 호버, 줌, 팬 등 웹 기반 차트에 최적
- **핵심 기능**:
  - Plotly Express(고수준 API) vs Graph Objects(저수준 API) 비교 및 선택 기준
  - 40+ 차트 유형 지원 (기본/통계/과학/금융/지도/3D/특수)
  - 서브플롯, 템플릿(plotly_white, plotly_dark 등), 커스터마이징
  - 인터랙티브 기능: 호버 툴팁, 범위 슬라이더, 버튼, 애니메이션
  - 내보내기: HTML(인터랙티브), PNG/PDF/SVG(정적, kaleido 필요)
  - Dash 대시보드 통합
- **과학 데이터**: scatter + trendline, heatmap, 3D surface
- **금융 데이터**: 캔들스틱, 시계열 + rangeslider

### references/ 파일 목록 및 핵심
| 파일 | 핵심 내용 |
|------|----------|
| `plotly-express.md` | 고수준 API 완전 가이드 |
| `graph-objects.md` | 저수준 API 세밀한 제어 가이드 |
| `chart-types.md` | 40+ 차트 유형 전체 카탈로그 + 예제 |
| `layouts-styling.md` | 서브플롯, 템플릿, 색상, 커스터마이징 |
| `export-interactivity.md` | 내보내기 옵션 + 인터랙티브 기능 |

### assets/ 및 scripts/ 파일
- **없음** (둘 다 없음)

### 현재 커버리지 분석
- **강점**: Plotly 핵심 기능과 API 비교가 잘 정리됨. 차트 유형이 매우 풍부
- **부족한 점**:
  - 학술 출판용 정적 이미지 내보내기 워크플로우가 약함
  - 학술 figure 스타일링(Nature/Science 규격 대응) 가이드 없음
  - Plotly 전용 스타일 프리셋/템플릿 파일 없음

---

## 4. seaborn

### SKILL.md 내용 요약
- **역할**: pandas 기반 통계적 시각화 라이브러리. matplotlib 위에 구축
- **핵심 기능**:
  - 설계 철학: Dataset-oriented, Semantic mapping, Statistical awareness, Aesthetic defaults
  - Function Interface(전통적) vs Objects Interface(모던, ggplot2 유사)
  - 카테고리별 플롯 함수:
    - Relational: scatterplot, lineplot, relplot
    - Distribution: histplot, kdeplot, ecdfplot, rugplot, displot, jointplot, pairplot
    - Categorical: stripplot, swarmplot, boxplot, violinplot, boxenplot, barplot, pointplot, countplot, catplot
    - Regression: regplot, lmplot, residplot
    - Matrix: heatmap, clustermap
  - Multi-Plot Grid: FacetGrid, PairGrid, JointGrid
  - Figure-level vs Axes-level 함수 구분 및 선택 기준
  - 색상 팔레트: qualitative, sequential, diverging
  - 테마/스타일: darkgrid, whitegrid, dark, white, ticks + 컨텍스트(paper, notebook, talk, poster)
- **데이터 구조**: Long-form(선호) vs Wide-form + 변환 방법

### references/ 파일 목록 및 핵심
| 파일 | 핵심 내용 |
|------|----------|
| `function_reference.md` | 모든 seaborn 함수의 파라미터 + 예제 |
| `objects_interface.md` | 모던 seaborn.objects API 상세 가이드 |
| `examples.md` | 일반적 사용 사례별 코드 패턴 |

### assets/ 및 scripts/ 파일
- **없음** (둘 다 없음)

### 현재 커버리지 분석
- **강점**: seaborn 함수 카탈로그가 매우 완전. 데이터 구조 요구사항과 인터페이스 선택 가이드가 실용적
- **부족한 점**:
  - 학술 출판 스타일 프리셋 없음 (scientific-visualization에 의존)
  - 도메인 특화 예제 부족 (예: 생물학 실험 데이터, 임상 데이터 시각화)

---

## 5. scientific-schematics

### SKILL.md 내용 요약
- **역할**: AI 기반 학술 다이어그램 생성. Nano Banana Pro로 생성 + Gemini 3 Pro로 품질 리뷰
- **핵심 기능**:
  - 자연어 설명으로 다이어그램 자동 생성 (코딩/템플릿/수동 그리기 불필요)
  - Smart Iterative Refinement: 품질이 임계값 미만일 때만 재생성
  - 문서 유형별 품질 임계값: journal(8.5), conference(8.0), thesis(8.0), grant(8.0), preprint(7.5), poster(7.0), presentation(6.5)
  - 평가 기준 5가지: Scientific Accuracy, Clarity/Readability, Label Quality, Layout/Composition, Professional Appearance (각 0-2점, 총 10점)
  - 지원 다이어그램 유형: CONSORT/PRISMA 플로우차트, 신경망 아키텍처, 생물학적 경로, 회로도, 시스템 아키텍처, 블록 다이어그램
- **프롬프트 엔지니어링**: 레이아웃, 정량적 세부사항, 시각적 스타일, 라벨, 색상 요구 등 포함 권장
- **품질 체크리스트**: 시각적 품질, 접근성, 타이포그래피, 출판 표준, 품질 검증, 문서화/버전 관리

### references/ 파일 목록 및 핵심
| 파일 | 핵심 내용 |
|------|----------|
| `best_practices.md` | 출판 표준, 접근성 가이드라인, 파일 포맷/해상도 요구사항 |
| `QUICK_REFERENCE.md` | 빠른 참조 가이드 |
| `README.md` | 레퍼런스 디렉토리 설명 |

### scripts/ 파일 목록
| 파일 | 핵심 내용 |
|------|----------|
| `generate_schematic.py` | CLI 진입점 |
| `generate_schematic_ai.py` | AI 생성 코어 로직 (ScientificSchematicGenerator) |
| `example_usage.sh` | 사용 예제 |

### assets/ 파일
- **없음**

### 현재 커버리지 분석
- **강점**: AI 기반 자동 생성 + 품질 리뷰 루프가 독창적. 문서 유형별 임계값 차등 적용
- **부족한 점**:
  - 코드 기반(matplotlib/tikz) 다이어그램 생성 옵션이 완전히 제거됨 (AI 전용)
  - 화학 반응 메커니즘 scheme 유형 미지원 (RDKit 등과의 연동 없음)
  - 수학/물리 다이어그램 (상태도, 페이즈 다이어그램 등) 구체적 가이드 부족

---

## 6. scientific-slides

### SKILL.md 내용 요약
- **역할**: 연구 발표용 슬라이드 덱 제작 (PowerPoint + LaTeX Beamer)
- **핵심 기능**:
  - Nano Banana Pro AI로 슬라이드 이미지 자동 생성 + Gemini 3 Pro 품질 리뷰
  - 두 가지 워크플로우: PDF 슬라이드(권장, AI 이미지 기반) / PPTX(python-pptx 기반)
  - LaTeX Beamer 대안 제공
  - 발표 구조 가이드: 서론 → 배경 → 방법 → 결과 → 토론 → 결론
  - 슬라이드 디자인 철학: 시각적으로 매력적, 연구 기반, 최소 텍스트, 전문적 디자인, 스토리 주도
  - 시각적 검증 워크플로우: PDF → 이미지 변환 → Gemini 3 Pro 리뷰
  - 타이밍 가이드라인 (발표 유형별)

### references/ 파일 목록 및 핵심
| 파일 | 핵심 내용 |
|------|----------|
| `presentation_structure.md` | 발표 구조 및 스토리텔링 가이드 |
| `slide_design_principles.md` | 슬라이드 디자인 원칙 (시각적 계층, 타이포그래피, 색상 이론) |
| `data_visualization_slides.md` | 발표용 데이터 시각화 가이드 (저널 figure와의 차이점 강조) |
| `beamer_guide.md` | LaTeX Beamer 완전 가이드 |
| `talk_types_guide.md` | 발표 유형별 가이드 (컨퍼런스, 세미나, 학위논문 방어 등) |
| `visual_review_workflow.md` | 시각적 검증 워크플로우 상세 |

### assets/ 파일 목록 및 핵심
| 파일 | 핵심 내용 |
|------|----------|
| `beamer_template_conference.tex` | 컨퍼런스 발표 Beamer 템플릿 |
| `beamer_template_defense.tex` | 학위논문 방어 Beamer 템플릿 |
| `beamer_template_seminar.tex` | 세미나 Beamer 템플릿 |
| `powerpoint_design_guide.md` | PowerPoint 디자인 가이드 |
| `timing_guidelines.md` | 발표 시간 관리 가이드 |

### scripts/ 파일 목록
| 파일 | 핵심 내용 |
|------|----------|
| `generate_slide_image.py` | AI 슬라이드 이미지 생성 (CLI) |
| `generate_slide_image_ai.py` | AI 생성 코어 로직 |
| `pdf_to_images.py` | PDF → 이미지 변환 |
| `slides_to_pdf.py` | 슬라이드 이미지 → PDF 결합 |
| `validate_presentation.py` | 발표 검증 |

### 현재 커버리지 분석
- **강점**: 발표 전 과정을 다루며, AI 생성 + 검증 루프가 있음. Beamer 템플릿이 실용적
- **부족한 점**:
  - Figure를 위한 슬라이드 최적화 가이드는 있지만, 저널 figure → 슬라이드 변환 자동화는 없음
  - 인터랙티브 슬라이드 (Reveal.js 등) 미지원
  - 발표 노트/스크립트 자동 생성 기능 없음

---

## 7. infographics

### SKILL.md 내용 요약
- **역할**: AI 기반 전문 인포그래픽 생성. Nano Banana Pro + Gemini 3 Pro 리뷰 + Perplexity Sonar 리서치
- **핵심 기능**:
  - 10가지 인포그래픽 유형: statistical, timeline, process, comparison, list, geographic, hierarchical, anatomical, resume, social
  - 8가지 산업 스타일: corporate, healthcare, technology, nature, education, marketing, finance, nonprofit
  - 3가지 색맹 안전 팔레트: Wong, IBM, Tol
  - 리서치 통합(`--research`): Perplexity Sonar로 최신 데이터/통계 자동 수집
  - 문서 유형별 품질 임계값: marketing(8.5), report(8.0), presentation(7.5), social(7.0), internal(7.0), draft(6.5)
  - Smart Iteration: 품질 미달 시에만 재생성
- **평가 기준 5가지**: Visual Hierarchy & Layout, Typography & Readability, Data Visualization, Color & Accessibility, Overall Impact

### references/ 파일 목록 및 핵심
| 파일 | 핵심 내용 |
|------|----------|
| `infographic_types.md` | 10+ 인포그래픽 유형별 확장 템플릿 |
| `design_principles.md` | 시각적 계층, 레이아웃 패턴, 타이포그래피, 색상 가이드 |
| `color_palettes.md` | 팔레트 전체 명세 |

### scripts/ 파일 목록
| 파일 | 핵심 내용 |
|------|----------|
| `generate_infographic.py` | CLI 진입점 |
| `generate_infographic_ai.py` | AI 생성 코어 로직 |

### assets/ 파일
- **없음**

### 현재 커버리지 분석
- **강점**: 다양한 인포그래픽 유형/스타일/팔레트 조합. 리서치 통합이 독창적
- **부족한 점**:
  - 학술 논문과의 직접적 연계가 약함 (마케팅/비즈니스 중심)
  - Graphical Abstract 용도로 활용 가능하지만 명시적 가이드 없음

---

## 8. 기타 시각화 관련 스킬

### 8-1. generate-image
- **역할**: AI 모델(FLUX, Nano Banana 2)을 사용한 범용 이미지 생성/편집
- **시각화 연관성**: 프레젠테이션 시각 자산, 컨셉 아트 등. 기술 다이어그램은 scientific-schematics로 위임
- **학술 figure 관련성**: 낮음 (범용 이미지 생성 용도)

### 8-2. markdown-mermaid-writing
- **역할**: Markdown + Mermaid 다이어그램을 기본 문서 표준으로 확립
- **시각화 연관성**: 텍스트 기반 다이어그램 (flowchart, sequence, gantt, ER 등 24가지 유형)
- **핵심 철학**: 이미지보다 텍스트 다이어그램을 source of truth로 사용 (Git 버전 관리, 토큰 효율성)
- **학술 figure 관련성**: 중간 (워크플로우 다이어그램, 방법론 흐름도 등에 유용하지만 출판 품질 이미지로의 변환이 약함)

### 8-3. networkx
- **역할**: 복잡한 네트워크/그래프 분석 및 시각화
- **시각화 연관성**: 네트워크 토폴로지 시각화, matplotlib 연동
- **학술 figure 관련성**: 높음 (소셜 네트워크, 생물학적 네트워크, PPI 네트워크 등의 시각화)

### 8-4. exploratory-data-analysis
- **역할**: 200+ 과학 파일 포맷에 대한 탐색적 데이터 분석
- **시각화 연관성**: 시각화 추천 기능 포함, 하지만 직접적 figure 생성보다는 분석 중심
- **학술 figure 관련성**: 낮음 (분석 결과 리포트 생성이 주 목적)

### 8-5. shap
- **역할**: ML 모델 해석/설명 (SHAP 값)
- **시각화 연관성**: waterfall, beeswarm, bar, scatter, force, heatmap 등 SHAP 특화 플롯
- **학술 figure 관련성**: 높음 (ML 논문의 핵심 figure 유형)

### 8-6. latex-posters
- **역할**: LaTeX(beamerposter, tikzposter, baposter)로 연구 포스터 제작
- **시각화 연관성**: 포스터에 figure 통합, AI 시각 요소 생성
- **학술 figure 관련성**: 중간 (figure를 포함하는 매체이지만 figure 자체의 품질 가이드보다는 포스터 레이아웃 중심)

### 8-7. pptx-posters
- **역할**: HTML/CSS 기반 연구 포스터 → PDF/PPTX 내보내기
- **시각화 연관성**: latex-posters의 PPTX 대안
- **학술 figure 관련성**: 낮음 (latex-posters 사용이 기본 권장)

---

## 종합 Gap 분석

### 현재 있는 것

| 영역 | 스킬 | 커버리지 수준 |
|------|------|-------------|
| 출판용 데이터 플롯 | scientific-visualization | ★★★★★ (매우 상세) |
| matplotlib 기본기 | matplotlib | ★★★★★ (포괄적) |
| 통계적 시각화 | seaborn | ★★★★★ (함수 카탈로그 완전) |
| 인터랙티브 시각화 | plotly | ★★★★☆ (차트 유형 풍부) |
| 과학 다이어그램 | scientific-schematics | ★★★★☆ (AI 기반) |
| 발표 슬라이드 | scientific-slides | ★★★★☆ (AI + Beamer) |
| 인포그래픽 | infographics | ★★★★☆ (다양한 유형/스타일) |
| 텍스트 기반 다이어그램 | markdown-mermaid-writing | ★★★★☆ (24 유형) |
| 저널별 규격 | scientific-visualization | ★★★★☆ (Nature/Science/Cell/PLOS) |
| 색맹 접근성 | scientific-visualization + 각 스킬 | ★★★★★ (Okabe-Ito/Wong/Tol) |
| 연구 포스터 | latex-posters + pptx-posters | ★★★★☆ |
| 네트워크 시각화 | networkx | ★★★☆☆ (분석 중심) |
| ML 해석 시각화 | shap | ★★★☆☆ (SHAP 특화) |

### 빠져있는 것 (학술 figure/table/scheme 품질 가이드라인)

#### 1. Table(표) 작성 가이드라인 -- 완전 부재
- 학술 논문에서 Figure와 동등하게 중요한 Table 작성 표준이 전혀 없음
- 필요한 내용: 3선 표(three-line table) 규칙, 유효숫자, 단위 표기, 통계 주석, pandas → LaTeX 테이블 변환, 저널별 테이블 포맷

#### 2. Scheme(반응 메커니즘 등) 가이드라인 -- 완전 부재
- 화학 논문의 Scheme (반응식, 합성 경로, 메커니즘) 작성 가이드 없음
- RDKit, ChemDraw 등과의 연동 가이드 필요

#### 3. Graphical Abstract(TOC Graphics) 가이드라인 -- 완전 부재
- 대부분의 저널에서 요구하는 Graphical Abstract/TOC Graphics 제작 가이드 없음
- 저널별 사이즈/해상도 규격, 디자인 원칙, 예제 필요

#### 4. Figure Caption 작성 가이드 -- 부족
- scientific-visualization에서 간략히 언급되지만, 체계적인 figure caption 작성 모범사례가 없음
- 필요한 내용: 구조(제목 문장 → 상세 설명 → 통계 정보), 패널별 설명, 약어 정의, n 수/에러바 정의 등

#### 5. 도메인 특화 시각화 패턴 -- 부족
- **생물학**: volcano plot, MA plot, survival curve(Kaplan-Meier), circos plot, genome browser track
- **화학**: 스펙트럼(NMR, IR, MS) 시각화, 결정구조 시각화
- **의학**: forest plot(메타분석), ROC curve, Bland-Altman plot
- **물리학**: Feynman diagram, phase diagram
- **공학**: Bode plot, Nyquist plot, P&ID diagram

#### 6. Figure Panel 조합 전략 -- 약함
- Multi-panel figure에서 panel 조합 전략 (어떤 데이터를 어떤 패널에 배치할지)
- Figure 스토리텔링 가이드 (한 figure가 전달해야 할 narrative)

#### 7. 재현성(Reproducibility) 가이드 -- 약함
- Figure 생성 코드의 재현성 보장 방법
- 랜덤 시드 고정, 데이터 스냅샷, 환경 고정 등

#### 8. Supplementary Figure 가이드 -- 완전 부재
- 본문 figure와 보충 figure의 품질/포맷 차이
- 보충 figure numbering (Fig. S1, S2 등) 관례

#### 9. 저널 간 Figure 재사용/변환 -- 부재
- 한 저널에서 reject 후 다른 저널 규격으로 변환하는 워크플로우
- 스타일 일괄 변환 스크립트

### 추천 개선 포인트

#### 우선순위 높음 (High)
1. **학술 Table 작성 스킬 신설 또는 scientific-visualization에 통합**: 3선 표 규칙, pandas/LaTeX 테이블 생성, 저널별 규격
2. **Figure Caption 가이드**: 체계적인 caption 작성 모범사례 + 템플릿
3. **Graphical Abstract / TOC Graphics 가이드**: 저널별 규격, 디자인 원칙, 예제 워크플로우
4. **도메인 특화 시각화 레시피북**: volcano plot, Kaplan-Meier, forest plot, Manhattan plot 등 자주 쓰이는 학술 figure 유형별 완전한 코드 + 가이드

#### 우선순위 중간 (Medium)
5. **Scheme 다이어그램 가이드**: RDKit/ChemDraw 연동, 화학 반응식/메커니즘 시각화
6. **Supplementary Figure/Material 가이드**: 포맷, numbering, 품질 기준
7. **Figure 스토리텔링 가이드**: panel 조합 전략, narrative 구성, reviewer가 기대하는 figure 흐름
8. **저널 간 Figure 변환 유틸리티**: 스타일 일괄 변환 스크립트 (예: Nature → Science 규격 변경)

#### 우선순위 낮음 (Low)
9. **Figure 재현성 가이드**: 코드 + 데이터 + 환경 고정 방법론
10. **인터랙티브 Figure 출판 가이드**: Plotly 기반 인터랙티브 supplementary figure 가이드
11. **Animation/Video Figure 가이드**: 보충 자료로 제출하는 동영상 figure 가이드

---

## 스킬 간 관계도

```
scientific-visualization (메타 스킬, 출판용 조율)
    ├── matplotlib (저수준 커스터마이징)
    ├── seaborn (통계적 시각화)
    └── plotly (인터랙티브 → 정적 내보내기)

scientific-schematics (AI 다이어그램)
    └── generate-image (범용 AI 이미지)

scientific-slides (발표)
    ├── scientific-visualization (figure 활용)
    └── scientific-schematics (다이어그램 활용)

infographics (인포그래픽)
    └── scientific-schematics (기술 다이어그램은 위임)

latex-posters / pptx-posters (포스터)
    └── scientific-visualization (figure 활용)

markdown-mermaid-writing (텍스트 다이어그램)
    └── 독립적 (downstream으로 이미지 변환 가능)
```

---

## 결론

현재 시각화 스킬 생태계는 **데이터 플롯 생성**(matplotlib/seaborn/plotly)과 **출판 규격 준수**(scientific-visualization)에서 매우 강력하다. AI 기반 다이어그램(scientific-schematics) 및 슬라이드(scientific-slides) 생성도 독창적이다.

가장 큰 Gap은 **학술 논문의 3대 시각 요소 중 Table과 Scheme에 대한 가이드가 완전히 부재**하다는 점이다. 또한 Graphical Abstract, 도메인 특화 플롯 유형, Figure caption 작성법 등 **실제 논문 작성 시 빈번히 필요한 가이드**가 부족하다. 이 Gap들을 메우면 학술 시각화 워크플로우를 end-to-end로 완전히 커버할 수 있을 것이다.
