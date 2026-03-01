# 학술 Figure/Table/Scheme 가이드라인 수집

> 수집일: 2026-03-02
> 출처: GitHub (ScientificFigures, paper-tips-and-tricks, matplotlib_for_papers),
>       PLOS CompBiol Ten Simple Rules, Nature/ACS/Elsevier 공식 가이드라인

---

## 1. 해상도 및 포맷

### 1.1 해상도 (DPI) 요구사항

| 이미지 유형 | Nature | ACS | Elsevier | 권장 기본값 |
|-------------|--------|-----|----------|-------------|
| Line art (선화/그래프) | 최소 300, 권장 450+ | 1200 dpi | 1000 dpi | **1200 dpi** |
| Halftone (사진/현미경) | 최소 300, 권장 450 | 300 dpi | 300 dpi | **300 dpi** |
| Combination (선+사진) | 최소 300, 권장 450 | 600 dpi | 500 dpi | **600 dpi** |
| Grayscale art | 최소 300 | 600 dpi | 500 dpi | **600 dpi** |

**구현 규칙:**
- `matplotlib.savefig(dpi=300)` 최소값; 선 그래프는 `dpi=600` 이상 사용
- 벡터 포맷(PDF/SVG/EPS) 사용 시 DPI 설정은 래스터 요소에만 적용
- 최종 인쇄 크기에서 DPI를 계산: `pixels = DPI x print_size_inches`

### 1.2 파일 형식 우선순위

**벡터 형식 (최우선):**
1. **PDF** - pdfLaTeX 호환, 가장 범용적. matplotlib `savefig('fig.pdf')`
2. **EPS** - 전통 LaTeX 호환. Nature/ACS/Elsevier 모두 허용
3. **SVG** - 후편집(Inkscape/Illustrator) 용이. 웹 호환

**래스터 형식 (벡터 불가 시):**
4. **TIFF** - 무손실 압축. Elsevier 선호. LZW 압축 권장
5. **PNG** - 무손실. 웹 및 프리프린트용
6. **JPEG** - **피할 것**. 손실 압축으로 텍스트/선 품질 저하. 사진만 불가피 시 최고 품질(quality=95+)

**절대 사용 금지:**
- GIF (색상 제한 256색)
- BMP (비압축, 파일 과대)
- 래스터화된 스크린캡처

**구현 규칙:**
```python
# matplotlib 기본 저장 설정
SAVE_FORMATS = {
    'publication': 'pdf',    # 논문 제출용
    'presentation': 'png',   # 발표용
    'editing': 'svg',        # 후편집용
    'fallback_raster': 'tiff'  # 래스터 필수 시
}
# JPEG 자동 차단: 파일명이 .jpg/.jpeg이면 경고 후 .pdf로 변환 제안
```

### 1.3 파일 크기 제한

- Nature: 개별 파일 **10 MB** 이하
- Elsevier: 개별 파일 **10 MB** 이하, 다수 figure 시 각 **7 MB** 이하
- ACS: 명시적 제한 없으나 10 MB 이하 권장

**구현 규칙:**
- 데이터 포인트 과다 시 래스터화: `ax.set_rasterization_zorder(1)` 또는 `rasterized=True`
- 텍스트/축은 벡터 유지, 데이터만 래스터화하여 파일 크기 절감

---

## 2. 색상 접근성

### 2.1 색맹 안전 팔레트

**Okabe-Ito 팔레트 (Nature 권장, 8색):**
```python
OKABE_ITO = {
    'orange':    '#E69F00',
    'sky_blue':  '#56B4E9',
    'green':     '#009E73',
    'yellow':    '#F0E442',
    'blue':      '#0072B2',
    'vermilion': '#D55E00',
    'purple':    '#CC79A7',
    'black':     '#000000',
}
```

**Viridis 계열 (연속 데이터용):**
- `viridis` - 기본 권장 (perceptually uniform)
- `magma` - 어두운 배경 강조
- `plasma` - 따뜻한 톤 강조
- `inferno` - 고대비
- `cividis` - 색각이상 최적화

### 2.2 색상 사용 규칙

**금지 조합:**
- 빨강 + 초록 단독 사용 (적녹색맹 8% 남성에게 구별 불가)
- Rainbow/Jet 컬러맵 (비균일 인지, 왜곡 유발)

**안전한 대체 조합:**
- 주황 + 파랑
- 보라 + 노랑
- 형광현미경: 마젠타/초록/파랑 또는 마젠타/노랑/시안

**히트맵 규칙:**
- 양극 데이터: 두 보색 끝점 + 흰색/검정 중앙
- 순차 데이터: 단일 색조 연속 (viridis 계열)
- 발산 데이터: 중앙 기준 대칭 2색 (예: 파랑-흰색-빨강)

**구현 규칙:**
```python
# 기본 색상 순환 설정
DEFAULT_COLOR_CYCLE = [
    '#0072B2',  # blue
    '#D55E00',  # vermilion
    '#009E73',  # green
    '#E69F00',  # orange
    '#CC79A7',  # purple
    '#56B4E9',  # sky blue
    '#F0E442',  # yellow
]

# 색상 외 구별 수단 병행 필수
DISTINGUISH_METHODS = ['color', 'linestyle', 'marker', 'hatch_pattern']
# 3개 이상 범주: 색상 + 마커 또는 선 스타일 병행
# 5개 이상 범주: 색상 + 마커 + 선 스타일 모두 병행
```

### 2.3 검증 도구

- **Color Oracle** (데스크톱): 실시간 색맹 시뮬레이션
- **Coblis** (웹): 이미지 업로드 색맹 검증
- **Viz Palette** (웹): HEX 코드 입력 후 접근성 확인
- 그레이스케일 변환 테스트: `fig.savefig('test_gray.pdf', ...) + 변환 확인`

**구현 규칙:**
```python
# 자동 그레이스케일 검증 함수
import numpy as np
from PIL import Image

def check_grayscale_distinguishability(image_path, n_colors):
    """그레이스케일 변환 후 색상 구별 가능 여부 확인"""
    img = Image.open(image_path).convert('L')  # 그레이스케일 변환
    arr = np.array(img)
    unique_grays = len(np.unique(arr))
    # 최소 구별 명도 차이: 256/n_colors 이상
    return unique_grays >= n_colors * 3  # 안전 마진 포함
```

---

## 3. 타이포그래피

### 3.1 폰트 패밀리

| 용도 | Nature | ACS | Elsevier | 권장 |
|------|--------|-----|----------|------|
| 기본 텍스트 | Arial, Helvetica | Helvetica, Arial | - | **Arial/Helvetica** |
| 아미노산 서열 | Courier | - | - | **Courier New** |
| 그리스 문자 | Symbol | - | - | **Symbol** |
| 수학 기호 | - | - | - | **Computer Modern** (LaTeX) |

**구현 규칙:**
```python
import matplotlib as mpl

FONT_CONFIG = {
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7,          # 기본 텍스트
    'mathtext.fontset': 'dejavusans',  # 수학 기호
}
mpl.rcParams.update(FONT_CONFIG)
```

### 3.2 폰트 크기

| 요소 | Nature | ACS | Elsevier | 최소 | 권장 |
|------|--------|-----|----------|------|------|
| 일반 텍스트 | 5-7 pt | 최소 4.5 pt | 7 pt | **6 pt** | **7-8 pt** |
| 패널 레이블 | 8 pt bold | - | - | **8 pt** | **8 pt bold** |
| 축 레이블 | 5-7 pt | 최소 4.5 pt | 7 pt | **6 pt** | **7 pt** |
| 축 눈금 | 5-7 pt | 최소 4.5 pt | 6 pt | **5 pt** | **6 pt** |
| 범례 | 5-7 pt | 최소 4.5 pt | 6 pt | **5 pt** | **6 pt** |
| 첨자/위첨자 | - | - | 최소 6 pt | **5 pt** | **6 pt** |

**핵심 원칙:**
- 최종 인쇄 크기 기준으로 폰트 크기 설정 (축소 후 크기 확인)
- Figure를 스크립트에서 최종 크기로 생성 (LaTeX에서 리사이징 금지)
- 모든 figure의 폰트 크기 일관성 유지

**구현 규칙:**
```python
# 최종 인쇄 크기에서의 폰트 크기 계산
def calculate_font_size(target_pt, figure_width_mm, column_width_mm):
    """Figure가 최종 컬럼 폭으로 축소될 때 적정 폰트 크기 역산"""
    scale_factor = column_width_mm / figure_width_mm
    return target_pt / scale_factor

# 예: 180mm로 만든 figure가 89mm 컬럼에 축소 시
# 7pt 목표 -> 실제 설정: 7 / (89/180) = 14.2pt
```

### 3.3 텍스트 서식 규칙

- **패널 레이블**: 소문자 볼드 (a, b, c...) - Nature 표준
- **축 레이블**: 단위는 괄호 안에 표기, 예: `"Temperature (°C)"`
- **축 눈금**: 소수점 자릿수 일관 유지
- **범례**: 데이터와 겹치지 않는 위치에 배치
- **텍스트 임베딩**: TrueType 2 또는 42 사용 (Type 3 금지)
- **텍스트 아웃라인화 금지**: 편집 가능 상태 유지

---

## 4. Figure 크기

### 4.1 저널별 Figure 폭 사양

| 저널 | 단일 컬럼 | 1.5 컬럼 | 이중 컬럼 | 최대 높이 |
|------|-----------|----------|-----------|-----------|
| **Nature** | 89 mm (3.5 in) | 120-136 mm | 183 mm (7.2 in) | 170 mm (6.7 in) |
| **ACS** | 84 mm (3.33 in) | - | 175 mm (7.0 in) | 233 mm (9.17 in) |
| **Elsevier** | 90 mm (3.54 in) | 140 mm (5.51 in) | 190 mm (7.48 in) | - |

**구현 규칙:**
```python
# 저널별 Figure 크기 사전 (mm)
JOURNAL_FIGURE_SIZES = {
    'nature': {
        'single_col': 89,
        'one_half_col': 120,
        'double_col': 183,
        'max_height': 170,
        'max_page': (183, 247),
    },
    'acs': {
        'single_col': 84,
        'double_col': 175,
        'max_height': 233,  # 캡션 포함
    },
    'elsevier': {
        'single_col': 90,
        'one_half_col': 140,
        'double_col': 190,
        'min_width': 30,
    },
}

def mm_to_inches(mm):
    return mm / 25.4

def get_figure_size(journal, layout='single_col', aspect_ratio=0.75):
    """저널 규격에 맞는 figure 크기 (inches) 반환"""
    w_mm = JOURNAL_FIGURE_SIZES[journal][layout]
    w_in = mm_to_inches(w_mm)
    h_in = w_in * aspect_ratio
    return (w_in, h_in)
```

### 4.2 황금비/표준 비율

- **기본 비율**: 4:3 (가로:세로) 또는 황금비(1:0.618)
- **넓은 데이터**: 16:9 (시계열 등)
- **정사각형**: 1:1 (히트맵, 상관행렬)
- 비율을 스크립트에서 명시적으로 설정하고, LaTeX에서 리사이징하지 않을 것

---

## 5. 테이블

### 5.1 Booktabs 스타일 (필수)

```latex
\usepackage{booktabs}

\begin{table}[t]  % 테이블은 상단 배치 권장
  \caption{캡션은 테이블 위에 배치}  % 핵심 규칙
  \label{tab:example}
  \centering
  \begin{tabular}{lcc}  % l=왼쪽, c=가운데, r=오른쪽
    \toprule
    항목 & 값 1 & 값 2 \\
    \midrule
    데이터 A & 1.23 & 4.56 \\
    데이터 B & 7.89 & 0.12 \\
    \addlinespace  % 그룹 간 간격
    데이터 C & 3.45 & 6.78 \\
    \bottomrule
  \end{tabular}
\end{table}
```

### 5.2 테이블 서식 규칙

**필수 규칙:**
- 세로선(vertical lines) 절대 사용 금지
- `\toprule`, `\midrule`, `\bottomrule`만 사용 (일반 `\hline` 금지)
- 그룹 구분: `\cmidrule{2-3}` (부분 가로선) 또는 `\addlinespace` (간격)
- 캡션은 테이블 **위에** 배치 (Figure 캡션은 **아래에**)

**정렬 규칙:**
- 텍스트 열: 왼쪽 정렬 (`l`)
- 숫자 열: 소수점 정렬 (`S` from siunitx) 또는 오른쪽 정렬 (`r`)
- 헤더: 문장 형식 대문자 (첫 글자만 대문자)
- 단위: 헤더에 괄호로 표기, 예: `"Mass (kg)"`

**숫자 서식:**
- 유효숫자 일관 유지
- 불확실도 표기: `1.23 ± 0.05` 또는 `1.23(5)`
- 볼드/색상으로 최적값 강조 가능 (단, 일관적으로)

**구현 규칙:**
```python
# pandas DataFrame -> booktabs LaTeX 변환
def df_to_booktabs(df, caption, label):
    """DataFrame을 booktabs 스타일 LaTeX 테이블로 변환"""
    latex = df.to_latex(
        index=False,
        escape=True,
        column_format='l' + 'c' * (len(df.columns) - 1),
        caption=caption,
        label=label,
        position='t',
    )
    # \hline -> booktabs 변환
    latex = latex.replace('\\hline', '')
    # 헤더 후 \midrule 삽입 등 후처리
    return latex
```

---

## 6. 스키마/다이어그램

### 6.1 패널 배치

**레이블 규칙:**
- 소문자 볼드: **a**, **b**, **c** (Nature 표준)
- 패널 레이블 크기: 8 pt bold
- 위치: 각 패널 좌상단, 데이터와 겹치지 않도록
- 패널 간 일관된 간격 유지

**레이아웃 원칙:**
- 읽기 순서: 좌→우, 상→하
- 관련 패널은 인접 배치
- 패널 간 간격: 2-5 mm

**구현 규칙:**
```python
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def create_multi_panel(n_rows, n_cols, journal='nature', layout='double_col'):
    """저널 규격 멀티패널 Figure 생성"""
    sizes = JOURNAL_FIGURE_SIZES[journal]
    w_mm = sizes[layout]
    w_in = mm_to_inches(w_mm)
    h_in = w_in * (n_rows / n_cols) * 0.8  # 종횡비 조정

    fig = plt.figure(figsize=(w_in, h_in))
    gs = gridspec.GridSpec(n_rows, n_cols, figure=fig,
                           hspace=0.3, wspace=0.3)

    axes = []
    labels = 'abcdefghijklmnopqrstuvwxyz'
    for i in range(n_rows):
        for j in range(n_cols):
            idx = i * n_cols + j
            ax = fig.add_subplot(gs[i, j])
            # 패널 레이블 추가
            ax.text(-0.15, 1.05, labels[idx],
                    transform=ax.transAxes,
                    fontsize=8, fontweight='bold',
                    va='bottom', ha='right')
            axes.append(ax)

    return fig, axes
```

### 6.2 화살표 및 주석

**화살표 규칙:**
- 화학 반응 화살표: 일관된 길이와 스타일
- 플로우차트 화살표: 직각 또는 직선, 곡선 최소화
- 화살촉: 데이터 겹침 방지, 크기 일관 유지

**주석 규칙:**
- `ax.annotate()` 사용, 직접 텍스트보다 선호
- 지시선(leader line)은 얇게 (0.5-1 pt)
- 배경색이 복잡한 경우 텍스트에 반투명 배경 박스 추가

**구현 규칙:**
```python
# 주석 스타일 사전
ANNOTATION_STYLE = {
    'fontsize': 7,
    'fontweight': 'normal',
    'arrowprops': dict(
        arrowstyle='->',
        color='black',
        lw=0.75,
        connectionstyle='arc3,rad=0.1',
    ),
    'bbox': dict(
        boxstyle='round,pad=0.3',
        facecolor='white',
        edgecolor='none',
        alpha=0.8,
    ),
}
```

### 6.3 스케일바

- 확대 배율 대신 **스케일바** 사용 (Nature 필수)
- 스케일바는 별도 편집 가능 레이어에 배치
- 이미지 위에 텍스트 직접 오버레이 금지 (복잡한 배경)

---

## 7. 저널별 요구사항 요약

### 7.1 Nature 계열

| 항목 | 사양 |
|------|------|
| 포맷 | PDF/EPS (벡터), TIFF (래스터) |
| 색상 모드 | RGB 권장 |
| 최소 DPI | 300 (사진), 450 권장 |
| 폰트 | Arial/Helvetica, 5-7 pt |
| 패널 레이블 | 8 pt bold 소문자 |
| 단일 컬럼 | 89 mm |
| 이중 컬럼 | 183 mm |
| 최대 높이 | 170 mm |
| 파일 크기 | 10 MB 이하 |
| 특이사항 | 스케일바 필수, 접근성 팔레트 필수, 텍스트 아웃라인화 금지 |

### 7.2 ACS 계열

| 항목 | 사양 |
|------|------|
| 포맷 | TIFF, EPS, PDF |
| 색상 모드 | RGB 또는 CMYK |
| Line art DPI | 1200 dpi |
| Halftone DPI | 300 dpi |
| 폰트 | Helvetica/Arial, 최소 4.5 pt (인쇄 후) |
| 단일 컬럼 | 84 mm (3.33 in) |
| 이중 컬럼 | 175 mm (7.0 in) |
| 최대 높이 | 233 mm (캡션 포함) |
| 선 두께 | 최소 0.5 pt |
| 특이사항 | 캡션 12 pt/줄 높이 허용, 유사 figure 간 일관된 해상도 |

### 7.3 Elsevier 계열

| 항목 | 사양 |
|------|------|
| 포맷 | EPS (벡터 선호), TIFF, PDF, JPEG |
| 색상 모드 | RGB 또는 CMYK |
| Line art DPI | 1000 dpi |
| Combination DPI | 500 dpi |
| Halftone DPI | 300 dpi |
| 폰트 | 7 pt 일반, 6 pt 첨자 |
| 단일 컬럼 | 90 mm |
| 1.5 컬럼 | 140 mm |
| 이중 컬럼 | 190 mm |
| 최소 폭 | 30 mm |
| 파일 크기 | 10 MB 이하 (다수 시 7 MB) |

### 7.4 공통 패턴 (3대 출판사 교집합)

```python
COMMON_REQUIREMENTS = {
    'min_dpi_photo': 300,
    'min_dpi_line': 1000,       # 안전한 공통 최소값
    'preferred_vector': ['pdf', 'eps'],
    'preferred_raster': ['tiff'],
    'color_mode': 'RGB',        # 3곳 모두 RGB 허용
    'min_font_pt': 6,           # 안전한 공통 최소값
    'recommended_font_pt': 7,
    'font_family': ['Arial', 'Helvetica'],
    'panel_label_pt': 8,
    'panel_label_style': 'bold lowercase',
    'max_file_mb': 10,
    'min_line_width_pt': 0.5,
    'caption_figure': 'below',  # Figure 캡션은 아래
    'caption_table': 'above',   # Table 캡션은 위
}
```

---

## 8. 흔한 실수들 (안티패턴)

### 8.1 해상도/포맷 실수

| 안티패턴 | 문제점 | 해결책 |
|----------|--------|--------|
| 스크린캡처를 figure로 사용 | 72 dpi, 래스터화 | 원본 데이터에서 벡터 재생성 |
| JPEG로 그래프 저장 | 텍스트/선 주변 아티팩트 | PDF/SVG 사용 |
| LaTeX에서 figure 리사이징 | 폰트 크기 불일치 | 스크립트에서 최종 크기로 생성 |
| PowerPoint에서 figure 생성 | 낮은 해상도, 비표준 포맷 | matplotlib/R 등 전용 도구 사용 |
| 래스터 figure를 확대 | 픽셀화 | 벡터 포맷으로 원본 생성 |

**구현 규칙 (자동 검증):**
```python
def validate_resolution(image_path, target_width_mm, image_type='line'):
    """이미지 해상도가 인쇄 요구사항을 충족하는지 검증"""
    from PIL import Image
    img = Image.open(image_path)
    width_px = img.size[0]
    width_in = target_width_mm / 25.4
    actual_dpi = width_px / width_in

    min_dpi = {'line': 1000, 'combination': 500, 'halftone': 300}
    required = min_dpi.get(image_type, 300)

    if actual_dpi < required:
        return False, f"DPI {actual_dpi:.0f} < 요구 {required} dpi"
    return True, f"DPI {actual_dpi:.0f} OK (>= {required})"
```

### 8.2 색상 실수

| 안티패턴 | 문제점 | 해결책 |
|----------|--------|--------|
| Rainbow/Jet 컬러맵 | 비균일 인지, 밝기 불연속 | viridis/magma/cividis 사용 |
| 빨강+초록만으로 구분 | 8% 남성 구별 불가 | 파랑+주황 또는 마커 병행 |
| 6개 이상 색상 사용 | 사전주의적(pre-attentive) 인지 한계 초과 | 5-7색 이하, 그룹화 활용 |
| 컬러로만 범주 구분 | 흑백 인쇄 시 정보 손실 | 마커/선 스타일/해칭 병행 |
| 의미 역전 색상 | 빨강=이득, 초록=손실 (문화적 혼동) | 관례적 색상-의미 매핑 준수 |

### 8.3 차트 유형 실수

| 안티패턴 | 문제점 | 해결책 |
|----------|--------|--------|
| 파이 차트 | 각도 비교 인간 인지 부정확 | 막대 그래프 사용 |
| 3D 막대/원 차트 | 왜곡된 크기 인지 | 2D 시각화 고수 |
| 이중 y축 | 스케일 조작으로 상관관계 왜곡 가능 | 별도 패널 분리 |
| 원 크기로 값 인코딩 | 면적 vs 지름 혼동 (6.5x 시각 왜곡) | 면적 비례 설정 또는 막대 사용 |
| 잘린 y축 | 차이 과장 | 0부터 시작 또는 명시적 축 절단 표시 |

**구현 규칙 (경고 시스템):**
```python
CHART_WARNINGS = {
    'pie': "파이 차트 대신 막대 그래프를 권장합니다 (각도 인지 부정확)",
    '3d_bar': "3D 막대 차트는 크기 왜곡을 유발합니다. 2D 사용을 권장합니다",
    'dual_yaxis': "이중 y축은 오해의 소지가 있습니다. 별도 패널을 권장합니다",
    'rainbow_cmap': "Rainbow/Jet 컬러맵은 비균일합니다. viridis를 권장합니다",
    'truncated_axis': "y축이 0에서 시작하지 않으면 차이가 과장될 수 있습니다",
}
```

### 8.4 레이아웃/타이포그래피 실수

| 안티패턴 | 문제점 | 해결책 |
|----------|--------|--------|
| 기본 설정(defaults) 그대로 사용 | "기본은 어디에나 적당하지만 최적은 아님" | 모든 요소 수동 조정 |
| 축 레이블 없음 | 단위/변수명 추정 불가 | 항상 축 레이블 + 단위 표기 |
| 범례가 데이터 위 | 데이터 가림 | figure 밖 또는 빈 공간에 배치 |
| 캡션 생략 | 독립적 해석 불가 | 읽는 방법 + 핵심 정보 포함 캡션 |
| Figure 간 폰트 불일치 | 비전문적 인상 | 통일된 스타일 시트 사용 |
| 불필요한 격자선 | 시각적 잡음 | 격자선 제거 또는 매우 연하게 |
| 과도한 장식 | "차트정크" | 데이터/잉크 비율 최대화 |
| 테이블에 세로선 | 시각적 잡음 | booktabs 스타일 사용 |

### 8.5 데이터 표현 실수

| 안티패턴 | 문제점 | 해결책 |
|----------|--------|--------|
| 평균±표준편차만 표시 | 정규분포 가정 위반 가능 | 중앙값+사분위수 또는 바이올린 플롯 |
| n수 미표기 | 통계적 신뢰도 판단 불가 | 항상 n 표시 |
| 모든 데이터를 하나의 figure에 | 과부하 | 여러 focused figure로 분리 |
| 시계열 순서 변경 | 시간적 논리 파괴 | 시간순 유지 |
| 개별 데이터 포인트 숨김 | 분포/이상치 은폐 | 박스플롯 위에 개별 점 오버레이 |

---

## 9. Ten Simple Rules 요약 (PLOS CompBiol)

| 규칙 | 핵심 내용 | 구현 지침 |
|------|-----------|-----------|
| 1. 청중을 알라 | 발표 대상에 따라 복잡도 조절 | `audience` 파라미터로 스타일 분기 |
| 2. 메시지를 먼저 정의하라 | Figure당 하나의 명확한 메시지 | 캡션 첫 문장에 메시지 명시 |
| 3. 매체에 맞춰라 | 발표/논문/포스터 별 스타일 분리 | `medium` 파라미터: 'paper'/'poster'/'slide' |
| 4. 캡션은 필수 | 읽는 방법 + 수치 정보 포함 | 캡션 없는 figure 경고 |
| 5. 기본값을 믿지 마라 | 모든 요소 수동 조정 | 커스텀 스타일 시트 적용 |
| 6. 색상을 효과적으로 | 강조 대상만 색상, 나머지 회색 | 기본 회색, 강조 색상 전략 |
| 7. 독자를 오도하지 마라 | 정직한 축 스케일, 적절한 차트 유형 | 축 범위 자동 검증 |
| 8. 차트정크 제거 | 불필요한 격자/배경/장식 제거 | `despine()` + 최소 격자 |
| 9. 메시지 > 미적 요소 | 명확성이 최우선 | 가독성 우선 레이아웃 |
| 10. 적절한 도구 사용 | matplotlib, R, Inkscape 등 목적별 도구 | 도구 추천 시스템 |

---

## 10. matplotlib 기본 스타일 설정 (통합)

```python
"""
학술 논문용 matplotlib 통합 스타일 설정
모든 저널 공통 요구사항의 교집합 기반
"""
import matplotlib as mpl

PUBLICATION_STYLE = {
    # 폰트
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 7,
    'axes.labelsize': 7,
    'axes.titlesize': 8,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'legend.fontsize': 6,

    # 선/마커
    'lines.linewidth': 1.0,
    'lines.markersize': 4,
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.minor.width': 0.4,
    'ytick.minor.width': 0.4,
    'xtick.major.size': 3,
    'ytick.major.size': 3,
    'xtick.minor.size': 1.5,
    'ytick.minor.size': 1.5,

    # 색상 (Okabe-Ito 기반)
    'axes.prop_cycle': mpl.cycler('color', [
        '#0072B2', '#D55E00', '#009E73', '#E69F00',
        '#CC79A7', '#56B4E9', '#F0E442',
    ]),

    # 레이아웃
    'figure.dpi': 150,         # 화면 표시용
    'savefig.dpi': 300,        # 저장 시 최소 DPI
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.02,
    'figure.constrained_layout.use': True,

    # 격자/프레임
    'axes.grid': False,        # 기본 격자 비활성화
    'axes.spines.top': False,  # 상단 프레임 제거
    'axes.spines.right': False,  # 우측 프레임 제거

    # 범례
    'legend.frameon': False,
    'legend.loc': 'best',

    # 저장
    'savefig.format': 'pdf',
    'pdf.fonttype': 42,        # TrueType (Type 3 방지)
    'ps.fonttype': 42,
}

def apply_publication_style():
    """학술 논문 스타일 일괄 적용"""
    mpl.rcParams.update(PUBLICATION_STYLE)
```

---

## 11. 출처

### GitHub 리포지토리
- [nrokh/ScientificFigures](https://github.com/nrokh/ScientificFigures) - Figure 평가 루브릭 및 가이드
- [Wookai/paper-tips-and-tricks](https://github.com/Wookai/paper-tips-and-tricks) - 논문 작성 팁 (Figure/Table 포함)
- [jbmouret/matplotlib_for_papers](https://github.com/jbmouret/matplotlib_for_papers) - matplotlib 논문용 설정

### 학술 논문
- [Ten Simple Rules for Better Figures (PLOS CompBiol)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003833)
- [Examining data visualization pitfalls in scientific publications (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8556474/)
- [Choosing color palettes for scientific figures (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7040535/)

### 저널 공식 가이드라인
- [Nature - Preparing figures specifications](https://research-figure-guide.nature.com/figures/preparing-figures-our-specifications/)
- [Nature - Final submission](https://www.nature.com/nature/for-authors/final-submission)
- [ACS - Preparing Manuscript Graphics](https://pubs.acs.org/page/4authors/submission/graphics_prep.html)
- [Elsevier - Artwork sizing](https://www.elsevier.com/about/policies-and-standards/author/artwork-and-media-instructions/artwork-sizing)
- [Elsevier - Artwork overview](https://www.elsevier.com/about/policies-and-standards/author/artwork-and-media-instructions/artwork-overview)

### 색상 접근성
- [NKI - Guidelines for color blind friendly figures](https://www.nki.nl/about-us/responsible-research/guidelines-color-blind-friendly-figures)
- [ASCB - How to make figures accessible to color-blind readers](https://www.ascb.org/diversity-equity-and-inclusion/how-to-make-scientific-figures-accessible-to-readers-with-color-blindness/)
