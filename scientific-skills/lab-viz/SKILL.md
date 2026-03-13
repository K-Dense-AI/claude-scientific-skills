---
name: lab-viz
description: >
  Unified lab data visualization skill (Haiku-powered for routine plots).
  Modes: general plots, HPLC chromatograms + deconvolution, kinetic progress
  curves, surface/contour plots, timecourse uncertainty bands, BO diagnostics,
  publication diagrams. Auto-selects mode from context. Okabe-Ito palette,
  Arial 8pt, despine, 300 DPI.
license: MIT license
metadata:
    skill-author: lab-internal
---

# Lab-Viz (실험 데이터 시각화 전용 스킬)

## 스킬 선택 기준 (시각화 라우팅)

| 요청 유형 | 사용 스킬 |
|----------|----------|
| 논문/publication figure, 저널 규격, 멀티패널 | **scientific-visualization** (메타스킬) |
| HPLC/크로마토그램/키네틱/BO/MC/surface plot | **lab-viz** (이 스킬) |
| rcParams 세밀 제어, 커스텀 axes 조작 | matplotlib/seaborn 직접 호출 |

> 이 스킬은 **실험실에서 생성된 데이터의 루틴 시각화 전용**이다.
> 논문 제출용 publication-quality figure가 필요하면 `scientific-visualization`을 사용할 것.

실험 데이터 시각화를 위한 통합 스킬. Haiku 에이전트로 저비용 처리.

## 모드 판별

| 진입 조건 | 모드 |
|----------|------|
| scatter / bar / line / 추세선 / 일반 | **A: 범용 시각화** |
| HPLC / 크로마토그램 / 피크 / 보유시간 | **B: 크로마토그램** (Haiku) |
| deconvolution / 피크 겹침 / Gaussian fit | **C: Deconvolution** (Haiku) |
| progress curve / kinetic / NADH / 흡광도 / 시계열 | **D: 키네틱 Progress Curve** (Haiku) |
| surface plot / heatmap / 농도 조합 / 2D 공간 / contour | **E: 표면 플롯** (Haiku) |
| Monte Carlo / MC bands / 불확실도 / 민감도 워터폴 | **F: 불확실도 시각화** (Haiku) |
| BO / Pareto / GP / 베이지안 최적화 결과 | **G: BO 진단 플롯** (Haiku) |
| 시스템 다이어그램 / 흐름도 / 반응 네트워크 | **H: 교육/논문용 다이어그램** (Haiku) |

> 모드 B~H: **반드시 Agent tool (model: haiku)으로 서브에이전트 실행**
> 모드 A: CEO가 직접 처리

---

## 공통 스타일 규칙 (전 모드 적용)

```python
OKABE_ITO = ['#E69F00','#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7','#000000']

import matplotlib as mpl
mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica'],
    'font.size': 8, 'axes.labelsize': 9,
    'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'axes.linewidth': 0.8,
})
# despine (top/right spine 제거), 300 DPI, bbox_inches='tight'
# 저장 후 start "" "<절대경로>" 로 자동 오픈
# 절대경로 사용: os.path.dirname(os.path.abspath(__file__))
```

---

## 모드 A: 범용 시각화 — 추세선 자동 선택

| 데이터 패턴 | 추세선 유형 |
|------------|-----------|
| 단조 증가/감소 | 선형 회귀 (y = ax + b) |
| 포화/정체 곡선 | 비선형 (Michaelis-Menten 등) |
| 최적점 존재 | 2차 다항식 |
| 지수적 변화 | 지수 피팅 |
| 포인트 < 3 또는 R² < 0.5 | 추세선 미적용 |

필수 출력: 방정식, R², 주요 지점(최대값/포화점/변곡점).

---

## 모드 B: HPLC 크로마토그램 (Haiku 에이전트)

**트리거**: "HPLC", "크로마토그램", "피크", "보유시간", "chromatogram", "peak annotation"

```
Agent(subagent_type="general-purpose", model="haiku", prompt=<아래>)
```

```
입력: {data_file}  (CSV 또는 Excel, time/intensity 컬럼)
출력: {output_dir}
피크 파일: {peaks_file} (선택, RT/name/height 포함)

## 처리 방식
- X축: 보유시간(min), Y축: 흡광도/강도
- 피크 주석: 수직 점선 + 화합물명 레이블 (fontsize=6)
- 피크 면적 음영 (있으면)
- 베이스라인 표시 (있으면)

## 참조 코드
# C:\Users\Jahyun\PeakPicker\src\peakpicker\infrastructure\exporters\plot_exporter.py
# ChromatogramPlotExporter.export_chromatogram()

출력: {basename}_chromatogram.png
파일 경로 + 1줄 요약만 보고
```

---

## 모드 C: Deconvolution 시각화 (Haiku 에이전트)

**트리거**: "deconvolution", "피크 겹침", "Gaussian fit", "피크 분리", "R²"

```
Agent(subagent_type="general-purpose", model="haiku", prompt=<아래>)
```

```
입력: {data_file}  (time/intensity + 피크 목록)
출력: {output_dir}
유형: single | batch | summary

### single — 개별 피크 2단계 시각화
- 위: 원본 + Gaussian 피크 오버레이
- 아래: 잔차 (residuals)

### batch — 전체 크로마토그램 deconvolution 오버레이

### summary — 4패널 통계 요약
- 히스토그램: 피크 개수, R², RMSE, shoulder 분포

## 참조 코드
# C:\Users\Jahyun\PeakPicker\src\deconvolution_visualizer.py
# DeconvolutionVisualizer.plot_single_deconvolution()

출력: {basename}_deconvolution_{유형}.png
```

---

## 모드 D: 키네틱 Progress Curve (Haiku 에이전트)

**트리거**: "progress curve", "kinetic data 시각화", "NADH 그래프", "조건별 그래프", "fitted curve"

```
Agent(subagent_type="general-purpose", model="haiku", prompt=<아래>)
```

```
입력: {input_file}  (long-format xlsx 또는 raw Excel)
출력: {output_dir}
모드: raw | fitted | comparison

### raw
- long-format: condition_id / replicate 컬럼 사용
- raw Excel: 헤더에서 substrate 농도 파싱
- 파악 불가 시: condition_0, condition_1, ... 자동 부여
- 조건별 subplot, 리플리케이트별 색상
- X: 시간(min), Y: NADH(mM) 또는 흡광도

### fitted
- long-format + JSON 파라미터 파일
- ic=[A0, B0, P0, Q0] (Q0 반드시 포함)
- sys.path.insert(0, r"C:\Users\Jahyun\Kinetic-modeling-and-optimization\KineticModeling\udh_characterization")
- from GDH_parameter_estimation import GDHKineticModel
- 데이터(점) + fitted curve(선) 오버레이

### comparison
- 2개 이상 모델/파일 (예: irreversible vs reversible)
- 동일 조건 overlay, 2행 배치

## 스타일
COLORS = ['#E69F00','#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7']
조건당 ~2인치, 최대 3열

출력: {basename}_progress_curves.png
파일 경로 + 1줄 요약만 보고, start "" "<경로>" 자동 오픈
```

---

## 모드 E: 표면 플롯 / Heatmap (Haiku 에이전트)

**트리거**: "surface plot", "heatmap", "농도 조합", "2D 최적화", "contour", "yield/titer 표면"

```
Agent(subagent_type="general-purpose", model="haiku", prompt=<아래>)
```

```
입력: {data_file}  (X, Y 조건 컬럼 + Z 결과 컬럼)
출력: {output_dir}
레이아웃: 1x3 | 2x2 | single
Z 컬럼: yield, titer, TTN 등 (복수 가능)

## 처리 방식
from scipy.interpolate import griddata
# 산점 → 정규 그리드 보간 (30×30)
xi = np.linspace(X.min(), X.max(), 30)
yi = np.linspace(Y.min(), Y.max(), 30)
zi = griddata((X, Y), Z, (xi[None,:], yi[:,None]), method='cubic')
# contourf + colorbar (viridis/coolwarm), 실험점 scatter 오버레이
# 레이아웃: 1x3(yield/ee/titer) 또는 2x2(yield/ee/titer/TTN)

## 참조 코드
# C:\Users\Jahyun\Kinetic-modeling-and-optimization\KineticModeling\l_ribose_cascade\model_simulation\surface_plot_utils.py
# combined_surface_plots_2x2(), combined_contour_plots_2x2()

출력: {basename}_surface.png
```

---

## 모드 F: 불확실도 시각화 / Monte Carlo (Haiku 에이전트)

**트리거**: "Monte Carlo", "MC band", "불확실도", "민감도", "워터폴", "confidence band", "Spearman"

```
Agent(subagent_type="general-purpose", model="haiku", prompt=<아래>)
```

```
입력: {mc_results_file}
출력: {output_dir}
유형: timecourse_band | waterfall | distribution | spearman

### timecourse_band
- 중앙값 + 5th/95th percentile 음영 (ax.fill_between)
- species별 subplot (D-Galactose, Galactitol, Tagatose, Formate)

### waterfall (민감도)
- 파라미터별 기여도 가로 바, 양/음 색상 구분, 내림차순 정렬

### distribution
- 최종 결과 분포 히스토그램 + KDE, 중앙값/5th/95th 수직선

### spearman
- 파라미터-결과 Spearman 상관 heatmap

## 참조 코드
# C:\Users\Jahyun\Kinetic-modeling-and-optimization\BayesianOptimization\tagatose_wholecell\utils\viz_utils.py
# C:\Users\Jahyun\biosteam-tagatose\biosteam\plots\plots.py (plot_montecarlo, plot_spearman)

출력: {basename}_uncertainty.png
```

---

## 모드 G: BO 진단 플롯 (Haiku 에이전트)

**트리거**: "BO", "베이지안 최적화", "Pareto", "GP", "acquisition function", "최적화 히스토리"

```
Agent(subagent_type="general-purpose", model="haiku", prompt=<아래>)
```

```
입력: {bo_results_xlsx}
출력: {output_dir}

### Pareto Front
- yield_pct vs titer_mM scatter
- in_silico(점) vs experiment(★) 구분, Pareto 프론트 라인

### BO History
- 라운드별 최대 yield 추이

### GP Landscape (있으면)
- GP 예측 평균 + 불확실도 오버레이 (acquisition function)
- 실험 데이터 포인트 오버레이

## 참조 코드
# C:\Users\Jahyun\Kinetic-modeling-and-optimization\BayesianOptimization\visualization\bo_viz.py
# C:\Users\Jahyun\Kinetic-modeling-and-optimization\BayesianOptimization\visualization\bayesian_optimization_diagram.py

출력: {basename}_bo_diagnostics.png
```

---

## 모드 H: 교육/논문용 다이어그램 (Haiku 에이전트)

**트리거**: "시스템 다이어그램", "반응 네트워크", "흐름도", "cascade 다이어그램", "ODE 모델 그림"

```
Agent(subagent_type="general-purpose", model="haiku", prompt=<아래>)
```

```
입력: {description}  (그릴 시스템/흐름 설명)
출력: {output_dir}

## 처리 방식
- matplotlib patches + FancyArrowPatch 사용
- 효소/반응물 박스: FancyBboxPatch (roundbox)
- 반응 화살표: 실선(→), 억제: 점선(⊣)
- 색상: Okabe-Ito + 회색 배경 박스

## 참조 코드
# C:\Users\Jahyun\Kinetic-modeling-and-optimization\docs\kinetic_modeling_diagram.py  (4-효소 캐스케이드)
# C:\Users\Jahyun\Kinetic-modeling-and-optimization\BayesianOptimization\visualization\bayesian_optimization_diagram.py  (GP 업데이트 흐름)

출력: {basename}_diagram.png (SVG도 선택적 저장)
```

---

## 레포별 참조 파일

| 용도 | 경로 |
|------|------|
| 크로마토그램 exporter | `C:\Users\Jahyun\PeakPicker\src\peakpicker\infrastructure\exporters\plot_exporter.py` |
| Deconvolution 시각화 | `C:\Users\Jahyun\PeakPicker\src\deconvolution_visualizer.py` |
| 한글 폰트 설정 | `C:\Users\Jahyun\PeakPicker\src\solid\utils\plot_utils.py` |
| 키네틱 ODE 모델 | `C:\Users\Jahyun\Kinetic-modeling-and-optimization\KineticModeling\udh_characterization\GDH_parameter_estimation.py` |
| 표면 플롯 유틸 (1190줄) | `C:\Users\Jahyun\Kinetic-modeling-and-optimization\KineticModeling\l_ribose_cascade\model_simulation\surface_plot_utils.py` |
| MC 불확실도 플롯 | `C:\Users\Jahyun\Kinetic-modeling-and-optimization\BayesianOptimization\tagatose_wholecell\utils\viz_utils.py` |
| BO 시각화 | `C:\Users\Jahyun\Kinetic-modeling-and-optimization\BayesianOptimization\visualization\bo_viz.py` |
| Biosteam 플롯 모음 | `C:\Users\Jahyun\biosteam-tagatose\biosteam\plots\plots.py` |
| 기존 분석 스크립트 | `C:\Users\Jahyun\Kinetic-modeling-and-optimization\analyses\` |

---

## 스킬 간 연동

- 결과 → experiment-hub 양상분석(M6)/비교(M7)
- 결과 → manuscript-writer Results/Figure
- 키네틱 raw data 없으면 → kinetic-data-prep 먼저 실행
- HPLC 피크 정량 → PeakPicker `quantify_peaks.py` 먼저 실행

## 주의사항

- ODE fitted 모드: ic에 Q0 반드시 포함 (0으로 두면 fitting 불일치)
- 추세선 외삽은 신뢰도 낮음을 명시
- 한글 레이블 필요 시 `Malgun Gothic` 폰트 사용 (plot_utils.setup_korean_font() 참조)

사용자 요청: $ARGUMENTS
