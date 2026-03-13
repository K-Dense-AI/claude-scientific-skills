---
name: lab-eda
description: >
  Lab-specific EDA skill (Haiku-powered). Complements the global
  exploratory-data-analysis skill with domain-specific analysis for
  Jahyun's data types: kinetic time-series, HPLC chromatography,
  protein sequence/clustering, and general tabular lab data.
  Generates structured markdown reports + quick summary figures.
license: MIT license
metadata:
    skill-author: lab-internal
---

# Lab-EDA (실험실 특화 탐색적 데이터 분석)

## lab-eda vs exploratory-data-analysis 분기 기준

> **핵심 판단**: 데이터가 Jahyun 실험실에서 직접 생성된 것인가?

| 판단 기준 | 사용 스킬 | 예시 |
|----------|----------|------|
| 실험실 생성 데이터 (키네틱/HPLC/클러스터/겔/감도분석) | **lab-eda** (이 스킬) | kinetic xlsx, HPLC csv, SDS-PAGE 이미지, Sobol 분석 |
| 일반 실험 Excel/CSV (조건-결과 테이블) | **lab-eda** (이 스킬) | 실험 결과 정리 엑셀, 배치 비교 CSV |
| 비표준 과학 포맷 (200+ 확장자) | **exploratory-data-analysis** (전역) | .pdb, .fasta, .fastq, .mzML, .nd2, .czi, .hdf5, .sdf |
| EDA 결과 → 시각화 필요 | **lab-viz** 연계 | EDA 후 그래프 요청 |

**간단 규칙**: 파일이 `.xlsx`, `.csv`, `.tsv`, 이미지(겔)이고 실험실 데이터면 → **lab-eda**.
파일 확장자가 과학 도메인 전용 포맷이면 → **exploratory-data-analysis**.

실험실 데이터 유형별 빠른 EDA. Haiku 에이전트로 저비용 처리.

## 모드 판별

| 진입 조건 | 모드 |
|----------|------|
| kinetic / NADH / progress curve / xlsx long-format | **A: 키네틱 데이터 EDA** (Haiku) |
| HPLC / 크로마토 / peak / Chemstation / Agilent | **B: HPLC 데이터 EDA** (Haiku) |
| sequence / clustering / PDB / foldseek / cluster | **C: 시퀀스/클러스터 EDA** (Haiku) |
| 일반 Excel/CSV / 실험 결과 테이블 | **D: 범용 테이블 EDA** (Haiku) |
| 여러 파일 일괄 / batch 비교 | **E: 배치 비교 EDA** (Haiku) |
| 겔 이미지 / SDS-PAGE / 웨스턴 / 밴드 | **F: 겔 이미지 정량 EDA** (Haiku) |
| 감도 분석 / Sobol / 파라미터 영향도 | **G: 감도 분석 EDA** (Haiku) |

---

## 공통 출력 형식

모든 모드의 기본 출력:
```
{output_dir}/{basename}_eda_report.md   ← 구조화된 마크다운 리포트
{output_dir}/{basename}_eda_summary.png ← 핵심 요약 figure (선택)
```
CEO에게: 파일 경로 + 3줄 요약만 보고.

---

## 모드 A: 키네틱 데이터 EDA (Haiku 에이전트)

**트리거**: "kinetic data 확인", "데이터 품질 확인", "노이즈 분석", "리플리케이트 변동성"

```
Agent(subagent_type="general-purpose", model="haiku", prompt=<아래>)
```

```
입력: {input_file}  (long-format xlsx 또는 raw Excel)
출력: {output_dir}

## 분석 항목

### 1. 구조 파악
- 시트명, 컬럼 목록, shape
- 조건 수 (condition_id 유니크) / 리플리케이트 수 / 시간 포인트 수
- 기질 농도 조합 테이블 (A0, B0, Q0 분포)

### 2. 데이터 품질
- 결측값 (NaN) 비율
- 음수값 (NADH < 0) 발생 여부 및 위치
- 시간축 단조증가 여부 확인
- t=0 보정 여부 (Q0 ≠ 0 이면 미보정 데이터)

### 3. 리플리케이트 변동성
- 조건별 replicate 간 std / mean (CV%) 계산
- CV > 20% 인 조건 플래그
- 이상치 replicate 탐지 (z-score > 2.5)

### 4. 반응 특성
- 각 조건의 최대 NADH (plateau 추정치)
- plateau 도달 여부 (마지막 5 포인트 std < 0.01 mM)
- 초기 반응 속도 (t=0~5min 기울기)

### 5. 시각화 (빠른 overview)
- 조건별 mean ± std 시계열 (1개 figure, 전체 조건 오버레이)
- 조건 × 지표 히트맵 (CV%, plateau %)

출력:
- {basename}_eda_report.md
- {basename}_eda_overview.png
파일 경로 + 3줄 요약만 보고

## 참조
# C:\Users\Jahyun\Kinetic-modeling-and-optimization\KineticModeling\gdh_characterization\analysis\GDH_noise_analysis.py (noise 분석 패턴)
# C:\Users\Jahyun\Kinetic-modeling-and-optimization\KineticModeling\gdh_characterization\analysis\rogdh_input_prep_and_fit.py (구조 파악)
```

---

## 모드 B: HPLC 데이터 EDA (Haiku 에이전트)

**트리거**: "HPLC 데이터 확인", "피크 품질", "크로마토 탐색", "batch HPLC 요약"

```
Agent(subagent_type="general-purpose", model="haiku", prompt=<아래>)
```

```
입력: {data_dir_or_file}  (Chemstation 내보내기 CSV 또는 디렉토리)
출력: {output_dir}

## 분석 항목

### 1. 파일 구조 파악
- 샘플 수, 검출기 채널, 분석 시간 범위
- 보유시간 범위 및 피크 수 분포

### 2. 크로마토그램 품질
- 베이스라인 안정성 (초기/말기 신호 표준편차)
- S/N ratio 추정 (피크 높이 / 베이스라인 노이즈)
- 피크 대칭성 (asymmetry factor = 1.0 기준)

### 3. 피크 요약 통계
- 각 화합물 RT 분포 (mean ± std)
- 면적/높이 CV% (배치간 재현성)
- 결측 피크 (미검출 샘플 수)

### 4. 이상치 탐지
- 면적 z-score > 3.0 플래그
- RT 이동 > 0.1 min 샘플 플래그
- 낮은 S/N (< 3) 피크 플래그

### 5. 정량 요약
- 각 화합물 농도 분포 (히스토그램)
- 조건별 mean ± SD 테이블

출력:
- {basename}_hplc_eda_report.md
- {basename}_hplc_eda_summary.png
파일 경로 + 3줄 요약만 보고

## 참조
# C:\Users\Jahyun\PeakPicker\scripts\hplc_analyzer_enhanced.py
# C:\Users\Jahyun\PeakPicker\src\peakpicker\ (피크 감지 모듈)
```

---

## 모드 C: 시퀀스/클러스터 EDA (Haiku 에이전트)

**트리거**: "클러스터 분석 확인", "시퀀스 분포", "noise 구조물", "cluster 품질"

```
Agent(subagent_type="general-purpose", model="haiku", prompt=<아래>)
```

```
입력: {result_file}  (클러스터링 결과 CSV/Excel 또는 PDB 디렉토리)
출력: {output_dir}

## 분석 항목

### 1. 클러스터 분포
- 전체 구조물 수 / 클러스터 수
- 클러스터 크기 분포 (히스토그램)
- noise 비율 (label=-1 또는 별도 noise 컬럼)
- 싱글턴 클러스터 수

### 2. 시퀀스 길이 분포 (있으면)
- mean ± std, Q1-Q3
- 이상치 (길이 기준 상하위 5%)

### 3. 클러스터 품질 지표
- 클러스터 내 평균 거리 (있으면)
- 대표 시퀀스 대비 분산

### 4. 이상 구조물 탐지
- noise 구조물 목록
- 클러스터 크기 = 1 목록

출력:
- {basename}_cluster_eda_report.md
파일 경로 + 3줄 요약만 보고

## 참조
# C:\Users\Jahyun\UDH_Clustering\analyze_noise_detailed.py
# C:\Users\Jahyun\UDH_Clustering\plot_cluster_analysis.py
```

---

## 모드 D: 범용 테이블 EDA (Haiku 에이전트)

**트리거**: "데이터 확인해줘", "엑셀 탐색", "CSV EDA", "컬럼 분포"

```
Agent(subagent_type="general-purpose", model="haiku", prompt=<아래>)
```

```
입력: {file}  (CSV, Excel, TSV)
출력: {output_dir}

## 분석 항목
1. 기본 구조: shape, dtypes, 컬럼명
2. 결측값: 컬럼별 null 수 / 비율
3. 요약 통계: describe() (수치형), value_counts() (범주형)
4. 상관관계: 수치 컬럼 corr() heatmap
5. 이상치: IQR 방법으로 각 컬럼 이상치 수
6. 중복행: 완전 중복 / 키 컬럼 중복

출력:
- {basename}_eda_report.md (마크다운 테이블 형식)
- {basename}_eda_corr.png (상관 heatmap, 컬럼 5개+ 일 때)
파일 경로 + 3줄 요약만 보고
```

---

## 모드 E: 배치 비교 EDA (Haiku 에이전트)

**트리거**: "여러 파일 비교", "배치 간 비교", "실험 간 일관성", "batch 요약"

```
Agent(subagent_type="general-purpose", model="haiku", prompt=<아래>)
```

```
입력: {file_list}  (같은 포맷 여러 파일)
출력: {output_dir}

## 분석 항목
1. 파일별 기본 통계 요약 테이블 (행/열 수, 키 지표 mean/std)
2. 배치 간 주요 지표 CV%
3. 이상 배치 플래그 (전체 평균 ± 2σ 벗어난 배치)
4. 공통 컬럼 확인 (파일 간 구조 일관성)

출력:
- {basename}_batch_comparison.md
- {basename}_batch_comparison.png
파일 경로 + 3줄 요약만 보고
```

---

---

## 모드 F: 겔 이미지 정량 EDA (Haiku 에이전트)

**트리거**: "겔 이미지", "SDS-PAGE", "웨스턴 블롯", "밴드 정량", "발현 수준"

```
Agent(subagent_type="general-purpose", model="haiku", prompt=<아래>)
```

```
입력: {image_file}  (PNG, TIFF, JPG — 겔 이미지)
출력: {output_dir}

## 처리 방식 (ImageJ 스타일)
from PIL import Image; import numpy as np; from scipy.signal import find_peaks

1. 레인 감지: 수직 강도 프로필 → 피크 위치 = 레인 중심
2. 강도 프로필: 레인별 가로 합 (또는 ROI 평균)
3. 밴드 감지: find_peaks(profile, prominence=threshold)
4. 적분: np.trapz(band_region)
5. 상대 발현도 = band_intensity / total_intensity (또는 loading control)

## 출력
- 레인 수, 밴드 수 요약
- 밴드별 상대 강도 테이블 (Excel)
- 레인별 강도 프로필 오버레이 figure

## 참조
# C:\Users\Jahyun\Kinetic-modeling-and-optimization\GelAnalysis\gel_analysis.py
```

---

## 모드 G: 감도 분석 EDA (Haiku 에이전트)

**트리거**: "감도 분석", "Sobol", "파라미터 영향도", "sensitivity", "토네이도 플롯"

```
Agent(subagent_type="general-purpose", model="haiku", prompt=<아래>)
```

```
입력: {params_file}  (파라미터 정의 + ODE 모델 또는 결과 매트릭스)
출력: {output_dir}

## 처리 방식
from SALib.sample import saltelli
from SALib.analyze import sobol

1. 파라미터 공간 정의 (bounds, names)
2. Saltelli 샘플링 (N×(2D+2) 점)
3. 모델 평가 (ODE 또는 대리 모델)
4. Sobol 지수 계산: S1 (1차), ST (전체 효과)
5. 토네이도 플롯: S1 기준 정렬, 오차 막대 포함

## 출력
- 파라미터별 S1/ST 표
- 토네이도 플롯 figure
- 상위 3 파라미터 해석 요약

## 참조
# C:\Users\Jahyun\Kinetic-modeling-and-optimization\BayesianOptimization\tagatose_wholecell\analysis\sensitivity_analysis.py
```

---

## 레포별 참조 파일

| 용도 | 경로 |
|------|------|
| 키네틱 noise 분석 | `C:\Users\Jahyun\Kinetic-modeling-and-optimization\KineticModeling\gdh_characterization\analysis\GDH_noise_analysis.py` |
| 키네틱 데이터 구조 파악 | `C:\Users\Jahyun\Kinetic-modeling-and-optimization\KineticModeling\gdh_characterization\analysis\rogdh_input_prep_and_fit.py` |
| HPLC 분석 파이프라인 | `C:\Users\Jahyun\PeakPicker\scripts\hplc_analyzer_enhanced.py` |
| 클러스터 noise 분석 | `C:\Users\Jahyun\UDH_Clustering\analyze_noise_detailed.py` |
| 클러스터 플롯 | `C:\Users\Jahyun\UDH_Clustering\plot_cluster_analysis.py` |
| 겔 이미지 정량화 | `C:\Users\Jahyun\Kinetic-modeling-and-optimization\GelAnalysis\gel_analysis.py` |
| 감도 분석 (Sobol) | `C:\Users\Jahyun\Kinetic-modeling-and-optimization\BayesianOptimization\tagatose_wholecell\analysis\sensitivity_analysis.py` |
| GP 시계열 불확실성 | `C:\Users\Jahyun\Kinetic-modeling-and-optimization\BayesianOptimization\tagatose_wholecell\analysis\260226_ipa_freeze_gp_analysis.py` |

---

## 스킬 간 연동

```
lab-eda  ──► lab-viz      (EDA 결과 → 시각화)
         ──► experiment-hub (이상치/패턴 발견 → 양상분석 M6)
         ──► kinetic-data-prep (raw 데이터면 → 전처리 먼저)
         ──► exploratory-data-analysis (비표준 포맷이면 → 전역 스킬)
```

## 주의사항

- 대용량 파일(>10MB xlsx): openpyxl read_only 모드 사용
- HPLC CSV: Chemstation 내보내기 포맷 (헤더 행 구조 가변)
- 결과 저장 경로: 해당 레포의 `analyses/` 폴더 또는 데이터와 동일 디렉토리 (MEMORY.md 분석 스크립트 저장 위치 규칙 참조)

사용자 요청: $ARGUMENTS
