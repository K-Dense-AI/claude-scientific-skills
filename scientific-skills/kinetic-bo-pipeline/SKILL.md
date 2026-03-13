---
name: kinetic-bo-pipeline
description: |
  효소 캐스케이드 Kinetic Fitting -> Global Parameter Optimization -> Bayesian Optimization
  통합 파이프라인 스킬. 개별 효소 피팅(progress-curve-fitting) + 전체 세포 ODE 보정 +
  BO MPSP 최적화(bo-mpsp-optimize) 3단계를 조율하는 메타 파이프라인.

  다음 상황에서 이 스킬을 사용:
  - 개별 효소 kinetics 피팅 후 whole-cell ODE에 통합하는 전체 흐름
  - Global Parameter Optimization 반복 (파라미터 고정/해제 결정 포함)
  - ODE 모델에 새 속도식(rate equation) 추가 후 BO 실행까지 연결
  - 여러 피팅 버전(v4a, v4b, ...) 비교 + 물리적 타당성 검증
  - Asana 태스크에 BO 결과 댓글 업데이트까지 포함한 전체 보고 흐름

  단위 작업만 필요한 경우:
  - 개별 효소 피팅만: progress-curve-fitting 스킬 사용
  - BO 실행만: bo-mpsp-optimize 스킬 사용
  - 원시 데이터 전처리만: kinetic-data-prep 스킬 사용
license: MIT license
metadata:
    skill-author: lab-internal
---

# Kinetic Fitting -> Global Optimization -> Bayesian Optimization Pipeline

효소 캐스케이드 전체 워크플로우: 개별 효소 kinetics -> whole-cell ODE -> BO MPSP 최적화.
2026-03-11 RoGDH 가역 모델 + v4e 피팅 + BO Phase 1 실제 사례 기반 설계.

## 대상 시스템

- **효소 캐스케이드**: RoGDH (galactitol -> D-tagatose), NOX, Formate dehydrogenase
- **ODE 모델**: `CalibratedCascade` (`calibrated_model.py`)
- **최적화 목표**: D-tagatose MPSP ($/kg) 최소화, yield >= 80% 제약

## 워크플로우 전체 흐름

```
[Step 1] 개별 효소 kinetics 피팅 (progress-curve-fitting 스킬)
   -> [Step 2] kinetic_params.py 업데이트
   -> [Step 3] ODE 모델 속도식 추가/수정 (calibrated_model.py)
   -> [Step 4] Global Parameter Optimization (반복)
         +-- 파라미터 고정/해제 결정 (의사결정 포인트 1)
         +-- 수렴 판정 + 물리적 타당성 검증 (의사결정 포인트 2)
   -> [Step 5] BO 실행 (bo-mpsp-optimize 스킬)
   -> [Step 6] 결과 비교 + Asana 보고
```

## 핵심 파일 위치

### 모델 파일

| 역할 | 파일 경로 |
|------|---------|
| 속도식 라이브러리 | `C:\Users\Jahyun\Kinetic-modeling-and-optimization\KineticModeling\l_ribose_cascade\parameter_estimation\rate_equations\RateEquationOrderedBiBi.py` |
| 전체 세포 파라미터 | `C:\Users\Jahyun\Kinetic-modeling-and-optimization\BayesianOptimization\tagatose_wholecell\models\kinetic_params.py` |
| ODE 모델 | `C:\Users\Jahyun\Kinetic-modeling-and-optimization\BayesianOptimization\tagatose_wholecell\models\calibrated_model.py` |
| 개별 효소 피팅 결과 | `C:\Users\Jahyun\Kinetic-modeling-and-optimization\analyses\RoGDH_fitting\results\RoGDH_model_reversible_fit.json` |

### 피팅 스크립트 네이밍 규칙

```
C:\Users\Jahyun\Kinetic-modeling-and-optimization\KineticModeling\l_ribose_cascade\cascade_fitting\fit_cascade_{YYMMDD}_{version}.py
예) fit_cascade_260308_v4e.py
```

### BO 스크립트 및 결과

```
C:\Users\Jahyun\Kinetic-modeling-and-optimization\BayesianOptimization\scripts\bo_phase1_flask_{version}.py
C:\Users\Jahyun\Kinetic-modeling-and-optimization\outputs\bo_phase1_flask_{version}\bo_flask_{version}_e8_final.csv
```

## Step 2: kinetic_params.py 업데이트

개별 효소 피팅 JSON 결과를 `kinetic_params.py`에 추가:

```python
# 예시: RoGDH 가역 모델 파라미터 추가
ROGDH_PARAMS_REV = {
    "kcat_gdh":   33.12,   # s-1  정반응 kcat
    "km_a_gdh":   0.093,   # mM   Km NAD+
    "km_b_gdh":   0.100,   # mM   Km galactitol
    "ki_a_gdh":   0.391,   # mM   Ki NAD+
    "ki_q_gdh":   0.00248, # mM   Ki NADH
    "kcat_r_gdh": 1.251,   # s-1  역반응 kcat
    "km_p_gdh":   0.0498,  # mM   Km D-tagatose
    "km_q_gdh":   0.00248, # mM   Km NADH (reverse)
    "ki_p_gdh":   0.0183,  # mM   Ki D-tagatose
    "keq_gdh":    0.0425,  # -    평형 상수
}
```

## Step 4: Global Parameter Optimization 반복 패턴

### 버전별 파라미터 구성 (v4 계열 실제 사례)

| 버전 | Free 파라미터 수 | km_b_gdh | futile | kLa | 물리적 타당성 |
|------|----------------|----------|--------|-----|------------|
| v4 | 6 | free (->0.009mM, 비물리적) | fixed=0 | fixed | 주의 |
| v4a | 5 | fixed=0.100 | fixed=0 | fixed | 실패 (alpha 폭주) |
| v4b | 6 | fixed=0.100 | fixed=0 | free (->0.0001) | 실패 |
| v4c | 6 | fixed=0.100 | free | fixed | 주의 |
| v4d | 7 | fixed=0.100 | free | free | OK |
| v4e | 8 | free (->0.0154) | free | free | 최우수 |

### 일반화된 반복 패턴

```
1. 베이스 버전 피팅 실행
2. 비물리적 파라미터 발견
   |
   +-- 해당 파라미터 고정 -> 다른 파라미터 폭주?
   |   +-- YES -> 누락된 메커니즘 파악 -> 새 free 파라미터 추가
   |   +-- NO  -> 고정값 유지하고 다음 단계
   |
   +-- 자유도 추가 (순차적으로)
       +-- warm-start로 재피팅 -> 수렴 확인
```

## 의사결정 포인트: 파라미터 물리성 판단 기준

| 파라미터 | 물리적 범위 | 비물리적 신호 |
|---------|-----------|------------|
| km_b_gdh (mM) | 0.01 ~ 10 | < 0.01 또는 > 50 |
| alpha_pntab | 1 ~ 5 | > 10 |
| qo2_basal (mmol O2/gDCW/h) | 0.05 ~ 0.30 | > 1.0 |
| effectiveness_factor | 0.15 ~ 1.0 | > 1.0 또는 < 0.1 |
| kLa_0rpm (h-1) | 0.005 ~ 0.5 | 경계값(0.0001 또는 상한) |
| vmax_futile_nadph (mM/s) | 0 ~ 0.01 | > 0.05 |

### 핵심 패턴

| 패턴 | 설명 | 해결책 |
|------|------|--------|
| 보상 효과 | 파라미터 A 고정 시 B가 비물리적으로 보상 | B 범위 내 새 메커니즘 파라미터 추가 |
| 경계값 수렴 | 파라미터가 상한/하한에 붙음 | 탐색 범위 재설정 또는 파라미터 해제 |
| warm-start | 이전 버전 파라미터를 초기값으로 사용 | 안정적 수렴, 새 자유도 추가 시 필수 |
| futile cycle | aerobic 조건 설명 시 NADPH oxidase 경로 필수 | vmax_futile_nadph를 free로 설정 |

## Step 5: BO 탐색 공간 (5D, Phase 1 Flask)

```python
SEARCH_SPACE = {
    "cell_conc_gL":     (5, 50),    # g DCW/L
    "dgal_mM":          (50, 200),  # mM
    "reaction_time_h":  (6, 72),    # h
    "formate_mM":       (200, 600), # mM
    "rpm":              (0, 250),   # rpm
}
```

**알고리즘**: Sobol 초기화 60점 -> 상위 20점 TEA -> GP-LCB 30회 추가 (총 50회)

**제약**: Feasible yield >= 80%, 목적함수 MPSP 최소화

## 입력/출력 명세

### 입력

| 입력 | 형식 | 설명 |
|------|------|------|
| 개별 효소 피팅 결과 | JSON | kcat, Km, Ki, keq 포함 |
| 기존 calibrated 파라미터 | JSON | warm-start용 베이스 버전 |
| 실험 데이터 | CSV/dict | RPM별 시계열 (농도 mM vs 시간 h) |
| BO 탐색 공간 정의 | dict | 5D bounds |
| Asana 태스크 GID | str | 결과 댓글 업데이트용 (선택) |

### 출력

| 출력 | 형식 | 설명 |
|------|------|------|
| 최적 피팅 파라미터 | JSON | `calibrated_v5_{date}_{version}_params.json` |
| BO 결과 | CSV | `bo_flask_{version}_e8_final.csv` |
| Feasible 최적 조건 | 텍스트/테이블 | yield >= 80% 중 MPSP 최소 조건 |
| 보고서 | Markdown | `{version}_result.md` |
| Asana 댓글 | Asana API | BO 결과 요약 + 실험 조건 비교 |

## 사용 예시

### 전체 파이프라인 (팀 오케스트레이터)

```python
mcp__team-orchestrator__start_project(
    description="""
    RoGDH 가역 모델 통합 + Global Fit v4e + BO Phase 1
    - Step 1: progress-curve-fitting 스킬로 RoGDH 가역 모델 파라미터 추정
    - Step 2-3: kinetic_params.py + calibrated_model.py 업데이트
    - Step 4: Global fit 반복 (v4d -> v4e warm-start)
    - Step 5: bo-mpsp-optimize 스킬로 BO 실행
    - Step 6: Asana 태스크 {GID}에 결과 댓글
    """
)
```

### Global Fit 스크립트 빠른 실행

```bash
# 새 버전 피팅 실행 (warm-start from v4d)
python C:\Users\Jahyun\Kinetic-modeling-and-optimization\KineticModeling\l_ribose_cascade\cascade_fitting\fit_cascade_260308_v4e.py

# BO 실행 (v4e 파라미터 기반)
python C:\Users\Jahyun\Kinetic-modeling-and-optimization\BayesianOptimization\scripts\bo_phase1_flask_v4e.py --n-initial 20 --n-calls 50 --ecoli-price 8
```

## 물리성 검증 체크리스트 (피팅 후 필수)

- [ ] 모든 파라미터가 물리적 범위 내 수렴 (위 기준표 참조)
- [ ] km_b_gdh < 0.01 mM이면 -> 실험적 재측정 또는 고정 검토
- [ ] alpha_pntab > 5이면 -> futile NADPH cycle 또는 새 메커니즘 추가
- [ ] kLa_0rpm이 경계값이면 -> DO 데이터 확보 또는 고정
- [ ] qo2_basal > 1.0 mmol/gDCW/h이면 -> 생물학적으로 비합리적

## 통합 스킬 체인

| 단계 | 사용 스킬 |
|------|---------|
| Step 0: 원시 데이터 전처리 | kinetic-data-prep |
| Step 1: 개별 효소 피팅 | progress-curve-fitting |
| Step 5: BO 실행 | bo-mpsp-optimize |
| 시각화 | lab-viz (모드 G: BO 진단) |
| 보고서 작성 | experiment-hub (모드 6: 패턴 분석) |

## 주의사항

- **파라미터 고정/해제 결정은 사람이 판단** — 생물학적 해석 필요
- ODE 속도식 추가 후 반드시 `sanity_check_ode()` 실행 (RPM=0/200 극단 조건)
- warm-start 없이 8-param 피팅 시 비수렴 위험 -> 반드시 lower-param 버전에서 시작
- BO 결과 채택 전 실험 타당성 검토 (BO 최적 조건이 실험 가능 범위인지 확인)
- Asana 댓글은 `html_text` 파라미터 사용 (html_notes 금지)
- fit_loss 개선이 5% 미만이면 버전 추가 대신 BO 진행 권장

## 실제 성과 (2026-03-11, v4e 기준)

| 파라미터 | v4e 값 | 물리적 타당성 |
|---------|-------|------------|
| km_b_gdh | 0.0154 mM | 허용 (하한 근접, 검증 권장) |
| alpha_pntab | 1.70 | OK |
| effectiveness_factor | 0.494 | OK |
| qo2_basal | 0.200 mmol O2/gDCW/h | OK |
| vmax_futile_nadph | 0.00197 mM/s | OK |
| fit_loss | 0.00492 | 4개 버전 중 최저 |

**BO 최적 조건 (Feasible, 최저 MPSP)**:
```
cell_conc_gL=46.5, dgal_mM=130, reaction_time_h=48, formate_mM=493, rpm=35
MPSP=$41.78/kg, yield=83.6%, titer=108 mM D-tagatose
```
