---
name: code-validator
description: 코드 검수 및 디버깅 통합 도구. (1) 레퍼런스 재현 테스트, (2) pytest 실행 및 출력 검증, (3) PyCharm 디버거 기반 에러 분석, (4) 과학적 assert 패턴 적용, (5) 수정→재테스트 루프 (최대 3회). P3 전담.
---

# Code Validator (코드 검수) — P3 전담

테스트 실행, 에러 분석, 수정, 재테스트를 하나의 루프로 관리하는 스킬.

## 핵심 원칙

**레퍼런스 먼저, 내 코드 다음.** 원본이 깨지지 않았음을 확인한 후에 내 수정분을 테스트한다.

## 실행 방식

**팀 생성 금지. Bash 직접 실행으로 처리한다.**

Step 1(레퍼런스 테스트)과 Step 2(pytest)가 독립적이면 Bash 2개를 **한 메시지에서 병렬 호출**한다.

## 워크플로우

```
P2 산출물 입력
    ↓
Step 1 + Step 2: Bash 병렬 실행
    ↓
Step 3: 에러 시 디버그 루프 (최대 3회)
    ↓
Step 4: 결과 기록 → P4 입력
```

---

## Step 1: 레퍼런스 재현 테스트

판정: ✅ 일치 → Step 2 / ⚠️ 미세 차이 → 기록 후 Step 2 / ❌ 불일치 → P2 회귀

---

## Step 2: 과학적 assert 패턴

```python
assert 0.0 <= conversion_rate <= 1.0       # 범위 검증
assert result == pytest.approx(expected, rel=0.05)  # 근사값
assert np.isfinite(array).all()             # NaN/Inf 없음
assert abs(mass_in - mass_out) / mass_in < 0.01  # 물질 수지
```

### Silent Failure 패턴 검출 (필수)

`dict.get()` 또는 None 반환이 데이터 누락을 조용히 숨기는 경우를 반드시 검출한다.

```python
# 위험: 없는 키 → 0 반환, 에러·경고 없음 → 차트/계산에서 버그 숨겨짐
def _unit_capex(info_dict, key):
    u = info_dict.get(key)
    if u is None:
        return 0.0  # ← silent failure

# 안전: 없는 키 → 경고 출력
import warnings
def _unit_capex(info_dict, key):
    u = info_dict.get(key)
    if u is None:
        warnings.warn(f"Unit key '{key}' not found in info dict — returning 0", stacklevel=2)
        return 0.0
```

또는 루프 진입 전 일괄 검증으로 조기 실패:

```python
missing = [k for k in unit_keys if k not in info_dict]
assert not missing, f"Missing unit keys: {missing}"
```

**적용 대상**: 플롯 스크립트, TEA 결과 추출, 비용 분해 함수 등 `info_dict.get()` 사용 부분 전체.

---

## Step 3: 디버그 루프 (최대 3회)

| 분류 | 설명 | 대응 |
|------|------|------|
| A. 코드 버그 | 타입 오류, 로직 실수 | P3 내 수정 |
| B. 데이터 문제 | 형식 불일치, 누락값 | P3 내 수정 |
| C. 의존성 문제 | 버전 충돌, API 변경 | P3 내 수정 |
| D. 환경 문제 | 경로, 권한, 메모리 | P3 내 수정 |
| E. 설계 문제 | 알고리즘 한계, 구조 결함 | **P2 회귀** |

3회 초과 시 미해결 에러를 로그에 기록하고 P5로 전달.

---

## Step 4: 결과 기록

```
레퍼런스 재현: ✅/⚠️/❌
전체 테스트: {passed}/{total}
디버그 루프: {N}회
```

커밋: `[P3] fix: {에러 설명} (수정 N/3)`

---

## 사전 확인 (모르면 물어보기)

작업 시작 전 다음 항목이 불명확하면 **추측하지 말고 사용자에게 질문한다**:
- 테스트 대상 파일/모듈이 어디인지
- 레퍼런스 기준값이 무엇인지
- 허용 오차 범위 (rel=0.05? 0.01?)
- 특정 테스트만 실행할지 전체 실행할지

## 완료 체크리스트 (끝나면 대조 보고)

작업 완료 시 사용자의 원래 요청을 항목별로 대조하여 보고한다:
```
✅ 레퍼런스 재현 테스트: 통과
✅ pytest 전체 실행: 70/70
✅ 에러 수정: 0건 (에러 없음)
⬜ [누락 항목 있으면 여기에 명시]
```

## 종료 조건

정상 → P4: 실행 결과 + 에러 로그
미해결 → P5: 미해결 에러 로그 + "분석 필요" 태그

## 스킬 간 연동

- code-implementer (P2): 코드 + 테스트를 입력으로 받음
- reference-surveyor (P1.5): reference_baseline을 재현 기준으로 사용
- research-discussion (P4): 실행 결과를 해석·비교에 전달
