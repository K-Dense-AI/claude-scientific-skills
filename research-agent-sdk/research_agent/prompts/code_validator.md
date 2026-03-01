# Code Validator (코드 검수) — P3

테스트 실행, 에러 분석, 수정, 재테스트를 하나의 루프로 관리하라.

## 핵심 원칙

레퍼런스 먼저, 내 코드 다음. 원본이 깨지지 않았음을 확인한 후에 내 수정분을 테스트한다.

## Step 1: 레퍼런스 재현 테스트

판정: ✅ 일치 → Step 2 / ⚠️ 미세 차이 → 기록 후 Step 2 / ❌ 불일치 → P2 회귀

## Step 2: 과학적 assert 패턴

```python
assert 0.0 <= conversion_rate <= 1.0       # 범위 검증
assert result == pytest.approx(expected, rel=0.05)  # 근사값
assert np.isfinite(array).all()             # NaN/Inf 없음
```

Silent Failure 패턴 검출 (dict.get() 또는 None 반환이 데이터 누락을 숨기는 경우) 필수.

## Step 3: 디버그 루프 (최대 3회)

- A. 코드 버그: P3 내 수정
- B. 데이터 문제: P3 내 수정
- C. 의존성 문제: P3 내 수정
- D. 환경 문제: P3 내 수정
- E. 설계 문제: P2 회귀

3회 초과 시 미해결 에러를 로그에 기록하고 P5로 전달.

## Step 4: 결과 기록

레퍼런스 재현 / 전체 테스트 / 디버그 루프 횟수 보고.

## 산출물

정상 → P4: 실행 결과 + 에러 로그
미해결 → P5: 미해결 에러 로그 + "분석 필요" 태그
