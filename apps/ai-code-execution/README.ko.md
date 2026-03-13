# MD Instruction Builder

Claude 웹/Code에 전달할 분석 지시 파일(`.md`)을 편집하고 다운로드하는 도구입니다.

---

## 어떻게 동작하나요?

```
1. 템플릿 편집  ──▶  2. .md 다운로드  ──▶  3. Claude 웹에 첨부  ──▶  4. Claude가 실행
```

코딩 없이 자연어로 작업을 설명하면 Claude가 코드를 작성하고 실행해 결과를 제공합니다.

---

## 빠른 시작

1. [앱 URL](https://claude-scientific-skills-7w7pxx6re8wlgwpbmrthdt.streamlit.app/) 접속
2. 모드 선택 (NO CODE / WITH CODE)
3. 템플릿 내용 편집
4. **📥 Download .md** 클릭
5. [claude.ai](https://claude.ai) → 📎 첨부 → `.md` + 데이터 파일 전송

---

## 두 가지 모드

### NO CODE 모드
코드 없이 자연어로 하고 싶은 것을 설명합니다.

```markdown
# Task: 월별 매출 분석

## Goal
sales_2024.xlsx 파일의 월별 총 매출을 계산하고 막대그래프를 그려주세요.

## Input Data
| 파일명 | 형식 | 설명 |
|--------|------|------|
| sales_2024.xlsx | Excel | 월별 매출 데이터 |

## Output
- 텍스트: 월별 합계
- 이미지: monthly_sales.png
```

### WITH CODE 모드
기존 코드를 붙여넣고 수정·완성을 요청합니다.

```markdown
# Task: 데이터 정제

## The Code
```python
import pandas as pd
df = pd.read_csv("raw_data.csv")
# TODO: 중복 제거
df.to_csv("cleaned_data.csv", index=False)
```

## What Needs to Be Done
- TODO 부분을 pandas 코드로 채워주세요
```

---

## Claude 웹에서 사용

1. [claude.ai](https://claude.ai) → 새 대화
2. 📎 → `.md` 파일 + 데이터 파일 첨부
3. 전송 (메시지 없이도 가능)
4. Claude가 분석 후 결과 제공

---

## 파일 구조

```
apps/ai-code-execution/
├── app.py           # Streamlit 앱 (템플릿 편집 + 다운로드)
├── requirements.txt # streamlit만 필요
└── README.ko.md     # 이 파일
```
