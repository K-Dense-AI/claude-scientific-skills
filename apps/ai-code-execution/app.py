"""MD Instruction Builder — create .md files to use with Claude web."""

import streamlit as st
import streamlit.components.v1 as components


NO_CODE_TEMPLATE = """\
# Task: [작업 제목 / Task Title]

## Goal / 목표
<!-- 하고 싶은 것을 자연어로 설명 / Describe what you want in plain language -->

예시: "업로드한 엑셀 파일에서 월별 총 매출을 계산하고,
상위 3개 제품을 찾아 막대그래프를 그려주세요."

---

## Input Data / 입력 데이터
<!-- 첨부할 파일 목록 / List the files you will attach -->

| 파일명 / File Name | 형식 / Format | 설명 / Description |
|-------------------|--------------|-------------------|
| sales_2024.xlsx   | Excel        | 월별 매출 (Date, Product, Amount 컬럼) |

---

## Output / 출력 요청

### 텍스트 / Text
- [ ] 월별 매출 합계 출력
- [ ] 상위 3개 제품 출력

### 차트 / Charts
- [ ] 월별 매출 막대그래프 (monthly_sales.png)

### 파일 / Files
- [ ] 결과 요약 엑셀 (summary.xlsx)

---

## 특별 조건 / Special Requirements
<!-- 필터, 단위, 날짜 형식 등 / Filters, units, date format, etc. -->

- 수량 = 0인 행 제외
- 금액 단위: 원(KRW)

---

## 참고 / Notes
<!-- AI가 이해하는 데 도움이 되는 배경 설명 -->

예시: "Price 컬럼은 단가이며, 총액 = Quantity × Price 입니다."
"""

WITH_CODE_TEMPLATE = """\
# Task: [작업 제목 / Task Title]

## 코드 설명 / What This Code Does
<!-- 코드의 목적을 설명 / Describe what the code is supposed to do -->

예시: "CSV 파일에서 중복 행을 제거하고 cleaned_data.csv로 저장하는 스크립트입니다."

---

## Input Data / 입력 데이터

| 파일명 / File Name     | 형식 / Format | 설명 / Description |
|-----------------------|--------------|-------------------|
| customer_data.csv     | CSV          | 고객 데이터        |

---

## 코드 / The Code

```python
import pandas as pd

# 첨부 파일을 읽을 때 파일명만 사용 (경로는 Claude가 처리)
df = pd.read_csv("customer_data.csv")

# TODO: 중복 제거
# df = df.drop_duplicates()

# 결과 저장
df.to_csv("cleaned_data.csv", index=False)
print(f"완료. {len(df)}개 행 처리됨.")
```

---

## 수정 요청 / What Needs to Be Done
<!-- AI에게 무엇을 고치거나 완성할지 알려주기 -->

- "TODO 부분을 pandas 코드로 채워주세요"
- "Age < 0 인 행을 추가로 제거해주세요"

---

## 오류 메시지 / Error Messages (있으면 / if any)

```
ValueError: could not convert string to float: 'N/A'
  File "script.py", line 15, in <module>
```

---

## 기대 출력 / Expected Output

```
완료. 1,248개 행 처리됨.
중복 제거: 34행
```
"""

CLAUDE_GUIDE_HTML = """
<script src="https://cdn.jsdelivr.net/npm/mermaid/dist/mermaid.min.js"></script>
<script>mermaid.initialize({startOnLoad:true, theme:'neutral'});</script>
<div class="mermaid">
flowchart LR
    A([":pencil: 템플릿 편집\\nEdit template"]) --> B([":floppy_disk: .md 다운로드\\nDownload .md"])
    B --> C([":paperclip: Claude 웹에 파일 첨부\\nAttach to Claude web"])
    C --> D([":robot_face: Claude가 실행\\nClaude runs it"])
    D --> E([":white_check_mark: 결과 확인\\nGet results"])
</div>
"""


def sidebar_guide():
    with st.sidebar:
        st.header("사용법 / How to Use")
        st.markdown("""
**1단계 / Step 1**
왼쪽 탭에서 모드 선택 후 내용 편집
Select a mode tab and edit the template

**2단계 / Step 2**
**📥 Download .md** 버튼으로 파일 저장
Click **📥 Download .md** to save the file

**3단계 / Step 3**
[claude.ai](https://claude.ai) 또는 Claude Code 열기
Open [claude.ai](https://claude.ai) or Claude Code

**4단계 / Step 4**
`.md` 파일 + 데이터 파일을 함께 첨부하고 전송
Attach the `.md` file + data files and send
        """)

        st.divider()
        st.subheader("두 가지 모드 / Two Modes")
        st.markdown("""
**NO CODE 모드**
- 하고 싶은 것을 자연어로만 설명
- Claude가 코드 작성 + 실행
- Describe in plain language; Claude writes & runs code

**WITH CODE 모드**
- 기존 코드를 붙여넣고 수정 요청
- Claude가 코드 완성·수정 후 실행
- Paste existing code; Claude fixes & runs it
        """)

        st.divider()
        st.subheader("팁 / Tips")
        st.markdown("""
- 컬럼명을 정확히 입력하세요 (대소문자 구분)
  List column names exactly as in your data
- 출력 형식을 구체적으로 명시 (차트, CSV, 텍스트)
  Be specific about output format
- 데이터 파일은 `.md`와 함께 첨부
  Attach data files together with `.md`
        """)


def main():
    st.set_page_config(
        page_title="MD Instruction Builder",
        page_icon="📝",
        layout="wide",
    )
    sidebar_guide()

    st.title("📝 MD Instruction Builder")
    st.caption("Claude 웹/Code용 분석 지시 파일 생성기 · Build instruction files for Claude web/Code")

    st.markdown("#### 실행 흐름 / Workflow")
    components.html(CLAUDE_GUIDE_HTML, height=130)

    st.divider()

    tab_no_code, tab_with_code = st.tabs(
        ["🗒️ NO CODE — AI가 코드 생성", "💻 WITH CODE — 코드 있음"]
    )

    with tab_no_code:
        st.markdown(
            "**코딩 없이 자연어로 설명** → Claude가 코드 작성 + 실행  \n"
            "_Describe in plain language; Claude writes and runs the code._"
        )
        no_code_text = st.text_area(
            "instruction.md 내용 편집 / Edit instruction.md",
            value=NO_CODE_TEMPLATE,
            height=500,
            key="no_code_editor",
        )
        st.download_button(
            label="📥 Download instruction_no_code.md",
            data=no_code_text.encode("utf-8"),
            file_name="instruction_no_code.md",
            mime="text/markdown",
            key="dl_no_code",
        )

    with tab_with_code:
        st.markdown(
            "**기존 코드를 붙여넣고 수정·완성 요청** → Claude가 실행  \n"
            "_Paste your code and describe what to fix; Claude completes and runs it._"
        )
        with_code_text = st.text_area(
            "instruction.md 내용 편집 / Edit instruction.md",
            value=WITH_CODE_TEMPLATE,
            height=500,
            key="with_code_editor",
        )
        st.download_button(
            label="📥 Download instruction_with_code.md",
            data=with_code_text.encode("utf-8"),
            file_name="instruction_with_code.md",
            mime="text/markdown",
            key="dl_with_code",
        )

    st.divider()
    st.markdown("""
#### Claude 웹에서 사용하는 방법 / How to use with Claude web

1. 위에서 `.md` 파일 다운로드
2. [claude.ai](https://claude.ai) → 새 대화 시작
3. 📎 첨부 버튼 → `.md` 파일 + 데이터 파일(CSV, Excel 등) 함께 첨부
4. 메시지 없이 전송 (또는 "이 파일을 분석해줘" 추가)
5. Claude가 코드를 작성하고 분석 환경에서 실행 후 결과 제공

> **Claude Code 사용 시**: `claude` 명령어 실행 → `.md` 파일 경로를 채팅에 입력
    """)


if __name__ == "__main__":
    main()
