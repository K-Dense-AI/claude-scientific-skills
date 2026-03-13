"""Streamlit UI for AI-Powered Code Execution System."""

import streamlit as st
import streamlit.components.v1 as components
from config import AI_BACKEND, MAX_RETRIES, MAX_MD_TOKENS
from ai_codegen import generate_code
from executor import run_in_sandbox


MERMAID_FLOW = """
flowchart TD
    A([":page_facing_up: Upload instruction.md\\n지시 파일 업로드"]) --> B{Has Python code?\\n파이썬 코드 있음?}
    B -- YES / 있음 --> C["Extract code\\n코드 추출"]
    B -- NO / 없음 --> D["AI generates code\\nAI 코드 생성 (Groq)"]
    C --> E[":rocket: Run in E2B Sandbox\\nE2B 샌드박스 실행"]
    D --> E
    E --> F{Success?\\n성공?}
    F -- YES / 성공 --> G([":white_check_mark: Show results\\n결과 출력"])
    F -- NO / 실패 --> H{Retry < 3\\n재시도 가능?}
    H -- YES --> D
    H -- NO --> I([":x: Show error\\n오류 표시"])
"""


def render_mermaid(diagram: str):
    html = f"""
    <script src="https://cdn.jsdelivr.net/npm/mermaid/dist/mermaid.min.js"></script>
    <script>mermaid.initialize({{startOnLoad:true, theme:'neutral'}});</script>
    <div class="mermaid">{diagram}</div>
    """
    components.html(html, height=420)


def detect_mode(md_text: str) -> str:
    if "```python" in md_text:
        return "with_code"
    return "no_code"


def estimate_tokens(text: str) -> int:
    return len(text) // 4


def sidebar_guide():
    with st.sidebar:
        st.header("How to use / 사용법")

        st.markdown("""
**Step 1 / 1단계**
Write an instruction `.md` file describing your task.
작업 내용을 `.md` 파일에 작성하세요.

**Step 2 / 2단계**
Upload the `.md` file (+ optional data files).
`.md` 파일과 데이터 파일(선택)을 업로드하세요.

**Step 3 / 3단계**
Click **▶ Run / 실행** and wait for results.
실행 버튼을 누르고 결과를 기다리세요.
        """)

        st.divider()
        st.subheader("Two modes / 두 가지 모드")
        st.markdown("""
**NO CODE mode / 코드 없음 모드**
- Describe what you want in plain language
- AI writes the Python code automatically
- 하고 싶은 것을 자연어로 설명
- AI가 파이썬 코드를 자동 생성

**WITH CODE mode / 코드 있음 모드**
- Include a ` ```python ``` ` block in your `.md`
- AI fixes/completes it if needed
- `.md` 파일에 파이썬 코드 블록 포함
- AI가 필요 시 수정·보완
        """)

        st.divider()
        st.subheader("Tips / 팁")
        st.markdown("""
- Be specific about **output format** (chart, CSV, text...)
  출력 형식을 구체적으로 명시하세요
- List **column names** exactly as in your data
  데이터의 컬럼명을 정확히 작성하세요
- Output files are saved to `/home/user/output/`
  출력 파일 경로: `/home/user/output/`
        """)

        st.divider()
        st.subheader("Flow / 실행 흐름")
        render_mermaid(MERMAID_FLOW)


def main():
    st.set_page_config(page_title="AI Code Execution", page_icon="▶", layout="wide")
    sidebar_guide()

    st.title("AI-Powered Code Execution / AI 코드 실행기")
    st.caption(f"Backend: **{AI_BACKEND}** | Max retries / 최대 재시도: {MAX_RETRIES}")

    # --- File uploads ---
    col1, col2 = st.columns(2)
    with col1:
        md_file = st.file_uploader(
            "Instruction file / 지시 파일 (.md)",
            type=["md", "txt"],
            key="md_upload",
            help="Upload a .md file describing the task. / 작업 내용을 담은 .md 파일을 업로드하세요.",
        )
    with col2:
        data_files = st.file_uploader(
            "Data files / 데이터 파일 (optional / 선택)",
            accept_multiple_files=True,
            key="data_upload",
            help="CSV, Excel, JSON, images, etc. / CSV, 엑셀, JSON, 이미지 등",
        )

    if not md_file:
        st.info("Upload an instruction .md file to get started. / 지시 파일(.md)을 업로드하면 시작됩니다.")
        return

    # --- Read instruction ---
    md_text = md_file.read().decode("utf-8", errors="replace")
    mode = detect_mode(md_text)

    token_est = estimate_tokens(md_text)
    if token_est > MAX_MD_TOKENS:
        st.warning(f"Instruction file is very large (~{token_est:,} tokens). Processing may be slow. / 파일이 매우 큽니다. 처리가 느릴 수 있습니다.")

    mode_label = "WITH CODE / 코드 있음" if mode == "with_code" else "NO CODE / 코드 없음 (AI가 생성)"
    st.markdown(f"**Mode / 모드:** `{mode_label}`")

    with st.expander("Preview instruction / 지시 내용 미리보기", expanded=False):
        st.markdown(md_text)

    # --- Prepare data files ---
    data_file_list: list[tuple[str, bytes]] = []
    data_filenames: list[str] = []
    for f in data_files:
        content = f.read()
        data_file_list.append((f.name, content))
        data_filenames.append(f.name)

    if data_filenames:
        st.markdown(f"**Data files / 데이터 파일:** {', '.join(data_filenames)}")

    # --- Run button ---
    if not st.button("▶ Run / 실행", type="primary"):
        return

    # --- Generate & Execute with retry ---
    error_feedback = None
    for attempt in range(1, MAX_RETRIES + 1):
        with st.status(f"Attempt {attempt}/{MAX_RETRIES}: Generating code... / 코드 생성 중...", expanded=True):
            try:
                code = generate_code(md_text, data_filenames, error_feedback=error_feedback)
            except Exception as e:
                st.error(f"AI code generation failed / AI 코드 생성 실패: {e}")
                return
            st.code(code, language="python")

        with st.status(f"Attempt {attempt}/{MAX_RETRIES}: Executing in sandbox... / 샌드박스 실행 중...", expanded=True):
            result = run_in_sandbox(code, data_file_list)

        if result.success:
            break

        st.warning(f"Attempt {attempt} failed / {attempt}회 실패: {result.error}")
        error_feedback = result.error

    # --- Display results ---
    st.divider()

    if not result.success:
        st.error(f"Execution failed after {MAX_RETRIES} attempts. / {MAX_RETRIES}회 시도 후 실패.\n지시 파일을 확인해주세요.")
        st.text(result.error)
        return

    st.success("Execution completed successfully! / 실행 완료!")

    if result.stdout:
        st.subheader("Output / 출력 결과")
        st.text(result.stdout)

    if result.output_files:
        st.subheader("Output Files / 출력 파일")
        for fname, fbytes in result.output_files:
            ext = fname.rsplit(".", 1)[-1].lower() if "." in fname else ""

            if ext in ("png", "jpg", "jpeg", "gif", "svg", "webp"):
                st.image(fbytes, caption=fname)
            elif ext in ("csv", "txt", "json", "log"):
                st.text(fbytes.decode("utf-8", errors="replace"))

            st.download_button(
                label=f"Download / 다운로드: {fname}",
                data=fbytes,
                file_name=fname,
                key=f"dl_{fname}",
            )


if __name__ == "__main__":
    main()
