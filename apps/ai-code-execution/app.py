"""Streamlit UI for AI-Powered Code Execution System."""

import streamlit as st
from config import AI_BACKEND, MAX_RETRIES, MAX_MD_TOKENS
from ai_codegen import generate_code
from executor import run_in_sandbox


def detect_mode(md_text: str) -> str:
    if "```python" in md_text:
        return "with_code"
    return "no_code"


def estimate_tokens(text: str) -> int:
    return len(text) // 4


def main():
    st.set_page_config(page_title="AI Code Execution", page_icon="▶", layout="wide")

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
