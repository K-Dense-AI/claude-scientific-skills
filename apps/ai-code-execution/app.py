"""Streamlit UI for AI-Powered Code Execution System."""

import io
import streamlit as st
from config import AI_BACKEND, MAX_RETRIES, MAX_MD_TOKENS
from ai_codegen import generate_code
from executor import run_in_sandbox


def detect_mode(md_text: str) -> str:
    """Detect instruction mode from .md content.

    Returns 'with_code' if python code blocks exist, else 'no_code'.
    """
    if "```python" in md_text:
        return "with_code"
    return "no_code"


def estimate_tokens(text: str) -> int:
    """Rough token estimate (1 token ~ 4 chars)."""
    return len(text) // 4


def main():
    st.set_page_config(page_title="AI Code Execution", page_icon="▶", layout="wide")
    st.title("AI-Powered Code Execution")
    st.caption(f"Backend: **{AI_BACKEND}** | Max retries: {MAX_RETRIES}")

    # --- File uploads ---
    col1, col2 = st.columns(2)
    with col1:
        md_file = st.file_uploader(
            "Upload instruction (.md)", type=["md", "txt"], key="md_upload"
        )
    with col2:
        data_files = st.file_uploader(
            "Upload data files (optional)",
            accept_multiple_files=True,
            key="data_upload",
        )

    if not md_file:
        st.info("Upload an instruction .md file to get started.")
        return

    # --- Read instruction ---
    md_text = md_file.read().decode("utf-8", errors="replace")
    mode = detect_mode(md_text)

    token_est = estimate_tokens(md_text)
    if token_est > MAX_MD_TOKENS:
        st.warning(f"Instruction file is very large (~{token_est:,} tokens). Processing may be slow.")

    st.markdown(f"**Mode detected:** `{'WITH CODE' if mode == 'with_code' else 'NO CODE'}`")

    with st.expander("Preview instruction", expanded=False):
        st.markdown(md_text)

    # --- Prepare data files ---
    data_file_list: list[tuple[str, bytes]] = []
    data_filenames: list[str] = []
    for f in data_files:
        content = f.read()
        data_file_list.append((f.name, content))
        data_filenames.append(f.name)

    if data_filenames:
        st.markdown(f"**Data files:** {', '.join(data_filenames)}")

    # --- Run button ---
    if not st.button("Run", type="primary"):
        return

    # --- Generate & Execute with retry ---
    error_feedback = None
    for attempt in range(1, MAX_RETRIES + 1):
        # Step 1: Generate code
        with st.status(f"Attempt {attempt}/{MAX_RETRIES}: Generating code...", expanded=True):
            try:
                code = generate_code(md_text, data_filenames, error_feedback=error_feedback)
            except Exception as e:
                st.error(f"AI code generation failed: {e}")
                return

            st.code(code, language="python")

        # Step 2: Execute in sandbox
        with st.status(f"Attempt {attempt}/{MAX_RETRIES}: Executing in sandbox...", expanded=True):
            result = run_in_sandbox(code, data_file_list)

        if result.success:
            break

        # Execution failed — prepare retry
        st.warning(f"Attempt {attempt} failed: {result.error}")
        error_feedback = result.error

    # --- Display results ---
    st.divider()

    if not result.success:
        st.error(f"Execution failed after {MAX_RETRIES} attempts.\nPlease review your instruction file.")
        st.text(result.error)
        return

    st.success("Execution completed successfully!")

    # Text output
    if result.stdout:
        st.subheader("Output")
        st.text(result.stdout)

    # Output files
    if result.output_files:
        st.subheader("Output Files")
        for fname, fbytes in result.output_files:
            ext = fname.rsplit(".", 1)[-1].lower() if "." in fname else ""

            # Display images inline
            if ext in ("png", "jpg", "jpeg", "gif", "svg", "webp"):
                st.image(fbytes, caption=fname)

            # Display CSV/text inline
            elif ext in ("csv", "txt", "json", "log"):
                st.text(fbytes.decode("utf-8", errors="replace"))

            # Download button for all files
            st.download_button(
                label=f"Download {fname}",
                data=fbytes,
                file_name=fname,
                key=f"dl_{fname}",
            )


if __name__ == "__main__":
    main()
