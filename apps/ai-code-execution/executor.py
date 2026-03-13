"""E2B sandbox executor — uploads data, runs code, collects outputs."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from e2b_code_interpreter import Sandbox

from config import SANDBOX_TIMEOUT, DATA_DIR, OUTPUT_DIR


@dataclass
class ExecutionResult:
    stdout: str = ""
    output_files: list[tuple[str, bytes]] = field(default_factory=list)
    success: bool = True
    error: str = ""


def run_in_sandbox(code: str,
                   data_files: list[tuple[str, bytes]] | None = None) -> ExecutionResult:
    """Execute Python code in an E2B sandbox.

    Args:
        code: Python code string to execute.
        data_files: List of (filename, file_bytes) tuples to upload.

    Returns:
        ExecutionResult with stdout, output files, success status, and error.
    """
    api_key = os.environ.get("E2B_API_KEY")
    if not api_key:
        return ExecutionResult(
            success=False,
            error="E2B_API_KEY not set. Add it to .env file.",
        )

    sbx = Sandbox.create(api_key=api_key, timeout=SANDBOX_TIMEOUT)
    try:
        # Create directories
        sbx.commands.run(f"mkdir -p {DATA_DIR} {OUTPUT_DIR}")

        # Upload data files
        if data_files:
            for filename, file_bytes in data_files:
                sbx.files.write(f"{DATA_DIR}{filename}", file_bytes)

        # Execute code
        execution = sbx.run_code(code)

        # Check for errors
        if execution.error:
            return ExecutionResult(
                stdout="".join(execution.logs.stdout),
                success=False,
                error=f"{execution.error.name}: {execution.error.value}\n{execution.error.traceback}",
            )

        # Collect stdout
        stdout = "".join(execution.logs.stdout)

        # Collect output files
        output_files = []
        try:
            file_list = sbx.files.list(OUTPUT_DIR)
            for f_info in file_list:
                if f_info.is_dir:
                    continue
                content = sbx.files.read(f"{OUTPUT_DIR}{f_info.name}", format="bytes")
                output_files.append((f_info.name, content))
        except Exception:
            pass  # No output files is OK

        return ExecutionResult(
            stdout=stdout,
            output_files=output_files,
            success=True,
        )

    except Exception as e:
        return ExecutionResult(success=False, error=str(e))
    finally:
        sbx.kill()
