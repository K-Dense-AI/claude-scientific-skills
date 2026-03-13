"""AI code generation backend — supports Gemini, Claude, and OpenAI."""

import os
import re
from config import AI_BACKEND, MODEL_DEFAULTS, API_KEY_VARS, DATA_DIR, OUTPUT_DIR


SYSTEM_PROMPT = f"""You are a Python code generator. Read the instruction below and write
Python code that fulfills it exactly.

Rules:
- Data files are available at: {DATA_DIR}{{filename}}
- Save all output files to: {OUTPUT_DIR}
- Print a summary of results to stdout
- Use only standard libraries + pandas, matplotlib, numpy, openpyxl, pillow, scipy, scikit-learn, seaborn
- Return ONLY the Python code block, no explanation"""


def _build_prompt(instruction_text: str, data_filenames: list[str],
                  error_feedback: str | None = None) -> str:
    """Build the user prompt for AI."""
    prompt = f"[INSTRUCTION]\n{instruction_text}\n\n"
    prompt += f"[DATA FILES AVAILABLE]\n{', '.join(data_filenames) if data_filenames else 'None'}\n"
    if error_feedback:
        prompt += f"\n[PREVIOUS ERROR]\nThis code failed with error:\n{error_feedback}\nFix it and return corrected code.\n"
    return prompt


def _extract_code(response_text: str) -> str:
    """Extract Python code from AI response (handles ```python blocks or raw code)."""
    match = re.search(r"```python\s*\n(.*?)```", response_text, re.DOTALL)
    if match:
        return match.group(1).strip()
    match = re.search(r"```\s*\n(.*?)```", response_text, re.DOTALL)
    if match:
        return match.group(1).strip()
    return response_text.strip()


def generate_code(instruction_text: str, data_filenames: list[str],
                  error_feedback: str | None = None,
                  backend: str | None = None) -> str:
    """Generate Python code from instruction using the configured AI backend.

    Args:
        instruction_text: Full .md content.
        data_filenames: List of uploaded file names.
        error_feedback: Previous error message for retry.
        backend: Override AI_BACKEND setting.

    Returns:
        Python code string ready to execute.
    """
    backend = backend or AI_BACKEND
    user_prompt = _build_prompt(instruction_text, data_filenames, error_feedback)

    if backend == "gemini":
        return _generate_gemini(user_prompt)
    elif backend == "claude":
        return _generate_claude(user_prompt)
    elif backend == "openai":
        return _generate_openai(user_prompt)
    elif backend == "groq":
        return _generate_groq(user_prompt)
    else:
        raise ValueError(f"Unknown AI backend: {backend}")


def _generate_gemini(user_prompt: str) -> str:
    import google.generativeai as genai

    api_key = os.environ.get(API_KEY_VARS["gemini"])
    if not api_key:
        raise ValueError("GEMINI_API_KEY not set. Add it to .env file.")

    genai.configure(api_key=api_key)
    model = genai.GenerativeModel(MODEL_DEFAULTS["gemini"])
    response = model.generate_content(
        f"{SYSTEM_PROMPT}\n\n{user_prompt}",
    )
    return _extract_code(response.text)


def _generate_claude(user_prompt: str) -> str:
    import anthropic

    api_key = os.environ.get(API_KEY_VARS["claude"])
    if not api_key:
        raise ValueError("ANTHROPIC_API_KEY not set. Add it to .env file.")

    client = anthropic.Anthropic(api_key=api_key)
    message = client.messages.create(
        model=MODEL_DEFAULTS["claude"],
        max_tokens=4096,
        system=SYSTEM_PROMPT,
        messages=[{"role": "user", "content": user_prompt}],
    )
    return _extract_code(message.content[0].text)


def _generate_groq(user_prompt: str) -> str:
    from openai import OpenAI

    api_key = os.environ.get(API_KEY_VARS["groq"])
    if not api_key:
        raise ValueError("GROQ_API_KEY not set. Add it to .env file.")

    client = OpenAI(api_key=api_key, base_url="https://api.groq.com/openai/v1")
    response = client.chat.completions.create(
        model=MODEL_DEFAULTS["groq"],
        messages=[
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": user_prompt},
        ],
    )
    return _extract_code(response.choices[0].message.content)


def _generate_openai(user_prompt: str) -> str:
    from openai import OpenAI

    api_key = os.environ.get(API_KEY_VARS["openai"])
    if not api_key:
        raise ValueError("OPENAI_API_KEY not set. Add it to .env file.")

    client = OpenAI(api_key=api_key)
    response = client.chat.completions.create(
        model=MODEL_DEFAULTS["openai"],
        messages=[
            {"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": user_prompt},
        ],
    )
    return _extract_code(response.choices[0].message.content)
