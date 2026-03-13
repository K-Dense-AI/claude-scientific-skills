"""Quick smoke test: executor → codegen → full pipeline."""
import sys, os, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# .env 로드 (config.py가 이미 하지만 직접 실행 시 경로 보장)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dotenv import load_dotenv
load_dotenv(os.path.join(os.path.dirname(__file__), ".env"))

from executor import run_in_sandbox
from ai_codegen import generate_code

# ── TEST 1: Executor 직접 테스트 ──────────────────────────────
print("=" * 50)
print("[TEST 1] E2B Executor — hello world")
result = run_in_sandbox("print('Hello from E2B!')\nprint(1 + 2)")
if result.success:
    print(f"[OK] stdout: {result.stdout!r}")
else:
    print(f"[FAIL] {result.error}")

# ── TEST 2: AI Codegen 테스트 (Claude) ────────────────────────
print()
print("=" * 50)
print("[TEST 2] AI Codegen (Claude) — 간단한 지시")
instruction = "Print the sum of numbers 1 to 10."
try:
    code = generate_code(instruction, [])
    print(f"[OK] Generated code:\n{code}")
except Exception as e:
    print(f"[FAIL] {e}")

# ── TEST 3: 전체 파이프라인 (AI → E2B) ───────────────────────
print()
print("=" * 50)
print("[TEST 3] Full pipeline — AI generates code, E2B runs it")
if 'code' in dir() and code:
    result2 = run_in_sandbox(code)
    if result2.success:
        print(f"[OK] stdout: {result2.stdout!r}")
    else:
        print(f"[FAIL] {result2.error}")
else:
    print("[SKIP] codegen failed, skipping")
