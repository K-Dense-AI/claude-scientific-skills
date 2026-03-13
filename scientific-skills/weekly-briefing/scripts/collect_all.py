"""주간 브리핑 데이터 수집 - 전체 스크립트 병렬 실행

한글 깨짐 해결: PYTHONIOENCODING=utf-8 환경변수 설정으로 서브프로세스 출력 UTF-8 강제.
"""
import subprocess
import sys
import io
import os
import time
from pathlib import Path

# Windows cp949 터미널 한글 깨짐 방지
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

PYTHON = r"C:\Users\Jahyun\anaconda3\python.exe"
SCRIPTS_DIR = Path(__file__).parent

SCRIPTS = [
    ("Asana",    SCRIPTS_DIR / "collect_asana.py"),
    ("GitHub",   SCRIPTS_DIR / "collect_github.py"),
    ("OneDrive", SCRIPTS_DIR / "collect_onedrive.py"),
    ("Orders",   SCRIPTS_DIR / "collect_orders.py"),
    ("Reagents", SCRIPTS_DIR / "collect_reagents.py"),
    ("GDrive",   SCRIPTS_DIR / "collect_gdrive.py"),
    ("Calendar", SCRIPTS_DIR / "collect_calendar.py"),
]

# UTF-8 강제 환경변수 (서브프로세스 한글 깨짐 방지)
ENV = {**os.environ, "PYTHONIOENCODING": "utf-8"}


def main():
    print("=== Weekly Briefing Collection Start ===")
    start = time.time()

    procs = []
    for name, script in SCRIPTS:
        p = subprocess.Popen(
            [PYTHON, str(script)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=ENV,
        )
        procs.append((name, p))
        print(f"[{name}] starting...")

    results = {}
    for name, p in procs:
        stdout, stderr = p.communicate()
        out = stdout.decode("utf-8", errors="replace") if stdout else ""
        err = stderr.decode("utf-8", errors="replace") if stderr else ""
        if p.returncode == 0:
            print(f"[{name}] OK\n{out.strip()}")
            results[name] = "ok"
        else:
            print(f"[{name}] FAILED\n{err.strip()[:300]}")
            results[name] = f"error: {err.strip()[:200]}"

    elapsed = round(time.time() - start, 1)
    print(f"\n=== Done ({elapsed}s) ===")
    ok = sum(1 for v in results.values() if v == "ok")
    print(f"Success: {ok}/{len(SCRIPTS)}")

    output_files = [
        "briefing_asana.json", "briefing_github.json", "briefing_onedrive.json",
        "briefing_orders.json", "briefing_reagents.json",
        "briefing_gdrive.json", "briefing_calendar.json",
    ]
    print("\nOutput files:")
    temp = Path(r"C:\Users\Jahyun\.claude\briefing_temp")
    for f in output_files:
        fp = temp / f
        print(f"  {'OK' if fp.exists() else 'MISSING'} {fp}")

    # Calendar auth 안내
    cal = results.get("Calendar", "")
    if "auth_required" in cal or ("WARN" in out and "reauth" in out):
        print("\n[Calendar] 최초 1회 인증 필요:")
        print(f"  {PYTHON} {SCRIPTS_DIR / 'reauth_calendar.py'}")


if __name__ == "__main__":
    main()
