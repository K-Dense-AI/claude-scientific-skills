"""주간 브리핑 데이터 수집 - 3개 스크립트 병렬 실행"""
import subprocess
import sys
import time
from pathlib import Path

PYTHON = r"C:\Users\Jahyun\anaconda3\python.exe"
SCRIPTS_DIR = Path(__file__).parent

SCRIPTS = [
    ("Asana",    SCRIPTS_DIR / "collect_asana.py"),
    ("GitHub",   SCRIPTS_DIR / "collect_github.py"),
    ("OneDrive", SCRIPTS_DIR / "collect_onedrive.py"),
    ("Orders",   SCRIPTS_DIR / "collect_orders.py"),
    ("Reagents", SCRIPTS_DIR / "collect_reagents.py"),
]

def main():
    print("=== Weekly Briefing Collection Start ===")
    start = time.time()

    procs = []
    for name, script in SCRIPTS:
        p = subprocess.Popen(
            [PYTHON, str(script)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        procs.append((name, p))
        print(f"[{name}] starting...")

    results = {}
    for name, p in procs:
        stdout, stderr = p.communicate()
        out = stdout.decode("utf-8", errors="replace") if stdout else ""
        err = stderr.decode("utf-8", errors="replace") if stderr else ""
        # strip replacement chars before printing to cp949 terminal
        out_safe = out.encode("ascii", errors="replace").decode("ascii")
        err_safe = err.encode("ascii", errors="replace").decode("ascii")
        if p.returncode == 0:
            print(f"[{name}] OK\n{out_safe.strip()}")
            results[name] = "ok"
        else:
            print(f"[{name}] FAILED\n{err_safe.strip()}")
            results[name] = f"error: {err_safe.strip()[:200]}"

    elapsed = round(time.time() - start, 1)
    print(f"\n=== Done ({elapsed}s) ===")
    ok = sum(1 for v in results.values() if v == "ok")
    print(f"Success: {ok}/{len(SCRIPTS)}")

    # Calendar is handled by CEO MCP
    print("\n[Calendar] CEO calls gcal_list_events MCP directly")
    print("\nOutput files:")
    temp = Path(r"C:\Users\Jahyun\lab-analyses\temp")
    for f in ["briefing_asana.json", "briefing_github.json", "briefing_onedrive.json", "briefing_orders.json", "briefing_reagents.json"]:
        fp = temp / f
        print(f"  {'OK' if fp.exists() else 'MISSING'} {fp}")

if __name__ == "__main__":
    main()
