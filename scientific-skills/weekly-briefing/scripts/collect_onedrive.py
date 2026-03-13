"""OneDrive 로컬 폴더 최근 변경 파일 수집 - 주간 브리핑용"""
import json
import sys
import io
import os
from datetime import datetime, timedelta
from pathlib import Path

# Windows cp949 환경에서 한글 출력 시 깨짐 방지
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

TEMP_DIR = Path(r"C:\Users\Jahyun\.claude\briefing_temp")
OUTPUT_FILE = TEMP_DIR / "briefing_onedrive.json"

SCAN_FOLDERS = {
    "호서대_바탕화면": r"C:\Users\Jahyun\OneDrive - 호서대학교\바탕 화면 [Labtop]",
    "호서대_논문": r"C:\Users\Jahyun\OneDrive - 호서대학교\논문",
    "고려대_저장소": r"C:\Users\Jahyun\OneDrive - 고려대학교\저장소",
}

SKIP_EXTS = {".tmp", ".lnk", ".ini", ".db", ".dat", ".log"}
SKIP_DIRS = {"$RECYCLE.BIN", ".git", "__pycache__", ".DS_Store"}

def collect():
    cutoff = datetime.now() - timedelta(days=7)
    results = {"folders": {}, "collected_at": str(datetime.now().date())}

    for label, folder_path in SCAN_FOLDERS.items():
        root = Path(folder_path)
        if not root.exists():
            results["folders"][label] = {"error": "폴더 없음"}
            continue

        recent_files = []
        try:
            for dirpath, dirnames, filenames in os.walk(root):
                # 불필요 폴더 스킵
                dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
                for fname in filenames:
                    ext = Path(fname).suffix.lower()
                    if ext in SKIP_EXTS:
                        continue
                    fpath = Path(dirpath) / fname
                    try:
                        mtime = datetime.fromtimestamp(fpath.stat().st_mtime)
                        if mtime >= cutoff:
                            rel = str(fpath.relative_to(root))
                            recent_files.append({
                                "name": fname,
                                "path": rel,
                                "modified": mtime.strftime("%Y-%m-%d %H:%M"),
                                "size_kb": round(fpath.stat().st_size / 1024, 1)
                            })
                    except (OSError, PermissionError):
                        continue
        except (OSError, PermissionError) as e:
            results["folders"][label] = {"error": str(e)}
            continue

        # 최신순 정렬
        recent_files.sort(key=lambda x: x["modified"], reverse=True)
        results["folders"][label] = {
            "path": folder_path,
            "recent_files": recent_files[:20],
            "total_count": len(recent_files)
        }

    TEMP_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_FILE.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    total = sum(f.get("total_count", 0) for f in results["folders"].values() if isinstance(f, dict))
    print(f"[OneDrive] done: {total} files changed in last 7 days")
    print(f"saved: {OUTPUT_FILE}")

if __name__ == "__main__":
    collect()
