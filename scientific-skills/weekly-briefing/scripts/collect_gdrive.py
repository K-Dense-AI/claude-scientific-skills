"""Google Drive 최근 변경 파일 수집 - weekly briefing"""
import json
import sys
import os
from pathlib import Path
from datetime import datetime, timedelta

sys.path.insert(0, str(Path(__file__).parent))
from google_auth import get_access_token
import requests

OUTPUT = Path(r"C:\Users\Jahyun\.claude\briefing_temp\briefing_gdrive.json")


def collect(days: int = 7, page_size: int = 20) -> dict:
    token = get_access_token()
    since = (datetime.now() - timedelta(days=days)).strftime("%Y-%m-%dT%H:%M:%S")

    resp = requests.get(
        "https://www.googleapis.com/drive/v3/files",
        headers={"Authorization": f"Bearer {token}"},
        params={
            "q": f"modifiedTime > '{since}'",
            "orderBy": "modifiedTime desc",
            "pageSize": page_size,
            "fields": "files(id,name,mimeType,modifiedTime,webViewLink,lastModifyingUser/displayName)",
        },
    )
    resp.raise_for_status()
    files = resp.json().get("files", [])

    result = {
        "collected_at": datetime.now().isoformat(),
        "period_days": days,
        "files_changed": len(files),
        "files": [
            {
                "name": f["name"],
                "type": f.get("mimeType", ""),
                "modified": f.get("modifiedTime", ""),
                "link": f.get("webViewLink", ""),
                "modifier": f.get("lastModifyingUser", {}).get("displayName", ""),
            }
            for f in files
        ],
    }

    OUTPUT.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return result


if __name__ == "__main__":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    r = collect()
    print(f"files_changed={r['files_changed']}")
    for f in r["files"][:5]:
        print(f"  {f['name']} ({f['modified'][:10]})")
