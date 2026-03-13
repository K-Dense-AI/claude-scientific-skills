"""Asana 태스크 수집 - 주간 브리핑용"""
import json
import sys
import io
import os
from datetime import datetime, timedelta, timezone
from pathlib import Path

# Windows cp949 환경에서 한글 출력 시 깨짐 방지
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

sys.path.insert(0, r"C:\Users\Jahyun\claude-scientific-skills\scientific-skills\asana-extended-api\scripts")
from asana_api import AsanaAPI

TEMP_DIR = Path(r"C:\Users\Jahyun\.claude\briefing_temp")
OUTPUT_FILE = TEMP_DIR / "briefing_asana.json"

# secrets.json에서 토큰 로드
secrets = json.loads(Path("~/.claude/secrets.json").expanduser().read_text())
api = AsanaAPI(pat=secrets["ASANA_PAT"])

PROJECTS = {
    "팀원관리": "1213473127047670",
}

def collect():
    today = datetime.now(timezone.utc).date()
    week_end = today + timedelta(days=7)
    week_start = today - timedelta(days=7)

    results = {"overdue": [], "due_this_week": [], "new_this_week": [], "collected_at": str(today)}

    for proj_name, proj_gid in PROJECTS.items():
        try:
            tasks = api.get_tasks(
                project_gid=proj_gid,
                opt_fields="name,assignee.name,due_on,completed,created_at,permalink_url"
            )
            for t in tasks:
                if t.get("completed"):
                    continue
                due = t.get("due_on")
                created = t.get("created_at", "")[:10] if t.get("created_at") else None
                name = t.get("name", "")
                assignee = t.get("assignee", {}).get("name", "미배정") if t.get("assignee") else "미배정"
                url = t.get("permalink_url") or f"https://app.asana.com/0/0/{t['gid']}/f"
                item = {"name": name, "assignee": assignee, "due_on": due, "url": url, "project": proj_name}

                if due:
                    due_date = datetime.strptime(due, "%Y-%m-%d").date()
                    if due_date < today:
                        item["days_overdue"] = (today - due_date).days
                        results["overdue"].append(item)
                    elif due_date <= week_end:
                        item["days_until"] = (due_date - today).days
                        results["due_this_week"].append(item)

                if created and created >= str(week_start):
                    results["new_this_week"].append(item)

        except Exception as e:
            results[f"error_{proj_name}"] = str(e)

    TEMP_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_FILE.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[Asana] done: overdue={len(results['overdue'])}, due_this_week={len(results['due_this_week'])}, new={len(results['new_this_week'])}")
    print(f"saved: {OUTPUT_FILE}")

if __name__ == "__main__":
    collect()
