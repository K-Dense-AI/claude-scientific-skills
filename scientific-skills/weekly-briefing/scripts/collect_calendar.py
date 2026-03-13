"""Google Calendar 이번 주 일정 수집 - weekly briefing.

사전 조건: python reauth_calendar.py 를 한 번 실행하여 ~/.claude/gcal_token.json 생성 필요.
"""
import json
import sys
import io
from pathlib import Path
from datetime import datetime, timedelta

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sys.path.insert(0, str(Path(__file__).parent))

import requests

OUTPUT = Path(r"C:\Users\Jahyun\.claude\briefing_temp\briefing_calendar.json")


def get_week_range():
    today = datetime.now()
    monday = today - timedelta(days=today.weekday())
    sunday = monday + timedelta(days=6)
    return (
        monday.strftime("%Y-%m-%dT00:00:00+09:00"),
        sunday.strftime("%Y-%m-%dT23:59:59+09:00"),
    )


def collect(calendar_id: str = "primary") -> dict:
    try:
        from google_auth import get_calendar_token
        token = get_calendar_token()
    except FileNotFoundError as e:
        result = {
            "collected_at": datetime.now().isoformat(),
            "status": "auth_required",
            "message": str(e),
            "events_count": 0,
            "events": [],
        }
        OUTPUT.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"[WARN] {e}")
        return result

    time_min, time_max = get_week_range()

    resp = requests.get(
        f"https://www.googleapis.com/calendar/v3/calendars/{calendar_id}/events",
        headers={"Authorization": f"Bearer {token}"},
        params={
            "timeMin": time_min,
            "timeMax": time_max,
            "singleEvents": True,
            "orderBy": "startTime",
        },
    )
    resp.raise_for_status()
    items = resp.json().get("items", [])

    events = []
    for e in items:
        start = e.get("start", {})
        start_str = start.get("dateTime", start.get("date", ""))
        events.append({
            "summary": e.get("summary", "(제목 없음)"),
            "start": start_str,
            "location": e.get("location", ""),
            "description": e.get("description", "")[:200],
        })

    result = {
        "collected_at": datetime.now().isoformat(),
        "week_range": {"start": time_min, "end": time_max},
        "status": "ok",
        "events_count": len(events),
        "events": events,
    }
    OUTPUT.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return result


if __name__ == "__main__":
    r = collect()
    if r["status"] == "ok":
        print(f"events_count={r['events_count']}")
        for e in r["events"]:
            print(f"  {e['start'][:10]} {e['summary']}")
    else:
        print(f"[AUTH REQUIRED] {r['message']}")
        print("실행: python reauth_calendar.py")
