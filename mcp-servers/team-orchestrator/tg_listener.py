#!/usr/bin/env python3
"""
tg_listener.py
Telegram getUpdates 폴링으로 사용자 답장을 받아 session_bridge inbox에 저장.

동작 방식:
  1. Telegram Bot API getUpdates 폴링 (long-polling, timeout=30)
  2. 메시지가 reply_to_message면 → reply 대상 message_id로 session_id 조회
  3. session_id 있으면 → session_bridge.add_inbox(session_id, text)
  4. 일반 메시지(reply 아님)면 → 최근 세션에 fallback 저장

실행:
    python tg_listener.py
"""
import json
import sys
import time
import urllib.request
import urllib.parse
from pathlib import Path

# Windows cp949 환경 UTF-8 출력
if sys.stdout and hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

SECRETS = json.loads((Path.home() / ".claude" / "secrets.json").read_text(encoding="utf-8"))
BOT_TOKEN = SECRETS["TELEGRAM_BOT_TOKEN"]
ALLOWED_IDS = set(SECRETS.get("TELEGRAM_ALLOWED_USER_IDS", []))

BASE_URL = f"https://api.telegram.org/bot{BOT_TOKEN}"

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from session_bridge import get_session_by_tg_msg_id, add_inbox, _get_conn
import state as st


def _get_updates(offset: int, timeout: int = 30) -> list | None:
    url = f"{BASE_URL}/getUpdates?offset={offset}&timeout={timeout}&allowed_updates=%5B%22message%22%5D"
    try:
        resp = urllib.request.urlopen(url, timeout=timeout + 5)
        data = json.loads(resp.read())
        return data.get("result", [])
    except Exception as e:
        print(f"[listener] getUpdates 실패: {e}", flush=True)
        return None


def _get_latest_session_id() -> str | None:
    """inbox에 보낼 세션이 없을 때 가장 최근 active 세션 반환."""
    try:
        conn = _get_conn()
        row = conn.execute(
            "SELECT session_id FROM sessions ORDER BY created_at DESC LIMIT 1"
        ).fetchone()
        conn.close()
        return row["session_id"] if row else None
    except Exception:
        return None


def _get_latest_running_team_id() -> str | None:
    """현재 running 상태인 가장 최근 팀(리드 우선) 반환."""
    try:
        st.init_db()
        teams = st.list_teams()
        running = [t for t in teams if t.status == "running"]
        if not running:
            return None
        # 리드(워커가 아닌 팀) 우선
        for t in running:
            if not hasattr(t, 'lead_id') or not t.lead_id:
                return t.team_id
        return running[0].team_id
    except Exception:
        return None


def _send_reply(chat_id: int, text: str, reply_to: int | None = None):
    """사용자에게 확인 메시지 전송."""
    payload = {"chat_id": chat_id, "text": text}
    if reply_to:
        payload["reply_to_message_id"] = reply_to
    data = json.dumps(payload).encode()
    req = urllib.request.Request(
        f"{BASE_URL}/sendMessage", data=data, method="POST",
        headers={"Content-Type": "application/json"}
    )
    try:
        urllib.request.urlopen(req, timeout=10)
    except Exception as e:
        print(f"[listener] sendMessage 실패: {e}", flush=True)


def handle_update(update: dict):
    msg = update.get("message")
    if not msg:
        return

    from_id = msg.get("from", {}).get("id")
    chat_id = msg.get("chat", {}).get("id")
    text = msg.get("text", "").strip()
    msg_id = msg.get("message_id")

    # 허용된 사용자만 처리
    if ALLOWED_IDS and from_id not in ALLOWED_IDS:
        print(f"[listener] 미허용 사용자 무시: {from_id}", flush=True)
        return

    if not text:
        return

    session_id = None

    # reply_to_message 있으면 → 해당 msg_id로 세션 조회
    reply_to_msg = msg.get("reply_to_message")
    if reply_to_msg:
        tg_msg_id = reply_to_msg.get("message_id")
        session_info = get_session_by_tg_msg_id(tg_msg_id)
        if session_info:
            session_id = session_info["session_id"]
            print(f"[listener] reply → session {session_id[:8]}...", flush=True)
        else:
            print(f"[listener] reply msg_id={tg_msg_id} 세션 없음 → fallback", flush=True)

    # fallback: 최근 세션
    if not session_id:
        session_id = _get_latest_session_id()
        if session_id:
            print(f"[listener] fallback → session {session_id[:8]}...", flush=True)

    if session_id:
        add_inbox(session_id, text)
        # 답장 대상 메시지에 team_id가 있으면 → state의 CEO 피드백 큐로도 전달
        team_id = None
        if reply_to_msg:
            tg_msg_id = reply_to_msg.get("message_id")
            session_info = get_session_by_tg_msg_id(tg_msg_id)
            if session_info:
                team_id = session_info.get("team_id")
        if not team_id:
            # fallback: 가장 최근 running 팀에 전달
            team_id = _get_latest_running_team_id()
        if team_id:
            try:
                st.init_db()
                st.push_message(team_id, "ceo", f"[사용자 텔레그램 답장] {text}")
                print(f"[listener] CEO 피드백 큐 전달: team={team_id}", flush=True)
            except Exception as e:
                print(f"[listener] CEO 피드백 전달 실패: {e}", flush=True)
        print(f"[listener] inbox 저장: '{text[:50]}'", flush=True)
        _send_reply(chat_id, f"✅ 전달됨: {text[:80]}", reply_to=msg_id)
    else:
        # 연결된 세션 없어도 running 팀이 있으면 직접 전달 시도
        team_id = _get_latest_running_team_id()
        if team_id:
            try:
                st.init_db()
                st.push_message(team_id, "ceo", f"[사용자 텔레그램 메시지] {text}")
                print(f"[listener] 세션 없음, running 팀에 직접 전달: {team_id}", flush=True)
                _send_reply(chat_id, f"✅ 실행 중 팀에 전달됨: {text[:80]}", reply_to=msg_id)
            except Exception as e:
                print(f"[listener] running 팀 전달 실패: {e}", flush=True)
                _send_reply(chat_id, "⚠️ 연결된 Claude 세션이 없습니다.", reply_to=msg_id)
        else:
            print(f"[listener] 저장할 세션 없음. 메시지 무시: '{text[:50]}'", flush=True)
            _send_reply(chat_id, "⚠️ 연결된 Claude 세션이 없습니다.", reply_to=msg_id)


def run():
    print("[listener] Telegram 폴링 시작...", flush=True)
    offset = 0
    retry_delay = 5
    while True:
        updates = _get_updates(offset)
        if updates is None:
            # 409 등 오류 → 대기 후 재시도
            print(f"[listener] {retry_delay}초 후 재시도...", flush=True)
            time.sleep(retry_delay)
            retry_delay = min(retry_delay * 2, 60)
            continue
        retry_delay = 5  # 성공 시 초기화
        for upd in updates:
            offset = upd["update_id"] + 1
            try:
                handle_update(upd)
            except Exception as e:
                print(f"[listener] handle_update 오류: {e}", flush=True)


if __name__ == "__main__":
    run()
