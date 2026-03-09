#!/usr/bin/env python3
"""
작업 완료 시 Telegram 메시지 전송 스크립트.
CEO 세션에서 에이전트 완료 후 사용자에게 알림.

사용법:
    python telegram_notify.py "메시지 내용"
    python telegram_notify.py --agent quick-XXXXX-1  # 에이전트 완료 대기 후 알림
"""
import io
import json
import sys
import time
import sqlite3
import urllib.request
from pathlib import Path

# Windows cp949 터미널에서 유니코드 출력 가능하도록
if sys.stdout and hasattr(sys.stdout, 'reconfigure'):
    try:
        sys.stdout.reconfigure(encoding='utf-8', errors='replace')
    except Exception:
        pass

SECRETS = json.loads((Path.home() / ".claude" / "secrets.json").read_text(encoding="utf-8"))
BOT_TOKEN = SECRETS["TELEGRAM_BOT_TOKEN"]
USER_IDS = SECRETS.get("TELEGRAM_ALLOWED_USER_IDS", [])
DB = Path(__file__).parent / "orchestrator.db"


def _format_table_for_telegram(text: str) -> str:
    """
    마크다운 테이블을 텔레그램 코드블록으로 감싸기.
    텔레그램이 마크다운 테이블을 HTML로 렌더링하면 복사할 때 망가지므로,
    코드블록으로 감싸면 순수 텍스트로 유지됨.
    """
    if "| " in text and " |" in text:  # 마크다운 테이블 감지
        lines = text.split("\n")
        result = []
        in_table = False
        for line in lines:
            if "| " in line and " |" in line:
                if not in_table:
                    result.append("```")
                    in_table = True
                result.append(line)
            else:
                if in_table:
                    result.append("```")
                    in_table = False
                result.append(line)
        if in_table:
            result.append("```")
        return "\n".join(result)
    return text


TG_MAX_LEN = 4000  # 텔레그램 4096 한도, 여유분 확보


def _split_message(text: str, max_len: int = TG_MAX_LEN) -> list[str]:
    """긴 메시지를 max_len 이하 청크로 분할. 줄 단위로 자름."""
    if len(text) <= max_len:
        return [text]
    chunks = []
    while text:
        if len(text) <= max_len:
            chunks.append(text)
            break
        # 줄바꿈 기준으로 자르기
        cut = text.rfind("\n", 0, max_len)
        if cut <= 0:
            cut = max_len
        chunks.append(text[:cut])
        text = text[cut:].lstrip("\n")
    # 파트 번호 추가
    if len(chunks) > 1:
        total = len(chunks)
        chunks = [f"[{i+1}/{total}]\n{c}" for i, c in enumerate(chunks)]
    return chunks


def send(message: str) -> int | None:
    """모든 허용 사용자에게 텔레그램 메시지 발송. 4096자 초과 시 자동 분할."""
    formatted_msg = _format_table_for_telegram(message)
    chunks = _split_message(formatted_msg)
    first_msg_id = None
    for uid in USER_IDS:
        for chunk in chunks:
            url = f"https://api.telegram.org/bot{BOT_TOKEN}/sendMessage"
            data = json.dumps({"chat_id": uid, "text": chunk}).encode()
            req = urllib.request.Request(url, data=data, method="POST",
                                         headers={"Content-Type": "application/json"})
            try:
                resp = urllib.request.urlopen(req, timeout=10)
                if first_msg_id is None:
                    result = json.loads(resp.read())
                    first_msg_id = result.get("result", {}).get("message_id")
            except Exception as e:
                print(f"[notify] 전송 실패 uid={uid}: {e}")
    return first_msg_id


def send_with_session(message: str, session_id: str, cwd: str,
                      team_id: str = None) -> int | None:
    """메시지 전송 후 tg_msg_id -> session_id 매핑을 session_bridge DB에 저장.
    team_id가 있으면 답장→팀 전달 경로도 함께 저장."""
    import sys as _sys
    _sys.path.insert(0, str(Path(__file__).parent))
    from session_bridge import save_session

    formatted_msg = _format_table_for_telegram(message)
    chunks = _split_message(formatted_msg)
    first_msg_id = None
    for uid in USER_IDS:
        for chunk in chunks:
            url = f"https://api.telegram.org/bot{BOT_TOKEN}/sendMessage"
            data = json.dumps({"chat_id": uid, "text": chunk}).encode()
            req = urllib.request.Request(url, data=data, method="POST",
                                         headers={"Content-Type": "application/json"})
            try:
                resp = urllib.request.urlopen(req, timeout=10)
                result = json.loads(resp.read())
                msg_id = result["result"]["message_id"]
                if first_msg_id is None:
                    first_msg_id = msg_id
                # 각 청크의 msg_id를 세션에 매핑 (어느 청크에 답장해도 추적 가능)
                save_session(session_id, cwd, msg_id, team_id=team_id)
            except Exception as e:
                print(f"[notify] 전송 실패 uid={uid}: {e}")
    return first_msg_id


def poll_and_notify(agent_id: str, interval: int = 20, timeout: int = 600):
    """에이전트 완료까지 폴링 후 결과 알림."""
    print(f"[poll] {agent_id} 모니터링 시작 (간격 {interval}초, 최대 {timeout}초)")
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            if DB.exists() and DB.stat().st_size > 0:
                conn = sqlite3.connect(str(DB))
                row = conn.execute(
                    "SELECT name, status, output FROM teams WHERE team_id=?", (agent_id,)
                ).fetchone()
                conn.close()
                if row and row[1] in ("done", "failed"):
                    name, status, output = row
                    emoji = "✅" if status == "done" else "❌"
                    summary = (output or "").strip()[:800]
                    if len((output or "").strip()) > 800:
                        summary += "\n...(잘림)"
                    msg = f"{emoji} 에이전트 완료: {name} [{status}]\n---\n{summary}"
                    formatted_msg = _format_table_for_telegram(msg)
                    send(formatted_msg)
                    print(f"[notify] 텔레그램 전송 완료: {agent_id} [{status}]")
                    sys.exit(0)
        except sqlite3.OperationalError:
            pass
        time.sleep(interval)
    formatted_msg = _format_table_for_telegram(f"⏰ 에이전트 타임아웃: {agent_id} ({timeout}초 초과)")
    send(formatted_msg)
    sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("사용법: telegram_notify.py '메시지' 또는 --agent <agent_id>")
        sys.exit(1)

    if sys.argv[1] == "--agent" and len(sys.argv) > 2:
        poll_and_notify(sys.argv[2])
    elif sys.argv[1] == "--session" and len(sys.argv) > 4:
        # --session <session_id> <cwd> <message...>
        session_id = sys.argv[2]
        cwd = sys.argv[3]
        msg = " ".join(sys.argv[4:])
        msg_id = send_with_session(msg, session_id, cwd)
        print(f"[notify] 세션 전송: {session_id[:8]}... msg_id={msg_id}")
    else:
        msg = " ".join(sys.argv[1:])
        send(msg)
        print(f"[notify] 전송: {msg}")
