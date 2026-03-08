"""
Telegram Bot Token을 secrets.json에 추가하는 헬퍼 스크립트.
실행: python setup_token.py
"""
import json
from pathlib import Path

SECRETS_FILE = Path.home() / ".claude" / "secrets.json"

secrets = json.loads(SECRETS_FILE.read_text(encoding="utf-8"))

print("=== Claude Telegram Bot 설정 ===\n")

# Bot Token
token = input("Telegram Bot Token을 입력하세요 (BotFather에서 발급): ").strip()
if token:
    secrets["TELEGRAM_BOT_TOKEN"] = token

# 허용 유저 ID (선택)
print("\n접근 제한 설정 (선택사항):")
print("  본인 Telegram User ID만 허용하려면 ID를 입력하세요.")
print("  (Telegram에서 @userinfobot 에게 메시지 보내면 ID 확인 가능)")
print("  빈칸이면 모든 사용자 허용")
user_ids_input = input("허용 User ID (콤마 구분, 예: 123456789,987654321): ").strip()
if user_ids_input:
    try:
        ids = [int(x.strip()) for x in user_ids_input.split(",") if x.strip()]
        secrets["TELEGRAM_ALLOWED_USER_IDS"] = ids
        print(f"허용 ID 설정: {ids}")
    except ValueError:
        print("숫자 형식 오류 — 허용 ID 미설정 (모두 허용)")

SECRETS_FILE.write_text(json.dumps(secrets, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"\nsecrets.json 업데이트 완료: {SECRETS_FILE}")
print("\n이제 start.bat (또는 python bot.py) 으로 봇을 실행하세요!")
