#!/usr/bin/env python3
"""
WebSocket 서버 + ngrok 터널 자동 실행
Notion embed용 공개 URL을 자동으로 생성
"""
import subprocess
import time
import sys
from pathlib import Path

from pyngrok import ngrok

# ngrok 토큰 설정
NGROK_TOKENS = [
    "TTNX2VY77S",
    "P6H7CEWZHP",
    "Z2G68SNM24",
    "CUEK9ZA3TY",
    "Q9HAY9MQVU",
    "WAZA7TNYHY",
    "CKP972EVYP",
    "6Y9C842SUX",
    "V3MBAWTMK9",
    "JR37WV4B7W",
]

def setup_ngrok():
    """ngrok 토큰 설정"""
    token = NGROK_TOKENS[0]  # 첫 번째 토큰 사용
    try:
        ngrok.set_auth_token(token)
        print(f"✅ ngrok 토큰 설정 완료: {token[:10]}...")
        return True
    except Exception as e:
        print(f"❌ ngrok 토큰 설정 실패: {e}")
        return False

def main():
    print("\n" + "="*60)
    print("  Team Orchestrator + ngrok 터널")
    print("="*60 + "\n")

    # 1. WebSocket 서버 시작 (백그라운드)
    print("[1/2] WebSocket 서버 시작 중 (포트 8765)...")
    script_dir = Path(__file__).parent

    try:
        ws_proc = subprocess.Popen(
            [sys.executable, str(script_dir / "ws_server.py")],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=str(script_dir)
        )
        print("✅ WebSocket 서버 시작됨 (PID: {})".format(ws_proc.pid))
    except Exception as e:
        print(f"❌ WebSocket 서버 실행 실패: {e}")
        return

    # 2. ngrok 토큰 설정
    time.sleep(1)
    print("\n[2/3] ngrok 토큰 설정 중...")
    if not setup_ngrok():
        ws_proc.terminate()
        return

    # 3. ngrok 터널 시작
    time.sleep(1)
    print("\n[3/3] ngrok 터널 시작 중...")
    try:
        public_url = ngrok.connect(8765, "http")
        print("\n" + "="*60)
        print(f"🎉 공개 URL 생성됨!")
        print("="*60)
        print(f"\n📱 Notion embed URL:")
        print(f"   {public_url}\n")
        print("📋 Notion 페이지에서:")
        print("   1. + 블록 추가 → Embed")
        print(f"   2. URL 입력: {public_url}")
        print("   3. 엔터 → 완성! ✨\n")
        print("💡 종료하려면 Ctrl+C 누르세요")
        print("="*60 + "\n")

        # ngrok 터널 유지
        ngrok_process = ngrok.get_ngrok_process()
        ngrok_process.proc.wait()

    except KeyboardInterrupt:
        print("\n\n⛔ 종료 중...")
        ngrok.kill()
        ws_proc.terminate()
        print("✅ 서버와 터널이 종료되었습니다.")
    except Exception as e:
        print(f"❌ ngrok 터널 실패: {e}")
        ws_proc.terminate()

if __name__ == "__main__":
    main()
