#!/usr/bin/env python3
"""
자동 시작 — 내부에서 직접 실행 (subprocess 방식)
"""
import subprocess
import sys
import os
import time
import threading

def run_ws_server():
    """WebSocket 서버 실행"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    subprocess.Popen(
        [sys.executable, os.path.join(script_dir, "ws_server.py")],
        cwd=script_dir,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )


def run_ngrok():
    """ngrok 실행"""
    time.sleep(2)
    try:
        from pyngrok import ngrok
        ngrok.set_auth_token("TTNX2VY77S")
        public_url = ngrok.connect(8765, "http")

        print("\n" + "="*70)
        print("🎉 성공! Notion Embed URL 생성됨")
        print("="*70)
        print(f"\n📱 복사할 URL:\n   {public_url}\n")
        print("📋 Notion 페이지 (이미 열려있음):")
        print("   https://www.notion.so/Orchestrator-Dashboard-31cf91aca96f8107819bef1d4b05900a")
        print("\n✅ URL을 복사하고 Notion에서 [+ 블록 추가] → [Embed]로 추가하세요!")
        print("="*70 + "\n")

        # ngrok 유지
        ngrok_proc = ngrok.get_ngrok_process()
        ngrok_proc.proc.wait()
    except Exception as e:
        print(f"❌ 오류: {e}")

if __name__ == "__main__":
    print("\n⏳ WebSocket 서버와 ngrok 시작 중...\n")

    # 백그라운드에서 실행
    thread1 = threading.Thread(target=run_ws_server, daemon=True)
    thread2 = threading.Thread(target=run_ngrok, daemon=True)

    thread1.start()
    thread2.start()

    # 프로세스 유지
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        print("\n⛔ 종료됨")
        sys.exit(0)
