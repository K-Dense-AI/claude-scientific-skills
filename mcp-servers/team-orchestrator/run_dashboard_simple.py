#!/usr/bin/env python3
"""
간단한 대시보드 실행 — VS Code에서 직접 실행 가능
python run_dashboard_simple.py 또는 VS Code Run 버튼 클릭
"""
import subprocess
import time
import sys
import os

# ngrok 토큰 (실제 토큰)
NGROK_TOKEN = "3AeEdW3mvUq4eSv8pn8XX25ppWn_zGyaccgCnPHuuixPGd56"

def main():
    print("\n" + "="*60)
    print("  🚀 Team Orchestrator + ngrok 시작")
    print("="*60 + "\n")

    # 1. WebSocket 서버 시작
    print("⏳ [1/2] WebSocket 서버 시작 (포트 8765)...")
    script_dir = os.path.dirname(os.path.abspath(__file__))

    try:
        ws_proc = subprocess.Popen(
            [sys.executable, os.path.join(script_dir, "ws_server.py")],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=script_dir
        )
        print("✅ WebSocket 서버 시작됨\n")
    except Exception as e:
        print(f"❌ 실패: {e}\n")
        return

    # 2. ngrok 설정
    time.sleep(2)
    print("⏳ [2/2] ngrok 터널 시작 중...")

    try:
        from pyngrok import ngrok

        # 토큰 설정
        ngrok.set_auth_token(NGROK_TOKEN)

        # 연결
        public_url = ngrok.connect(8765, "http")

        print("\n" + "="*60)
        print("🎉 완료! 공개 URL이 생성되었습니다")
        print("="*60)
        print(f"\n📱 Notion embed URL:")
        print(f"\n   {public_url}\n")
        print("📋 Notion 페이지 설정:")
        print("   1️⃣  Notion 페이지 열기")
        print("   2️⃣  위의 URL 복사")
        print("   3️⃣  [+ 블록 추가] → [Embed]")
        print("   4️⃣  URL 붙여넣기 → 엔터")
        print("   5️⃣  실시간 대시보드가 Notion 안에 나타남! ✨\n")
        print("💡 팁: 이 터미널을 계속 열어두세요 (종료하면 URL 사라짐)")
        print("="*60 + "\n")

        # ngrok 프로세스 유지
        try:
            ngrok_proc = ngrok.get_ngrok_process()
            ngrok_proc.proc.wait()
        except KeyboardInterrupt:
            print("\n⛔ 종료 중...")
            ngrok.kill()
            ws_proc.terminate()
            print("✅ 서버가 종료되었습니다.")

    except ImportError:
        print("❌ pyngrok 라이브러리가 없습니다.")
        print("설치: pip install pyngrok\n")
        ws_proc.terminate()
    except Exception as e:
        print(f"❌ 오류: {e}\n")
        ws_proc.terminate()

if __name__ == "__main__":
    main()
