@echo off
REM WebSocket 서버 백그라운드 실행 (Windows)

setlocal enabledelayedexpansion

REM 현재 디렉토리로 이동
cd /d "%~dp0"

echo.
echo 📦 의존성 설치 중...
pip install -q -r requirements.txt

echo.
echo 🚀 WebSocket 서버 시작...
echo    📍 WS: ws://localhost:8765/ws
echo    📍 HTTP: http://localhost:8765/status
echo    📍 Health: http://localhost:8765/health
echo.

REM WebSocket 서버를 새 CMD 창에서 실행
start "Team Orchestrator WebSocket" cmd /k python ws_server.py

timeout /t 2

echo ✅ WebSocket 서버 시작됨
echo.
echo 🌐 대시보드 열기:
echo    → 브라우저에서 다음 주소 입력: file:///%CD%\notion_dashboard.html
echo       (또는 Ctrl+Alt+O 누르기)
echo.
echo 💡 팁: dashboard.html 파일을 Notion Embed 블록에 드래그앤드롭하세요
echo.
pause
