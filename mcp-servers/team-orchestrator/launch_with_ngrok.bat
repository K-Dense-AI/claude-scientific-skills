@echo off
REM WebSocket 서버 + ngrok 동시 실행 (Windows)

echo.
echo ========================================
echo   Team Orchestrator + ngrok 터널링
echo ========================================
echo.

cd /d "%~dp0"

REM 1. WebSocket 서버를 백그라운드로 시작
echo [1/2] WebSocket 서버 시작 (포트 8765)...
start /B python ws_server.py

REM 2. ngrok 터널링 시작
timeout /t 2 /nobreak
echo.
echo [2/2] ngrok 터널 시작 중...
echo.
ngrok http 8765 --log=stdout

REM ngrok 종료 시 Python도 함께 종료
taskkill /F /IM python.exe >nul 2>&1
echo.
echo WebSocket 서버와 ngrok이 종료되었습니다.
pause
