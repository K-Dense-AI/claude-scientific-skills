#!/bin/bash
# WebSocket 서버 + ngrok 동시 실행 (macOS/Linux)

echo ""
echo "========================================"
echo "  Team Orchestrator + ngrok 터널링"
echo "========================================"
echo ""

cd "$(dirname "$0")"

# 1. WebSocket 서버를 백그라운드로 시작
echo "[1/2] WebSocket 서버 시작 (포트 8765)..."
python ws_server.py &
WS_PID=$!

# 2. ngrok 터널링 시작
sleep 2
echo ""
echo "[2/2] ngrok 터널 시작 중..."
echo ""
ngrok http 8765 --log=stdout

# ngrok 종료 시 Python도 함께 종료
kill $WS_PID 2>/dev/null
echo ""
echo "WebSocket 서버와 ngrok이 종료되었습니다."
