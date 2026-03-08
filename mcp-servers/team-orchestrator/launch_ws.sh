#!/bin/bash
# WebSocket 서버 백그라운드 실행 스크립트

set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

# 의존성 설치
echo "📦 의존성 설치 중..."
pip install -q -r requirements.txt

# WebSocket 서버 실행 (백그라운드)
echo ""
echo "🚀 WebSocket 서버 시작..."
echo "   📍 WS: ws://localhost:8765/ws"
echo "   📍 HTTP: http://localhost:8765/status"
echo "   📍 Health: http://localhost:8765/health"
echo ""

python3 ws_server.py &
WS_PID=$!

echo "✅ WebSocket 서버 시작 (PID: $WS_PID)"
echo ""
echo "🌐 Notion 대시보드 열기:"
echo "   → file://$(pwd)/notion_dashboard.html"
echo ""
echo "📊 Notion 페이지에 Embed하기:"
echo "   → <iframe src=\"file://$(pwd)/notion_dashboard.html\" width=\"100%\" height=\"600\"></iframe>"
echo ""
echo "⏹️  종료하려면: kill $WS_PID"
echo ""

# SIGTERM 수신 시 정리
trap "echo '🛑 WebSocket 서버 종료...' && kill $WS_PID 2>/dev/null || true; exit 0" SIGTERM SIGINT

# 백그라운드 프로세스 대기
wait $WS_PID
