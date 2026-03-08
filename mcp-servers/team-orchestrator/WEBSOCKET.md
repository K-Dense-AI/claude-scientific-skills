# Team Orchestrator WebSocket 서버

실시간 에이전트 상태 모니터링을 위한 WebSocket 서버입니다.

## 🚀 빠른 시작

### Windows
```bash
.\launch_ws.bat
```

### macOS / Linux
```bash
bash launch_ws.sh
```

## 📡 WebSocket 엔드포인트

### 기본 연결
```
ws://localhost:8765/ws
```

### 프로젝트별 필터링
```
ws://localhost:8765/ws?project_id=abc123
```

## 📊 메시지 형식

### 수신 (Broadcast)

```json
{
  "timestamp": 1700000000.123,
  "project_id": "abc123",
  "status_counts": {
    "running": 2,
    "pending": 1,
    "done": 3,
    "failed": 0,
    "paused": 0
  },
  "teams": [
    {
      "team_id": "abc123-1",
      "name": "research-team",
      "type": "research",
      "status": "running",
      "elapsed_seconds": 325.4,
      "task": "Literature search for..."
    },
    ...
  ]
}
```

## 🌐 대시보드 사용

### 1. 로컬 HTML 파일 열기
```bash
# 브라우저에서 열기
open notion_dashboard.html    # macOS
start notion_dashboard.html   # Windows
xdg-open notion_dashboard.html # Linux
```

### 2. Notion 페이지에 Embed하기

1. Notion 페이지에서 `/embed` 또는 `/code` 블록 추가
2. 다음 HTML 복사:

```html
<!DOCTYPE html>
<html>
<head>
    <style>
        body { margin: 0; }
        iframe { width: 100%; height: 600px; border: none; }
    </style>
</head>
<body>
    <iframe src="file:///path/to/notion_dashboard.html"></iframe>
</body>
</html>
```

3. 또는 직접 파일 경로 사용:
   - Windows: `file:///C:/Users/Jahyun/claude-scientific-skills/mcp-servers/team-orchestrator/notion_dashboard.html`
   - macOS/Linux: `file:///path/to/notion_dashboard.html`

### 3. 프로젝트별 필터링
대시보드에서 `?project_id=abc123` 파라미터 사용:

```
file:///path/to/notion_dashboard.html?project_id=abc123
```

## 🔍 HTTP 폴링 대체

WebSocket을 사용할 수 없으면 HTTP 폴링으로 대체:

### 전체 상태
```bash
curl http://localhost:8765/status
```

### 프로젝트별 상태
```bash
curl "http://localhost:8765/status?project_id=abc123"
```

### 헬스 체크
```bash
curl http://localhost:8765/health
```

## 📈 성능 지표

- **업데이트 주기**: 100ms (10 Hz)
- **지연 시간**: <100ms
- **대역폭**: ~1KB/메시지
- **동시 연결**: 무제한 (제한 없음)

## 🔧 개발 & 디버깅

### WebSocket 클라이언트 테스트
```bash
python3 -m websockets ws://localhost:8765/ws
```

### 서버 로그 확인
```bash
tail -f ws_server.log
```

### HTTP API 테스트
```bash
# 상태 조회
curl -s http://localhost:8765/status | python3 -m json.tool

# 헬스 체크
curl -s http://localhost:8765/health | python3 -m json.tool
```

## ⚙️ 설정

### 포트 변경
`ws_server.py`의 다음 부분 수정:

```python
if __name__ == "__main__":
    port = 9000  # 변경
    ...
```

### 업데이트 주기 변경
```python
await asyncio.sleep(0.05)  # 50ms로 변경
```

### 호스트 변경 (원격 접속)
```python
uvicorn.run(
    app,
    host="0.0.0.0",  # 모든 IP에서 수락
    port=port,
)
```

## 🚨 트러블슈팅

### "Address already in use" 에러
```bash
# 포트 8765 사용 중인 프로세스 찾기
lsof -i :8765  # macOS/Linux
netstat -ano | findstr :8765  # Windows

# 프로세스 종료
kill -9 <PID>  # macOS/Linux
taskkill /PID <PID> /F  # Windows
```

### WebSocket 연결 안 됨
1. 방화벽 확인: `firewall-cmd --list-all`
2. 포트 열려있는지 확인: `nc -zv localhost 8765`
3. 서버 재시작

### Notion에서 iframe 로드 안 됨
1. 브라우저 콘솔 확인 (F12)
2. CORS 설정 확인
3. 파일 경로 절대경로인지 확인

## 📚 API 레퍼런스

### WebSocket 클라이언트 예제 (JavaScript)

```javascript
const ws = new WebSocket('ws://localhost:8765/ws?project_id=abc123');

ws.onopen = () => {
    console.log('연결됨');
};

ws.onmessage = (event) => {
    const data = JSON.parse(event.data);
    console.log('상태 업데이트:', data);

    // 타이머 표시
    data.teams.forEach(team => {
        if (team.status === 'running') {
            const mins = Math.floor(team.elapsed_seconds / 60);
            const secs = Math.round(team.elapsed_seconds % 60);
            console.log(`🔄 ${team.name}: ${mins}m ${secs}s`);
        }
    });
};

ws.onerror = (error) => {
    console.error('에러:', error);
};

ws.onclose = () => {
    console.log('연결 해제. 3초 후 재연결...');
    setTimeout(() => location.reload(), 3000);
};
```

### Python 클라이언트 예제

```python
import asyncio
import json
import websockets

async def monitor():
    async with websockets.connect('ws://localhost:8765/ws?project_id=abc123') as ws:
        while True:
            data = json.loads(await ws.recv())
            print(f"🔄 실행: {data['status_counts']['running']}")
            print(f"✅ 완료: {data['status_counts']['done']}")

            for team in data['teams']:
                if team['status'] == 'running':
                    elapsed = team['elapsed_seconds']
                    print(f"  ⏱️ {team['name']}: {elapsed:.1f}s")

asyncio.run(monitor())
```

## 📄 라이센스

MIT

## 🤝 기여

이슈 또는 PR 환영합니다!
