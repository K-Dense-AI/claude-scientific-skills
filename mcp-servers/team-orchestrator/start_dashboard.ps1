# WebSocket 서버 + ngrok 자동 실행 (PowerShell)
# 관리자 권한 필요 없음

Write-Host "`n============================================================" -ForegroundColor Cyan
Write-Host "   Team Orchestrator + ngrok 터널" -ForegroundColor Cyan
Write-Host "============================================================`n" -ForegroundColor Cyan

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$pythonExe = "python"  # conda 환경에서 활성화된 python 사용

# 1. WebSocket 서버 시작
Write-Host "[1/2] WebSocket 서버 시작 중 (포트 8765)..." -ForegroundColor Yellow
$wsProcess = Start-Process -FilePath $pythonExe -ArgumentList "ws_server.py" `
    -WorkingDirectory $scriptDir `
    -NoNewWindow -PassThru

Write-Host "✅ WebSocket 서버 시작됨 (PID: $($wsProcess.Id))" -ForegroundColor Green
Start-Sleep -Seconds 2

# 2. ngrok 토큰 설정 및 터널 시작
Write-Host "`n[2/2] ngrok 설정 및 터널 시작 중..." -ForegroundColor Yellow

$ngrokToken = "TTNX2VY77S"
$pythonCode = @"
import sys
import time
from pyngrok import ngrok

try:
    ngrok.set_auth_token('$ngrokToken')
    print('✅ ngrok 토큰 설정 완료')

    time.sleep(1)

    public_url = ngrok.connect(8765, 'http')
    print('\n' + '='*60)
    print('🎉 공개 URL 생성됨!')
    print('='*60)
    print(f'\n📱 Notion embed URL:')
    print(f'   {public_url}\n')
    print('📋 Notion 페이지에서:')
    print('   1. + 블록 추가 → Embed')
    print(f'   2. URL 입력: {public_url}')
    print('   3. 엔터 → 완성! ✨\n')
    print('💡 종료하려면 Ctrl+C 누르세요')
    print('='*60 + '\n')

    # ngrok 프로세스 유지
    ngrok_proc = ngrok.get_ngrok_process()
    ngrok_proc.proc.wait()

except KeyboardInterrupt:
    print('\n⛔ 종료 중...')
    ngrok.kill()
    sys.exit(0)
except Exception as e:
    print(f'❌ 오류: {e}')
    sys.exit(1)
"@

# Python 인라인 코드 실행
$pythonProcess = Start-Process -FilePath $pythonExe -ArgumentList "-c", $pythonCode `
    -NoNewWindow -Wait -PassThru

# 정리
$wsProcess.Kill()
Write-Host "`n✅ 서버와 터널이 종료되었습니다." -ForegroundColor Green
