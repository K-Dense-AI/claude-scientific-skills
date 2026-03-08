@echo off
setlocal enabledelayedexpansion

REM 스크립트 디렉토리로 이동
cd /d "%~dp0"

REM 로그 파일 설정
set LOG_FILE=%~dp0bot_startup.log
echo [%date% %time%] Telegram Bot 시작 시도 >> "%LOG_FILE%"

REM Python 경로 확인
set PYTHON_PATH=C:\Users\Jahyun\anaconda3\pythonw.exe
if not exist "%PYTHON_PATH%" (
    echo [%date% %time%] 에러: Python을 찾을 수 없음 - %PYTHON_PATH% >> "%LOG_FILE%"
    exit /b 1
)
echo [%date% %time%] Python 경로 확인됨: %PYTHON_PATH% >> "%LOG_FILE%"

REM bot.py 확인
if not exist "bot.py" (
    echo [%date% %time%] 에러: bot.py를 찾을 수 없음 >> "%LOG_FILE%"
    exit /b 1
)
echo [%date% %time%] bot.py 확인됨 >> "%LOG_FILE%"

REM 프로세스 시작
start "" "%PYTHON_PATH%" bot.py
echo [%date% %time%] bot.py 프로세스 시작됨 >> "%LOG_FILE%"

timeout /t 2 /nobreak >nul
endlocal
