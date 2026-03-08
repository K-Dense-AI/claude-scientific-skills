"""Google OAuth token 관리 - drive.readonly 토큰 자동 갱신"""
import json
import time
import requests
from pathlib import Path

TOKEN_PATH = Path("~/.claude/gdrive_token.json").expanduser()
CREDS_PATH = Path("~/.claude/gdrive_credentials.json").expanduser()

def get_access_token():
    token = json.loads(TOKEN_PATH.read_text())
    expiry = token.get("expiry_date", 0) / 1000  # ms -> s

    if time.time() < expiry - 60:
        return token["access_token"]

    # 만료됐으면 refresh
    creds_raw = json.loads(CREDS_PATH.read_text())
    creds = creds_raw.get("installed", creds_raw)  # handle nested structure
    resp = requests.post("https://oauth2.googleapis.com/token", data={
        "client_id": creds["client_id"],
        "client_secret": creds["client_secret"],
        "refresh_token": token["refresh_token"],
        "grant_type": "refresh_token",
    })
    resp.raise_for_status()
    new = resp.json()

    token["access_token"] = new["access_token"]
    token["expiry_date"] = int((time.time() + new.get("expires_in", 3600)) * 1000)
    TOKEN_PATH.write_text(json.dumps(token, indent=2))
    return token["access_token"]

def download_sheet_xlsx(file_id: str) -> bytes:
    """Google Sheets를 XLSX로 다운로드 (drive.readonly 스코프 사용)"""
    token = get_access_token()
    url = f"https://www.googleapis.com/drive/v3/files/{file_id}/export"
    resp = requests.get(url, params={"mimeType": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"},
                        headers={"Authorization": f"Bearer {token}"})
    resp.raise_for_status()
    return resp.content
