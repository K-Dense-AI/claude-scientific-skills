#!/bin/bash
# install.sh — claude-scientific-skills 설치 스크립트
# 사용법: bash install.sh
# 다른 컴퓨터에서 git clone 후 이 스크립트를 실행하세요.

set -e

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLAUDE_DIR="${HOME}/.claude"

echo "=== claude-scientific-skills 설치 ==="
echo "레포: $REPO_DIR"
echo "Claude 설정 디렉토리: $CLAUDE_DIR"
echo ""

# 1) ~/.claude 디렉토리 생성
mkdir -p "$CLAUDE_DIR/hooks"
mkdir -p "$CLAUDE_DIR/skills"

# 2) 훅 파일 복사
echo "[1/4] 훅 파일 복사 중..."
cp "$REPO_DIR/hooks/stop_hook.sh"              "$CLAUDE_DIR/hooks/"
cp "$REPO_DIR/hooks/prompt_hook.sh"            "$CLAUDE_DIR/hooks/"
cp "$REPO_DIR/hooks/skill-activation-prompt.sh" "$CLAUDE_DIR/hooks/"
cp "$REPO_DIR/hooks/skill-activation-prompt.py" "$CLAUDE_DIR/hooks/"
chmod +x "$CLAUDE_DIR/hooks/"*.sh
echo "    완료: $CLAUDE_DIR/hooks/"

# 3) 스킬 파일 복사 (scientific-skills/)
echo "[2/4] 스킬 파일 복사 중..."
if [ -d "$REPO_DIR/scientific-skills" ]; then
    cp -r "$REPO_DIR/scientific-skills/"* "$CLAUDE_DIR/skills/" 2>/dev/null || true
    echo "    완료: $CLAUDE_DIR/skills/"
else
    echo "    건너뜀: scientific-skills/ 디렉토리 없음"
fi

# 4) skill-rules.json 복사
if [ -f "$REPO_DIR/scientific-skills/skill-rules.json" ]; then
    cp "$REPO_DIR/scientific-skills/skill-rules.json" "$CLAUDE_DIR/skills/"
    echo "    skill-rules.json 복사 완료"
fi

# 5) settings.json에 훅 등록 (기존 설정 보존)
echo "[3/4] settings.json 훅 등록 중..."
SETTINGS="$CLAUDE_DIR/settings.json"

if [ ! -f "$SETTINGS" ]; then
    echo '{}' > "$SETTINGS"
fi

# Python으로 settings.json 업데이트 (hooks 섹션만 추가/갱신)
python3 - "$SETTINGS" "$CLAUDE_DIR" << 'PYEOF'
import json, sys, os

settings_path = sys.argv[1]
claude_dir = sys.argv[2]

# 경로를 bash 스타일로 변환 (Windows Git Bash용)
def to_bash_path(p):
    # C:\Users\foo\.claude -> /c/Users/foo/.claude
    p = p.replace('\\', '/')
    if len(p) >= 2 and p[1] == ':':
        drive = p[0].lower()
        p = '/' + drive + p[2:]
    return p

hooks_dir = to_bash_path(claude_dir) + "/hooks"

with open(settings_path) as f:
    settings = json.load(f)

settings["hooks"] = {
    "Stop": [{
        "hooks": [{
            "type": "command",
            "command": f"bash '{hooks_dir}/stop_hook.sh'"
        }]
    }],
    "UserPromptSubmit": [{
        "hooks": [
            {
                "type": "command",
                "command": f"bash '{hooks_dir}/prompt_hook.sh'"
            },
            {
                "type": "command",
                "command": f"bash '{hooks_dir}/skill-activation-prompt.sh'"
            }
        ]
    }]
}

with open(settings_path, 'w', encoding='utf-8') as f:
    json.dump(settings, f, indent=2, ensure_ascii=False)

print(f"    settings.json 업데이트 완료: {settings_path}")
PYEOF

# 6) ccusage 설치 확인
echo "[4/4] ccusage 설치 확인 중..."
if command -v ccusage &>/dev/null; then
    echo "    ccusage 이미 설치됨: $(ccusage --version 2>/dev/null || echo '버전 확인 불가')"
else
    echo "    ccusage 미설치. 설치 중..."
    npm install -g ccusage 2>/dev/null || echo "    경고: npm 설치 실패. 수동으로 'npm install -g ccusage' 실행하세요."
fi

echo ""
echo "=== 설치 완료 ==="
echo "Claude Code를 재시작하면 토큰 사용량 보고가 활성화됩니다."
echo ""
echo "추가 설정 필요:"
echo "  - secrets.json: ASANA_PAT 토큰 직접 입력 필요"
echo "    echo '{\"ASANA_PAT\": \"your_token\"}' > $CLAUDE_DIR/secrets.json"
echo "  - CLAUDE.md: 필요 시 $CLAUDE_DIR/CLAUDE.md 수동 복사"
