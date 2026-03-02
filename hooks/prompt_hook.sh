#!/bin/bash
# UserPromptSubmit hook: 사용량 보고 + 진행 중인 작업 인터럽트 감지

PENDING='/c/Users/Jahyun/.claude/pending_usage.txt'
TODOS_DIR='/c/Users/Jahyun/.claude/todos'

# 1) 이전 세션 사용량 보고 - 첫 줄에 그대로 출력
if [ -f "$PENDING" ] && [ -s "$PENDING" ]; then
    USAGE=$(cat "$PENDING")
    rm -f "$PENDING"
    echo "=== [사용량 보고] 응답 맨 첫 줄에 아래 텍스트를 그대로 복사해서 출력하세요 (수정 금지) ==="
    echo "$USAGE"
    echo "==="
fi

# 2) 진행 중인 작업 감지 (최근 2분 내 수정된 todos 파일 확인)
if [ -d "$TODOS_DIR" ]; then
    INTERRUPT_MSG=$(python3 -c "
import json, os, time, glob

todos_dir = '/c/Users/Jahyun/.claude/todos'
now = time.time()
two_min_ago = now - 120

# 최근 2분 내 수정된 파일 중 in_progress 항목 확인
recent_files = [
    f for f in glob.glob(todos_dir + '/*.json')
    if os.path.getmtime(f) > two_min_ago
]

if not recent_files:
    exit(0)

# 가장 최근 수정된 파일
latest = max(recent_files, key=os.path.getmtime)

try:
    with open(latest) as f:
        todos = json.load(f)
    active = [t for t in todos if isinstance(t, dict) and t.get('status') == 'in_progress']
    if active:
        print('⚠️ [인터럽트 감지] 진행 중인 작업이 있습니다:')
        for t in active:
            print(f'  - {t[\"content\"]} ({t.get(\"activeForm\", \"\")})')
        print()
        print('[처리 지침] ~/.claude/CLAUDE.md의 중간 인터럽트 처리 절차를 따르세요.')
        print('에이전트 추가 여부와 기존 작업 재개 여부를 판단하고 사용자에게 먼저 알리세요.')
except Exception:
    pass
" 2>/dev/null)

    if [ -n "$INTERRUPT_MSG" ]; then
        echo "$INTERRUPT_MSG"
    fi
fi

exit 0
