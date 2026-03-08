"""
Claude Telegram Bot
claude CLI를 subprocess로 호출 → Claude Code 구독 내에서 실행.
MCP 도구(Asana, Notion, GitHub 등) 자동 사용 가능.

실행: python bot.py
"""

import asyncio
import json
import logging
import os
import subprocess
import tempfile
from pathlib import Path

from telegram import Update
from telegram.ext import (
    Application,
    CommandHandler,
    MessageHandler,
    filters,
    ContextTypes,
)

from tools import agent_list, agent_get_output

# session_bridge 경로 추가
import sys as _sys
_sys.path.insert(0, str(Path.home() / "claude-scientific-skills" / "mcp-servers" / "team-orchestrator"))
try:
    from session_bridge import get_session_by_tg_msg_id, add_inbox
    _BRIDGE_AVAILABLE = True
except ImportError:
    _BRIDGE_AVAILABLE = False
    logger_placeholder = None

# ── 설정 ───────────────────────────────────────────────────────────────────────

SECRETS_FILE = Path.home() / ".claude" / "secrets.json"
secrets = json.loads(SECRETS_FILE.read_text(encoding="utf-8"))

BOT_TOKEN: str = secrets["TELEGRAM_BOT_TOKEN"]
ALLOWED_USER_IDS: list[int] = secrets.get("TELEGRAM_ALLOWED_USER_IDS", [])

# 사용자별 대화 히스토리: {user_id: [{"role": ..., "content": ...}]}
conversations: dict[int, list] = {}

logging.basicConfig(
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    level=logging.INFO,
)
logger = logging.getLogger("claude-bot")

SYSTEM_PROMPT = """당신은 Jahyun의 개인 AI 어시스턴트입니다. Telegram을 통해 대화합니다.

MCP 도구(Asana, Notion, GitHub, Google Calendar 등)와 파일 시스템을 자유롭게 사용하세요.
항상 한국어로 답변하고, Telegram 마크다운/HTML 없이 일반 텍스트로 답변하세요.
도구 실행 결과는 핵심만 요약하고 URL을 포함하세요."""


# ── claude CLI subprocess를 사용한 호출 (구독 사용, API 크레딧 불필요) ─────────────

async def ask_claude(user_id: int, user_message: str) -> str:
    """claude -p CLI로 쿼리. Claude Code 구독 내에서 실행."""
    history = conversations.get(user_id, [])

    # 이전 대화 컨텍스트 포함
    parts = [SYSTEM_PROMPT, ""]
    if history:
        parts.append("=== 이전 대화 ===")
        for msg in history[-6:]:
            role = "사용자" if msg["role"] == "user" else "Claude"
            parts.append(f"{role}: {msg['content'][:800]}")
        parts.append("=== 대화 끝 ===\n")
    parts.append(f"사용자: {user_message}")
    full_prompt = "\n".join(parts)

    try:
        loop = asyncio.get_event_loop()

        def _run():
            # CLAUDECODE 제거 (중첩 세션 차단 우회) + BOT_CALL 플래그 (hook 발동 억제)
            clean_env = {k: v for k, v in os.environ.items() if k != "CLAUDECODE"}
            clean_env["CLAUDE_BOT_CALL"] = "1"
            return subprocess.run(
                ["claude", "-p", full_prompt, "--output-format", "text"],
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=180,
                cwd=str(Path.home()),
                env=clean_env,
                shell=True,
            )

        result = await asyncio.wait_for(
            loop.run_in_executor(None, _run),
            timeout=200,
        )
        response = (result.stdout or "").strip()
        if not response:
            stderr = (result.stderr or "").strip()
            response = f"(응답 없음){f' — {stderr[:200]}' if stderr else ''}"

    except (asyncio.TimeoutError, subprocess.TimeoutExpired):
        response = "응답 시간 초과 (3분). 더 짧은 질문으로 나눠서 요청해 주세요."
    except Exception as e:
        response = f"[오류] {type(e).__name__}: {str(e)[:300]}"

    history.append({"role": "user", "content": user_message})
    history.append({"role": "assistant", "content": response})
    conversations[user_id] = history[-40:]

    return response


# ── Telegram 핸들러 ────────────────────────────────────────────────────────────

def is_allowed(user_id: int) -> bool:
    return not ALLOWED_USER_IDS or user_id in ALLOWED_USER_IDS


async def cmd_start(update: Update, context: ContextTypes.DEFAULT_TYPE):
    if not is_allowed(update.effective_user.id):
        return
    await update.message.reply_text(
        "안녕하세요! Claude AI 어시스턴트입니다.\n\n"
        "Asana, Notion, GitHub, Google Calendar, 파일 작업,\n"
        "에이전트 실행 등을 지원합니다.\n\n"
        "명령어:\n"
        "/clear - 대화 기록 초기화\n"
        "/agents - 에이전트 목록 조회\n"
        "/help - 사용 예시"
    )


async def cmd_clear(update: Update, context: ContextTypes.DEFAULT_TYPE):
    if not is_allowed(update.effective_user.id):
        return
    conversations.pop(update.effective_user.id, None)
    await update.message.reply_text("대화 기록을 초기화했습니다.")


async def cmd_agents(update: Update, context: ContextTypes.DEFAULT_TYPE):
    if not is_allowed(update.effective_user.id):
        return
    await context.bot.send_chat_action(update.effective_chat.id, "typing")
    try:
        agents = agent_list()
        if not agents or (len(agents) == 1 and "error" in agents[0]):
            await update.message.reply_text("실행된 에이전트가 없습니다.")
            return
        icons = {"running": "[실행중]", "done": "[완료]", "failed": "[실패]", "pending": "[대기]"}
        lines = ["최근 에이전트 목록:\n"]
        for a in agents:
            lines.append(f"{icons.get(a['status'], '[?]')} {a['name']} ({a['type']})\n   ID: {a['agent_id']}")
        await update.message.reply_text("\n".join(lines))
    except Exception as e:
        await update.message.reply_text(f"조회 오류: {e}")


async def cmd_help(update: Update, context: ContextTypes.DEFAULT_TYPE):
    if not is_allowed(update.effective_user.id):
        return
    await update.message.reply_text(
        "사용 예시:\n\n"
        "Asana:\n"
        "  오늘 논문 리뷰 태스크 만들어줘\n"
        "  내 Asana 할일 보여줘\n\n"
        "Notion:\n"
        "  Notion에서 실험 프로토콜 찾아줘\n\n"
        "GitHub:\n"
        "  jahyunlee00299/PeakPicker 이슈 목록 보여줘\n\n"
        "구글 캘린더:\n"
        "  오늘 일정 알려줘\n"
        "  내일 오후 2시에 실험 미팅 추가해줘\n\n"
        "파일:\n"
        "  ~/lab-analyses 폴더 보여줘\n\n"
        "기타:\n"
        "  /clear - 대화 초기화\n"
        "  /agents - 에이전트 목록"
    )


async def handle_message(update: Update, context: ContextTypes.DEFAULT_TYPE):
    if not update.message or not update.message.text:
        return

    user_id = update.effective_user.id
    if not is_allowed(user_id):
        await update.message.reply_text("접근 권한이 없습니다.")
        return

    # ── 답장(reply) 메시지 → 세션 inbox 라우팅 ─────────────────────────────
    if _BRIDGE_AVAILABLE and update.message.reply_to_message:
        reply_msg_id = update.message.reply_to_message.message_id
        session = get_session_by_tg_msg_id(reply_msg_id)
        if session:
            add_inbox(session["session_id"], update.message.text)
            cwd = session.get("cwd", "알 수 없음")
            sid_short = session["session_id"][:8]
            await update.message.reply_text(
                f"📨 [{cwd}] 세션에 전달 완료\n"
                f"다음 질문 입력 시 Claude가 자동으로 확인합니다. (세션 {sid_short}...)"
            )
            return
    # ────────────────────────────────────────────────────────────────────────

    is_reply = bool(update.message.reply_to_message)
    status_text = "💬 대화 이어서 응답 생성 중..." if is_reply else "🔄 새 세션 열림 · 응답 생성 중..."

    status_msg = await update.message.reply_text("✅ 수신확인")
    await asyncio.sleep(0.4)
    await status_msg.edit_text(status_text)
    await context.bot.send_chat_action(update.effective_chat.id, "typing")

    FOOTER = "\n\n↩️ 이 메시지에 답장하면 대화가 이어집니다"

    try:
        response = await ask_claude(user_id, update.message.text)

        first_chunk = response[:4096 - len(FOOTER)]
        await status_msg.edit_text(first_chunk + FOOTER)
        if len(response) > 4096 - len(FOOTER):
            for i in range(4096 - len(FOOTER), len(response), 4000):
                await update.message.reply_text(response[i : i + 4000])
                await asyncio.sleep(0.3)

    except Exception as e:
        logger.exception("메시지 처리 오류")
        await status_msg.edit_text(f"❌ 오류 발생: {e}")


# ── 메인 ───────────────────────────────────────────────────────────────────────

PID_FILE = Path(__file__).parent / "bot.pid"


def _kill_old_instance():
    """이전 봇 인스턴스를 PID 파일로 찾아 종료."""
    if not PID_FILE.exists():
        return
    try:
        old_pid = int(PID_FILE.read_text().strip())
        import ctypes
        PROCESS_TERMINATE = 1
        handle = ctypes.windll.kernel32.OpenProcess(PROCESS_TERMINATE, False, old_pid)
        if handle:
            ctypes.windll.kernel32.TerminateProcess(handle, 0)
            ctypes.windll.kernel32.CloseHandle(handle)
            logger.info(f"이전 봇 인스턴스 종료 (PID {old_pid})")
            import time as _t; _t.sleep(1)
    except Exception as e:
        logger.info(f"이전 인스턴스 종료 시도: {e}")
    PID_FILE.unlink(missing_ok=True)


def _write_pid():
    PID_FILE.write_text(str(os.getpid()))


def main():
    _kill_old_instance()
    _write_pid()
    app = Application.builder().token(BOT_TOKEN).build()

    app.add_handler(CommandHandler("start", cmd_start))
    app.add_handler(CommandHandler("clear", cmd_clear))
    app.add_handler(CommandHandler("agents", cmd_agents))
    app.add_handler(CommandHandler("help", cmd_help))
    app.add_handler(MessageHandler(filters.TEXT & ~filters.COMMAND, handle_message))

    logger.info("Claude Telegram Bot 시작 (claude CLI 모드)...")
    app.run_polling(drop_pending_updates=True)


if __name__ == "__main__":
    main()
