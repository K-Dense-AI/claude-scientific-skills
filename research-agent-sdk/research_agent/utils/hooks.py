"""Pre/Post tool use hooks for logging and safety."""

from __future__ import annotations

from datetime import datetime

BLOCKED_COMMANDS = {"rm -rf /", "rm -rf ~", "rm -rf *", ":(){ :|:& };:"}


async def pre_tool_hook(input_data: dict, tool_use_id: str, context: dict) -> dict:
    """Log tool invocations and block dangerous commands."""
    tool_name = input_data.get("tool_name", "unknown")
    timestamp = datetime.now().strftime("%H:%M:%S")
    print(f"  [{timestamp}] -> {tool_name}", flush=True)

    # Block dangerous bash commands
    if tool_name == "Bash":
        command = input_data.get("tool_input", {}).get("command", "")
        for blocked in BLOCKED_COMMANDS:
            if blocked in command:
                return {
                    "hookSpecificOutput": {
                        "permissionDecision": "deny",
                        "permissionDecisionReason": f"Blocked dangerous command: {blocked}",
                    }
                }

    return {}


async def post_tool_hook(input_data: dict, tool_use_id: str, context: dict) -> dict:
    """Post-execution logging (no-op for now)."""
    return {}
