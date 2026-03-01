"""Handle streamed messages from the Agent SDK."""

from __future__ import annotations


def handle_message(msg: object) -> None:
    """Process and print a message from the SDK.

    Handles AssistantMessage (text blocks) and ResultMessage (cost info).
    Uses duck typing to avoid hard dependency on SDK types at import time.
    """
    msg_type = type(msg).__name__

    if msg_type == "AssistantMessage":
        for block in getattr(msg, "content", []):
            block_type = type(block).__name__
            if block_type == "TextBlock":
                print(getattr(block, "text", ""), end="", flush=True)

    elif msg_type == "ResultMessage":
        result = getattr(msg, "result", None)
        if result:
            print(f"\n{result}")
        cost = getattr(msg, "total_cost_usd", None)
        if cost:
            print(f"[Cost: ${cost:.4f}]")
