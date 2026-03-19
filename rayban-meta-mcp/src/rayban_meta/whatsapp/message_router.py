"""Central message router – dispatches WhatsApp messages to LLM + tools.

Flow:
  1. Extract message type (text / audio / image)
  2. Audio → transcribe → treat as text
  3. Image → download + store + send to vision LLM
  4. Text → classify intent → dispatch to tool or conversation
  5. Record activity for RAG learning
  6. Send response back via WhatsApp
"""

from __future__ import annotations

import logging
from datetime import datetime

from rayban_meta.config import settings
from rayban_meta.llm.router import LLMRouter
from rayban_meta.llm.transcription import transcribe
from rayban_meta.memory.sqlite_store import SQLiteStore
from rayban_meta.tools.base import ToolResult
from rayban_meta.tools.registry import (
    INTENT_CATEGORIES,
    ToolRegistry,
    detect_save_command,
)
from rayban_meta.whatsapp.api import WhatsAppClient

logger = logging.getLogger(__name__)

# Singletons (initialized on first call)
_wa_client: WhatsAppClient | None = None
_llm_router: LLMRouter | None = None
_tool_registry: ToolRegistry | None = None


def _get_wa_client() -> WhatsAppClient:
    global _wa_client
    if _wa_client is None:
        _wa_client = WhatsAppClient()
    return _wa_client


def _get_llm_router() -> LLMRouter:
    global _llm_router
    if _llm_router is None:
        _llm_router = LLMRouter()
    return _llm_router


def _get_tool_registry() -> ToolRegistry:
    global _tool_registry
    if _tool_registry is None:
        from rayban_meta.tools.calendar import CalendarTool
        from rayban_meta.tools.knowledge_base import AdaptiveKnowledgeTool
        from rayban_meta.tools.notes import SaveNoteTool, SearchKnowledgeTool
        from rayban_meta.tools.vision import VisionTool
        from rayban_meta.tools.web_search import WebSearchTool

        _tool_registry = ToolRegistry()
        _tool_registry.register(WebSearchTool())
        _tool_registry.register(CalendarTool())
        _tool_registry.register(SaveNoteTool())
        _tool_registry.register(SearchKnowledgeTool())
        _tool_registry.register(VisionTool())
        _tool_registry.register(AdaptiveKnowledgeTool())
    return _tool_registry


SYSTEM_PROMPT = """You are a concise AI assistant running on Ray-Ban Meta smart glasses via WhatsApp.
Keep answers SHORT (1-3 sentences max) because the glasses display truncates long text.
Be direct and helpful. The user is speaking to you through voice commands.
Current time: {now}

Available tools: {tools}
User's frequent actions: {frequent}
"""


async def handle_message(app, raw_message: dict) -> None:
    """Process a single incoming WhatsApp message."""
    store: SQLiteStore = app.state.store
    wa = _get_wa_client()
    llm_router = _get_llm_router()
    registry = _get_tool_registry()

    msg_type = raw_message.get("type", "text")
    phone = raw_message.get("from", "")
    msg_id = raw_message.get("id", "")

    # Mark as read
    await wa.mark_as_read(msg_id)

    llm = llm_router.get()
    response_text = ""

    try:
        if msg_type == "audio":
            response_text = await _handle_audio(raw_message, phone, store, llm, registry, wa)
        elif msg_type == "image":
            response_text = await _handle_image(raw_message, phone, store, llm, wa)
        elif msg_type == "text":
            text = raw_message.get("text", {}).get("body", "")
            response_text = await _handle_text(text, phone, store, llm, registry)
        else:
            response_text = f"I received a {msg_type} message but I can only handle text, voice, and images for now."

    except Exception:
        logger.exception("Error handling message")
        response_text = "Sorry, something went wrong. Please try again."

    # Send response
    if response_text:
        await wa.send_text(phone, response_text)


async def _handle_audio(raw_message: dict, phone: str, store, llm, registry, wa) -> str:
    """Download audio → transcribe → process as text."""
    media_id = raw_message.get("audio", {}).get("id", "")
    if not media_id:
        return "Couldn't process the voice message."

    audio_bytes = await wa.download_media(media_id)
    text = await transcribe(audio_bytes)
    logger.info("Transcribed audio from %s: %s", phone, text[:100])

    # Store the transcription
    await store.add_message(phone, "user", f"[voice] {text}", media_id=media_id, media_type="audio")

    response = await _handle_text(text, phone, store, llm, registry)

    await store.add_message(phone, "assistant", response)
    return response


async def _handle_image(raw_message: dict, phone: str, store, llm, wa) -> str:
    """Download image → analyze with vision LLM."""
    media_id = raw_message.get("image", {}).get("id", "")
    caption = raw_message.get("image", {}).get("caption", "")
    if not media_id:
        return "Couldn't process the image."

    image_bytes = await wa.download_media(media_id)
    mime_type = raw_message.get("image", {}).get("mime_type", "image/jpeg")

    # Store image
    await store.store_media(media_id, image_bytes, mime_type, phone)
    await store.add_message(phone, "user", f"[photo] {caption or 'Image sent'}", media_id=media_id, media_type="image")

    # Analyze with vision
    prompt = caption or "Describe what you see in this image. Be concise."
    history = await _build_context(phone, store)
    response = await llm.complete_with_vision(
        messages=history + [{"role": "user", "content": prompt}],
        images=[image_bytes],
        system="You are a concise assistant on smart glasses. Describe images in 1-2 sentences.",
        max_tokens=150,
    )

    await store.add_message(phone, "assistant", response)
    return response


async def _handle_text(text: str, phone: str, store, llm, registry) -> str:
    """Classify intent → dispatch to tool or general conversation."""
    if not text.strip():
        return "I didn't catch that. Could you repeat?"

    # Check for explicit save command first
    is_save, folder = detect_save_command(text)
    if is_save:
        tool = registry.get("save")
        if tool:
            result = await tool.execute(text, store, llm)
            await store.record_activity(phone, "save", text)
            await store.add_message(phone, "user", text, intent="save")
            await store.add_message(phone, "assistant", result.text)
            return result.text

    # Classify intent using LLM
    # First check user patterns for faster classification
    frequent = await store.get_frequent_intents(phone, limit=5)
    intent = await llm.classify(text, INTENT_CATEGORIES)
    logger.info("Intent for '%s': %s", text[:50], intent)

    # Record activity for learning
    await store.record_activity(phone, intent, text)
    await store.add_message(phone, "user", text, intent=intent)

    # Dispatch to tool
    tool = registry.match_intent(intent)
    if tool and intent != "conversation":
        result = await tool.execute(text, store, llm)

        # Auto-save to folder if tool requests it
        if result.save_to_folder:
            await store.add_knowledge(result.text, source=tool.name, folder=result.save_to_folder)

        await store.add_message(phone, "assistant", result.text)
        return result.text

    # General conversation
    history = await _build_context(phone, store)
    frequent_str = ", ".join(f"{f['intent']}({f['count']}x)" for f in frequent) if frequent else "none yet"
    tools_str = registry.tools_description()

    system = SYSTEM_PROMPT.format(
        now=datetime.now().strftime("%Y-%m-%d %H:%M"),
        tools=tools_str,
        frequent=frequent_str,
    )

    response = await llm.complete(
        messages=history + [{"role": "user", "content": text}],
        system=system,
        max_tokens=200,
    )
    await store.add_message(phone, "assistant", response)
    return response


async def _build_context(phone: str, store) -> list[dict]:
    """Build conversation context from recent messages."""
    history = await store.get_history(phone, limit=10)
    return [{"role": m.role, "content": m.content} for m in history]
