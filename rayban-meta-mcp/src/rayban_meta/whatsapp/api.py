"""WhatsApp Cloud API client – send messages, download media."""

from __future__ import annotations

import asyncio
import logging

import httpx

from rayban_meta.config import settings
from rayban_meta.utils.chunking import chunk_for_glasses

logger = logging.getLogger(__name__)

GRAPH_URL = "https://graph.facebook.com/v20.0"


class WhatsAppClient:
    """Async wrapper around the WhatsApp Cloud API."""

    def __init__(self) -> None:
        self._http = httpx.AsyncClient(
            timeout=30.0,
            headers={"Authorization": f"Bearer {settings.whatsapp_auth_token}"},
        )

    # ── Sending ──────────────────────────────────────────────

    async def send_text(self, to: str, text: str) -> None:
        """Send a text message, auto-chunking for the glasses display."""
        chunks = chunk_for_glasses(text)
        for i, chunk in enumerate(chunks):
            await self._send_raw_text(to, chunk)
            if i < len(chunks) - 1:
                await asyncio.sleep(2)

    async def _send_raw_text(self, to: str, text: str) -> dict:
        url = f"{GRAPH_URL}/{settings.whatsapp_phone_number_id}/messages"
        payload = {
            "messaging_product": "whatsapp",
            "recipient_type": "individual",
            "to": to,
            "type": "text",
            "text": {"preview_url": False, "body": text},
        }
        resp = await self._http.post(url, json=payload)
        resp.raise_for_status()
        return resp.json()

    async def send_audio(self, to: str, audio_url: str) -> dict:
        """Send an audio message (for TTS responses)."""
        url = f"{GRAPH_URL}/{settings.whatsapp_phone_number_id}/messages"
        payload = {
            "messaging_product": "whatsapp",
            "to": to,
            "type": "audio",
            "audio": {"link": audio_url},
        }
        resp = await self._http.post(url, json=payload)
        resp.raise_for_status()
        return resp.json()

    # ── Media download ───────────────────────────────────────

    async def download_media(self, media_id: str) -> bytes:
        """Two-step media download: get URL, then fetch binary."""
        # Step 1: resolve download URL
        meta_resp = await self._http.get(f"{GRAPH_URL}/{media_id}")
        meta_resp.raise_for_status()
        download_url = meta_resp.json()["url"]

        # Step 2: download binary
        data_resp = await self._http.get(download_url)
        data_resp.raise_for_status()
        return data_resp.content

    # ── Read receipts ────────────────────────────────────────

    async def mark_as_read(self, message_id: str) -> None:
        url = f"{GRAPH_URL}/{settings.whatsapp_phone_number_id}/messages"
        payload = {
            "messaging_product": "whatsapp",
            "status": "read",
            "message_id": message_id,
        }
        try:
            await self._http.post(url, json=payload)
        except Exception:
            logger.warning("Failed to mark message %s as read", message_id)

    async def close(self) -> None:
        await self._http.aclose()
