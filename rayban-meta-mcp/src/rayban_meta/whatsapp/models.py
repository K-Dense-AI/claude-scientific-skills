"""Pydantic models for WhatsApp Cloud API webhook payloads."""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field


class WhatsAppMedia(BaseModel):
    mime_type: str = ""
    sha256: str = ""
    id: str = ""


class WhatsAppText(BaseModel):
    body: str = ""


class WhatsAppMessage(BaseModel):
    from_: str = Field("", alias="from")
    id: str = ""
    timestamp: str = ""
    type: str = "text"
    text: WhatsAppText | None = None
    image: WhatsAppMedia | None = None
    audio: WhatsAppMedia | None = None
    video: WhatsAppMedia | None = None

    model_config = {"populate_by_name": True}


class WhatsAppContact(BaseModel):
    wa_id: str = ""
    profile: dict[str, Any] = {}


class WebhookMetadata(BaseModel):
    display_phone_number: str = ""
    phone_number_id: str = ""


class WebhookValue(BaseModel):
    messaging_product: str = ""
    metadata: WebhookMetadata = WebhookMetadata()
    contacts: list[WhatsAppContact] = []
    messages: list[WhatsAppMessage] = []
    statuses: list[dict[str, Any]] = []


class WebhookChange(BaseModel):
    value: WebhookValue = WebhookValue()
    field: str = ""


class WebhookEntry(BaseModel):
    id: str = ""
    changes: list[WebhookChange] = []


class WebhookPayload(BaseModel):
    object: str = ""
    entry: list[WebhookEntry] = []
