#!/usr/bin/env python3
"""
Publish or schedule content to social media platforms via Composio MCP.

Reads a content file (JSON) specifying the platform, caption, media, and
schedule time, then constructs the appropriate Composio MCP tool call
for publishing.

Usage:
    python publish_content.py \
        --content-file post.json \
        --platform instagram \
        --schedule "2026-04-01T10:00:00-04:00" \
        --output publish_log.json

    python publish_content.py \
        --content-file post.json \
        --platform linkedin \
        --schedule now \
        --output publish_log.json
"""

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

# --- Composio MCP tool mappings ---
COMPOSIO_TOOLS = {
    "instagram": {
        "tool_name": "instagram-automation",
        "capabilities": ["feed_post", "reel", "story", "carousel"],
        "required_fields": ["caption", "media_paths"],
        "optional_fields": ["alt_text", "first_comment", "location_tag", "cover_image"],
        "media_requirements": {
            "feed_post": {"formats": ["jpg", "png"], "max_size_mb": 8, "aspect_ratios": ["1:1", "4:5", "1.91:1"]},
            "reel": {"formats": ["mp4", "mov"], "max_size_mb": 250, "max_duration_sec": 90},
            "story": {"formats": ["jpg", "png", "mp4"], "max_size_mb": 30, "aspect_ratio": "9:16"},
            "carousel": {"formats": ["jpg", "png", "mp4"], "max_slides": 20, "max_size_mb": 8},
        },
    },
    "facebook": {
        "tool_name": "facebook",
        "capabilities": ["page_post", "story", "reel", "link_share"],
        "required_fields": ["caption"],
        "optional_fields": ["media_paths", "link_url", "link_title", "link_description"],
        "media_requirements": {
            "page_post": {"formats": ["jpg", "png", "gif"], "max_size_mb": 10},
            "reel": {"formats": ["mp4", "mov"], "max_size_mb": 1000, "max_duration_sec": 90},
            "story": {"formats": ["jpg", "png", "mp4"], "max_size_mb": 30},
        },
    },
    "linkedin": {
        "tool_name": "linkedin-automation",
        "capabilities": ["text_post", "image_post", "document_share", "article"],
        "required_fields": ["caption"],
        "optional_fields": ["media_paths", "document_path", "article_url", "article_title"],
        "media_requirements": {
            "image_post": {"formats": ["jpg", "png", "gif"], "max_size_mb": 10},
            "document_share": {"formats": ["pdf", "pptx", "docx"], "max_size_mb": 100},
        },
    },
    "tiktok": {
        "tool_name": "tiktok-automation",
        "capabilities": ["video_upload"],
        "required_fields": ["caption", "media_paths"],
        "optional_fields": ["cover_image", "sound_id", "allow_comments", "allow_duet", "allow_stitch"],
        "media_requirements": {
            "video_upload": {"formats": ["mp4", "mov", "webm"], "max_size_mb": 287, "max_duration_sec": 600},
        },
    },
}


def load_content(path: str) -> dict:
    """Load content from JSON file."""
    with open(path, "r") as f:
        return json.load(f)


def validate_content(content: dict, platform: str) -> list:
    """Validate content against platform requirements. Returns list of issues."""
    issues = []
    tool_config = COMPOSIO_TOOLS.get(platform)

    if not tool_config:
        issues.append(f"Unsupported platform: {platform}")
        return issues

    # Check required fields
    for field in tool_config["required_fields"]:
        if field not in content or not content[field]:
            issues.append(f"Missing required field: {field}")

    # Validate post type
    post_type = content.get("post_type", "feed_post")
    if post_type not in tool_config["capabilities"]:
        issues.append(
            f"Post type '{post_type}' not supported on {platform}. "
            f"Supported: {', '.join(tool_config['capabilities'])}"
        )

    # Validate media files exist (if provided)
    media_paths = content.get("media_paths", [])
    for media_path in media_paths:
        p = Path(media_path)
        if not p.exists():
            issues.append(f"Media file not found: {media_path}")
        else:
            # Check file extension
            ext = p.suffix.lower().lstrip(".")
            media_reqs = tool_config["media_requirements"].get(post_type, {})
            allowed_formats = media_reqs.get("formats", [])
            if allowed_formats and ext not in allowed_formats:
                issues.append(
                    f"File format '{ext}' not supported for {post_type} on {platform}. "
                    f"Allowed: {', '.join(allowed_formats)}"
                )

    # Validate caption length
    caption = content.get("caption", "")
    platform_caption_limits = {
        "instagram": 2200,
        "facebook": 63206,
        "linkedin": 3000,
        "tiktok": 4000,
    }
    limit = platform_caption_limits.get(platform, 5000)
    if len(caption) > limit:
        issues.append(f"Caption exceeds {platform} limit of {limit} characters (got {len(caption)}).")

    return issues


def build_composio_payload(content: dict, platform: str, schedule_time: str) -> dict:
    """Build the payload for the Composio MCP tool call."""
    tool_config = COMPOSIO_TOOLS[platform]

    payload = {
        "tool": tool_config["tool_name"],
        "action": "publish" if schedule_time == "now" else "schedule",
        "platform": platform,
        "post_type": content.get("post_type", "feed_post"),
        "caption": content.get("caption", ""),
    }

    # Add media
    if content.get("media_paths"):
        payload["media_paths"] = content["media_paths"]

    # Add optional fields
    for field in tool_config["optional_fields"]:
        if field in content and content[field]:
            payload[field] = content[field]

    # Add schedule time
    if schedule_time != "now":
        payload["schedule_time"] = schedule_time

    return payload


def create_publish_log(
    content: dict,
    platform: str,
    schedule_time: str,
    payload: dict,
    status: str,
    error: str = None,
) -> dict:
    """Create a publish log entry."""
    return {
        "post_id": f"{platform}_{datetime.now(timezone.utc).strftime('%Y%m%d%H%M%S')}",
        "platform": platform,
        "post_type": content.get("post_type", "feed_post"),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "schedule_time": schedule_time,
        "status": status,
        "caption_preview": content.get("caption", "")[:100] + ("..." if len(content.get("caption", "")) > 100 else ""),
        "media_count": len(content.get("media_paths", [])),
        "composio_tool": COMPOSIO_TOOLS.get(platform, {}).get("tool_name", "unknown"),
        "composio_payload": payload,
        "response_data": None,
        "error_message": error,
    }


def main():
    parser = argparse.ArgumentParser(description="Publish or schedule content via Composio MCP.")
    parser.add_argument("--content-file", required=True, help="Path to content JSON file.")
    parser.add_argument("--platform", required=True,
                        choices=["instagram", "facebook", "linkedin", "tiktok"],
                        help="Target platform.")
    parser.add_argument("--schedule", required=True,
                        help="Schedule time as ISO datetime (e.g., '2026-04-01T10:00:00-04:00') or 'now' for immediate publish.")
    parser.add_argument("--output", required=True, help="Output log file path (JSON).")
    parser.add_argument("--dry-run", action="store_true",
                        help="Validate and build payload without actually publishing.")

    args = parser.parse_args()

    # Load content
    content = load_content(args.content_file)

    # Override platform if specified
    content["platform"] = args.platform

    # Validate
    issues = validate_content(content, args.platform)
    if issues:
        print("Validation issues found:", file=sys.stderr)
        for issue in issues:
            print(f"  - {issue}", file=sys.stderr)
        if not args.dry_run:
            log = create_publish_log(content, args.platform, args.schedule, {}, "failed", "; ".join(issues))
            with open(args.output, "w") as f:
                json.dump(log, f, indent=2)
            print(f"\nFailed. Log written to {args.output}", file=sys.stderr)
            sys.exit(1)

    # Build payload
    payload = build_composio_payload(content, args.platform, args.schedule)

    if args.dry_run:
        print("DRY RUN - No content will be published.\n")
        print("Composio MCP payload:")
        print(json.dumps(payload, indent=2))

        log = create_publish_log(content, args.platform, args.schedule, payload, "dry_run")
    else:
        # In production, this is where the Composio MCP tool call would be made.
        # The script constructs the payload and documents the integration point.
        #
        # Integration pattern:
        #   1. The payload is passed to the appropriate Composio MCP tool.
        #   2. For Instagram: use composio.instagram-automation with action=publish or schedule.
        #   3. For Facebook: use composio.facebook with action=publish or schedule.
        #   4. For LinkedIn: use composio.linkedin-automation with action=publish or schedule.
        #   5. For TikTok: use composio.tiktok-automation with action=publish or schedule.
        #
        # The MCP tool returns a response with post_id, status, and platform-specific data.
        # That response is captured in the publish log.

        print(f"Publishing to {args.platform}...")
        print(f"Tool: {COMPOSIO_TOOLS[args.platform]['tool_name']}")
        print(f"Schedule: {'Immediate' if args.schedule == 'now' else args.schedule}")

        # Placeholder for actual MCP integration
        print("\n[Composio MCP integration point]")
        print("Payload constructed. Pass to Composio MCP tool for execution.")
        print(json.dumps(payload, indent=2))

        status = "scheduled" if args.schedule != "now" else "pending_publish"
        log = create_publish_log(content, args.platform, args.schedule, payload, status)

    # Write log
    output_path = Path(args.output)

    # Append to existing log if it exists
    existing_logs = []
    if output_path.exists():
        try:
            with open(output_path, "r") as f:
                existing = json.load(f)
                if isinstance(existing, list):
                    existing_logs = existing
                else:
                    existing_logs = [existing]
        except (json.JSONDecodeError, KeyError):
            existing_logs = []

    existing_logs.append(log)

    with open(output_path, "w") as f:
        json.dump(existing_logs if len(existing_logs) > 1 else log, f, indent=2)

    print(f"\nPublish log written to {args.output}")
    print(f"Status: {log['status']}")


if __name__ == "__main__":
    main()
