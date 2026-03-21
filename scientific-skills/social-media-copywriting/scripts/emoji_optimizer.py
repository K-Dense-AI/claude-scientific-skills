#!/usr/bin/env python3
"""
Add strategic emoji placement to social media captions based on
platform norms and brand voice.

Usage:
    python emoji_optimizer.py \
        --input "Check out our new collection" \
        --platform instagram \
        --brand-voice playful \
        --output optimized.json
"""

import argparse
import json
import re
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Emoji databases by category and sentiment
# ---------------------------------------------------------------------------
EMOJI_BY_CATEGORY = {
    "positive": ["✨", "🎉", "🙌", "💯", "🔥", "⭐", "🌟", "💪", "👏", "🎊"],
    "love": ["❤️", "💕", "🥰", "💖", "💗", "😍", "🤍", "💛", "💚", "💙"],
    "pointing": ["👇", "👆", "👉", "➡️", "⬇️", "📌", "🔗"],
    "celebration": ["🎉", "🥳", "🎊", "🍾", "🏆", "🎯", "💐", "🎈"],
    "business": ["📈", "💼", "🎯", "✅", "📊", "💡", "🔑", "📋"],
    "education": ["📚", "💡", "🧠", "📝", "✏️", "🎓", "📖", "💭"],
    "food": ["🍕", "🍳", "☕", "🍽️", "🥗", "🧁", "🍷", "🥑"],
    "nature": ["🌿", "🌸", "🌻", "🌊", "🌅", "🍃", "☀️", "🌙"],
    "tech": ["💻", "📱", "⚡", "🤖", "🔧", "⚙️", "🛠️", "📡"],
    "health": ["💚", "🧘", "🏃", "💪", "🥗", "🧠", "😌", "🌱"],
    "money": ["💰", "💵", "📈", "🤑", "💎", "🏦", "📊", "💲"],
    "warning": ["⚠️", "🚨", "❗", "‼️", "🔴", "⛔"],
    "checkmarks": ["✅", "☑️", "✔️", "💯", "👍"],
    "arrows": ["➡️", "⬇️", "👇", "📍", "🔜"],
    "time": ["⏰", "🕐", "⏳", "📅", "🗓️"],
    "creative": ["🎨", "🖌️", "📸", "🎬", "🎵", "🎹"],
}

# Brand voice to emoji density and style mapping
VOICE_PROFILES = {
    "corporate": {
        "density": "low",
        "max_per_post": 2,
        "allowed_categories": ["checkmarks", "arrows", "business"],
        "no_faces": True,
        "style": "Minimal, professional symbols only",
    },
    "professional": {
        "density": "low-medium",
        "max_per_post": 3,
        "allowed_categories": ["checkmarks", "arrows", "business", "education", "pointing"],
        "no_faces": True,
        "style": "Subtle emphasis, no face emoji",
    },
    "friendly": {
        "density": "medium",
        "max_per_post": 5,
        "allowed_categories": [
            "positive", "love", "pointing", "celebration",
            "checkmarks", "education", "nature",
        ],
        "no_faces": False,
        "style": "Warm and relatable, occasional face emoji",
    },
    "playful": {
        "density": "high",
        "max_per_post": 8,
        "allowed_categories": [
            "positive", "love", "pointing", "celebration",
            "creative", "food", "nature",
        ],
        "no_faces": False,
        "style": "Expressive, creative, fun combinations",
    },
    "gen-z": {
        "density": "high",
        "max_per_post": 10,
        "allowed_categories": list(EMOJI_BY_CATEGORY.keys()),
        "no_faces": False,
        "style": "Ironic, layered, trend-specific, sometimes absurdist",
    },
}

# Platform-specific emoji norms
PLATFORM_EMOJI_NORMS = {
    "instagram": {
        "hook_emoji": True,
        "bullet_emoji": True,
        "cta_emoji": True,
        "line_break_emoji": True,
        "density_modifier": 1.0,
        "notes": "Liberal emoji use, especially in hooks and bullet lists",
    },
    "linkedin": {
        "hook_emoji": False,
        "bullet_emoji": True,
        "cta_emoji": False,
        "line_break_emoji": False,
        "density_modifier": 0.4,
        "notes": "Minimal emoji, checkmarks and arrows for structure",
    },
    "tiktok": {
        "hook_emoji": True,
        "bullet_emoji": False,
        "cta_emoji": True,
        "line_break_emoji": False,
        "density_modifier": 0.7,
        "notes": "Moderate emoji, trend-dependent",
    },
    "facebook": {
        "hook_emoji": True,
        "bullet_emoji": True,
        "cta_emoji": True,
        "line_break_emoji": False,
        "density_modifier": 0.8,
        "notes": "Moderate use, similar to Instagram but slightly less",
    },
    "twitter": {
        "hook_emoji": True,
        "bullet_emoji": False,
        "cta_emoji": False,
        "line_break_emoji": False,
        "density_modifier": 0.5,
        "notes": "Strategic placement, 1-2 per tweet max",
    },
}

# Keyword to emoji mapping for content-aware placement
KEYWORD_EMOJI_MAP = {
    # Positive sentiment
    "new": "✨",
    "launch": "🚀",
    "announce": "📣",
    "excited": "🎉",
    "love": "❤️",
    "amazing": "🤩",
    "incredible": "🔥",
    "beautiful": "✨",
    "perfect": "💯",
    "best": "🏆",
    # Actions
    "save": "📌",
    "share": "📤",
    "click": "👆",
    "follow": "👉",
    "comment": "💬",
    "download": "⬇️",
    "subscribe": "🔔",
    "join": "🤝",
    # Topics
    "tip": "💡",
    "secret": "🤫",
    "free": "🎁",
    "sale": "🏷️",
    "limited": "⏰",
    "exclusive": "👑",
    "results": "📈",
    "money": "💰",
    "learn": "📚",
    "grow": "🌱",
    "health": "💚",
    "food": "🍽️",
    "coffee": "☕",
    "morning": "🌅",
    "night": "🌙",
    "travel": "✈️",
    "home": "🏡",
    "work": "💼",
    "team": "👥",
    "win": "🏆",
    "goal": "🎯",
    "idea": "💡",
    "warning": "⚠️",
    "important": "❗",
    "question": "❓",
    "check": "✅",
    "step": "📍",
    # Emotions
    "happy": "😊",
    "sad": "😢",
    "funny": "😂",
    "surprised": "😮",
    "grateful": "🙏",
    "proud": "💪",
    "inspired": "✨",
}


def count_existing_emoji(text: str) -> int:
    """Count emoji already present in the text."""
    emoji_pattern = re.compile(
        "["
        "\U0001F600-\U0001F64F"
        "\U0001F300-\U0001F5FF"
        "\U0001F680-\U0001F6FF"
        "\U0001F1E0-\U0001F1FF"
        "\U00002702-\U000027B0"
        "\U000024C2-\U0001F251"
        "\U0001F900-\U0001F9FF"
        "\U0001FA00-\U0001FA6F"
        "\U0001FA70-\U0001FAFF"
        "\U00002600-\U000026FF"
        "\U00002700-\U000027BF"
        "]+",
        flags=re.UNICODE,
    )
    return len(emoji_pattern.findall(text))


def detect_keywords(text: str) -> list[tuple[str, str]]:
    """Detect keywords in text and return matching emoji suggestions."""
    text_lower = text.lower()
    matches = []
    for keyword, emoji in KEYWORD_EMOJI_MAP.items():
        if keyword in text_lower:
            matches.append((keyword, emoji))
    return matches


def get_allowed_emoji(brand_voice: str, platform: str) -> list[str]:
    """Get the list of emoji allowed by brand voice and platform norms."""
    profile = VOICE_PROFILES.get(brand_voice, VOICE_PROFILES["friendly"])
    allowed_cats = profile["allowed_categories"]

    all_emoji = []
    for cat in allowed_cats:
        all_emoji.extend(EMOJI_BY_CATEGORY.get(cat, []))

    return list(set(all_emoji))


def calculate_max_emoji(brand_voice: str, platform: str) -> int:
    """Calculate maximum emoji count based on voice and platform."""
    profile = VOICE_PROFILES.get(brand_voice, VOICE_PROFILES["friendly"])
    norms = PLATFORM_EMOJI_NORMS.get(platform, PLATFORM_EMOJI_NORMS["instagram"])

    base_max = profile["max_per_post"]
    modifier = norms["density_modifier"]

    return max(1, int(base_max * modifier))


def suggest_placements(
    text: str,
    platform: str,
    brand_voice: str,
) -> dict:
    """
    Analyze text and suggest strategic emoji placements.

    Returns a structured suggestion with positions and reasoning.
    """
    profile = VOICE_PROFILES.get(brand_voice, VOICE_PROFILES["friendly"])
    norms = PLATFORM_EMOJI_NORMS.get(platform, PLATFORM_EMOJI_NORMS["instagram"])
    max_emoji = calculate_max_emoji(brand_voice, platform)
    existing_count = count_existing_emoji(text)
    remaining_budget = max(0, max_emoji - existing_count)

    # Detect keywords for content-aware emoji
    keyword_matches = detect_keywords(text)

    # Split text into lines for positional analysis
    lines = text.strip().split("\n")

    suggestions = []

    # 1. Hook emoji (first line)
    if norms["hook_emoji"] and remaining_budget > 0 and lines:
        first_line = lines[0]
        hook_emoji_candidates = [
            emoji for kw, emoji in keyword_matches
            if kw.lower() in first_line.lower()
        ]
        if hook_emoji_candidates:
            suggestions.append({
                "position": "hook_end",
                "line": 0,
                "emoji": hook_emoji_candidates[0],
                "reasoning": f"Reinforces hook keyword",
            })
            remaining_budget -= 1

    # 2. Bullet point emoji
    if norms["bullet_emoji"] and remaining_budget > 0:
        for i, line in enumerate(lines):
            if line.strip().startswith(("-", "*", "•")) and remaining_budget > 0:
                line_matches = [
                    emoji for kw, emoji in keyword_matches
                    if kw.lower() in line.lower()
                ]
                if line_matches:
                    suggestions.append({
                        "position": "bullet_prefix",
                        "line": i,
                        "emoji": line_matches[0],
                        "reasoning": "Replace bullet with themed emoji",
                    })
                    remaining_budget -= 1

    # 3. CTA emoji
    cta_keywords = ["link", "click", "follow", "save", "share", "comment", "dm", "join"]
    if norms["cta_emoji"] and remaining_budget > 0:
        for i, line in enumerate(lines):
            line_lower = line.lower()
            if any(kw in line_lower for kw in cta_keywords):
                cta_emoji_match = [
                    emoji for kw, emoji in keyword_matches
                    if kw.lower() in line_lower
                ]
                emoji_choice = cta_emoji_match[0] if cta_emoji_match else "👇"
                suggestions.append({
                    "position": "cta_prefix",
                    "line": i,
                    "emoji": emoji_choice,
                    "reasoning": "Draw attention to CTA",
                })
                remaining_budget -= 1
                break

    # 4. Fill remaining budget with keyword-matched emoji
    used_lines = {s["line"] for s in suggestions}
    for kw, emoji in keyword_matches:
        if remaining_budget <= 0:
            break
        for i, line in enumerate(lines):
            if i not in used_lines and kw.lower() in line.lower():
                suggestions.append({
                    "position": "keyword_accent",
                    "line": i,
                    "emoji": emoji,
                    "reasoning": f"Accent for keyword '{kw}'",
                })
                used_lines.add(i)
                remaining_budget -= 1
                break

    return {
        "existing_emoji_count": existing_count,
        "max_emoji_allowed": max_emoji,
        "emoji_budget_remaining": max(0, max_emoji - existing_count - len(suggestions)),
        "suggestions": suggestions,
    }


def apply_suggestions(text: str, suggestions: list[dict]) -> str:
    """Apply emoji suggestions to the text (best-effort)."""
    lines = text.strip().split("\n")

    # Sort suggestions by line number (reverse to avoid index shifting)
    sorted_suggestions = sorted(suggestions, key=lambda s: s["line"], reverse=True)

    for suggestion in sorted_suggestions:
        line_idx = suggestion["line"]
        if line_idx >= len(lines):
            continue

        emoji = suggestion["emoji"]
        position = suggestion["position"]

        if position in ("hook_end", "keyword_accent"):
            # Add emoji at end of line
            lines[line_idx] = f"{lines[line_idx].rstrip()} {emoji}"
        elif position == "bullet_prefix":
            # Replace bullet character with emoji
            line = lines[line_idx]
            line = re.sub(r"^(\s*)[-*•]\s*", rf"\1{emoji} ", line)
            lines[line_idx] = line
        elif position == "cta_prefix":
            # Add emoji before CTA line
            lines[line_idx] = f"{emoji} {lines[line_idx].lstrip()}"

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Add strategic emoji placement to social media captions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Caption text (string) or path to a text/JSON file",
    )
    parser.add_argument(
        "--platform",
        required=True,
        choices=list(PLATFORM_EMOJI_NORMS.keys()),
        help="Target social media platform",
    )
    parser.add_argument(
        "--brand-voice",
        default="friendly",
        choices=list(VOICE_PROFILES.keys()),
        help="Brand voice profile (default: friendly)",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Apply suggestions directly to the text (default: suggest only)",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output file path (default: stdout)",
    )

    args = parser.parse_args()

    # Check for .env
    env_path = Path(".env")
    if env_path.exists():
        print("Note: .env file detected in working directory", file=sys.stderr)

    # Load input text
    input_text = args.input
    input_path = Path(args.input)
    if input_path.exists():
        with open(input_path, "r", encoding="utf-8") as f:
            content = f.read()
            # Try parsing as JSON
            try:
                data = json.loads(content)
                # Extract text from caption structure
                if "caption" in data:
                    parts = []
                    for section in data["caption"].values():
                        if isinstance(section, dict) and "text" in section:
                            parts.append(section["text"])
                    input_text = "\n\n".join(parts)
                else:
                    input_text = content
            except json.JSONDecodeError:
                input_text = content

    # Get suggestions
    placement_data = suggest_placements(input_text, args.platform, args.brand_voice)

    # Build output
    output = {
        "platform": args.platform,
        "brand_voice": args.brand_voice,
        "voice_profile": VOICE_PROFILES.get(args.brand_voice, {}),
        "platform_norms": PLATFORM_EMOJI_NORMS.get(args.platform, {}),
        "original_text": input_text,
        "analysis": placement_data,
    }

    if args.apply:
        optimized_text = apply_suggestions(input_text, placement_data["suggestions"])
        output["optimized_text"] = optimized_text

    output_json = json.dumps(output, indent=2, ensure_ascii=False)

    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(output_json)
        print(f"Emoji analysis written to {args.output}", file=sys.stderr)
    else:
        print(output_json)


if __name__ == "__main__":
    main()
