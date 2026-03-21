#!/usr/bin/env python3
"""Generate a structured content brief for a specific platform and format.

Produces a JSON brief with hook, angle, key messages, visual direction,
caption notes, CTA, hashtags, and target emotion -- ready for handoff
to visual-generator and copywriting skills.
"""

import argparse
import json
import uuid
from datetime import datetime, timezone
from pathlib import Path


PLATFORM_FORMATS = {
    "instagram": ["post", "carousel", "reel", "story", "quiz"],
    "tiktok": ["video", "carousel", "story", "live"],
    "linkedin": ["text", "carousel", "video", "poll", "article", "newsletter"],
    "twitter": ["tweet", "thread", "poll"],
    "youtube": ["video", "short", "live", "community"],
    "facebook": ["post", "video", "reel", "story", "carousel", "live"],
}

CONTENT_PILLARS = [
    "authority",
    "relatability",
    "community",
    "aspiration",
    "entertainment",
]

TARGET_EMOTIONS = [
    "awe",
    "curiosity",
    "inspiration",
    "urgency",
    "humor",
    "empathy",
    "surprise",
    "confidence",
    "belonging",
    "determination",
]

HOOK_TYPES = [
    "curiosity_gap",
    "controversy",
    "storytelling",
    "educational",
    "social_proof",
    "fomo",
    "pattern_interrupt",
]

MIX_CATEGORIES = ["value", "curated", "promotional"]


def select_hook_type(format_type: str) -> str:
    """Select the most effective hook type for a given format."""
    format_hook_map = {
        "carousel": ["curiosity_gap", "educational", "storytelling"],
        "reel": ["pattern_interrupt", "storytelling", "controversy"],
        "video": ["pattern_interrupt", "storytelling", "curiosity_gap"],
        "short": ["pattern_interrupt", "curiosity_gap", "fomo"],
        "story": ["pattern_interrupt", "fomo", "social_proof"],
        "post": ["curiosity_gap", "controversy", "educational"],
        "text": ["storytelling", "controversy", "educational"],
        "thread": ["curiosity_gap", "educational", "storytelling"],
        "tweet": ["controversy", "curiosity_gap", "pattern_interrupt"],
        "poll": ["controversy", "curiosity_gap", "social_proof"],
        "quiz": ["curiosity_gap", "educational", "pattern_interrupt"],
        "live": ["fomo", "social_proof", "curiosity_gap"],
        "article": ["curiosity_gap", "educational", "storytelling"],
        "newsletter": ["curiosity_gap", "educational", "storytelling"],
        "community": ["social_proof", "curiosity_gap", "pattern_interrupt"],
    }
    options = format_hook_map.get(format_type, HOOK_TYPES[:3])
    return options[0]


def select_emotion(topic: str, hook_type: str) -> str:
    """Select target emotion based on hook type."""
    emotion_map = {
        "curiosity_gap": "curiosity",
        "controversy": "surprise",
        "storytelling": "empathy",
        "educational": "confidence",
        "social_proof": "belonging",
        "fomo": "urgency",
        "pattern_interrupt": "surprise",
    }
    return emotion_map.get(hook_type, "curiosity")


def generate_brief(
    brand_dna: str,
    topic: str,
    platform: str,
    format_type: str,
) -> dict:
    """Generate a structured content brief."""
    hook_type = select_hook_type(format_type)
    emotion = select_emotion(topic, hook_type)
    pillar = CONTENT_PILLARS[hash(topic) % len(CONTENT_PILLARS)]
    mix_cat = MIX_CATEGORIES[0]  # Default to value content

    brief = {
        "brief_id": str(uuid.uuid4()),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "brand_context": {
            "brand_dna": brand_dna,
            "content_pillar": pillar,
            "target_audience": f"Derived from brand DNA: {brand_dna}",
        },
        "content_spec": {
            "platform": platform,
            "format": format_type,
            "hook": {
                "type": hook_type,
                "template": f"[{hook_type.replace('_', ' ').title()} hook about: {topic}]",
                "guidance": (
                    f"Write a {hook_type.replace('_', ' ')} hook that connects "
                    f"'{topic}' to the brand positioning: '{brand_dna}'. "
                    f"First 5-8 words must earn the rest of the sentence."
                ),
            },
            "angle": (
                f"Approach '{topic}' through the lens of {pillar} content. "
                f"Use {hook_type.replace('_', ' ')} psychology to open, "
                f"then deliver actionable value tied to '{brand_dna}'."
            ),
            "key_messages": [
                f"[Primary insight about {topic} relevant to target audience]",
                f"[Supporting data point or proof element]",
                f"[Actionable takeaway the audience can apply immediately]",
            ],
            "visual_direction": (
                f"Platform: {platform} | Format: {format_type} | "
                f"Mood: {emotion} | "
                f"Brand alignment: {brand_dna}. "
                f"See references/platform_formats.md for exact dimensions."
            ),
            "caption_notes": (
                f"Open with {hook_type.replace('_', ' ')} hook. "
                f"Target emotion: {emotion}. "
                f"End with CTA appropriate for {platform}. "
                f"Match brand voice derived from: {brand_dna}."
            ),
            "cta": f"[Platform-appropriate CTA for {platform} {format_type}]",
            "target_emotion": emotion,
        },
        "distribution": {
            "hashtags": {
                "high_volume": [f"[3-5 high-volume hashtags for {topic}]"],
                "medium": [f"[5-8 medium-volume hashtags for {topic}]"],
                "niche": [f"[3-5 niche hashtags for {topic}]"],
                "branded": [f"[1-2 branded hashtags from {brand_dna}]"],
            },
            "posting_time": f"[Optimal posting window for {platform}]",
            "cross_post_adaptations": {
                p: f"Adapt to {p} native format"
                for p in PLATFORM_FORMATS
                if p != platform
            },
        },
        "metadata": {
            "content_mix_category": mix_cat,
            "estimated_engagement": "medium",
            "trend_alignment": None,
            "hook_type": hook_type,
            "content_pillar": pillar,
        },
    }

    return brief


def main():
    parser = argparse.ArgumentParser(
        description="Generate a structured content brief for social media.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example:\n"
            "  python generate_content_brief.py \\\n"
            '    --brand-dna "SaaS productivity tool for remote teams" \\\n'
            '    --topic "async communication best practices" \\\n'
            "    --platform instagram \\\n"
            "    --format carousel \\\n"
            "    --output brief.json"
        ),
    )
    parser.add_argument(
        "--brand-dna",
        required=True,
        help="Brand positioning statement or description",
    )
    parser.add_argument(
        "--topic",
        required=True,
        help="Content topic or subject",
    )
    parser.add_argument(
        "--platform",
        required=True,
        choices=list(PLATFORM_FORMATS.keys()),
        help="Target social media platform",
    )
    parser.add_argument(
        "--format",
        dest="format_type",
        required=True,
        help="Content format (platform-dependent: post, carousel, reel, story, etc.)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output JSON file path (prints to stdout if omitted)",
    )

    args = parser.parse_args()

    valid_formats = PLATFORM_FORMATS.get(args.platform, [])
    if args.format_type not in valid_formats:
        parser.error(
            f"Format '{args.format_type}' is not valid for {args.platform}. "
            f"Valid formats: {', '.join(valid_formats)}"
        )

    brief = generate_brief(
        brand_dna=args.brand_dna,
        topic=args.topic,
        platform=args.platform,
        format_type=args.format_type,
    )

    output_json = json.dumps(brief, indent=2)

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output_json)
        print(f"Brief written to {args.output}")
    else:
        print(output_json)


if __name__ == "__main__":
    main()
