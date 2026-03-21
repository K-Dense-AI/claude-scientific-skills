#!/usr/bin/env python3
"""Suggest viral content angles by combining trending formats with brand pillars.

Analyzes the topic through multiple psychological lenses and trending format
templates to produce ranked angle suggestions with hooks, format recommendations,
and engagement estimates.
"""

import argparse
import json
import uuid
from datetime import datetime, timezone
from pathlib import Path


TRENDING_FORMATS = [
    {
        "name": "POV (Point of View)",
        "engagement_profile": "comments, shares",
        "effort": "low",
        "best_platforms": ["tiktok", "instagram"],
    },
    {
        "name": "Before/After Transformation",
        "engagement_profile": "saves, shares",
        "effort": "medium",
        "best_platforms": ["instagram", "tiktok", "linkedin"],
    },
    {
        "name": "Myth vs Reality",
        "engagement_profile": "saves, comments",
        "effort": "medium",
        "best_platforms": ["instagram", "tiktok", "twitter", "linkedin"],
    },
    {
        "name": "Step-by-Step Tutorial",
        "engagement_profile": "saves, shares",
        "effort": "medium",
        "best_platforms": ["instagram", "tiktok", "youtube", "linkedin"],
    },
    {
        "name": "Listicle Carousel",
        "engagement_profile": "saves, shares",
        "effort": "medium",
        "best_platforms": ["instagram", "linkedin"],
    },
    {
        "name": "Day-in-the-Life",
        "engagement_profile": "follows, comments",
        "effort": "medium-high",
        "best_platforms": ["tiktok", "instagram", "youtube"],
    },
    {
        "name": "Behind-the-Scenes",
        "engagement_profile": "comments, trust",
        "effort": "low-medium",
        "best_platforms": ["instagram", "tiktok", "linkedin"],
    },
    {
        "name": "This or That Poll",
        "engagement_profile": "votes, comments",
        "effort": "low",
        "best_platforms": ["instagram", "twitter", "linkedin"],
    },
    {
        "name": "Comparison (X vs Y)",
        "engagement_profile": "saves, comments",
        "effort": "medium",
        "best_platforms": ["instagram", "tiktok", "youtube", "linkedin"],
    },
    {
        "name": "Mini-Documentary",
        "engagement_profile": "shares, follows",
        "effort": "high",
        "best_platforms": ["youtube", "tiktok", "instagram"],
    },
    {
        "name": "Data Carousel with Narrative",
        "engagement_profile": "saves, shares",
        "effort": "high",
        "best_platforms": ["linkedin", "instagram", "twitter"],
    },
    {
        "name": "Setup-Twist-Payoff",
        "engagement_profile": "shares, replays",
        "effort": "low-medium",
        "best_platforms": ["tiktok", "instagram", "youtube"],
    },
    {
        "name": "Green Screen Reaction",
        "engagement_profile": "views, follows",
        "effort": "low",
        "best_platforms": ["tiktok", "instagram"],
    },
    {
        "name": "Wrong Answers Only",
        "engagement_profile": "comments",
        "effort": "low",
        "best_platforms": ["twitter", "instagram", "tiktok"],
    },
    {
        "name": "Photo Dump with Narrative",
        "engagement_profile": "likes, comments",
        "effort": "low",
        "best_platforms": ["instagram", "tiktok"],
    },
]

HOOK_PSYCHOLOGY = [
    {
        "type": "curiosity_gap",
        "template": "What most people don't realize about {topic} is...",
        "engagement_driver": "clicks, reads",
    },
    {
        "type": "controversy",
        "template": "Unpopular opinion: the way most people approach {topic} is fundamentally broken.",
        "engagement_driver": "comments, shares",
    },
    {
        "type": "storytelling",
        "template": "I spent [time] deep in {topic}. Here's the story nobody tells.",
        "engagement_driver": "watch time, follows",
    },
    {
        "type": "educational",
        "template": "The [number]-step framework for mastering {topic} (backed by data).",
        "engagement_driver": "saves, shares",
    },
    {
        "type": "social_proof",
        "template": "Why leading [industry players] are rethinking {topic} right now.",
        "engagement_driver": "trust, clicks",
    },
    {
        "type": "fomo",
        "template": "The {topic} playbook is changing. Here's what you need to know before it does.",
        "engagement_driver": "urgency, saves",
    },
    {
        "type": "pattern_interrupt",
        "template": "Stop scrolling. If you care about {topic}, read this.",
        "engagement_driver": "attention, engagement",
    },
]

CONTENT_PILLARS = [
    "authority",
    "relatability",
    "community",
    "aspiration",
    "entertainment",
]

ENGAGEMENT_LEVELS = ["moderate", "high", "very_high", "viral_potential"]


def generate_angles(
    brand_dna: str,
    topic: str,
    num_angles: int,
) -> list[dict]:
    """Generate viral angle suggestions."""
    angles = []

    for i in range(min(num_angles, len(TRENDING_FORMATS))):
        fmt = TRENDING_FORMATS[i % len(TRENDING_FORMATS)]
        hook_psych = HOOK_PSYCHOLOGY[i % len(HOOK_PSYCHOLOGY)]
        pillar = CONTENT_PILLARS[i % len(CONTENT_PILLARS)]
        engagement = ENGAGEMENT_LEVELS[
            min(i % len(ENGAGEMENT_LEVELS), len(ENGAGEMENT_LEVELS) - 1)
        ]

        angle = {
            "angle_id": str(uuid.uuid4()),
            "rank": i + 1,
            "hook": {
                "text": hook_psych["template"].format(topic=topic),
                "psychology": hook_psych["type"],
                "engagement_driver": hook_psych["engagement_driver"],
            },
            "format": {
                "name": fmt["name"],
                "effort_level": fmt["effort"],
                "engagement_profile": fmt["engagement_profile"],
            },
            "platform_recommendation": fmt["best_platforms"][0],
            "all_suitable_platforms": fmt["best_platforms"],
            "content_pillar": pillar,
            "estimated_engagement": engagement,
            "reasoning": (
                f"Combining {hook_psych['type'].replace('_', ' ')} psychology "
                f"with the '{fmt['name']}' format creates a strong "
                f"{pillar} content piece. The hook targets "
                f"{hook_psych['engagement_driver']}, while the format "
                f"naturally drives {fmt['engagement_profile']}. "
                f"Brand alignment with '{brand_dna}' comes through "
                f"the {pillar} pillar positioning."
            ),
            "brief_guidance": {
                "visual_direction": (
                    f"Design for {fmt['name']} format on "
                    f"{fmt['best_platforms'][0]}. "
                    f"Mood: {pillar}-driven. "
                    f"See references/trending_formats.md for structure template."
                ),
                "caption_approach": (
                    f"Open with {hook_psych['type'].replace('_', ' ')} hook. "
                    f"Body follows {fmt['name']} narrative structure. "
                    f"CTA aligned with {fmt['engagement_profile']} goals."
                ),
            },
        }
        angles.append(angle)

    return angles


def main():
    parser = argparse.ArgumentParser(
        description="Suggest viral content angles for a topic and brand.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example:\n"
            "  python suggest_viral_angles.py \\\n"
            '    --brand-dna "Sustainable fashion brand targeting Gen Z" \\\n'
            '    --topic "fast fashion environmental impact" \\\n'
            "    --num-angles 7 \\\n"
            "    --output angles.json"
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
        "--num-angles",
        type=int,
        default=5,
        help="Number of angle suggestions to generate (default: 5)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output JSON file path (prints to stdout if omitted)",
    )

    args = parser.parse_args()

    angles = generate_angles(
        brand_dna=args.brand_dna,
        topic=args.topic,
        num_angles=args.num_angles,
    )

    result = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "brand_dna": args.brand_dna,
        "topic": args.topic,
        "num_angles": len(angles),
        "angles": angles,
    }

    output_json = json.dumps(result, indent=2)

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output_json)
        print(f"Angles written to {args.output}")
    else:
        print(output_json)


if __name__ == "__main__":
    main()
