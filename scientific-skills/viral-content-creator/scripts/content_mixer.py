#!/usr/bin/env python3
"""Generate a balanced content calendar following the 60/30/10 rule.

Produces a structured content calendar with pillar rotation,
platform-specific format assignments, and content mix category
distribution across the specified duration.
"""

import argparse
import json
import uuid
import math
from datetime import datetime, timezone, timedelta
from pathlib import Path


MIX_RULE = {
    "value": {
        "ratio": 0.60,
        "description": "Educate, entertain, or inspire the audience",
        "examples": [
            "tutorials",
            "tips",
            "behind-the-scenes",
            "storytelling",
            "how-tos",
            "frameworks",
            "industry insights",
        ],
    },
    "curated": {
        "ratio": 0.30,
        "description": "Build community and show industry awareness",
        "examples": [
            "industry news commentary",
            "UGC reposts",
            "collaborations",
            "polls",
            "Q&A",
            "challenges",
            "community spotlights",
        ],
    },
    "promotional": {
        "ratio": 0.10,
        "description": "Drive conversions and revenue",
        "examples": [
            "product features",
            "launches",
            "offers",
            "testimonials",
            "case studies",
            "demos",
            "CTAs",
        ],
    },
}

CONTENT_PILLARS = [
    {
        "name": "authority",
        "description": "Demonstrate expertise",
        "best_formats": ["carousel", "thread", "video", "article"],
        "best_mix_categories": ["value"],
    },
    {
        "name": "relatability",
        "description": "Show the human side",
        "best_formats": ["reel", "story", "post", "video"],
        "best_mix_categories": ["value", "curated"],
    },
    {
        "name": "community",
        "description": "Foster belonging and participation",
        "best_formats": ["poll", "story", "post", "quiz"],
        "best_mix_categories": ["curated"],
    },
    {
        "name": "aspiration",
        "description": "Paint the future and showcase transformations",
        "best_formats": ["carousel", "reel", "video", "post"],
        "best_mix_categories": ["value", "promotional"],
    },
    {
        "name": "entertainment",
        "description": "Create joy and shareable moments",
        "best_formats": ["reel", "video", "post", "story"],
        "best_mix_categories": ["value", "curated"],
    },
]

PLATFORM_FORMAT_WEIGHTS = {
    "instagram": {
        "carousel": 0.30,
        "reel": 0.35,
        "post": 0.15,
        "story": 0.15,
        "quiz": 0.05,
    },
    "tiktok": {
        "video": 0.70,
        "carousel": 0.15,
        "story": 0.10,
        "live": 0.05,
    },
    "linkedin": {
        "text": 0.30,
        "carousel": 0.30,
        "video": 0.15,
        "poll": 0.15,
        "article": 0.10,
    },
    "twitter": {
        "tweet": 0.40,
        "thread": 0.35,
        "poll": 0.15,
        "tweet": 0.10,
    },
    "youtube": {
        "video": 0.40,
        "short": 0.45,
        "live": 0.10,
        "community": 0.05,
    },
    "facebook": {
        "post": 0.25,
        "video": 0.25,
        "reel": 0.25,
        "story": 0.15,
        "carousel": 0.10,
    },
}

POSTING_WINDOWS = {
    "instagram": ["09:00-11:00", "12:00-13:00", "17:00-19:00"],
    "tiktok": ["07:00-09:00", "12:00-15:00", "19:00-22:00"],
    "linkedin": ["07:30-08:30", "12:00-13:00", "17:00-18:00"],
    "twitter": ["08:00-10:00", "12:00-13:00", "17:00-18:00"],
    "youtube": ["14:00-16:00", "17:00-19:00"],
    "facebook": ["09:00-11:00", "13:00-15:00", "18:00-20:00"],
}


def distribute_mix_categories(total_posts: int) -> dict[str, int]:
    """Distribute posts across mix categories per the 60/30/10 rule."""
    counts = {}
    remaining = total_posts
    for category, config in MIX_RULE.items():
        count = max(1, round(total_posts * config["ratio"]))
        counts[category] = min(count, remaining)
        remaining -= counts[category]

    # Distribute remainder to value
    if remaining > 0:
        counts["value"] = counts.get("value", 0) + remaining

    return counts


def select_format(platform: str, pillar: dict) -> str:
    """Select a content format based on platform weights and pillar fit."""
    formats = PLATFORM_FORMAT_WEIGHTS.get(platform, {"post": 1.0})
    pillar_formats = set(pillar["best_formats"])

    # Prefer formats that match both platform weight and pillar
    for fmt in sorted(formats, key=formats.get, reverse=True):
        if fmt in pillar_formats:
            return fmt

    # Fallback to highest-weighted platform format
    return max(formats, key=formats.get)


def generate_calendar(
    brand_dna: str,
    duration_days: int,
    platforms: list[str],
    posts_per_week: int,
) -> dict:
    """Generate a content calendar with balanced mix and pillar rotation."""
    total_weeks = math.ceil(duration_days / 7)
    total_posts_per_platform = total_weeks * posts_per_week
    start_date = datetime.now(timezone.utc).date()

    calendar = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "brand_dna": brand_dna,
        "duration_days": duration_days,
        "start_date": start_date.isoformat(),
        "end_date": (start_date + timedelta(days=duration_days)).isoformat(),
        "platforms": platforms,
        "posts_per_week_per_platform": posts_per_week,
        "total_posts": total_posts_per_platform * len(platforms),
        "mix_distribution": {
            cat: {
                "target_ratio": config["ratio"],
                "description": config["description"],
            }
            for cat, config in MIX_RULE.items()
        },
        "content_pillars": [
            {"name": p["name"], "description": p["description"]}
            for p in CONTENT_PILLARS
        ],
        "platform_calendars": {},
    }

    for platform in platforms:
        mix_counts = distribute_mix_categories(total_posts_per_platform)
        windows = POSTING_WINDOWS.get(platform, ["10:00-12:00"])

        posts = []
        post_index = 0
        mix_pool = []
        for cat, count in mix_counts.items():
            mix_pool.extend([cat] * count)

        for week in range(total_weeks):
            week_start = start_date + timedelta(weeks=week)

            for day_offset in range(min(posts_per_week, 7)):
                if post_index >= total_posts_per_platform:
                    break

                post_date = week_start + timedelta(
                    days=day_offset * (7 // posts_per_week)
                )
                if (post_date - start_date).days >= duration_days:
                    break

                pillar = CONTENT_PILLARS[post_index % len(CONTENT_PILLARS)]
                mix_category = (
                    mix_pool[post_index]
                    if post_index < len(mix_pool)
                    else "value"
                )
                fmt = select_format(platform, pillar)
                window = windows[post_index % len(windows)]

                post = {
                    "post_id": str(uuid.uuid4()),
                    "date": post_date.isoformat(),
                    "day_of_week": post_date.strftime("%A"),
                    "posting_window": window,
                    "platform": platform,
                    "format": fmt,
                    "content_pillar": pillar["name"],
                    "mix_category": mix_category,
                    "topic_guidance": (
                        f"[{mix_category.title()} content] "
                        f"{pillar['description']} through {fmt} format. "
                        f"Topic aligned with: {brand_dna}"
                    ),
                    "hook_type_suggestion": _suggest_hook_type(
                        mix_category, pillar["name"]
                    ),
                    "brief_status": "pending",
                }
                posts.append(post)
                post_index += 1

        # Compute actual distribution
        actual_mix = {}
        for cat in MIX_RULE:
            count = sum(1 for p in posts if p["mix_category"] == cat)
            actual_mix[cat] = {
                "count": count,
                "actual_ratio": round(count / len(posts), 2) if posts else 0,
                "target_ratio": MIX_RULE[cat]["ratio"],
            }

        actual_pillars = {}
        for pillar in CONTENT_PILLARS:
            count = sum(
                1 for p in posts if p["content_pillar"] == pillar["name"]
            )
            actual_pillars[pillar["name"]] = count

        calendar["platform_calendars"][platform] = {
            "total_posts": len(posts),
            "mix_distribution": actual_mix,
            "pillar_distribution": actual_pillars,
            "posts": posts,
        }

    return calendar


def _suggest_hook_type(mix_category: str, pillar: str) -> str:
    """Suggest a hook type based on mix category and pillar."""
    suggestions = {
        ("value", "authority"): "educational",
        ("value", "relatability"): "storytelling",
        ("value", "community"): "pattern_interrupt",
        ("value", "aspiration"): "curiosity_gap",
        ("value", "entertainment"): "pattern_interrupt",
        ("curated", "authority"): "social_proof",
        ("curated", "relatability"): "storytelling",
        ("curated", "community"): "pattern_interrupt",
        ("curated", "aspiration"): "social_proof",
        ("curated", "entertainment"): "controversy",
        ("promotional", "authority"): "social_proof",
        ("promotional", "relatability"): "storytelling",
        ("promotional", "community"): "fomo",
        ("promotional", "aspiration"): "fomo",
        ("promotional", "entertainment"): "curiosity_gap",
    }
    return suggestions.get((mix_category, pillar), "curiosity_gap")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate a balanced content calendar following the "
            "60/30/10 rule with pillar rotation."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example:\n"
            "  python content_mixer.py \\\n"
            '    --brand-dna "B2B data analytics platform" \\\n'
            "    --duration-days 30 \\\n"
            "    --platforms instagram,linkedin,tiktok \\\n"
            "    --posts-per-week 5 \\\n"
            "    --output calendar.json"
        ),
    )
    parser.add_argument(
        "--brand-dna",
        required=True,
        help="Brand positioning statement or description",
    )
    parser.add_argument(
        "--duration-days",
        type=int,
        required=True,
        help="Calendar duration in days",
    )
    parser.add_argument(
        "--platforms",
        required=True,
        help="Comma-separated list of platforms (e.g., instagram,linkedin,tiktok)",
    )
    parser.add_argument(
        "--posts-per-week",
        type=int,
        default=5,
        help="Number of posts per week per platform (default: 5)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output JSON file path (prints to stdout if omitted)",
    )

    args = parser.parse_args()

    platforms = [p.strip().lower() for p in args.platforms.split(",")]
    valid_platforms = set(PLATFORM_FORMAT_WEIGHTS.keys())
    for p in platforms:
        if p not in valid_platforms:
            parser.error(
                f"Unknown platform '{p}'. "
                f"Valid platforms: {', '.join(sorted(valid_platforms))}"
            )

    calendar = generate_calendar(
        brand_dna=args.brand_dna,
        duration_days=args.duration_days,
        platforms=platforms,
        posts_per_week=args.posts_per_week,
    )

    output_json = json.dumps(calendar, indent=2)

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output_json)
        print(f"Content calendar written to {args.output}")
    else:
        print(output_json)


if __name__ == "__main__":
    main()
