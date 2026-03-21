#!/usr/bin/env python3
"""Generate a tiered hashtag strategy for social media content.

Produces a structured hashtag set organized into four tiers:
high-volume (reach), medium (discoverability), niche (targeted),
and branded -- optimized for the specified platform.
"""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


PLATFORM_HASHTAG_LIMITS = {
    "instagram": {"max": 30, "recommended": 12, "placement": "caption or first comment"},
    "tiktok": {"max": 10, "recommended": 4, "placement": "caption"},
    "linkedin": {"max": 30, "recommended": 4, "placement": "end of post"},
    "twitter": {"max": 10, "recommended": 2, "placement": "inline or end of tweet"},
    "youtube": {"max": 60, "recommended": 5, "placement": "description, first 3 visible above title"},
    "facebook": {"max": 30, "recommended": 3, "placement": "end of post"},
}

TIER_RATIOS = {
    "high_volume": 0.25,
    "medium": 0.35,
    "niche": 0.30,
    "branded": 0.10,
}

TIER_DESCRIPTIONS = {
    "high_volume": (
        "Broad, popular hashtags with millions of posts. "
        "Provide maximum exposure but high competition. "
        "Your content appears briefly in the feed before being buried. "
        "Use for reach, not discoverability."
    ),
    "medium": (
        "Moderately popular hashtags with 100K-1M posts. "
        "Balance between reach and staying visible longer. "
        "Best for discoverability -- your content can rank in 'Top' tabs. "
        "Core of the hashtag strategy."
    ),
    "niche": (
        "Specific, targeted hashtags with 10K-100K posts. "
        "Highly relevant audience with lower competition. "
        "Content stays visible longest here. "
        "Drives the most qualified engagement."
    ),
    "branded": (
        "Brand-specific or campaign-specific hashtags. "
        "Build a searchable content library. "
        "Track campaign performance. "
        "Encourage UGC with a memorable branded tag."
    ),
}


def calculate_tier_counts(total: int) -> dict[str, int]:
    """Distribute hashtags across tiers based on ratios."""
    counts = {}
    remaining = total
    for tier, ratio in TIER_RATIOS.items():
        count = max(1, round(total * ratio))
        counts[tier] = min(count, remaining)
        remaining -= counts[tier]
        if remaining <= 0:
            break

    # Distribute any remainder to medium tier
    if remaining > 0:
        counts["medium"] = counts.get("medium", 0) + remaining

    return counts


def generate_hashtag_templates(
    topic: str,
    platform: str,
    brand_dna: str,
    num_hashtags: int,
) -> dict:
    """Generate a tiered hashtag strategy."""
    platform_config = PLATFORM_HASHTAG_LIMITS.get(platform, PLATFORM_HASHTAG_LIMITS["instagram"])
    effective_count = min(num_hashtags, platform_config["max"])
    tier_counts = calculate_tier_counts(effective_count)

    # Generate template hashtags per tier
    topic_words = topic.lower().replace("-", " ").split()
    primary_keyword = topic_words[0] if topic_words else "content"

    tiers = {}
    for tier_name, count in tier_counts.items():
        tiers[tier_name] = {
            "count": count,
            "description": TIER_DESCRIPTIONS[tier_name],
            "hashtags": [
                f"[{tier_name} hashtag {i+1} for '{topic}']"
                for i in range(count)
            ],
            "guidance": _tier_guidance(tier_name, primary_keyword, brand_dna),
        }

    strategy = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "topic": topic,
        "platform": platform,
        "brand_dna": brand_dna,
        "platform_config": {
            "max_allowed": platform_config["max"],
            "recommended_count": platform_config["recommended"],
            "placement": platform_config["placement"],
        },
        "total_hashtags": effective_count,
        "tier_distribution": tier_counts,
        "tiers": tiers,
        "usage_guidelines": {
            "rotation": (
                "Rotate hashtag sets every 3-5 posts to avoid shadowban risk "
                "and test different combinations. Keep branded hashtags constant."
            ),
            "research_method": (
                f"Search each hashtag on {platform} before using. "
                "Verify: (1) post volume is in the target tier range, "
                "(2) content in the hashtag matches your brand, "
                "(3) the hashtag isn't associated with spam or banned content."
            ),
            "tracking": (
                "Track which hashtag sets drive the most profile visits "
                "and follows (not just impressions). Use platform analytics "
                "to identify top-performing combinations."
            ),
            "platform_specific": _platform_notes(platform),
        },
    }

    return strategy


def _tier_guidance(tier: str, keyword: str, brand_dna: str) -> str:
    """Generate tier-specific selection guidance."""
    guidance = {
        "high_volume": (
            f"Search for broad hashtags related to '{keyword}' with 1M+ posts. "
            f"Examples: industry-wide terms, popular category tags, "
            f"general audience descriptors. Use sparingly -- these provide "
            f"reach but your content competes with millions of posts."
        ),
        "medium": (
            f"Search for hashtags related to '{keyword}' with 100K-1M posts. "
            f"These should be more specific than high-volume tags but still "
            f"have active communities. Look for topic-specific tags, "
            f"format-specific tags, and audience-segment tags."
        ),
        "niche": (
            f"Search for highly specific hashtags related to '{keyword}' "
            f"with 10K-100K posts. These should match your exact topic, "
            f"audience segment, or content angle. Your content can stay "
            f"visible in these feeds for days."
        ),
        "branded": (
            f"Create or use hashtags specific to the brand: '{brand_dna}'. "
            f"Include a general brand hashtag and a campaign-specific one. "
            f"Keep them short, memorable, and unique to your brand."
        ),
    }
    return guidance.get(tier, "")


def _platform_notes(platform: str) -> str:
    """Platform-specific hashtag best practices."""
    notes = {
        "instagram": (
            "Place hashtags in the caption (first comment placement has "
            "reduced effectiveness since 2023). Mix hashtag sizes. "
            "Avoid banned or spammy hashtags. Use the 'Recent' tab to "
            "verify your content appears under each hashtag."
        ),
        "tiktok": (
            "Use fewer, more targeted hashtags (3-5). TikTok's algorithm "
            "relies more on content analysis than hashtags. Include one "
            "trending hashtag if relevant. #FYP and #ForYou have minimal "
            "proven effect -- use topic-specific tags instead."
        ),
        "linkedin": (
            "Use 3-5 hashtags maximum. Place at the end of the post. "
            "LinkedIn hashtags are more about categorization than "
            "discoverability. Follow your own hashtags to engage with "
            "the community using them."
        ),
        "twitter": (
            "Use 1-2 hashtags maximum. Inline placement reads more "
            "naturally than appended hashtags. Hashtags in tweets reduce "
            "readability if overused. For threads, hashtag only the "
            "first tweet."
        ),
        "youtube": (
            "Place the most important 3 hashtags in the first line of "
            "the description (they appear above the title). Use remaining "
            "hashtags in the description body. Hashtags in titles can "
            "look spammy -- use description placement instead."
        ),
        "facebook": (
            "Hashtags have limited discoverability on Facebook. Use 1-3 "
            "if any. They work better in Groups than on Pages. Focus on "
            "branded hashtags for campaign tracking rather than "
            "discoverability."
        ),
    }
    return notes.get(platform, "No platform-specific notes available.")


def main():
    parser = argparse.ArgumentParser(
        description="Generate a tiered hashtag strategy for social media content.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example:\n"
            "  python hashtag_optimizer.py \\\n"
            '    --topic "home workout routines" \\\n'
            "    --platform tiktok \\\n"
            '    --brand-dna "Fitness app for busy professionals" \\\n'
            "    --num-hashtags 25 \\\n"
            "    --output hashtags.json"
        ),
    )
    parser.add_argument(
        "--topic",
        required=True,
        help="Content topic or subject",
    )
    parser.add_argument(
        "--platform",
        required=True,
        choices=list(PLATFORM_HASHTAG_LIMITS.keys()),
        help="Target social media platform",
    )
    parser.add_argument(
        "--brand-dna",
        required=True,
        help="Brand positioning statement or description",
    )
    parser.add_argument(
        "--num-hashtags",
        type=int,
        default=15,
        help="Total number of hashtags to generate (default: 15)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output JSON file path (prints to stdout if omitted)",
    )

    args = parser.parse_args()

    strategy = generate_hashtag_templates(
        topic=args.topic,
        platform=args.platform,
        brand_dna=args.brand_dna,
        num_hashtags=args.num_hashtags,
    )

    output_json = json.dumps(strategy, indent=2)

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output_json)
        print(f"Hashtag strategy written to {args.output}")
    else:
        print(output_json)


if __name__ == "__main__":
    main()
