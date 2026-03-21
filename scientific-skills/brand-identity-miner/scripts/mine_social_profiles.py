#!/usr/bin/env python3
"""
mine_social_profiles.py - Analyze social media profiles for brand identity signals.

Extracts bio information, content themes, visual patterns, hashtag strategy,
and engagement indicators from social media profiles. Supports Instagram,
Facebook, LinkedIn, and TikTok.

This script primarily serves as a structural template and data schema
definition. Actual social media scraping typically requires Chrome automation
MCP or platform-specific APIs due to authentication requirements and
anti-scraping measures.

Usage:
    python mine_social_profiles.py --profiles "instagram:@acme" "linkedin:/company/acme" --output social_analysis.json
    python mine_social_profiles.py --input raw_social_data.json --output social_analysis.json

Requirements:
    pip install requests beautifulsoup4

Environment Variables:
    INSTAGRAM_ACCESS_TOKEN - Optional. For Instagram Graph API access.
    LINKEDIN_ACCESS_TOKEN  - Optional. For LinkedIn API access.
"""

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def load_env(env_path: str = ".env") -> None:
    """Load environment variables from .env file if it exists."""
    env_file = Path(env_path)
    if env_file.exists():
        with open(env_file) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#") and "=" in line:
                    key, _, value = line.partition("=")
                    os.environ.setdefault(key.strip(), value.strip().strip("\"'"))


PLATFORM_URLS = {
    "instagram": "https://www.instagram.com/{handle}/",
    "facebook": "https://www.facebook.com/{handle}/",
    "linkedin": "https://www.linkedin.com/{handle}",
    "tiktok": "https://www.tiktok.com/@{handle}",
    "x_twitter": "https://x.com/{handle}",
    "youtube": "https://www.youtube.com/{handle}",
}

EMPTY_PROFILE_ANALYSIS = {
    "platform": "",
    "handle": "",
    "url": "",
    "fetched_at": "",
    "bio": {
        "text": "",
        "links": [],
        "hashtags": [],
        "mentions": [],
    },
    "metrics": {
        "followers": None,
        "following": None,
        "total_posts": None,
    },
    "content_analysis": {
        "recent_posts_analyzed": 0,
        "content_themes": [],
        "content_mix": {
            "educational": 0,
            "promotional": 0,
            "entertainment": 0,
            "inspirational": 0,
            "community": 0,
            "behind_the_scenes": 0,
            "news_commentary": 0,
            "values_driven": 0,
        },
        "formats_used": [],
        "posting_frequency": "",
        "avg_caption_length": None,
        "common_ctas": [],
    },
    "visual_patterns": {
        "dominant_colors": [],
        "filter_style": "",
        "grid_layout_pattern": "",
        "visual_consistency_score": None,
    },
    "hashtag_strategy": {
        "branded_hashtags": [],
        "community_hashtags": [],
        "industry_hashtags": [],
        "avg_hashtags_per_post": None,
        "top_hashtags": [],
    },
    "engagement": {
        "avg_likes": None,
        "avg_comments": None,
        "avg_shares": None,
        "engagement_rate": None,
        "top_performing_post": None,
        "low_performing_post": None,
    },
    "audience_signals": {
        "comment_sentiment": "",
        "comment_themes": [],
        "demographic_signals": [],
    },
}


def parse_profile_arg(profile_str: str) -> tuple[str, str]:
    """Parse a profile argument in format 'platform:handle'.

    Args:
        profile_str: String in format "instagram:@handle" or "linkedin:/company/name".

    Returns:
        Tuple of (platform, handle).

    Raises:
        ValueError: If the format is invalid or platform is unsupported.
    """
    if ":" not in profile_str:
        raise ValueError(
            f"Invalid profile format: '{profile_str}'. "
            f"Expected 'platform:handle' (e.g., 'instagram:@acme')."
        )

    platform, _, handle = profile_str.partition(":")
    platform = platform.strip().lower()
    handle = handle.strip().lstrip("@")

    supported = list(PLATFORM_URLS.keys())
    if platform not in supported:
        raise ValueError(
            f"Unsupported platform: '{platform}'. Supported: {', '.join(supported)}"
        )

    return platform, handle


def build_profile_url(platform: str, handle: str) -> str:
    """Build the public profile URL for a given platform and handle."""
    template = PLATFORM_URLS.get(platform, "")
    return template.format(handle=handle)


def create_empty_analysis(platform: str, handle: str) -> dict[str, Any]:
    """Create a blank analysis structure for a profile.

    The structure defines the expected output schema. In practice,
    Chrome automation MCP or platform APIs fill in these fields.
    """
    analysis = json.loads(json.dumps(EMPTY_PROFILE_ANALYSIS))
    analysis["platform"] = platform
    analysis["handle"] = handle
    analysis["url"] = build_profile_url(platform, handle)
    analysis["fetched_at"] = datetime.now(timezone.utc).isoformat()
    return analysis


def analyze_bio_text(bio_text: str) -> dict[str, Any]:
    """Analyze a social media bio for brand signals.

    Args:
        bio_text: The raw bio/about text from the profile.

    Returns:
        Dictionary with extracted hashtags, mentions, links, and text.
    """
    import re

    result: dict[str, Any] = {
        "text": bio_text,
        "hashtags": [],
        "mentions": [],
        "links": [],
        "emojis_used": False,
        "line_count": len(bio_text.strip().split("\n")) if bio_text else 0,
        "character_count": len(bio_text) if bio_text else 0,
    }

    if not bio_text:
        return result

    # Extract hashtags
    result["hashtags"] = re.findall(r"#(\w+)", bio_text)

    # Extract mentions
    result["mentions"] = re.findall(r"@(\w+)", bio_text)

    # Extract URLs
    url_pattern = r"https?://[^\s<>\"')\]]+"
    result["links"] = re.findall(url_pattern, bio_text)

    # Detect emoji usage (basic check for common emoji ranges)
    emoji_pattern = re.compile(
        "["
        "\U0001f600-\U0001f64f"  # emoticons
        "\U0001f300-\U0001f5ff"  # symbols & pictographs
        "\U0001f680-\U0001f6ff"  # transport & map
        "\U0001f1e0-\U0001f1ff"  # flags
        "\U00002702-\U000027b0"
        "\U000024c2-\U0001f251"
        "]+",
        flags=re.UNICODE,
    )
    result["emojis_used"] = bool(emoji_pattern.search(bio_text))

    return result


def classify_post_content(caption: str) -> str:
    """Classify a post caption into a content theme category.

    Args:
        caption: The post caption text.

    Returns:
        One of: educational, promotional, entertainment, inspirational,
        community, behind_the_scenes, news_commentary, values_driven.
    """
    caption_lower = caption.lower() if caption else ""

    # Simple keyword-based classification
    # In production, use an LLM for more accurate classification
    educational_keywords = [
        "how to", "tips", "guide", "learn", "tutorial", "step by step",
        "did you know", "here's how", "pro tip", "explained",
    ]
    promotional_keywords = [
        "shop now", "buy", "sale", "discount", "limited time", "new arrival",
        "launch", "available now", "order", "get yours",
    ]
    entertainment_keywords = [
        "meme", "funny", "lol", "relatable", "mood", "vibes",
        "who else", "tag someone", "caption this",
    ]
    inspirational_keywords = [
        "dream", "believe", "motivation", "inspire", "success story",
        "never give up", "journey", "milestone", "achievement",
    ]
    community_keywords = [
        "thank you", "community", "together", "your story", "share with us",
        "poll", "question", "what do you think", "tell us",
    ]
    bts_keywords = [
        "behind the scenes", "making of", "process", "our team",
        "day in the life", "how we", "sneak peek", "work in progress",
    ]
    news_keywords = [
        "breaking", "announcement", "update", "industry", "trend",
        "report", "study shows", "according to", "news",
    ]
    values_keywords = [
        "sustainability", "impact", "giving back", "cause", "mission",
        "environment", "community support", "diversity", "inclusion",
    ]

    scores = {
        "educational": sum(1 for kw in educational_keywords if kw in caption_lower),
        "promotional": sum(1 for kw in promotional_keywords if kw in caption_lower),
        "entertainment": sum(1 for kw in entertainment_keywords if kw in caption_lower),
        "inspirational": sum(1 for kw in inspirational_keywords if kw in caption_lower),
        "community": sum(1 for kw in community_keywords if kw in caption_lower),
        "behind_the_scenes": sum(1 for kw in bts_keywords if kw in caption_lower),
        "news_commentary": sum(1 for kw in news_keywords if kw in caption_lower),
        "values_driven": sum(1 for kw in values_keywords if kw in caption_lower),
    }

    if max(scores.values()) == 0:
        return "promotional"  # Default if no keywords match

    return max(scores, key=scores.get)


def analyze_posts(posts: list[dict[str, Any]]) -> dict[str, Any]:
    """Analyze a list of post objects for content patterns.

    Args:
        posts: List of post dicts with at minimum 'caption' and optional
               'likes', 'comments', 'shares', 'type', 'hashtags' fields.

    Returns:
        Content analysis dictionary matching the schema.
    """
    if not posts:
        return {
            "recent_posts_analyzed": 0,
            "content_themes": [],
            "content_mix": {},
            "formats_used": [],
        }

    theme_counts: dict[str, int] = {}
    all_hashtags: list[str] = []
    caption_lengths: list[int] = []
    formats: set[str] = set()

    for post in posts:
        caption = post.get("caption", "")
        theme = classify_post_content(caption)
        theme_counts[theme] = theme_counts.get(theme, 0) + 1
        caption_lengths.append(len(caption))

        if post.get("hashtags"):
            all_hashtags.extend(post["hashtags"])

        if post.get("type"):
            formats.add(post["type"])

    total = len(posts)
    content_mix = {k: round(v / total * 100, 1) for k, v in theme_counts.items()}

    # Count hashtag frequency
    hashtag_freq: dict[str, int] = {}
    for tag in all_hashtags:
        hashtag_freq[tag] = hashtag_freq.get(tag, 0) + 1
    top_hashtags = sorted(hashtag_freq.items(), key=lambda x: x[1], reverse=True)[:20]

    return {
        "recent_posts_analyzed": total,
        "content_themes": sorted(theme_counts, key=theme_counts.get, reverse=True),
        "content_mix": content_mix,
        "formats_used": list(formats),
        "avg_caption_length": round(sum(caption_lengths) / len(caption_lengths)) if caption_lengths else None,
        "top_hashtags": [{"tag": tag, "count": count} for tag, count in top_hashtags],
    }


def merge_social_analyses(analyses: list[dict[str, Any]]) -> dict[str, Any]:
    """Merge individual platform analyses into a summary.

    Args:
        analyses: List of per-platform analysis dictionaries.

    Returns:
        Merged summary with cross-platform insights.
    """
    summary: dict[str, Any] = {
        "analyzed_at": datetime.now(timezone.utc).isoformat(),
        "platforms_analyzed": len(analyses),
        "profiles": analyses,
        "cross_platform_summary": {
            "total_followers": 0,
            "platforms": [],
            "consistent_bio_themes": [],
            "shared_hashtags": [],
            "dominant_content_type": "",
        },
    }

    all_themes: dict[str, int] = {}
    all_hashtags_cross: dict[str, int] = {}

    for a in analyses:
        platform = a.get("platform", "unknown")
        summary["cross_platform_summary"]["platforms"].append(platform)

        followers = a.get("metrics", {}).get("followers")
        if followers:
            summary["cross_platform_summary"]["total_followers"] += followers

        # Aggregate content themes
        for theme in a.get("content_analysis", {}).get("content_themes", []):
            all_themes[theme] = all_themes.get(theme, 0) + 1

        # Aggregate hashtags
        for ht in a.get("hashtag_strategy", {}).get("top_hashtags", []):
            tag = ht.get("tag", "") if isinstance(ht, dict) else ht
            if tag:
                all_hashtags_cross[tag] = all_hashtags_cross.get(tag, 0) + 1

    if all_themes:
        summary["cross_platform_summary"]["dominant_content_type"] = max(
            all_themes, key=all_themes.get
        )

    # Hashtags appearing on 2+ platforms
    summary["cross_platform_summary"]["shared_hashtags"] = [
        tag for tag, count in all_hashtags_cross.items() if count >= 2
    ]

    return summary


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Analyze social media profiles for brand identity signals.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            '  python mine_social_profiles.py --profiles "instagram:@acme" "linkedin:/company/acme" --output social.json\n'
            "  python mine_social_profiles.py --input raw_data.json --output social.json\n"
            "\n"
            "Profile format: platform:handle\n"
            "Supported platforms: instagram, facebook, linkedin, tiktok, x_twitter, youtube\n"
            "\n"
            "NOTE: Most social platforms require authentication for scraping.\n"
            "This script generates the analysis structure. For actual data\n"
            "collection, use Chrome automation MCP or platform APIs."
        ),
    )

    parser.add_argument(
        "--profiles",
        nargs="+",
        default=[],
        help='Social profiles in "platform:handle" format (e.g., "instagram:@acme")',
    )
    parser.add_argument(
        "--input",
        default="",
        help="Path to pre-collected raw social data JSON (alternative to --profiles)",
    )
    parser.add_argument(
        "--output",
        default="social_analysis.json",
        help="Output JSON file path (default: social_analysis.json)",
    )

    args = parser.parse_args()

    load_env()

    analyses = []

    if args.input:
        # Load pre-collected data
        input_path = Path(args.input)
        if not input_path.exists():
            print(f"Error: Input file not found: {args.input}", file=sys.stderr)
            sys.exit(1)

        with open(input_path) as f:
            raw_data = json.load(f)

        if isinstance(raw_data, list):
            analyses = raw_data
        elif isinstance(raw_data, dict) and "profiles" in raw_data:
            analyses = raw_data["profiles"]
        else:
            analyses = [raw_data]

        print(f"Loaded {len(analyses)} profile(s) from {args.input}", file=sys.stderr)

    elif args.profiles:
        for profile_str in args.profiles:
            try:
                platform, handle = parse_profile_arg(profile_str)
            except ValueError as e:
                print(f"Error: {e}", file=sys.stderr)
                sys.exit(1)

            print(
                f"Creating analysis template for {platform}:{handle}",
                file=sys.stderr,
            )
            analysis = create_empty_analysis(platform, handle)

            # Note about actual data collection
            analysis["_note"] = (
                "This is a template structure. Populate with actual data using "
                "Chrome automation MCP or platform APIs. Fields with null values "
                "need to be filled from real profile data."
            )

            analyses.append(analysis)

        print(
            "\nNOTE: Profile templates created. To populate with real data:",
            file=sys.stderr,
        )
        print(
            "  1. Use Chrome automation MCP to visit each profile URL",
            file=sys.stderr,
        )
        print(
            "  2. Extract bio, metrics, and recent posts",
            file=sys.stderr,
        )
        print(
            "  3. Re-run with --input to analyze the collected data",
            file=sys.stderr,
        )
    else:
        parser.print_help()
        sys.exit(1)

    # Merge into summary
    result = merge_social_analyses(analyses)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        json.dump(result, f, indent=2)

    print(f"\nAnalysis saved to: {output_path}", file=sys.stderr)
    print(
        f"Profiles processed: {result['platforms_analyzed']}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
