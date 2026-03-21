#!/usr/bin/env python3
"""
Plan A/B tests for social media content.

Generates test variants with a measurement plan and tracking structure
based on a base content piece and the variable to test.

Usage:
    python ab_test_planner.py \
        --variable hook \
        --base-content base_post.json \
        --num-variants 3 \
        --duration-days 7 \
        --output ab_test_plan.json
"""

import argparse
import json
import sys
from datetime import datetime, timedelta
from pathlib import Path

# --- Measurement frameworks per variable ---
MEASUREMENT_FRAMEWORKS = {
    "hook": {
        "primary_metrics": ["engagement_rate", "video_watch_time_avg_pct"],
        "secondary_metrics": ["saves", "comments"],
        "tertiary_metrics": ["shares"],
        "what_changes": "Opening line or first 3 seconds of video",
        "what_stays_constant": ["visual", "cta", "posting_time", "format", "hashtags"],
        "min_duration_days": 7,
        "min_impressions_per_variant": 1000,
        "min_engagements_per_variant": 30,
    },
    "visual": {
        "primary_metrics": ["impressions", "engagement_rate"],
        "secondary_metrics": ["profile_visits", "saves"],
        "tertiary_metrics": ["follows"],
        "what_changes": "Image style, thumbnail, colors, text overlay",
        "what_stays_constant": ["caption", "cta", "posting_time", "format", "hashtags"],
        "min_duration_days": 7,
        "min_impressions_per_variant": 1000,
        "min_engagements_per_variant": 30,
    },
    "cta": {
        "primary_metrics": ["target_action_rate", "click_through_rate"],
        "secondary_metrics": ["saves", "shares"],
        "tertiary_metrics": ["comments"],
        "what_changes": "Call-to-action text and placement",
        "what_stays_constant": ["hook", "visual", "posting_time", "format", "hashtags"],
        "min_duration_days": 7,
        "min_impressions_per_variant": 1000,
        "min_engagements_per_variant": 30,
    },
    "time": {
        "primary_metrics": ["reach", "impressions"],
        "secondary_metrics": ["engagement_velocity_1h"],
        "tertiary_metrics": ["engagement_rate"],
        "what_changes": "Hour and day of week",
        "what_stays_constant": ["caption", "visual", "cta", "format", "hashtags"],
        "min_duration_days": 14,
        "min_impressions_per_variant": 1000,
        "min_engagements_per_variant": 30,
    },
    "format": {
        "primary_metrics": ["reach", "engagement_rate"],
        "secondary_metrics": ["time_on_content", "saves"],
        "tertiary_metrics": ["shares"],
        "what_changes": "Content format (reel vs carousel vs static etc.)",
        "what_stays_constant": ["topic", "messaging", "posting_time", "hashtags"],
        "min_duration_days": 14,
        "min_impressions_per_variant": 1000,
        "min_engagements_per_variant": 30,
    },
    "caption_length": {
        "primary_metrics": ["engagement_rate", "comments"],
        "secondary_metrics": ["saves", "time_on_post"],
        "tertiary_metrics": ["shares"],
        "what_changes": "Caption length (short <50 words, medium 50-150, long 150+)",
        "what_stays_constant": ["hook_structure", "visual", "cta", "posting_time", "format"],
        "min_duration_days": 7,
        "min_impressions_per_variant": 1000,
        "min_engagements_per_variant": 30,
    },
    "hashtags": {
        "primary_metrics": ["reach_from_hashtags", "impressions"],
        "secondary_metrics": ["new_follower_rate"],
        "tertiary_metrics": ["engagement_rate"],
        "what_changes": "Hashtag set (niche vs broad vs mixed)",
        "what_stays_constant": ["caption", "visual", "cta", "posting_time", "format"],
        "min_duration_days": 14,
        "min_impressions_per_variant": 2000,
        "min_engagements_per_variant": 30,
    },
}

# --- Variant generation hints per variable ---
VARIANT_HINTS = {
    "hook": [
        "Question-based hook (e.g., 'Are you making this mistake?')",
        "Statement hook with bold claim (e.g., 'This changed everything')",
        "Numbered list hook (e.g., '5 things I wish I knew')",
        "Contrarian hook (e.g., 'Stop doing X')",
        "Story hook (e.g., 'Last week something happened...')",
    ],
    "visual": [
        "Clean graphic with bold text overlay",
        "Authentic photo (no text overlay)",
        "UGC-style / raw aesthetic",
        "Dark background with contrast colors",
        "Bright, colorful lifestyle image",
    ],
    "cta": [
        "Save this for later",
        "Share with someone who needs this",
        "Drop a comment with your answer",
        "Link in bio for the full guide",
        "DM me 'START' for the free template",
    ],
    "time": [
        "Early morning (7-8 AM)",
        "Mid-morning (9-10 AM)",
        "Lunch (12-1 PM)",
        "Afternoon (3-4 PM)",
        "Evening (7-8 PM)",
    ],
    "format": [
        "Short-form video (Reel/TikTok)",
        "Carousel (multi-slide)",
        "Static image post",
        "Text-only post",
        "Document/PDF slides",
    ],
    "caption_length": [
        "Short (under 50 words, punchy)",
        "Medium (50-150 words, balanced)",
        "Long (150+ words, detailed storytelling)",
    ],
    "hashtags": [
        "Niche-specific (low volume, high relevance)",
        "Broad/popular (high volume, lower relevance)",
        "Mixed (combination of niche and broad)",
        "Branded + community hashtags only",
        "No hashtags (control for organic reach)",
    ],
}

VARIANT_NAMES = ["CONTROL", "VAR-A", "VAR-B", "VAR-C", "VAR-D"]


def load_base_content(path: str) -> dict:
    """Load base content from JSON file."""
    with open(path, "r") as f:
        return json.load(f)


def generate_test_name(platform: str, variable: str) -> str:
    """Generate a standardized test name."""
    platform_code = {
        "instagram": "IG",
        "facebook": "FB",
        "linkedin": "LI",
        "tiktok": "TT",
    }.get(platform, platform.upper()[:2])

    date_str = datetime.now().strftime("%Y%m%d")
    return f"{platform_code}-{variable.upper()}-{date_str}-01"


def generate_variants(variable: str, base_content: dict, num_variants: int) -> list:
    """Generate test variants from base content."""
    hints = VARIANT_HINTS.get(variable, ["Variant description placeholder"])
    variants = []

    for i in range(num_variants):
        name = VARIANT_NAMES[i] if i < len(VARIANT_NAMES) else f"VAR-{chr(65 + i - 1)}"
        hint = hints[i % len(hints)]

        variant = {
            "variant_name": name,
            "description": hint if i > 0 else "Baseline (current approach)",
            "is_control": i == 0,
        }

        # Copy base content fields
        for key in ["platform", "caption", "visual_description", "cta", "posting_time", "format", "hashtags"]:
            if key in base_content:
                variant[key] = base_content[key]

        # Mark which field is varied
        variant["variable_tested"] = variable
        if i > 0:
            variant[f"{variable}_variant"] = hint

        variants.append(variant)

    return variants


def generate_tracking_structure(test_name: str, variants: list, framework: dict, duration_days: int) -> dict:
    """Generate a tracking spreadsheet structure."""
    start_date = datetime.now().date()
    end_date = start_date + timedelta(days=duration_days)

    rows = []
    for variant in variants:
        row = {
            "test_name": test_name,
            "variant": variant["variant_name"],
            "start_date": start_date.isoformat(),
            "end_date": end_date.isoformat(),
        }
        # Add metric columns (to be filled during test)
        for metric in framework["primary_metrics"] + framework["secondary_metrics"] + framework["tertiary_metrics"]:
            row[metric] = None
        row["impressions"] = None
        row["total_engagements"] = None
        row["notes"] = ""
        rows.append(row)

    return {
        "columns": list(rows[0].keys()) if rows else [],
        "rows": rows,
    }


def main():
    parser = argparse.ArgumentParser(description="Plan A/B tests for social media content.")
    parser.add_argument("--variable", required=True,
                        choices=["hook", "visual", "cta", "time", "format", "caption_length", "hashtags"],
                        help="Variable to test.")
    parser.add_argument("--base-content", required=True, help="Path to base content JSON file.")
    parser.add_argument("--num-variants", type=int, default=3,
                        help="Number of variants including control (default: 3, min: 2, max: 5).")
    parser.add_argument("--duration-days", type=int, default=7,
                        help="Test duration in days (default: 7).")
    parser.add_argument("--output", required=True, help="Output file path (JSON).")

    args = parser.parse_args()

    # Validate
    num_variants = max(2, min(5, args.num_variants))
    framework = MEASUREMENT_FRAMEWORKS.get(args.variable)
    if not framework:
        print(f"Error: Unknown variable '{args.variable}'", file=sys.stderr)
        sys.exit(1)

    # Enforce minimum duration
    min_duration = framework["min_duration_days"]
    if args.duration_days < min_duration:
        print(f"Warning: Minimum recommended duration for {args.variable} tests is {min_duration} days. "
              f"Adjusting from {args.duration_days} to {min_duration}.")
        duration_days = min_duration
    else:
        duration_days = args.duration_days

    # Load base content
    base_content = load_base_content(args.base_content)
    platform = base_content.get("platform", "instagram")

    # Generate test plan
    test_name = generate_test_name(platform, args.variable)
    variants = generate_variants(args.variable, base_content, num_variants)
    tracking = generate_tracking_structure(test_name, variants, framework, duration_days)

    start_date = datetime.now().date()
    end_date = start_date + timedelta(days=duration_days)

    output = {
        "test_name": test_name,
        "variable": args.variable,
        "platform": platform,
        "hypothesis": f"[Enter your hypothesis: which variant do you expect to win and why]",
        "start_date": start_date.isoformat(),
        "end_date": end_date.isoformat(),
        "duration_days": duration_days,
        "num_variants": num_variants,
        "measurement_framework": {
            "what_changes": framework["what_changes"],
            "what_stays_constant": framework["what_stays_constant"],
            "primary_metrics": framework["primary_metrics"],
            "secondary_metrics": framework["secondary_metrics"],
            "tertiary_metrics": framework["tertiary_metrics"],
            "min_impressions_per_variant": framework["min_impressions_per_variant"],
            "min_engagements_per_variant": framework["min_engagements_per_variant"],
        },
        "variants": variants,
        "tracking_structure": tracking,
        "rules": [
            "Test one variable at a time.",
            "Run all variants during the same time period.",
            f"Minimum test duration: {min_duration} days.",
            f"Minimum impressions per variant: {framework['min_impressions_per_variant']}.",
            "Do not declare a winner before the test end date.",
            "Document external factors (algorithm changes, trending topics, etc.).",
        ],
        "result_template": {
            "winner": "[CONTROL / VAR-A / VAR-B / ...]",
            "primary_metric_delta": "[e.g., +38% engagement rate]",
            "statistical_significance": "[p-value or 'directional']",
            "learning": "[What did you learn?]",
            "next_step": "[Roll out winner / re-test / etc.]",
            "external_factors": "[None noted / describe any]",
        },
        "generated_at": datetime.now().isoformat(),
    }

    with open(args.output, "w") as f:
        json.dump(output, f, indent=2)

    print(f"A/B test plan written to {args.output}")
    print(f"\nTest: {test_name}")
    print(f"Variable: {args.variable}")
    print(f"Platform: {platform}")
    print(f"Variants: {num_variants}")
    print(f"Duration: {duration_days} days ({start_date} to {end_date})")
    print(f"\nPrimary metrics: {', '.join(framework['primary_metrics'])}")
    print(f"\nVariants:")
    for v in variants:
        label = "CONTROL" if v["is_control"] else v["variant_name"]
        print(f"  {label}: {v['description']}")


if __name__ == "__main__":
    main()
