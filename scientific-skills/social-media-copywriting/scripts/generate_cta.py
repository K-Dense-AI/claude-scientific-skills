#!/usr/bin/env python3
"""
Select and customize CTAs from the library based on goal and platform.

Usage:
    python generate_cta.py \
        --goal engage \
        --platform instagram \
        --brand-dna brand_dna.json \
        --output cta.json
"""

import argparse
import json
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# CTA Library (embedded for standalone operation)
# ---------------------------------------------------------------------------
CTA_LIBRARY = {
    "engage": {
        "comment": [
            {
                "text": "Drop a {emoji} if you agree",
                "platforms": ["instagram", "tiktok"],
                "context": "Opinion posts, relatable content",
                "customizable": {"emoji": ["🔥", "💯", "🙌", "☕", "❤️"]},
            },
            {
                "text": "Tag someone who needs to see this",
                "platforms": ["instagram", "facebook"],
                "context": "Helpful tips, educational content",
            },
            {
                "text": "What's your take? Comment below",
                "platforms": ["linkedin", "facebook"],
                "context": "Controversial topics, industry debates",
            },
            {
                "text": "Agree or disagree?",
                "platforms": ["linkedin", "twitter"],
                "context": "Bold statements, hot takes",
            },
            {
                "text": "Fill in the blank: ___",
                "platforms": ["facebook", "linkedin"],
                "context": "Community engagement, fun prompts",
            },
            {
                "text": "Wrong answers only",
                "platforms": ["instagram", "tiktok"],
                "context": "Humor, community building",
            },
            {
                "text": "What would you add to this list?",
                "platforms": ["linkedin", "instagram"],
                "context": "Listicle posts, collaborative content",
            },
            {
                "text": "Type YES if you want the free guide",
                "platforms": ["instagram", "facebook"],
                "context": "Lead magnets, freebies",
            },
            {
                "text": "Who else feels this way?",
                "platforms": ["instagram", "facebook"],
                "context": "Relatable struggles and experiences",
            },
            {
                "text": "How many of these have you tried?",
                "platforms": ["instagram", "tiktok"],
                "context": "Listicles, checklists",
            },
        ],
        "save_share": [
            {
                "text": "Save this for later",
                "platforms": ["instagram"],
                "context": "Tutorials, checklists, guides",
            },
            {
                "text": "Bookmark this post",
                "platforms": ["twitter", "instagram"],
                "context": "Reference material, tips",
            },
            {
                "text": "Share this with someone who needs it",
                "platforms": ["instagram", "facebook"],
                "context": "Helpful advice, resources",
            },
            {
                "text": "Repost if you found this helpful",
                "platforms": ["linkedin", "twitter"],
                "context": "Frameworks, actionable advice",
            },
            {
                "text": "Pin this for later",
                "platforms": ["tiktok"],
                "context": "Tutorials, recipes, guides",
            },
            {
                "text": "Screenshot this",
                "platforms": ["instagram", "tiktok"],
                "context": "Quick tips, formulas",
            },
        ],
        "reaction": [
            {
                "text": "Double tap if this is you",
                "platforms": ["instagram"],
                "context": "Relatable memes, lifestyle content",
            },
            {
                "text": "Leave a {emoji} for part 2",
                "platforms": ["tiktok", "instagram"],
                "context": "Series content, cliffhangers",
                "customizable": {"emoji": ["🔥", "❤️", "👇", "✨"]},
            },
        ],
    },
    "traffic": [
        {
            "text": "Link in bio",
            "platforms": ["instagram", "tiktok"],
            "context": "Any post requiring a link",
        },
        {
            "text": "Click the link below",
            "platforms": ["facebook", "linkedin"],
            "context": "Posts with direct links",
        },
        {
            "text": "Head to {url} to learn more",
            "platforms": ["linkedin", "twitter"],
            "context": "Educational content, reports",
            "customizable": {"url": "brand.com"},
        },
        {
            "text": "Full article in the comments",
            "platforms": ["linkedin"],
            "context": "Blog promotion, thought leadership",
        },
        {
            "text": "DM us '{keyword}' for the link",
            "platforms": ["instagram", "tiktok"],
            "context": "Lead generation, exclusivity",
            "customizable": {"keyword": "INFO"},
        },
        {
            "text": "Download the free template",
            "platforms": ["linkedin", "instagram"],
            "context": "Lead magnets",
        },
        {
            "text": "Grab your free copy",
            "platforms": ["instagram", "linkedin", "facebook", "twitter"],
            "context": "Ebooks, guides, downloads",
        },
        {
            "text": "Try it free for 14 days",
            "platforms": ["linkedin", "twitter"],
            "context": "SaaS, subscription products",
        },
        {
            "text": "Register now (it's free)",
            "platforms": ["linkedin", "facebook"],
            "context": "Webinars, events",
        },
        {
            "text": "Watch the full video on YouTube",
            "platforms": ["instagram", "tiktok"],
            "context": "Cross-platform promotion",
        },
    ],
    "convert": [
        {
            "text": "Shop now",
            "platforms": ["instagram", "facebook"],
            "context": "Direct product promotion",
        },
        {
            "text": "Get yours today",
            "platforms": ["instagram", "facebook"],
            "context": "Product launches, restocks",
        },
        {
            "text": "Use code {code} for {discount}% off",
            "platforms": ["instagram", "facebook", "twitter", "tiktok"],
            "context": "Discount promotions",
            "customizable": {"code": "SAVE20", "discount": "20"},
        },
        {
            "text": "Limited spots available",
            "platforms": ["linkedin", "instagram"],
            "context": "Services, courses, workshops",
        },
        {
            "text": "Book your free consultation",
            "platforms": ["linkedin", "facebook"],
            "context": "Service businesses",
        },
        {
            "text": "Start your free trial",
            "platforms": ["linkedin", "twitter"],
            "context": "SaaS products",
        },
        {
            "text": "Add to cart before it's gone",
            "platforms": ["instagram", "facebook"],
            "context": "Scarcity-driven sales",
        },
        {
            "text": "Join the waitlist",
            "platforms": ["instagram", "linkedin", "facebook", "twitter"],
            "context": "Pre-launch, building demand",
        },
        {
            "text": "Get instant access",
            "platforms": ["instagram", "linkedin"],
            "context": "Digital products",
        },
        {
            "text": "Lock in this price",
            "platforms": ["instagram", "facebook"],
            "context": "Pre-launch pricing, early bird",
        },
    ],
    "community": [
        {
            "text": "Join our community",
            "platforms": ["facebook", "linkedin"],
            "context": "Group promotion",
        },
        {
            "text": "Share your story with #{hashtag}",
            "platforms": ["instagram", "facebook"],
            "context": "User-generated content campaigns",
            "customizable": {"hashtag": "BrandStories"},
        },
        {
            "text": "Let's connect",
            "platforms": ["linkedin"],
            "context": "Networking posts",
        },
        {
            "text": "Help us decide",
            "platforms": ["instagram", "twitter"],
            "context": "Polls, community input",
        },
        {
            "text": "We want to hear from you",
            "platforms": ["instagram", "linkedin", "facebook", "twitter"],
            "context": "Feedback requests",
        },
        {
            "text": "DM us anytime",
            "platforms": ["instagram", "tiktok"],
            "context": "Building trust, accessibility",
        },
        {
            "text": "What's your biggest challenge with {topic}?",
            "platforms": ["linkedin", "instagram"],
            "context": "Research, engagement, community building",
            "customizable": {"topic": "marketing"},
        },
        {
            "text": "Introduce yourself in the comments",
            "platforms": ["facebook", "linkedin"],
            "context": "New community posts, welcome threads",
        },
        {
            "text": "Tag your accountability partner",
            "platforms": ["instagram", "facebook"],
            "context": "Challenge content, goal-setting",
        },
    ],
}


def flatten_ctas(goal_data) -> list[dict]:
    """Flatten nested CTA structures into a flat list."""
    if isinstance(goal_data, list):
        return goal_data
    if isinstance(goal_data, dict):
        flat = []
        for subcategory in goal_data.values():
            if isinstance(subcategory, list):
                flat.extend(subcategory)
        return flat
    return []


def filter_by_platform(ctas: list[dict], platform: str) -> list[dict]:
    """Filter CTAs that are suitable for the given platform."""
    return [
        cta for cta in ctas
        if platform in cta.get("platforms", [])
    ]


def customize_cta(cta: dict, brand_dna: dict) -> dict:
    """Apply brand-specific customizations to a CTA."""
    result = dict(cta)
    text = result["text"]

    # Apply customizable fields
    customizable = result.get("customizable", {})
    for key, default_value in customizable.items():
        placeholder = "{" + key + "}"
        if placeholder in text:
            # Check brand_dna for overrides
            brand_value = brand_dna.get("cta_overrides", {}).get(key, default_value)
            if isinstance(brand_value, list):
                brand_value = brand_value[0]
            text = text.replace(placeholder, str(brand_value))

    result["text"] = text
    result["brand"] = brand_dna.get("brand_name", "Brand")
    return result


def load_brand_dna(path: str | None) -> dict:
    """Load brand DNA configuration."""
    if not path:
        return {"brand_name": "Brand"}
    file_path = Path(path)
    if not file_path.exists():
        print(f"Warning: Brand DNA file not found: {path}", file=sys.stderr)
        return {"brand_name": "Brand"}
    with open(file_path, "r", encoding="utf-8") as f:
        return json.load(f)


def main():
    parser = argparse.ArgumentParser(
        description="Select and customize CTAs from the library",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--goal",
        required=True,
        choices=["engage", "traffic", "convert", "community"],
        help="CTA goal category",
    )
    parser.add_argument(
        "--platform",
        required=True,
        choices=["instagram", "linkedin", "tiktok", "facebook", "twitter"],
        help="Target social media platform",
    )
    parser.add_argument(
        "--brand-dna",
        default=None,
        help="Path to brand_dna.json for customization",
    )
    parser.add_argument(
        "--count",
        type=int,
        default=5,
        help="Number of CTA suggestions to return (default: 5)",
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

    # Load brand DNA
    brand_dna = load_brand_dna(args.brand_dna)

    # Get CTAs for the goal
    goal_data = CTA_LIBRARY.get(args.goal, [])
    all_ctas = flatten_ctas(goal_data)

    # Filter by platform
    platform_ctas = filter_by_platform(all_ctas, args.platform)

    if not platform_ctas:
        # Fall back to all CTAs for this goal
        platform_ctas = all_ctas
        print(
            f"Note: No CTAs specifically tagged for {args.platform}, "
            f"showing all {args.goal} CTAs",
            file=sys.stderr,
        )

    # Customize and limit
    customized = [customize_cta(cta, brand_dna) for cta in platform_ctas]
    selected = customized[: args.count]

    output = {
        "goal": args.goal,
        "platform": args.platform,
        "count": len(selected),
        "ctas": selected,
    }

    output_json = json.dumps(output, indent=2, ensure_ascii=False)

    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(output_json)
        print(f"CTAs written to {args.output}", file=sys.stderr)
    else:
        print(output_json)


if __name__ == "__main__":
    main()
