#!/usr/bin/env python3
"""
Brand Brief Generator

Transforms a brand DNA JSON (from brand-identity-miner) into a
campaign-specific creative brief. Adds campaign goals, tone adjustments,
content focus areas, and platform-specific directives.

Usage:
    python brand_brief_generator.py \
        --brand-dna ./campaign_output/brand_dna.json \
        --campaign-type product_launch \
        --output ./campaign_output/campaign_brief.json
"""

import argparse
import json
import os
import sys
from datetime import datetime


# ---------------------------------------------------------------------------
# Campaign-Specific Adjustments
# ---------------------------------------------------------------------------

CAMPAIGN_TONE_ADJUSTMENTS = {
    "product_launch": {
        "tone_shift": "Add excitement and exclusivity. Build anticipation in early phases, shift to confident celebration at launch.",
        "vocabulary_additions": [
            "introducing",
            "unveiling",
            "first look",
            "behind the scenes",
            "exclusive",
            "launching",
            "available now",
        ],
        "cta_emphasis": "pre-order, sign up for early access, be the first",
        "visual_direction": "Tease with partial reveals, dramatic reveals at launch. High-energy, aspirational imagery.",
    },
    "brand_awareness": {
        "tone_shift": "Lead with value and education. Be helpful and authoritative without pushing sales.",
        "vocabulary_additions": [
            "discover",
            "learn",
            "explore",
            "did you know",
            "here is why",
            "the truth about",
        ],
        "cta_emphasis": "follow for more, save this, share with someone who needs this",
        "visual_direction": "Bright, shareable graphics. Infographics and educational carousels. Consistent brand aesthetic.",
    },
    "engagement_growth": {
        "tone_shift": "Conversational and inclusive. Ask questions, invite participation, celebrate community.",
        "vocabulary_additions": [
            "what do you think",
            "tell us",
            "your turn",
            "community",
            "together",
            "you asked",
            "we listened",
        ],
        "cta_emphasis": "comment below, tag a friend, share your story, vote now",
        "visual_direction": "UGC-style, authentic, less polished. Feature community members. Interactive formats.",
    },
    "lead_generation": {
        "tone_shift": "Authoritative and solution-oriented. Address pain points directly, position as the answer.",
        "vocabulary_additions": [
            "free guide",
            "download",
            "step-by-step",
            "proven",
            "results",
            "transform",
            "get started",
        ],
        "cta_emphasis": "download now, get your free copy, book a demo, start free trial",
        "visual_direction": "Professional, data-driven. Charts, results screenshots, before/after. Clean, trustworthy design.",
    },
    "event_promotion": {
        "tone_shift": "Energetic and urgent. Build FOMO before, capture energy during, extend value after.",
        "vocabulary_additions": [
            "join us",
            "do not miss",
            "limited spots",
            "live",
            "happening now",
            "countdown",
            "save your spot",
        ],
        "cta_emphasis": "register now, get tickets, save your spot, watch live",
        "visual_direction": "Event branding, speaker spotlights, countdown graphics. During: live/raw energy. After: highlight reels.",
    },
    "seasonal": {
        "tone_shift": "Festive and timely. Blend brand identity with seasonal warmth. Create urgency with limited-time framing.",
        "vocabulary_additions": [
            "limited time",
            "seasonal",
            "celebrate",
            "gift",
            "special",
            "this season",
            "only until",
        ],
        "cta_emphasis": "shop now, limited time offer, grab yours before they are gone",
        "visual_direction": "Seasonal colors and motifs integrated with brand palette. Warm, festive, but still on-brand.",
    },
}

PLATFORM_CONTENT_DIRECTIVES = {
    "instagram": {
        "formats": ["Reels (15-60s)", "Carousels (5-10 slides)", "Stories (daily)", "Feed posts"],
        "copy_guidelines": "Hook in first line. Use line breaks. 2,200 char limit. 20-30 hashtags. Emoji-friendly.",
        "visual_specs": {
            "feed": "1080x1080 or 1080x1350",
            "stories": "1080x1920",
            "reels": "1080x1920 (9:16)",
            "carousel": "1080x1080 per slide",
        },
    },
    "linkedin": {
        "formats": ["Text posts", "Document carousels (PDF)", "Articles", "Video"],
        "copy_guidelines": "Strong hook in first 2 lines (visible before 'see more'). Professional tone. 3,000 char limit. 3-5 hashtags.",
        "visual_specs": {
            "post_image": "1200x627",
            "document": "1080x1080 per page (PDF)",
            "video": "1920x1080 or 1080x1080",
        },
    },
    "tiktok": {
        "formats": ["Short videos (15-60s)", "Longer videos (1-3min)", "Duets", "Stitches"],
        "copy_guidelines": "Hook in first 3 seconds of video. Captions short and punchy. 2,200 char limit. Trending hashtags.",
        "visual_specs": {
            "video": "1080x1920 (9:16)",
        },
    },
    "facebook": {
        "formats": ["Feed posts", "Stories", "Reels", "Group posts", "Events"],
        "copy_guidelines": "Conversational tone. 63,206 char limit but keep concise. Links in post body work. 1-3 hashtags.",
        "visual_specs": {
            "feed": "1200x630",
            "stories": "1080x1920",
            "event_cover": "1920x1005",
        },
    },
    "twitter": {
        "formats": ["Tweets", "Threads (5-15 tweets)", "Polls", "Spaces"],
        "copy_guidelines": "280 char limit per tweet. Punchy, opinionated. Thread first tweet is the hook. 1-3 hashtags.",
        "visual_specs": {
            "image": "1200x675",
            "video": "1920x1080",
        },
    },
}


# ---------------------------------------------------------------------------
# Core Functions
# ---------------------------------------------------------------------------

def load_brand_dna(path: str) -> dict:
    """Load and validate brand DNA JSON."""
    if not os.path.exists(path):
        print(f"Error: Brand DNA file not found: {path}", file=sys.stderr)
        sys.exit(1)
    with open(path) as f:
        data = json.load(f)
    required_keys = ["brand_name", "visual_identity", "voice", "audience"]
    missing = [k for k in required_keys if k not in data]
    if missing:
        print(
            f"Warning: Brand DNA is missing keys: {', '.join(missing)}. "
            "Brief will use defaults for missing fields.",
            file=sys.stderr,
        )
    return data


def generate_content_pillars(brand_dna: dict, campaign_type: str) -> list[dict]:
    """Generate campaign-specific content pillars from brand DNA."""
    base_pillars = brand_dna.get("content_pillars", [])
    if not base_pillars:
        base_pillars = ["brand_story", "value_content", "community"]

    pillar_details = []
    for pillar in base_pillars[:5]:
        pillar_details.append({
            "name": pillar,
            "description": f"Content centered on {pillar}",
            "content_ratio": round(1.0 / len(base_pillars[:5]), 2),
            "campaign_angle": CAMPAIGN_TONE_ADJUSTMENTS[campaign_type]["tone_shift"],
        })
    return pillar_details


def generate_brief(
    brand_dna: dict,
    campaign_type: str,
    platforms: list[str] | None = None,
) -> dict:
    """Generate a complete campaign-specific creative brief."""
    adjustments = CAMPAIGN_TONE_ADJUSTMENTS.get(campaign_type, {})

    # Determine platforms from brand DNA audience or default
    if platforms is None:
        audience_platforms = brand_dna.get("audience", {}).get("platforms", [])
        platforms = audience_platforms if audience_platforms else ["instagram", "linkedin"]

    # Build platform directives
    platform_directives = {}
    for platform in platforms:
        if platform in PLATFORM_CONTENT_DIRECTIVES:
            platform_directives[platform] = PLATFORM_CONTENT_DIRECTIVES[platform]

    brief = {
        "brief_id": f"brief_{campaign_type}_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
        "generated_at": datetime.now().isoformat(),
        "campaign_type": campaign_type,
        "brand": {
            "name": brand_dna.get("brand_name", "Unknown Brand"),
            "tagline": brand_dna.get("tagline", ""),
            "industry": brand_dna.get("industry", ""),
            "unique_selling_points": brand_dna.get("unique_selling_points", []),
        },
        "visual_identity": brand_dna.get("visual_identity", {}),
        "voice": {
            **brand_dna.get("voice", {}),
            "campaign_tone_adjustment": adjustments.get("tone_shift", ""),
            "additional_vocabulary": adjustments.get("vocabulary_additions", []),
        },
        "audience": brand_dna.get("audience", {}),
        "campaign_specifics": {
            "cta_emphasis": adjustments.get("cta_emphasis", ""),
            "visual_direction": adjustments.get("visual_direction", ""),
            "content_pillars": generate_content_pillars(brand_dna, campaign_type),
            "content_mix": {
                "value": "60% — educational, entertaining, or inspiring content",
                "curated": "30% — industry news, UGC, partner content",
                "promotional": "10% — direct sales, offers, product features",
            },
        },
        "platform_directives": platform_directives,
        "brand_guardrails": {
            "do": brand_dna.get("voice", {}).get("do_list", []),
            "dont": brand_dna.get("voice", {}).get("dont_list", []),
            "required_elements": [
                "Brand colors in all visuals",
                "Consistent typography",
                "Logo placement per brand guidelines",
                "Tone alignment with brand voice",
            ],
        },
        "competitors": brand_dna.get("competitors", []),
    }

    return brief


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate a campaign-specific creative brief from brand DNA."
    )
    parser.add_argument(
        "--brand-dna",
        required=True,
        help="Path to brand_dna.json from brand-identity-miner",
    )
    parser.add_argument(
        "--campaign-type",
        required=True,
        choices=list(CAMPAIGN_TONE_ADJUSTMENTS.keys()),
        help="Type of campaign",
    )
    parser.add_argument(
        "--platforms",
        default=None,
        help="Comma-separated platforms (overrides brand DNA defaults)",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output path for the campaign brief JSON",
    )
    args = parser.parse_args()

    # Load brand DNA
    print(f"Loading brand DNA from: {args.brand_dna}")
    brand_dna = load_brand_dna(args.brand_dna)
    print(f"  Brand: {brand_dna.get('brand_name', 'Unknown')}")

    # Parse platforms
    platforms = None
    if args.platforms:
        platforms = [p.strip().lower() for p in args.platforms.split(",")]

    # Generate brief
    print(f"Generating {args.campaign_type} campaign brief...")
    brief = generate_brief(brand_dna, args.campaign_type, platforms)

    # Write output
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(brief, f, indent=2)
    print(f"Campaign brief written to: {args.output}")

    # Summary
    print(f"\n{'=' * 50}")
    print("CAMPAIGN BRIEF SUMMARY")
    print(f"{'=' * 50}")
    print(f"Brand:          {brief['brand']['name']}")
    print(f"Campaign:       {brief['campaign_type']}")
    print(f"Platforms:      {', '.join(brief['platform_directives'].keys())}")
    print(f"Tone shift:     {brief['voice']['campaign_tone_adjustment'][:80]}...")
    print(f"CTA emphasis:   {brief['campaign_specifics']['cta_emphasis']}")
    print(f"Content pillars: {len(brief['campaign_specifics']['content_pillars'])}")
    for pillar in brief["campaign_specifics"]["content_pillars"]:
        print(f"  - {pillar['name']} ({pillar['content_ratio']:.0%})")
    print(f"{'=' * 50}")


if __name__ == "__main__":
    main()
