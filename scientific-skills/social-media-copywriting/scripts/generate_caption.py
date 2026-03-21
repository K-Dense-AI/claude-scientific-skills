#!/usr/bin/env python3
"""
Generate platform-optimized social media captions.

Produces caption text, hashtag block, CTA, and emoji placements
while enforcing platform character limits.

Usage:
    python generate_caption.py \
        --brief '{"topic": "product launch", "key_message": "New organic skincare", "tone": "excited"}' \
        --platform instagram \
        --language en \
        --brand-dna brand_dna.json \
        --output caption.json
"""

import argparse
import json
import os
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Platform specifications
# ---------------------------------------------------------------------------
PLATFORM_SPECS = {
    "instagram": {
        "caption_limit": 2200,
        "visible_before_fold": 125,
        "optimal_length": {"min": 150, "max": 300},
        "max_hashtags": 30,
        "recommended_hashtags": {"min": 5, "max": 15},
        "link_in_caption": False,
        "emoji_density": "high",
        "tone": "personal, aspirational, visual-first",
    },
    "linkedin": {
        "caption_limit": 3000,
        "visible_before_fold": 140,
        "optimal_length": {"min": 150, "max": 300},
        "max_hashtags": 10,
        "recommended_hashtags": {"min": 3, "max": 5},
        "link_in_caption": True,
        "emoji_density": "low",
        "tone": "professional, insightful, thought-leadership",
    },
    "tiktok": {
        "caption_limit": 4000,
        "visible_before_fold": 55,
        "optimal_length": {"min": 50, "max": 150},
        "max_hashtags": 10,
        "recommended_hashtags": {"min": 3, "max": 5},
        "link_in_caption": False,
        "emoji_density": "medium",
        "tone": "casual, authentic, trend-aware",
    },
    "facebook": {
        "caption_limit": 63206,
        "visible_before_fold": 110,
        "optimal_length": {"min": 80, "max": 250},
        "max_hashtags": 5,
        "recommended_hashtags": {"min": 1, "max": 3},
        "link_in_caption": True,
        "emoji_density": "medium",
        "tone": "conversational, community-oriented",
    },
    "twitter": {
        "caption_limit": 280,
        "visible_before_fold": 280,
        "optimal_length": {"min": 71, "max": 100},
        "max_hashtags": 3,
        "recommended_hashtags": {"min": 1, "max": 2},
        "link_in_caption": True,
        "emoji_density": "low-medium",
        "tone": "witty, concise, opinionated",
    },
}

# ---------------------------------------------------------------------------
# Copywriting frameworks
# ---------------------------------------------------------------------------
FRAMEWORKS = {
    "product_launch": "AIDA",
    "educational": "PAS",
    "testimonial": "BAB",
    "behind_the_scenes": "BAB",
    "promotional": "AIDA",
    "storytelling": "BAB",
    "engagement": "PAS",
    "announcement": "AIDA",
    "how_to": "PAS",
    "question": "4U",
    "milestone": "BAB",
}


def load_json_file(path: str) -> dict:
    """Load a JSON file and return its contents."""
    file_path = Path(path)
    if not file_path.exists():
        print(f"Warning: File not found: {path}", file=sys.stderr)
        return {}
    with open(file_path, "r", encoding="utf-8") as f:
        return json.load(f)


def load_brand_dna(path: str | None) -> dict:
    """Load brand DNA configuration if provided."""
    if not path:
        return {
            "brand_name": "Brand",
            "voice": {
                "tone": ["friendly"],
                "formality": "casual-professional",
                "emoji_density": "medium",
            },
            "language": {"primary": "en", "ro_formality": "tu"},
            "hashtags": {"branded": [], "always_include": [], "never_use": []},
        }
    return load_json_file(path)


def parse_brief(brief_str: str) -> dict:
    """Parse the content brief from JSON string or file path."""
    brief_str = brief_str.strip()
    if brief_str.startswith("{"):
        return json.loads(brief_str)
    return load_json_file(brief_str)


def select_framework(topic: str) -> str:
    """Select the best copywriting framework based on content topic."""
    topic_lower = topic.lower().replace(" ", "_").replace("-", "_")
    for key, framework in FRAMEWORKS.items():
        if key in topic_lower:
            return framework
    return "AIDA"


def build_caption_structure(
    brief: dict,
    platform: str,
    language: str,
    brand_dna: dict,
) -> dict:
    """
    Build the caption structure with hook, body, CTA, and hashtags.

    Returns a dict describing each section and constraints for the AI
    to fill in during the generation step.
    """
    specs = PLATFORM_SPECS.get(platform, PLATFORM_SPECS["instagram"])
    framework = select_framework(brief.get("topic", "general"))

    # Determine language-specific labels
    if language == "ro":
        section_labels = {
            "hook": "HOOK (prima linie, vizibila inainte de pliere)",
            "body": "CONTINUT (livrarea valorii)",
            "cta": "INDEMN LA ACTIUNE",
            "hashtags": "HASHTAG-URI",
        }
    else:
        section_labels = {
            "hook": "HOOK (first line, visible before fold)",
            "body": "BODY (value delivery)",
            "cta": "CALL TO ACTION",
            "hashtags": "HASHTAGS",
        }

    # Build branded hashtag list
    branded = brand_dna.get("hashtags", {}).get("always_include", [])
    never_use = set(brand_dna.get("hashtags", {}).get("never_use", []))

    structure = {
        "platform": platform,
        "language": language,
        "framework": framework,
        "specs": {
            "caption_limit": specs["caption_limit"],
            "visible_before_fold": specs["visible_before_fold"],
            "optimal_length": specs["optimal_length"],
            "max_hashtags": specs["recommended_hashtags"]["max"],
            "emoji_density": specs["emoji_density"],
            "tone": specs["tone"],
        },
        "sections": section_labels,
        "brief": brief,
        "brand_voice": brand_dna.get("voice", {}),
        "branded_hashtags": branded,
        "banned_hashtags": list(never_use),
        "formality": (
            brand_dna.get("language", {}).get("ro_formality", "tu")
            if language == "ro"
            else "n/a"
        ),
    }

    return structure


def generate_caption_prompt(structure: dict) -> dict:
    """
    Produce the final caption scaffold that an AI or human can fill in.

    This function creates a structured output with placeholders and
    guidelines for each section, ready for content generation.
    """
    specs = structure["specs"]
    brief = structure["brief"]
    framework = structure["framework"]
    language = structure["language"]
    platform = structure["platform"]

    # Framework-specific section guidance
    framework_guidance = {
        "AIDA": {
            "hook": "Attention-grabbing statement that stops the scroll",
            "body_sections": [
                "Interest: Expand with benefits or surprising facts",
                "Desire: Social proof or emotional appeal",
            ],
        },
        "PAS": {
            "hook": "Name the specific pain point your audience feels",
            "body_sections": [
                "Agitation: Amplify the frustration, show consequences",
                "Solution: Present the offer as the natural resolution",
            ],
        },
        "BAB": {
            "hook": "Describe the painful or relatable 'before' state",
            "body_sections": [
                "After: Paint the desirable future outcome",
                "Bridge: Show how to get from before to after",
            ],
        },
        "4U": {
            "hook": "Ultra-specific, useful statement with urgency",
            "body_sections": [
                "Expand on uniqueness and specificity",
            ],
        },
    }

    guidance = framework_guidance.get(framework, framework_guidance["AIDA"])

    caption_scaffold = {
        "metadata": {
            "platform": platform,
            "language": language,
            "framework": framework,
            "character_limit": specs["caption_limit"],
            "optimal_length_range": specs["optimal_length"],
            "visible_before_fold": specs["visible_before_fold"],
            "emoji_density": specs["emoji_density"],
            "tone": specs["tone"],
        },
        "brief": brief,
        "caption": {
            "hook": {
                "max_characters": specs["visible_before_fold"],
                "guidance": guidance["hook"],
                "text": f"[WRITE HOOK HERE - max {specs['visible_before_fold']} chars]",
            },
            "body": {
                "sections": guidance["body_sections"],
                "text": "[WRITE BODY HERE]",
            },
            "cta": {
                "text": "[WRITE CTA HERE]",
                "guidance": f"Single clear action appropriate for {platform}",
            },
        },
        "hashtags": {
            "max_count": specs["max_hashtags"],
            "branded": structure["branded_hashtags"],
            "banned": structure["banned_hashtags"],
            "suggested": [],
        },
        "emoji_placements": {
            "density": specs["emoji_density"],
            "suggested_positions": ["hook_end", "bullet_points", "cta"],
        },
    }

    if language == "ro":
        caption_scaffold["localization"] = {
            "formality": structure["formality"],
            "note": "Write naturally in Romanian, not a translation from English",
        }

    return caption_scaffold


def validate_output(caption_data: dict) -> list[str]:
    """Validate the generated caption against platform constraints."""
    warnings = []
    meta = caption_data.get("metadata", {})
    caption = caption_data.get("caption", {})

    # Check hook length
    hook_text = caption.get("hook", {}).get("text", "")
    max_hook = meta.get("visible_before_fold", 125)
    if len(hook_text) > max_hook and not hook_text.startswith("["):
        warnings.append(
            f"Hook exceeds visible length: {len(hook_text)} > {max_hook} chars"
        )

    # Check hashtag count
    hashtags = caption_data.get("hashtags", {})
    max_tags = hashtags.get("max_count", 15)
    suggested = hashtags.get("suggested", [])
    branded = hashtags.get("branded", [])
    total_tags = len(suggested) + len(branded)
    if total_tags > max_tags:
        warnings.append(
            f"Too many hashtags: {total_tags} > {max_tags} recommended"
        )

    return warnings


def main():
    parser = argparse.ArgumentParser(
        description="Generate platform-optimized social media captions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--brief",
        required=True,
        help='Content brief as JSON string or path to JSON file. '
             'Example: \'{"topic": "product launch", "key_message": "...", "tone": "excited"}\'',
    )
    parser.add_argument(
        "--platform",
        required=True,
        choices=list(PLATFORM_SPECS.keys()),
        help="Target social media platform",
    )
    parser.add_argument(
        "--language",
        default="en",
        choices=["en", "ro"],
        help="Output language (default: en)",
    )
    parser.add_argument(
        "--brand-dna",
        default=None,
        help="Path to brand_dna.json for voice consistency",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output file path (default: stdout)",
    )

    args = parser.parse_args()

    # Check for .env if present
    env_path = Path(".env")
    if env_path.exists():
        print("Note: .env file detected in working directory", file=sys.stderr)

    # Load inputs
    brief = parse_brief(args.brief)
    brand_dna = load_brand_dna(args.brand_dna)

    # Build caption structure
    structure = build_caption_structure(
        brief=brief,
        platform=args.platform,
        language=args.language,
        brand_dna=brand_dna,
    )

    # Generate caption scaffold
    caption_data = generate_caption_prompt(structure)

    # Validate
    warnings = validate_output(caption_data)
    if warnings:
        caption_data["warnings"] = warnings
        for w in warnings:
            print(f"Warning: {w}", file=sys.stderr)

    # Output
    output_json = json.dumps(caption_data, indent=2, ensure_ascii=False)

    if args.output:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(output_json)
        print(f"Caption scaffold written to {args.output}", file=sys.stderr)
    else:
        print(output_json)


if __name__ == "__main__":
    main()
