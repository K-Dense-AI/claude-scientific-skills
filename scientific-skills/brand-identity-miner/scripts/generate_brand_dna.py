#!/usr/bin/env python3
"""
generate_brand_dna.py - Synthesize all analyses into a unified brand_dna.json.

Merges website analysis, social profile analysis, and competitor analysis
into a single brand DNA document following the schema defined in
references/brand_dna_schema.md.

This script handles data merging, conflict resolution, and structural
validation. Brand archetype classification and content pillar generation
are best handled by LLM analysis of the merged data.

Usage:
    python generate_brand_dna.py --website website_analysis.json --output brand_dna.json
    python generate_brand_dna.py --website website.json --social social.json --competitors competitors.json --output brand_dna.json
    python generate_brand_dna.py --input merged_raw.json --output brand_dna.json

Requirements:
    No external dependencies (stdlib only).
"""

import argparse
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "1.0.0"

VALID_ARCHETYPES = [
    "Innocent", "Explorer", "Sage", "Hero", "Outlaw", "Magician",
    "Regular Guy", "Lover", "Jester", "Caregiver", "Creator", "Ruler",
]

VALID_VOCABULARY_LEVELS = ["technical", "casual", "mixed"]
VALID_FORMALITY_LEVELS = ["formal", "semiformal", "casual"]
VALID_PHOTO_STYLES = [
    "lifestyle", "product-focused", "abstract", "people-centric",
    "illustration", "data-visualization", "stock-minimal", "ugc-style",
]


def load_json_file(path: str) -> dict[str, Any] | None:
    """Load a JSON file, returning None if it doesn't exist or is invalid."""
    file_path = Path(path)
    if not file_path.exists():
        print(f"Warning: File not found: {path}", file=sys.stderr)
        return None
    try:
        with open(file_path) as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError) as e:
        print(f"Warning: Failed to load {path}: {e}", file=sys.stderr)
        return None


def normalize_hex_color(color: str) -> str:
    """Normalize a hex color to uppercase 6-digit format."""
    color = color.strip()
    if not color.startswith("#"):
        return color

    hex_part = color[1:]
    if len(hex_part) == 3:
        hex_part = "".join(c * 2 for c in hex_part)

    return f"#{hex_part.upper()}"


def deduplicate_colors(colors: list[str]) -> list[str]:
    """Deduplicate color list, normalizing hex values."""
    seen: set[str] = set()
    result: list[str] = []
    for c in colors:
        normalized = normalize_hex_color(c)
        if normalized not in seen:
            seen.add(normalized)
            result.append(normalized)
    return result


def extract_visual_identity(
    website: dict[str, Any] | None,
    social: dict[str, Any] | None,
) -> dict[str, Any]:
    """Build the visual_identity section from available data.

    Prefers website data as the canonical source for colors and fonts.
    Social data supplements with observed visual patterns.
    """
    visual = {
        "primary_colors": [],
        "secondary_colors": [],
        "accent_colors": [],
        "fonts": {"heading": "", "body": "", "accent": ""},
        "logo_description": "",
        "photography_style": "",
        "visual_mood": "",
    }

    if website:
        vi = website.get("visual_identity", {})

        # Colors from website
        color_data = vi.get("colors", {})
        all_colors = color_data.get("all_colors", [])
        custom_props = color_data.get("custom_properties", {})

        # Try to use CSS custom properties for primary identification
        for prop_name, prop_value in custom_props.items():
            prop_lower = prop_name.lower()
            if "primary" in prop_lower or "brand" in prop_lower:
                if prop_value.startswith("#"):
                    visual["primary_colors"].append(prop_value)
            elif "secondary" in prop_lower:
                if prop_value.startswith("#"):
                    visual["secondary_colors"].append(prop_value)
            elif "accent" in prop_lower:
                if prop_value.startswith("#"):
                    visual["accent_colors"].append(prop_value)

        # If no custom properties, use the most frequent colors
        if not visual["primary_colors"] and all_colors:
            # First 3 colors as primary
            visual["primary_colors"] = all_colors[:3]
            # Next 3 as secondary
            visual["secondary_colors"] = all_colors[3:6]
            # Remaining as accent
            visual["accent_colors"] = all_colors[6:8]

        # Fonts
        font_data = vi.get("fonts", {})
        css_fonts = font_data.get("css_fonts", [])
        external_fonts = font_data.get("external_fonts", [])

        all_fonts = external_fonts + css_fonts  # Prefer external (explicit choices)
        if len(all_fonts) >= 2:
            visual["fonts"]["heading"] = all_fonts[0]
            visual["fonts"]["body"] = all_fonts[1]
            if len(all_fonts) >= 3:
                visual["fonts"]["accent"] = all_fonts[2]
        elif len(all_fonts) == 1:
            visual["fonts"]["heading"] = all_fonts[0]
            visual["fonts"]["body"] = all_fonts[0]

        # Logo
        logo_data = vi.get("logo", {})
        if logo_data.get("logo_alt"):
            visual["logo_description"] = logo_data["logo_alt"]
        elif logo_data.get("logo_src"):
            visual["logo_description"] = f"Logo image at: {logo_data['logo_src']}"

    # Supplement with social visual data
    if social:
        for profile in social.get("profiles", []):
            visual_patterns = profile.get("visual_patterns", {})
            if visual_patterns.get("dominant_colors") and not visual["primary_colors"]:
                visual["primary_colors"] = visual_patterns["dominant_colors"]

    # Deduplicate colors
    visual["primary_colors"] = deduplicate_colors(visual["primary_colors"])
    visual["secondary_colors"] = deduplicate_colors(visual["secondary_colors"])
    visual["accent_colors"] = deduplicate_colors(visual["accent_colors"])

    return visual


def extract_voice(
    website: dict[str, Any] | None,
    social: dict[str, Any] | None,
) -> dict[str, Any]:
    """Build the voice section from available content data.

    Analyzes text content from website and social profiles to infer
    voice characteristics. The actual tone classification is best
    performed by LLM analysis of the sample text.
    """
    voice = {
        "tone": [],
        "personality_traits": [],
        "vocabulary_level": "mixed",
        "formality": "semiformal",
        "voice_dimensions": {
            "formal_casual": 5,
            "serious_playful": 5,
            "technical_simple": 5,
            "authoritative_friendly": 5,
        },
        "do_say": [],
        "dont_say": [],
        "sample_phrases": [],
    }

    if website:
        content = website.get("content", {})

        # Collect sample phrases from hero text and headings
        for text in content.get("hero_text", [])[:5]:
            if text and len(text) > 10:
                voice["sample_phrases"].append(text)

        for heading in content.get("headings", [])[:5]:
            text = heading.get("text", "") if isinstance(heading, dict) else heading
            if text and len(text) > 10:
                voice["sample_phrases"].append(text)

        # Collect CTA phrases
        for cta in content.get("cta_text", [])[:5]:
            if cta:
                voice["do_say"].append(cta)

        # Basic formality detection from paragraphs
        paragraphs = content.get("paragraphs", [])
        if paragraphs:
            combined = " ".join(paragraphs[:10]).lower()

            # Check for contractions (casual indicator)
            contraction_count = sum(
                1 for w in ["don't", "we're", "you'll", "it's", "can't", "won't"]
                if w in combined
            )
            # Check for formal indicators
            formal_count = sum(
                1 for w in ["therefore", "furthermore", "hereby", "pursuant", "accordingly"]
                if w in combined
            )

            if contraction_count > 3:
                voice["formality"] = "casual"
                voice["voice_dimensions"]["formal_casual"] = 7
            elif formal_count > 2:
                voice["formality"] = "formal"
                voice["voice_dimensions"]["formal_casual"] = 2
            else:
                voice["formality"] = "semiformal"
                voice["voice_dimensions"]["formal_casual"] = 5

            # Check for technical vocabulary
            tech_count = sum(
                1 for w in ["api", "algorithm", "infrastructure", "integration", "scalable", "optimize"]
                if w in combined
            )
            if tech_count > 3:
                voice["vocabulary_level"] = "technical"
                voice["voice_dimensions"]["technical_simple"] = 3
            elif tech_count == 0:
                voice["vocabulary_level"] = "casual"
                voice["voice_dimensions"]["technical_simple"] = 8

    # Supplement with social voice data
    if social:
        for profile in social.get("profiles", []):
            bio = profile.get("bio", {})
            bio_text = bio.get("text", "") if isinstance(bio, dict) else ""
            if bio_text and len(bio_text) > 10:
                voice["sample_phrases"].append(bio_text)

    # Deduplicate sample phrases
    seen: set[str] = set()
    unique_phrases: list[str] = []
    for phrase in voice["sample_phrases"]:
        if phrase not in seen:
            seen.add(phrase)
            unique_phrases.append(phrase)
    voice["sample_phrases"] = unique_phrases[:10]

    return voice


def extract_target_audience(website: dict[str, Any] | None) -> dict[str, Any]:
    """Build the target_audience section.

    Returns a template structure. Actual audience analysis requires
    LLM interpretation of the content signals.
    """
    audience = {
        "demographics": {
            "age_range": "",
            "gender": "all",
            "location": "",
            "income_level": "",
            "education": "",
            "occupation": "",
        },
        "psychographics": {
            "interests": [],
            "values": [],
            "lifestyle": "",
            "media_consumption": [],
        },
        "pain_points": [],
        "aspirations": [],
    }

    # Extract audience signals from meta keywords if available
    if website:
        meta = website.get("meta", {})
        keywords = meta.get("meta_keywords", "")
        if keywords:
            audience["psychographics"]["interests"] = [
                k.strip() for k in keywords.split(",")[:10]
            ]

    return audience


def extract_competitors_section(
    competitor_data: dict[str, Any] | None,
) -> list[dict[str, Any]]:
    """Build the competitors array from competitor analysis data."""
    competitors = []

    if not competitor_data:
        return competitors

    for comp in competitor_data.get("competitors", []):
        entry = {
            "name": comp.get("name", ""),
            "url": comp.get("url", ""),
            "positioning": comp.get("positioning", {}).get("target_market", ""),
            "strengths": comp.get("swot", {}).get("strengths", []),
            "weaknesses": comp.get("swot", {}).get("weaknesses", []),
            "social_profiles": comp.get("social_profiles", {}),
            "estimated_sov": "",
        }

        # Try to get SOV from the share_of_voice section
        sov_data = competitor_data.get("share_of_voice", {})
        for platform_data in sov_data.get("platforms", {}).values():
            for breakdown in platform_data.get("breakdown", []):
                if breakdown.get("name") == entry["name"]:
                    entry["estimated_sov"] = f"{breakdown.get('share_pct', '')}%"
                    break

        competitors.append(entry)

    return competitors


def generate_brand_dna(
    website: dict[str, Any] | None,
    social: dict[str, Any] | None,
    competitor_data: dict[str, Any] | None,
) -> dict[str, Any]:
    """Generate the complete brand_dna.json from all input sources.

    Args:
        website: Website analysis JSON (from mine_website.py).
        social: Social profile analysis JSON (from mine_social_profiles.py).
        competitor_data: Competitor analysis JSON (from analyze_competitors.py).

    Returns:
        Complete brand DNA document following the schema.
    """
    # Extract company name from website data
    company_name = ""
    tagline = ""
    website_url = ""

    if website:
        meta = website.get("meta", {})
        company_name = (
            meta.get("og_title", "")
            or meta.get("page_title", "")
        ).split("|")[0].split("-")[0].strip()

        website_url = website.get("url", "")

        content = website.get("content", {})
        hero = content.get("hero_text", [])
        if hero:
            tagline = hero[0]

    # Build social profiles map
    social_profiles = {
        "instagram": "",
        "facebook": "",
        "linkedin": "",
        "tiktok": "",
        "x_twitter": "",
        "youtube": "",
    }
    if social:
        for profile in social.get("profiles", []):
            platform = profile.get("platform", "")
            handle = profile.get("handle", "")
            if platform and handle:
                social_profiles[platform] = handle

    brand_dna: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "generated_by": "brand-identity-miner",

        "company_name": company_name,
        "tagline": tagline,
        "website_url": website_url,
        "industry": "",  # Requires LLM classification
        "brand_archetype": "",  # Requires LLM classification

        "visual_identity": extract_visual_identity(website, social),
        "voice": extract_voice(website, social),

        "values": [],  # Requires LLM extraction from content
        "mission": "",  # Requires LLM extraction from content
        "unique_selling_propositions": [],  # Requires LLM extraction

        "target_audience": extract_target_audience(website),
        "competitors": extract_competitors_section(competitor_data),
        "differentiators": [],  # Requires LLM analysis

        "content_pillars": [],  # Requires LLM generation

        "social_profiles": social_profiles,

        "_generation_notes": {
            "sources_used": {
                "website": website is not None,
                "social": social is not None,
                "competitors": competitor_data is not None,
            },
            "fields_requiring_llm": [
                "industry",
                "brand_archetype",
                "values",
                "mission",
                "unique_selling_propositions",
                "target_audience (refinement)",
                "differentiators",
                "content_pillars",
                "voice.tone",
                "voice.personality_traits",
                "voice.do_say (refinement)",
                "voice.dont_say",
                "visual_identity.photography_style",
                "visual_identity.visual_mood",
            ],
            "instructions": (
                "This brand DNA was auto-generated from extracted data. "
                "Fields listed in 'fields_requiring_llm' need human or LLM "
                "review to be finalized. Run the brand DNA through an LLM "
                "prompt to classify archetypes, define voice characteristics, "
                "and generate content pillars based on the raw data."
            ),
        },
    }

    return brand_dna


def validate_brand_dna(dna: dict[str, Any]) -> list[str]:
    """Validate a brand DNA document against the schema.

    Returns:
        List of validation warnings (empty if valid).
    """
    warnings: list[str] = []

    # Check required top-level fields
    required_fields = [
        "schema_version", "company_name", "tagline", "website_url",
        "industry", "brand_archetype", "visual_identity", "voice",
        "values", "mission", "target_audience",
    ]
    for field in required_fields:
        if not dna.get(field):
            warnings.append(f"Missing or empty required field: {field}")

    # Validate archetype
    archetype = dna.get("brand_archetype", "")
    if archetype and archetype not in VALID_ARCHETYPES:
        warnings.append(
            f"Invalid brand_archetype: '{archetype}'. "
            f"Must be one of: {', '.join(VALID_ARCHETYPES)}"
        )

    # Validate vocabulary level
    vocab = dna.get("voice", {}).get("vocabulary_level", "")
    if vocab and vocab not in VALID_VOCABULARY_LEVELS:
        warnings.append(f"Invalid vocabulary_level: '{vocab}'")

    # Validate formality
    formality = dna.get("voice", {}).get("formality", "")
    if formality and formality not in VALID_FORMALITY_LEVELS:
        warnings.append(f"Invalid formality: '{formality}'")

    # Validate voice dimensions range
    dims = dna.get("voice", {}).get("voice_dimensions", {})
    for dim_name, dim_val in dims.items():
        if isinstance(dim_val, (int, float)) and not (1 <= dim_val <= 10):
            warnings.append(f"Voice dimension '{dim_name}' out of range: {dim_val} (must be 1-10)")

    # Validate color formats
    for color_key in ["primary_colors", "secondary_colors", "accent_colors"]:
        colors = dna.get("visual_identity", {}).get(color_key, [])
        for color in colors:
            if isinstance(color, str) and color.startswith("#"):
                hex_part = color[1:]
                if len(hex_part) not in (3, 6) or not all(
                    c in "0123456789abcdefABCDEF" for c in hex_part
                ):
                    warnings.append(f"Invalid hex color in {color_key}: {color}")

    # Check content pillars count
    pillars = dna.get("content_pillars", [])
    if pillars and (len(pillars) < 3 or len(pillars) > 5):
        warnings.append(
            f"content_pillars should have 3-5 entries, found {len(pillars)}"
        )

    # Check values count
    values = dna.get("values", [])
    if values and (len(values) < 3 or len(values) > 7):
        warnings.append(f"values should have 3-7 entries, found {len(values)}")

    return warnings


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Generate unified brand_dna.json from analysis inputs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python generate_brand_dna.py --website website.json --output brand_dna.json\n"
            "  python generate_brand_dna.py --website w.json --social s.json --competitors c.json -o dna.json\n"
            "  python generate_brand_dna.py --validate existing_brand_dna.json\n"
            "\n"
            "The generated brand DNA will have some fields left empty that\n"
            "require LLM analysis to fill (archetype, content pillars, etc.).\n"
            "These are listed in the _generation_notes.fields_requiring_llm array."
        ),
    )

    parser.add_argument(
        "--website",
        default="",
        help="Path to website analysis JSON (from mine_website.py)",
    )
    parser.add_argument(
        "--social",
        default="",
        help="Path to social profile analysis JSON (from mine_social_profiles.py)",
    )
    parser.add_argument(
        "--competitors",
        default="",
        help="Path to competitor analysis JSON (from analyze_competitors.py)",
    )
    parser.add_argument(
        "--input",
        default="",
        help="Path to a single merged input JSON (alternative to separate files)",
    )
    parser.add_argument(
        "--output", "-o",
        default="brand_dna.json",
        help="Output JSON file path (default: brand_dna.json)",
    )
    parser.add_argument(
        "--validate",
        default="",
        help="Validate an existing brand_dna.json file and report issues",
    )

    args = parser.parse_args()

    # Validation mode
    if args.validate:
        dna = load_json_file(args.validate)
        if not dna:
            print(f"Error: Could not load {args.validate}", file=sys.stderr)
            sys.exit(1)

        warnings = validate_brand_dna(dna)
        if warnings:
            print(f"Validation found {len(warnings)} issue(s):", file=sys.stderr)
            for w in warnings:
                print(f"  - {w}", file=sys.stderr)
            sys.exit(1)
        else:
            print("Validation passed. Brand DNA is valid.", file=sys.stderr)
            sys.exit(0)

    # Generation mode
    website = None
    social = None
    competitor_data = None

    if args.input:
        merged = load_json_file(args.input)
        if merged:
            website = merged.get("website")
            social = merged.get("social")
            competitor_data = merged.get("competitors")
            # If it's a flat structure, treat as website data
            if not any([website, social, competitor_data]):
                website = merged
    else:
        if args.website:
            website = load_json_file(args.website)
        if args.social:
            social = load_json_file(args.social)
        if args.competitors:
            competitor_data = load_json_file(args.competitors)

    if not any([website, social, competitor_data]):
        print(
            "Error: No input data provided. Specify --website, --social, "
            "--competitors, or --input.",
            file=sys.stderr,
        )
        parser.print_help()
        sys.exit(1)

    sources = []
    if website:
        sources.append("website")
    if social:
        sources.append("social")
    if competitor_data:
        sources.append("competitors")
    print(f"Generating brand DNA from: {', '.join(sources)}", file=sys.stderr)

    brand_dna = generate_brand_dna(website, social, competitor_data)

    # Validate
    warnings = validate_brand_dna(brand_dna)
    if warnings:
        print(f"\nValidation warnings ({len(warnings)}):", file=sys.stderr)
        for w in warnings:
            print(f"  - {w}", file=sys.stderr)

    # Write output
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        json.dump(brand_dna, f, indent=2)

    print(f"\nBrand DNA saved to: {output_path}", file=sys.stderr)

    llm_fields = brand_dna.get("_generation_notes", {}).get("fields_requiring_llm", [])
    if llm_fields:
        print(
            f"\nNOTE: {len(llm_fields)} fields require LLM analysis to complete.",
            file=sys.stderr,
        )
        print(
            "Run the brand DNA through an LLM to finalize archetype, "
            "content pillars, voice characteristics, etc.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
