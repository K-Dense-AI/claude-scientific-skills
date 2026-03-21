#!/usr/bin/env python3
"""
Generate a multi-slide carousel for social media with brand consistency.

Produces individual slide images or documents Gamma MCP carousel creation.
Maintains visual consistency across all slides through shared brand context.

Usage:
    python generate_carousel.py \
        --slides '[{"title":"Slide 1","body":"Content here"},{"title":"Slide 2","body":"More content"}]' \
        --platform instagram --style minimalist \
        --brand-dna brand_dna.json --output-dir output/carousel/

    python generate_carousel.py \
        --slides-file slides.json \
        --platform linkedin --style corporate \
        --brand-dna brand_dna.json --output-dir output/linkedin_carousel/
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

CAROUSEL_SPECS = {
    "instagram": {
        "square": (1080, 1080),
        "portrait": (1080, 1350),
        "max_slides": 20,
    },
    "linkedin": {
        "square": (1080, 1080),
        "portrait": (1080, 1350),
        "max_slides": 300,  # PDF pages; practical limit ~10-15
    },
    "facebook": {
        "landscape": (1200, 630),
        "max_slides": 10,
    },
    "twitter": {
        "landscape": (1600, 900),
        "max_slides": 4,
    },
}

STYLE_KEYWORDS = {
    "minimalist": "clean, minimal, whitespace, simple shapes, muted tones, flat design",
    "bold": "high contrast, saturated, strong shapes, graphic, impactful",
    "corporate": "professional, polished, structured, muted blues and grays",
    "playful": "bright colors, rounded shapes, whimsical, fun, bouncy",
    "elegant": "refined, gold accents, soft lighting, luxury feel, serif typography",
    "tech": "dark background, neon accents, sleek, futuristic, gradient mesh",
    "organic": "earth tones, natural textures, warm, handcrafted feel",
    "data": "clean charts, infographic style, structured grid, clear hierarchy",
}


def load_brand_dna(path: str) -> dict:
    """Load brand DNA configuration from JSON file."""
    brand_path = Path(path)
    if not brand_path.exists():
        print(f"Warning: brand_dna.json not found at {path}. Using defaults.")
        return {}
    with open(brand_path, "r") as f:
        return json.load(f)


def load_slides(slides_json: str | None, slides_file: str | None) -> list[dict]:
    """Load slide content from JSON string or file."""
    if slides_file:
        with open(slides_file, "r") as f:
            slides = json.load(f)
    elif slides_json:
        slides = json.loads(slides_json)
    else:
        print("Error: Provide either --slides or --slides-file.")
        sys.exit(1)

    if not isinstance(slides, list) or len(slides) == 0:
        print("Error: Slides must be a non-empty JSON array.")
        sys.exit(1)

    return slides


def get_carousel_dimensions(platform: str, aspect: str) -> tuple[int, int]:
    """Get dimensions for carousel slides on a given platform."""
    platform = platform.lower()
    if platform not in CAROUSEL_SPECS:
        print(f"Error: Platform '{platform}' not supported for carousels.")
        print(f"Available: {', '.join(CAROUSEL_SPECS.keys())}")
        sys.exit(1)

    specs = CAROUSEL_SPECS[platform]
    if aspect in specs:
        return specs[aspect]

    # Default to first available non-max_slides key
    for key, val in specs.items():
        if key != "max_slides" and isinstance(val, tuple):
            return val

    return (1080, 1080)


def build_slide_prompt(
    slide: dict,
    slide_index: int,
    total_slides: int,
    brand: dict,
    style: str,
) -> str:
    """Build a generation prompt for an individual carousel slide."""
    parts = []

    # Brand color prefix
    color_prefix = brand.get("color_prompt_prefix", "")
    if color_prefix:
        parts.append(color_prefix.strip())

    # Style keywords
    style_kw = STYLE_KEYWORDS.get(style, style)
    parts.append(f"Style: {style_kw}.")

    # Slide context
    slide_type = "title slide" if slide_index == 0 else "content slide"
    if slide_index == total_slides - 1:
        slide_type = "call-to-action slide"

    parts.append(f"Social media carousel {slide_type}, slide {slide_index + 1} of {total_slides}.")

    # Title and body
    title = slide.get("title", "")
    body = slide.get("body", "")
    if title:
        parts.append(f"Heading: '{title}'.")
    if body:
        parts.append(f"Supporting text concept: {body}.")

    # Visual direction if specified
    visual = slide.get("visual", "")
    if visual:
        parts.append(f"Visual element: {visual}.")

    # Brand suffix
    style_suffix = brand.get("style_prompt_suffix", "")
    if style_suffix:
        parts.append(style_suffix.strip())

    # Consistent layout instruction
    parts.append("Consistent padding and margins. Clean layout suitable for mobile viewing.")

    return " ".join(parts)


def generate_gamma_mcp_config(
    slides: list[dict],
    platform: str,
    style: str,
    brand: dict,
    width: int,
    height: int,
) -> dict:
    """
    Build a Gamma MCP configuration for carousel generation.
    Returns a dict that documents the MCP call parameters.
    """
    # Map dimensions to Gamma dimension codes
    if width == height:
        gamma_dimensions = "1x1"
    elif width < height:
        gamma_dimensions = "4x5"
    else:
        gamma_dimensions = "16x9"

    # Build inputText for Gamma
    input_lines = []
    for i, slide in enumerate(slides):
        title = slide.get("title", f"Slide {i + 1}")
        body = slide.get("body", "")
        input_lines.append(f"## {title}")
        if body:
            input_lines.append(body)
        input_lines.append("")

    input_text = "\n".join(input_lines)

    # Style-to-tone mapping
    tone_map = {
        "minimalist": "professional",
        "bold": "bold",
        "corporate": "professional",
        "playful": "casual",
        "elegant": "professional",
        "tech": "professional",
    }

    return {
        "tool": "Gamma MCP generate",
        "parameters": {
            "inputText": input_text,
            "format": "social",
            "numCards": len(slides),
            "cardOptions": {
                "dimensions": gamma_dimensions,
            },
            "textOptions": {
                "tone": tone_map.get(style, "professional"),
                "amount": "brief",
            },
        },
        "notes": [
            f"Target platform: {platform}",
            f"Pixel dimensions: {width}x{height} per slide",
            f"Style: {style}",
            "Use get_themes to find a matching theme if specific visual style needed.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(
        description="Generate a multi-slide carousel for social media."
    )
    parser.add_argument(
        "--slides",
        default=None,
        help='JSON array of slide objects, e.g., \'[{"title":"...","body":"..."}]\'.',
    )
    parser.add_argument(
        "--slides-file",
        default=None,
        help="Path to a JSON file containing slide array.",
    )
    parser.add_argument(
        "--platform",
        required=True,
        choices=list(CAROUSEL_SPECS.keys()),
        help="Target social media platform.",
    )
    parser.add_argument(
        "--style",
        default="minimalist",
        choices=list(STYLE_KEYWORDS.keys()),
        help="Visual style for the carousel.",
    )
    parser.add_argument(
        "--aspect",
        default="square",
        choices=["square", "portrait", "landscape"],
        help="Aspect ratio for slides.",
    )
    parser.add_argument(
        "--brand-dna",
        default="brand_dna.json",
        help="Path to brand_dna.json file.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory for generated slide images.",
    )
    parser.add_argument(
        "--method",
        choices=["prompts", "gamma"],
        default="prompts",
        help="Generation method: 'prompts' outputs per-slide prompts, 'gamma' outputs Gamma MCP config.",
    )
    parser.add_argument(
        "--generate-image-script",
        default=None,
        help="Path to the generate-image skill script for direct generation.",
    )
    parser.add_argument(
        "--base-seed",
        type=int,
        default=42,
        help="Base seed for consistent style across slides.",
    )

    args = parser.parse_args()

    # Load inputs
    brand = load_brand_dna(args.brand_dna)
    slides = load_slides(args.slides, args.slides_file)

    # Check slide count
    max_slides = CAROUSEL_SPECS[args.platform].get("max_slides", 20)
    if len(slides) > max_slides:
        print(f"Warning: {args.platform} supports max {max_slides} slides. "
              f"You have {len(slides)}. Truncating.")
        slides = slides[:max_slides]

    # Get dimensions
    width, height = get_carousel_dimensions(args.platform, args.aspect)
    print(f"Carousel: {len(slides)} slides at {width}x{height} for {args.platform}")

    # Ensure output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.method == "gamma":
        # Output Gamma MCP configuration
        gamma_config = generate_gamma_mcp_config(
            slides, args.platform, args.style, brand, width, height
        )
        config_path = output_dir / "gamma_config.json"
        with open(config_path, "w") as f:
            json.dump(gamma_config, f, indent=2)
        print(f"\nGamma MCP configuration saved to: {config_path}")
        print(json.dumps(gamma_config, indent=2))
    else:
        # Output per-slide generation prompts
        manifest = {
            "platform": args.platform,
            "style": args.style,
            "dimensions": {"width": width, "height": height},
            "base_seed": args.base_seed,
            "slides": [],
        }

        for i, slide in enumerate(slides):
            prompt = build_slide_prompt(slide, i, len(slides), brand, args.style)
            slide_output = str(output_dir / f"slide_{i + 1:02d}.png")
            seed = args.base_seed + i

            slide_entry = {
                "index": i + 1,
                "title": slide.get("title", f"Slide {i + 1}"),
                "prompt": prompt,
                "output": slide_output,
                "seed": seed,
                "width": width,
                "height": height,
            }
            manifest["slides"].append(slide_entry)

            print(f"\n--- Slide {i + 1}/{len(slides)} ---")
            print(f"Title: {slide.get('title', 'Untitled')}")
            print(f"Prompt: {prompt}")
            print(f"Output: {slide_output}")
            print(f"Seed: {seed}")

            # If generate-image script provided, run it
            if args.generate_image_script and Path(args.generate_image_script).exists():
                cmd = [
                    sys.executable,
                    args.generate_image_script,
                    "--prompt", prompt,
                    "--width", str(width),
                    "--height", str(height),
                    "--output", slide_output,
                    "--seed", str(seed),
                ]
                print(f"Generating slide {i + 1}...")
                result = subprocess.run(cmd, capture_output=True, text=True)
                if result.returncode != 0:
                    print(f"Error on slide {i + 1}: {result.stderr}")

        # Save manifest
        manifest_path = output_dir / "carousel_manifest.json"
        with open(manifest_path, "w") as f:
            json.dump(manifest, f, indent=2)
        print(f"\nManifest saved to: {manifest_path}")


if __name__ == "__main__":
    main()
