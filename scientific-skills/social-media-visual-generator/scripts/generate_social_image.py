#!/usr/bin/env python3
"""
Generate a platform-optimized social media image with brand consistency.

Constructs an enhanced prompt by prepending brand colors and style keywords
from brand_dna.json, selects correct dimensions based on platform and format,
then generates the image via the generate-image script or documents fal.ai MCP usage.

Usage:
    python generate_social_image.py \
        --prompt "A person working at a standing desk" \
        --platform instagram --format post \
        --brand-dna brand_dna.json --output output/ig_post.png

    python generate_social_image.py \
        --prompt "Sunset over a modern city skyline" \
        --platform youtube --format thumbnail \
        --brand-dna brand_dna.json --output output/yt_thumb.png \
        --aspect-variant landscape
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

# Platform dimension specs: (width, height)
# When multiple aspect ratios exist, the first is the default.
PLATFORM_SPECS = {
    "instagram": {
        "post": [(1080, 1080), (1080, 1350), (1080, 566)],
        "story": [(1080, 1920)],
        "reel": [(1080, 1920)],
        "carousel": [(1080, 1080), (1080, 1350)],
        "profile": [(320, 320)],
    },
    "facebook": {
        "post": [(1200, 630)],
        "story": [(1080, 1920)],
        "cover": [(820, 312)],
        "event": [(1920, 1005)],
        "ad": [(1200, 628)],
        "profile": [(360, 360)],
    },
    "linkedin": {
        "post": [(1200, 627), (1080, 1080)],
        "article": [(1200, 644)],
        "cover": [(1584, 396)],
        "profile": [(400, 400)],
        "carousel": [(1080, 1080), (1080, 1350)],
    },
    "tiktok": {
        "video": [(1080, 1920)],
        "profile": [(200, 200)],
        "cover": [(1080, 1920)],
    },
    "twitter": {
        "post": [(1600, 900)],
        "header": [(1500, 500)],
        "profile": [(400, 400)],
    },
    "youtube": {
        "thumbnail": [(1280, 720)],
        "channel_art": [(2560, 1440)],
        "short": [(1080, 1920)],
    },
}

ASPECT_VARIANT_MAP = {
    "square": 0,
    "portrait": 1,
    "landscape": 2,
}


def load_brand_dna(path: str) -> dict:
    """Load brand DNA configuration from JSON file."""
    brand_path = Path(path)
    if not brand_path.exists():
        print(f"Warning: brand_dna.json not found at {path}. Using defaults.")
        return {}
    with open(brand_path, "r") as f:
        return json.load(f)


def get_dimensions(platform: str, fmt: str, variant_index: int = 0) -> tuple[int, int]:
    """Get pixel dimensions for a platform and format combination."""
    platform = platform.lower()
    fmt = fmt.lower()

    if platform not in PLATFORM_SPECS:
        print(f"Error: Unknown platform '{platform}'.")
        print(f"Available platforms: {', '.join(PLATFORM_SPECS.keys())}")
        sys.exit(1)

    if fmt not in PLATFORM_SPECS[platform]:
        print(f"Error: Unknown format '{fmt}' for platform '{platform}'.")
        print(f"Available formats: {', '.join(PLATFORM_SPECS[platform].keys())}")
        sys.exit(1)

    variants = PLATFORM_SPECS[platform][fmt]
    idx = min(variant_index, len(variants) - 1)
    return variants[idx]


def build_enhanced_prompt(base_prompt: str, brand: dict) -> str:
    """Construct an enhanced prompt with brand context prepended and appended."""
    parts = []

    # Prepend brand color directive
    color_prefix = brand.get("color_prompt_prefix", "")
    if not color_prefix:
        colors = brand.get("colors", {})
        color_names = []
        for key in ["primary_name", "secondary_name", "accent_name"]:
            if key in colors:
                color_names.append(colors[key])
        if color_names:
            color_prefix = f"Color palette: {', '.join(color_names)}."

    if color_prefix:
        parts.append(color_prefix.strip())

    # Prepend style keywords
    style_suffix = brand.get("style_prompt_suffix", "")
    if not style_suffix:
        style_kw = brand.get("style_keywords", [])
        if style_kw:
            style_suffix = f"Style: {', '.join(style_kw)}."
    if style_suffix:
        parts.append(style_suffix.strip())

    # Core prompt
    parts.append(base_prompt.strip())

    # Append brand mood
    mood = brand.get("mood", "")
    if mood:
        parts.append(f"Brand mood: {mood}.")

    # Default: no text in image unless specified
    if "text" not in base_prompt.lower():
        parts.append("No text in image unless specified.")

    return " ".join(parts)


def build_negative_prompt(brand: dict) -> str:
    """Build negative prompt from brand DNA."""
    negatives = []
    color_neg = brand.get("color_negative_prompt", "")
    if color_neg:
        negatives.append(color_neg)
    style_neg = brand.get("style_negative", "")
    if style_neg:
        negatives.append(style_neg)
    return ", ".join(negatives) if negatives else ""


def generate_image(
    prompt: str,
    width: int,
    height: int,
    output: str,
    negative_prompt: str = "",
    generate_image_script: str | None = None,
) -> None:
    """
    Generate an image using the generate-image skill script.

    If the script path is not provided, prints the command that would be run
    and documents how to use fal.ai MCP as an alternative.
    """
    if generate_image_script and Path(generate_image_script).exists():
        cmd = [
            sys.executable,
            generate_image_script,
            "--prompt",
            prompt,
            "--width",
            str(width),
            "--height",
            str(height),
            "--output",
            output,
        ]
        if negative_prompt:
            cmd.extend(["--negative-prompt", negative_prompt])

        print(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"Image generated: {output}")
            if result.stdout:
                print(result.stdout)
        else:
            print(f"Error generating image: {result.stderr}")
            sys.exit(1)
    else:
        # Output the generation parameters for manual use or MCP integration
        generation_params = {
            "prompt": prompt,
            "negative_prompt": negative_prompt,
            "width": width,
            "height": height,
            "output": output,
        }
        print("Image generation parameters:")
        print(json.dumps(generation_params, indent=2))
        print()
        print("To generate this image, use one of:")
        print()
        print("1. generate-image skill:")
        print(f'   python generate_image.py --prompt "{prompt}" '
              f'--width {width} --height {height} --output {output}')
        print()
        print("2. fal.ai MCP (fal-generate tool):")
        print(f"   Use prompt above with image_size: {width}x{height}")
        print()
        print("3. Gamma MCP (for social cards):")
        print(f'   Use generate tool with format: "social", '
              f'dimensions matching {width}x{height}')


def main():
    parser = argparse.ArgumentParser(
        description="Generate a platform-optimized social media image with brand consistency."
    )
    parser.add_argument(
        "--prompt", required=True, help="Base image generation prompt."
    )
    parser.add_argument(
        "--platform",
        required=True,
        choices=list(PLATFORM_SPECS.keys()),
        help="Target social media platform.",
    )
    parser.add_argument(
        "--format",
        required=True,
        dest="fmt",
        help="Image format (post, story, reel, thumbnail, cover, etc.).",
    )
    parser.add_argument(
        "--brand-dna",
        default="brand_dna.json",
        help="Path to brand_dna.json file.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output file path for the generated image.",
    )
    parser.add_argument(
        "--aspect-variant",
        choices=["square", "portrait", "landscape"],
        default=None,
        help="Aspect ratio variant when multiple options exist (e.g., Instagram post).",
    )
    parser.add_argument(
        "--generate-image-script",
        default=None,
        help="Path to the generate-image skill script. If not provided, outputs params only.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility.",
    )

    args = parser.parse_args()

    # Load brand DNA
    brand = load_brand_dna(args.brand_dna)

    # Determine aspect variant index
    variant_index = 0
    if args.aspect_variant and args.aspect_variant in ASPECT_VARIANT_MAP:
        variant_index = ASPECT_VARIANT_MAP[args.aspect_variant]

    # Get dimensions
    width, height = get_dimensions(args.platform, args.fmt, variant_index)
    print(f"Target: {args.platform} {args.fmt} at {width}x{height}")

    # Build enhanced prompt
    enhanced_prompt = build_enhanced_prompt(args.prompt, brand)
    print(f"Enhanced prompt: {enhanced_prompt}")

    # Build negative prompt
    negative_prompt = build_negative_prompt(brand)
    if negative_prompt:
        print(f"Negative prompt: {negative_prompt}")

    # Ensure output directory exists
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Generate image
    generate_image(
        prompt=enhanced_prompt,
        width=width,
        height=height,
        output=args.output,
        negative_prompt=negative_prompt,
        generate_image_script=args.generate_image_script,
    )


if __name__ == "__main__":
    main()
