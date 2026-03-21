#!/usr/bin/env python3
"""
Generate images and videos using Google Gemini API for social media content.

Supports:
- Image generation via Gemini 2.0 Flash (Imagen 3)
- Image editing (inpainting, style transfer)
- Video generation via Veo 2 (when available)

Usage:
    # Generate an image
    python generate_with_gemini.py --prompt "Modern coffee shop interior, warm tones" \
        --type image --output output/coffee_shop.png

    # Generate with brand DNA context
    python generate_with_gemini.py --prompt "Product showcase on marble surface" \
        --type image --brand-dna brand_dna.json --platform instagram --format post \
        --output output/product.png

    # Edit an existing image
    python generate_with_gemini.py --prompt "Change the background to a sunset" \
        --type edit --input original.png --output output/edited.png

    # Generate a video storyboard frame sequence
    python generate_with_gemini.py --prompt "Smooth camera pan across luxury products" \
        --type video-frames --num-frames 4 --output output/frames/

Environment:
    GOOGLE_API_KEY: Your Gemini API key (also checks .env file)
"""

import argparse
import base64
import json
import os
import sys
from pathlib import Path

try:
    import requests
except ImportError:
    print("Error: 'requests' package required. Install with: pip install requests")
    sys.exit(1)


# ---------------------------------------------------------------------------
# Platform dimensions (duplicated from generate_social_image.py for standalone use)
# ---------------------------------------------------------------------------
PLATFORM_SPECS = {
    "instagram": {
        "post": (1080, 1080), "story": (1080, 1920), "reel": (1080, 1920),
        "carousel": (1080, 1080),
    },
    "facebook": {
        "post": (1200, 630), "story": (1080, 1920), "cover": (820, 312),
    },
    "linkedin": {
        "post": (1200, 627), "article": (1200, 644), "carousel": (1080, 1080),
    },
    "tiktok": {
        "video": (1080, 1920), "cover": (1080, 1920),
    },
    "twitter": {
        "post": (1600, 900),
    },
    "youtube": {
        "thumbnail": (1280, 720), "short": (1080, 1920),
    },
}


def load_env():
    """Load API key from environment or .env file."""
    key = os.environ.get("GOOGLE_API_KEY")
    if key:
        return key

    # Search for .env in current dir, script dir, and parent dirs
    search_dirs = [
        Path.cwd(),
        Path(__file__).parent,
        Path(__file__).parent.parent,
        Path.home(),
    ]
    for d in search_dirs:
        env_file = d / ".env"
        if env_file.exists():
            with open(env_file) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("GOOGLE_API_KEY=") and not line.startswith("#"):
                        key = line.split("=", 1)[1].strip().strip('"').strip("'")
                        if key:
                            return key
    return None


def load_brand_dna(path: str) -> dict:
    """Load brand DNA JSON."""
    p = Path(path)
    if not p.exists():
        print(f"Warning: brand_dna.json not found at {path}. Using defaults.")
        return {}
    with open(p) as f:
        return json.load(f)


def enhance_prompt(base_prompt: str, brand: dict, platform: str = "", fmt: str = "") -> str:
    """Enhance prompt with brand context."""
    parts = []

    # Brand colors
    visual = brand.get("visual_identity", {})
    colors = visual.get("primary_colors", [])
    if colors:
        parts.append(f"Color palette: {', '.join(colors)}.")

    # Photography style
    photo_style = visual.get("photography_style", "")
    if photo_style:
        parts.append(f"Photography style: {photo_style}.")

    # Visual mood
    mood = visual.get("visual_mood", "")
    if mood:
        parts.append(f"Mood: {mood}.")

    # Core prompt
    parts.append(base_prompt.strip())

    # Platform context
    if platform and fmt:
        parts.append(f"Optimized for {platform} {fmt}.")

    # Default: professional quality
    parts.append("Professional quality, high resolution, clean composition.")

    return " ".join(parts)


def generate_image_gemini(api_key: str, prompt: str, output_path: str,
                          width: int = 1024, height: int = 1024) -> bool:
    """Generate an image using Gemini API with Imagen 3."""
    url = f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.0-flash-exp:generateContent?key={api_key}"

    payload = {
        "contents": [
            {
                "parts": [
                    {
                        "text": f"Generate a high-quality image: {prompt}"
                    }
                ]
            }
        ],
        "generationConfig": {
            "responseModalities": ["TEXT", "IMAGE"],
        }
    }

    print(f"Generating image with Gemini...")
    print(f"Prompt: {prompt}")
    print(f"Target size: {width}x{height}")

    try:
        response = requests.post(url, json=payload, timeout=120)
        response.raise_for_status()
        data = response.json()

        # Extract image from response
        candidates = data.get("candidates", [])
        if not candidates:
            print("Error: No candidates in Gemini response.")
            print(f"Response: {json.dumps(data, indent=2)[:500]}")
            return False

        parts = candidates[0].get("content", {}).get("parts", [])
        for part in parts:
            if "inlineData" in part:
                image_data = part["inlineData"]["data"]
                mime_type = part["inlineData"].get("mimeType", "image/png")

                # Determine file extension
                ext = ".png"
                if "jpeg" in mime_type or "jpg" in mime_type:
                    ext = ".jpg"
                elif "webp" in mime_type:
                    ext = ".webp"

                # Ensure output has correct extension
                out = Path(output_path)
                if out.suffix not in [".png", ".jpg", ".jpeg", ".webp"]:
                    out = out.with_suffix(ext)

                out.parent.mkdir(parents=True, exist_ok=True)
                with open(out, "wb") as f:
                    f.write(base64.b64decode(image_data))

                print(f"Image saved to: {out}")
                return True

        # If no image found, print text response
        for part in parts:
            if "text" in part:
                print(f"Gemini response (text only): {part['text'][:300]}")

        print("Error: No image data in response.")
        return False

    except requests.exceptions.RequestException as e:
        print(f"Error calling Gemini API: {e}")
        return False


def edit_image_gemini(api_key: str, prompt: str, input_path: str,
                      output_path: str) -> bool:
    """Edit an existing image using Gemini API."""
    input_file = Path(input_path)
    if not input_file.exists():
        print(f"Error: Input image not found: {input_path}")
        return False

    # Read and encode input image
    with open(input_file, "rb") as f:
        image_bytes = f.read()
    image_b64 = base64.b64encode(image_bytes).decode("utf-8")

    # Determine MIME type
    suffix = input_file.suffix.lower()
    mime_map = {".png": "image/png", ".jpg": "image/jpeg", ".jpeg": "image/jpeg",
                ".webp": "image/webp"}
    mime_type = mime_map.get(suffix, "image/png")

    url = f"https://generativelanguage.googleapis.com/v1beta/models/gemini-2.0-flash-exp:generateContent?key={api_key}"

    payload = {
        "contents": [
            {
                "parts": [
                    {
                        "inlineData": {
                            "mimeType": mime_type,
                            "data": image_b64,
                        }
                    },
                    {
                        "text": f"Edit this image: {prompt}"
                    }
                ]
            }
        ],
        "generationConfig": {
            "responseModalities": ["TEXT", "IMAGE"],
        }
    }

    print(f"Editing image with Gemini...")
    print(f"Input: {input_path}")
    print(f"Edit prompt: {prompt}")

    try:
        response = requests.post(url, json=payload, timeout=120)
        response.raise_for_status()
        data = response.json()

        candidates = data.get("candidates", [])
        if not candidates:
            print("Error: No candidates in response.")
            return False

        parts = candidates[0].get("content", {}).get("parts", [])
        for part in parts:
            if "inlineData" in part:
                out_data = part["inlineData"]["data"]
                out = Path(output_path)
                out.parent.mkdir(parents=True, exist_ok=True)
                with open(out, "wb") as f:
                    f.write(base64.b64decode(out_data))
                print(f"Edited image saved to: {out}")
                return True

        print("Error: No image in edit response.")
        return False

    except requests.exceptions.RequestException as e:
        print(f"Error: {e}")
        return False


def generate_video_frames(api_key: str, prompt: str, num_frames: int,
                          output_dir: str, brand: dict = None) -> bool:
    """Generate a sequence of keyframes for video storyboard using Gemini."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Generating {num_frames} video storyboard frames...")

    frame_prompts = []
    for i in range(num_frames):
        progress = (i / max(num_frames - 1, 1)) * 100
        frame_prompt = (
            f"Frame {i + 1}/{num_frames} of a video sequence ({progress:.0f}% through): "
            f"{prompt}. "
            f"This frame shows the {'opening' if i == 0 else 'closing' if i == num_frames - 1 else 'middle'} "
            f"of the sequence. Cinematic composition, 16:9 aspect ratio."
        )
        frame_prompts.append(frame_prompt)

    success_count = 0
    for i, fp in enumerate(frame_prompts):
        output_path = out_dir / f"frame_{i + 1:03d}.png"
        print(f"\n--- Frame {i + 1}/{num_frames} ---")
        if generate_image_gemini(api_key, fp, str(output_path), 1920, 1080):
            success_count += 1

    print(f"\nGenerated {success_count}/{num_frames} frames in {out_dir}")

    # Save storyboard metadata
    meta = {
        "prompt": prompt,
        "num_frames": num_frames,
        "frames": [
            {"index": i + 1, "file": f"frame_{i + 1:03d}.png", "prompt": fp}
            for i, fp in enumerate(frame_prompts)
        ],
    }
    meta_path = out_dir / "storyboard.json"
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"Storyboard metadata: {meta_path}")

    return success_count > 0


def main():
    parser = argparse.ArgumentParser(
        description="Generate images and videos using Google Gemini API for social media."
    )
    parser.add_argument("--prompt", required=True, help="Generation prompt.")
    parser.add_argument(
        "--type", required=True, choices=["image", "edit", "video-frames"],
        help="Generation type: image, edit, or video-frames."
    )
    parser.add_argument("--output", required=True, help="Output file/directory path.")
    parser.add_argument("--input", default=None, help="Input image for editing.")
    parser.add_argument("--brand-dna", default=None, help="Path to brand_dna.json.")
    parser.add_argument(
        "--platform", default=None,
        choices=list(PLATFORM_SPECS.keys()),
        help="Target platform for dimension optimization."
    )
    parser.add_argument("--format", default=None, dest="fmt", help="Content format.")
    parser.add_argument(
        "--num-frames", type=int, default=4,
        help="Number of frames for video-frames type (default: 4)."
    )

    args = parser.parse_args()

    # Load API key
    api_key = load_env()
    if not api_key:
        print("Error: GOOGLE_API_KEY not found.")
        print("Set it via:")
        print("  1. Environment variable: export GOOGLE_API_KEY=your-key")
        print("  2. .env file: GOOGLE_API_KEY=your-key")
        print("  Get a key at: https://aistudio.google.com/apikey")
        sys.exit(1)

    # Load brand DNA
    brand = {}
    if args.brand_dna:
        brand = load_brand_dna(args.brand_dna)

    # Enhance prompt with brand context
    prompt = enhance_prompt(args.prompt, brand, args.platform or "", args.fmt or "")

    # Determine dimensions
    width, height = 1024, 1024
    if args.platform and args.fmt:
        key = args.platform.lower()
        fmt = args.fmt.lower()
        if key in PLATFORM_SPECS and fmt in PLATFORM_SPECS[key]:
            dims = PLATFORM_SPECS[key][fmt]
            if isinstance(dims, tuple):
                width, height = dims
            elif isinstance(dims, list):
                width, height = dims[0]

    # Execute
    if args.type == "image":
        success = generate_image_gemini(api_key, prompt, args.output, width, height)
    elif args.type == "edit":
        if not args.input:
            print("Error: --input required for edit type.")
            sys.exit(1)
        success = edit_image_gemini(api_key, prompt, args.input, args.output)
    elif args.type == "video-frames":
        success = generate_video_frames(api_key, prompt, args.num_frames, args.output, brand)
    else:
        print(f"Unknown type: {args.type}")
        sys.exit(1)

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
