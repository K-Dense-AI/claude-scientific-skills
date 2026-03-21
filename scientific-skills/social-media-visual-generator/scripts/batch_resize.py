#!/usr/bin/env python3
"""
Batch-resize a single image to all required platform dimensions.

Uses Pillow (PIL) for high-quality resizing with proper aspect ratio handling.
Supports smart cropping to maintain subject focus.

Usage:
    python batch_resize.py \
        --input hero.png --platforms instagram \
        --formats post,story,reel,profile --output-dir output/resized/

    python batch_resize.py \
        --input photo.jpg --platforms instagram,linkedin,twitter \
        --formats post --output-dir output/all_platforms/

    python batch_resize.py \
        --input banner.png --platforms youtube \
        --formats thumbnail,channel_art,short \
        --output-dir output/yt/ --mode fill

Dependencies:
    pip install Pillow
"""

import argparse
import sys
from pathlib import Path

try:
    from PIL import Image, ImageFilter
except ImportError:
    print("Error: Pillow is required. Install with: pip install Pillow")
    sys.exit(1)

# Complete platform dimension specs
PLATFORM_SPECS = {
    "instagram": {
        "post": (1080, 1080),
        "post_portrait": (1080, 1350),
        "post_landscape": (1080, 566),
        "story": (1080, 1920),
        "reel": (1080, 1920),
        "carousel": (1080, 1080),
        "carousel_portrait": (1080, 1350),
        "profile": (320, 320),
    },
    "facebook": {
        "post": (1200, 630),
        "story": (1080, 1920),
        "cover": (820, 312),
        "event": (1920, 1005),
        "ad": (1200, 628),
        "profile": (360, 360),
    },
    "linkedin": {
        "post": (1200, 627),
        "post_square": (1080, 1080),
        "article": (1200, 644),
        "cover": (1584, 396),
        "profile": (400, 400),
        "carousel": (1080, 1080),
    },
    "tiktok": {
        "video": (1080, 1920),
        "profile": (200, 200),
        "cover": (1080, 1920),
    },
    "twitter": {
        "post": (1600, 900),
        "header": (1500, 500),
        "profile": (400, 400),
    },
    "youtube": {
        "thumbnail": (1280, 720),
        "channel_art": (2560, 1440),
        "short": (1080, 1920),
    },
}


def resize_fit(img: Image.Image, target_w: int, target_h: int, bg_color: tuple = (255, 255, 255)) -> Image.Image:
    """
    Resize to fit within target dimensions, maintaining aspect ratio.
    Adds padding (letterbox/pillarbox) to fill remaining space.
    """
    img_ratio = img.width / img.height
    target_ratio = target_w / target_h

    if img_ratio > target_ratio:
        # Image is wider: fit to width, pad height
        new_w = target_w
        new_h = int(target_w / img_ratio)
    else:
        # Image is taller: fit to height, pad width
        new_h = target_h
        new_w = int(target_h * img_ratio)

    resized = img.resize((new_w, new_h), Image.LANCZOS)

    # Create background and paste centered
    result = Image.new("RGB", (target_w, target_h), bg_color)
    offset_x = (target_w - new_w) // 2
    offset_y = (target_h - new_h) // 2
    result.paste(resized, (offset_x, offset_y))

    return result


def resize_fill(img: Image.Image, target_w: int, target_h: int) -> Image.Image:
    """
    Resize to fill target dimensions, cropping excess.
    Centers the crop to maintain subject focus.
    """
    img_ratio = img.width / img.height
    target_ratio = target_w / target_h

    if img_ratio > target_ratio:
        # Image is wider: fit to height, crop width
        new_h = target_h
        new_w = int(target_h * img_ratio)
    else:
        # Image is taller: fit to width, crop height
        new_w = target_w
        new_h = int(target_w / img_ratio)

    resized = img.resize((new_w, new_h), Image.LANCZOS)

    # Center crop
    left = (new_w - target_w) // 2
    top = (new_h - target_h) // 2
    right = left + target_w
    bottom = top + target_h

    return resized.crop((left, top, right, bottom))


def resize_stretch(img: Image.Image, target_w: int, target_h: int) -> Image.Image:
    """
    Resize to exactly match target dimensions. May distort aspect ratio.
    Use only when exact dimensions are critical and distortion is acceptable.
    """
    return img.resize((target_w, target_h), Image.LANCZOS)


def resize_blur_fill(img: Image.Image, target_w: int, target_h: int) -> Image.Image:
    """
    Resize to fit within target, fill background with a blurred version
    of the original image. Creates a polished, content-aware background.
    """
    # Create blurred background at target size
    bg = resize_fill(img.copy(), target_w, target_h)
    bg = bg.filter(ImageFilter.GaussianBlur(radius=30))

    # Darken the blurred background slightly
    from PIL import ImageEnhance
    enhancer = ImageEnhance.Brightness(bg)
    bg = enhancer.enhance(0.5)

    # Resize original to fit
    img_ratio = img.width / img.height
    target_ratio = target_w / target_h

    if img_ratio > target_ratio:
        new_w = target_w
        new_h = int(target_w / img_ratio)
    else:
        new_h = target_h
        new_w = int(target_h * img_ratio)

    resized = img.resize((new_w, new_h), Image.LANCZOS)

    # Paste centered on blurred background
    offset_x = (target_w - new_w) // 2
    offset_y = (target_h - new_h) // 2
    bg.paste(resized, (offset_x, offset_y))

    return bg


RESIZE_MODES = {
    "fit": resize_fit,
    "fill": resize_fill,
    "stretch": resize_stretch,
    "blur_fill": resize_blur_fill,
}


def main():
    parser = argparse.ArgumentParser(
        description="Batch-resize a single image to all required platform dimensions."
    )
    parser.add_argument(
        "--input", required=True, help="Path to the source image."
    )
    parser.add_argument(
        "--platforms",
        required=True,
        help="Comma-separated list of platforms (e.g., instagram,linkedin,twitter). Use 'all' for all platforms.",
    )
    parser.add_argument(
        "--formats",
        required=True,
        help="Comma-separated list of formats (e.g., post,story,reel). Use 'all' for all formats on selected platforms.",
    )
    parser.add_argument(
        "--output-dir", required=True, help="Output directory for resized images."
    )
    parser.add_argument(
        "--mode",
        choices=list(RESIZE_MODES.keys()),
        default="fill",
        help="Resize mode: 'fit' (letterbox), 'fill' (crop), 'stretch', or 'blur_fill'. Default: fill.",
    )
    parser.add_argument(
        "--quality",
        type=int,
        default=95,
        help="JPEG quality (1-100). Default: 95.",
    )
    parser.add_argument(
        "--output-format",
        choices=["png", "jpg", "webp"],
        default="png",
        help="Output file format. Default: png.",
    )
    parser.add_argument(
        "--bg-color",
        default="255,255,255",
        help="Background color for 'fit' mode as R,G,B. Default: 255,255,255 (white).",
    )

    args = parser.parse_args()

    # Validate input
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)

    # Parse platforms
    if args.platforms.lower() == "all":
        platforms = list(PLATFORM_SPECS.keys())
    else:
        platforms = [p.strip().lower() for p in args.platforms.split(",")]
        for p in platforms:
            if p not in PLATFORM_SPECS:
                print(f"Error: Unknown platform '{p}'. Available: {', '.join(PLATFORM_SPECS.keys())}")
                sys.exit(1)

    # Parse background color
    try:
        bg_color = tuple(int(c.strip()) for c in args.bg_color.split(","))
    except ValueError:
        bg_color = (255, 255, 255)

    # Load source image
    print(f"Loading: {args.input}")
    img = Image.open(input_path)
    if img.mode == "RGBA":
        # Composite onto background for formats that don't support transparency
        bg = Image.new("RGB", img.size, bg_color)
        bg.paste(img, mask=img.split()[3])
        img_rgb = bg
    else:
        img_rgb = img.convert("RGB")

    print(f"Source: {img.width}x{img.height}")

    # Prepare output
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    resize_fn = RESIZE_MODES[args.mode]
    generated = 0

    for platform in platforms:
        specs = PLATFORM_SPECS[platform]

        if args.formats.lower() == "all":
            formats = list(specs.keys())
        else:
            formats = [f.strip().lower() for f in args.formats.split(",")]

        for fmt in formats:
            if fmt not in specs:
                print(f"  Skipping: {platform}/{fmt} (format not available)")
                continue

            target_w, target_h = specs[fmt]
            ext = args.output_format
            filename = f"{platform}_{fmt}_{target_w}x{target_h}.{ext}"
            output_path = output_dir / filename

            print(f"  Generating: {filename} ({target_w}x{target_h}, mode={args.mode})")

            if args.mode == "fit":
                result = resize_fn(img_rgb, target_w, target_h, bg_color)
            else:
                result = resize_fn(img_rgb, target_w, target_h)

            # Save
            save_kwargs = {}
            if ext == "jpg":
                save_kwargs["quality"] = args.quality
                save_kwargs["optimize"] = True
            elif ext == "webp":
                save_kwargs["quality"] = args.quality
            elif ext == "png":
                save_kwargs["optimize"] = True

            result.save(output_path, **save_kwargs)
            generated += 1

    print(f"\nGenerated {generated} resized images in: {args.output_dir}")


if __name__ == "__main__":
    main()
