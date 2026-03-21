#!/usr/bin/env python3
"""
Add brand elements to an existing image: color border, logo watermark, CTA text overlay.

Uses Pillow for all image manipulation. Reads brand configuration from brand_dna.json
to ensure consistent branding across all visual assets.

Usage:
    python add_brand_overlay.py \
        --input photo.png --brand-dna brand_dna.json \
        --logo logo.png --cta-text "Learn More" \
        --output output/branded.png

    python add_brand_overlay.py \
        --input hero.png --brand-dna brand_dna.json \
        --border-width 8 --logo logo.png --logo-position bottom-right \
        --logo-opacity 0.7 --output output/branded_hero.png

    python add_brand_overlay.py \
        --input slide.png --brand-dna brand_dna.json \
        --cta-text "Swipe for more" --cta-position bottom-center \
        --output output/slide_cta.png

Dependencies:
    pip install Pillow
"""

import argparse
import json
import sys
from pathlib import Path

try:
    from PIL import Image, ImageDraw, ImageFont
except ImportError:
    print("Error: Pillow is required. Install with: pip install Pillow")
    sys.exit(1)


def hex_to_rgb(hex_color: str) -> tuple[int, int, int]:
    """Convert hex color string to RGB tuple."""
    hex_color = hex_color.lstrip("#")
    return tuple(int(hex_color[i:i + 2], 16) for i in (0, 2, 4))


def load_brand_dna(path: str) -> dict:
    """Load brand DNA configuration from JSON file."""
    brand_path = Path(path)
    if not brand_path.exists():
        print(f"Warning: brand_dna.json not found at {path}. Using defaults.")
        return {}
    with open(brand_path, "r") as f:
        return json.load(f)


def get_brand_colors(brand: dict) -> dict:
    """Extract and convert brand colors from brand DNA."""
    colors = brand.get("colors", {})
    result = {
        "primary": hex_to_rgb(colors.get("primary", "#1A2744")),
        "secondary": hex_to_rgb(colors.get("secondary", "#FF6B6B")),
        "accent": hex_to_rgb(colors.get("accent", "#F5A623")),
        "neutral_light": hex_to_rgb(colors.get("neutral_light", "#F8F9FA")),
        "neutral_dark": hex_to_rgb(colors.get("neutral_dark", "#2D3436")),
        "background": hex_to_rgb(colors.get("background", "#FFFFFF")),
    }
    return result


def load_font(font_name: str | None, size: int) -> ImageFont.FreeTypeFont:
    """Load a font, falling back to default if not found."""
    if font_name:
        # Try common font paths
        font_paths = [
            f"/usr/share/fonts/truetype/{font_name}.ttf",
            f"/usr/share/fonts/{font_name}.ttf",
            f"/System/Library/Fonts/{font_name}.ttc",
            f"/Library/Fonts/{font_name}.ttf",
            f"/Library/Fonts/{font_name}.ttc",
            f"C:/Windows/Fonts/{font_name}.ttf",
            font_name,  # Try as direct path
        ]
        for fp in font_paths:
            try:
                return ImageFont.truetype(fp, size)
            except (OSError, IOError):
                continue

    # Try system defaults
    default_fonts = [
        "/System/Library/Fonts/Helvetica.ttc",
        "/System/Library/Fonts/SFNSText.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
        "C:/Windows/Fonts/arial.ttf",
    ]
    for fp in default_fonts:
        try:
            return ImageFont.truetype(fp, size)
        except (OSError, IOError):
            continue

    # Last resort: default bitmap font
    return ImageFont.load_default()


def add_border(img: Image.Image, border_width: int, color: tuple[int, int, int]) -> Image.Image:
    """Add a colored border around the image."""
    if border_width <= 0:
        return img

    new_w = img.width + 2 * border_width
    new_h = img.height + 2 * border_width

    bordered = Image.new("RGB", (new_w, new_h), color)
    bordered.paste(img, (border_width, border_width))

    return bordered


def add_inner_border(img: Image.Image, border_width: int, color: tuple[int, int, int]) -> Image.Image:
    """Add a colored border inside the image without changing dimensions."""
    if border_width <= 0:
        return img

    result = img.copy()
    draw = ImageDraw.Draw(result)

    # Draw border rectangles
    w, h = result.size
    for i in range(border_width):
        draw.rectangle([i, i, w - 1 - i, h - 1 - i], outline=color + (255,))

    return result


def add_logo(
    img: Image.Image,
    logo_path: str,
    position: str = "bottom-right",
    max_size_ratio: float = 0.15,
    opacity: float = 0.8,
    padding: int = 20,
) -> Image.Image:
    """Add a logo watermark to the image."""
    logo_file = Path(logo_path)
    if not logo_file.exists():
        print(f"Warning: Logo file not found: {logo_path}. Skipping logo.")
        return img

    logo = Image.open(logo_file).convert("RGBA")

    # Scale logo to max_size_ratio of image dimensions
    max_logo_w = int(img.width * max_size_ratio)
    max_logo_h = int(img.height * max_size_ratio)

    logo_ratio = logo.width / logo.height
    if logo.width / max_logo_w > logo.height / max_logo_h:
        new_w = max_logo_w
        new_h = int(max_logo_w / logo_ratio)
    else:
        new_h = max_logo_h
        new_w = int(max_logo_h * logo_ratio)

    logo = logo.resize((new_w, new_h), Image.LANCZOS)

    # Apply opacity
    if opacity < 1.0:
        alpha = logo.split()[3]
        alpha = alpha.point(lambda p: int(p * opacity))
        logo.putalpha(alpha)

    # Calculate position
    positions = {
        "top-left": (padding, padding),
        "top-center": ((img.width - new_w) // 2, padding),
        "top-right": (img.width - new_w - padding, padding),
        "center": ((img.width - new_w) // 2, (img.height - new_h) // 2),
        "bottom-left": (padding, img.height - new_h - padding),
        "bottom-center": ((img.width - new_w) // 2, img.height - new_h - padding),
        "bottom-right": (img.width - new_w - padding, img.height - new_h - padding),
    }

    pos = positions.get(position, positions["bottom-right"])

    # Composite logo onto image
    result = img.convert("RGBA")
    result.paste(logo, pos, logo)

    return result.convert("RGB")


def add_cta_text(
    img: Image.Image,
    text: str,
    position: str = "bottom-center",
    font_name: str | None = None,
    font_size: int = 0,
    text_color: tuple[int, int, int] = (255, 255, 255),
    bg_color: tuple[int, int, int, int] = (0, 0, 0, 180),
    padding_x: int = 30,
    padding_y: int = 15,
    margin: int = 40,
    border_radius: int = 8,
) -> Image.Image:
    """Add a CTA text overlay with a semi-transparent background."""
    result = img.convert("RGBA")
    draw = ImageDraw.Draw(result)

    # Auto-size font based on image width if not specified
    if font_size == 0:
        font_size = max(24, img.width // 25)

    font = load_font(font_name, font_size)

    # Calculate text size
    bbox = draw.textbbox((0, 0), text, font=font)
    text_w = bbox[2] - bbox[0]
    text_h = bbox[3] - bbox[1]

    # Background box dimensions
    box_w = text_w + 2 * padding_x
    box_h = text_h + 2 * padding_y

    # Calculate position
    positions = {
        "top-left": (margin, margin),
        "top-center": ((img.width - box_w) // 2, margin),
        "top-right": (img.width - box_w - margin, margin),
        "center": ((img.width - box_w) // 2, (img.height - box_h) // 2),
        "bottom-left": (margin, img.height - box_h - margin),
        "bottom-center": ((img.width - box_w) // 2, img.height - box_h - margin),
        "bottom-right": (img.width - box_w - margin, img.height - box_h - margin),
    }

    box_pos = positions.get(position, positions["bottom-center"])

    # Draw background box with rounded corners
    overlay = Image.new("RGBA", result.size, (0, 0, 0, 0))
    overlay_draw = ImageDraw.Draw(overlay)

    box_rect = [
        box_pos[0],
        box_pos[1],
        box_pos[0] + box_w,
        box_pos[1] + box_h,
    ]
    overlay_draw.rounded_rectangle(box_rect, radius=border_radius, fill=bg_color)

    result = Image.alpha_composite(result, overlay)

    # Draw text
    draw = ImageDraw.Draw(result)
    text_pos = (box_pos[0] + padding_x, box_pos[1] + padding_y)
    draw.text(text_pos, text, fill=text_color + (255,), font=font)

    return result.convert("RGB")


def add_gradient_bar(
    img: Image.Image,
    position: str = "bottom",
    height: int = 80,
    color: tuple[int, int, int] = (0, 0, 0),
    opacity: float = 0.6,
) -> Image.Image:
    """Add a gradient bar (useful for text readability on photo backgrounds)."""
    result = img.convert("RGBA")

    gradient = Image.new("RGBA", (img.width, height), (0, 0, 0, 0))
    for y in range(height):
        if position == "bottom":
            alpha = int(255 * opacity * (y / height))
        else:
            alpha = int(255 * opacity * (1 - y / height))

        for x in range(img.width):
            gradient.putpixel((x, y), color + (alpha,))

    if position == "bottom":
        paste_y = img.height - height
    else:
        paste_y = 0

    result.paste(gradient, (0, paste_y), gradient)

    return result.convert("RGB")


def main():
    parser = argparse.ArgumentParser(
        description="Add brand elements to an existing image."
    )
    parser.add_argument(
        "--input", required=True, help="Path to the source image."
    )
    parser.add_argument(
        "--brand-dna",
        default="brand_dna.json",
        help="Path to brand_dna.json file.",
    )
    parser.add_argument(
        "--output", required=True, help="Output path for the branded image."
    )

    # Border options
    parser.add_argument(
        "--border-width",
        type=int,
        default=0,
        help="Border width in pixels. 0 = no border. Default: 0.",
    )
    parser.add_argument(
        "--border-color",
        default="primary",
        help="Border color: 'primary', 'secondary', 'accent', or hex code. Default: primary.",
    )
    parser.add_argument(
        "--inner-border",
        action="store_true",
        help="Draw border inside the image (no dimension change) instead of outside.",
    )

    # Logo options
    parser.add_argument(
        "--logo", default=None, help="Path to logo image file."
    )
    parser.add_argument(
        "--logo-position",
        default="bottom-right",
        choices=["top-left", "top-center", "top-right", "center", "bottom-left", "bottom-center", "bottom-right"],
        help="Logo position. Default: bottom-right.",
    )
    parser.add_argument(
        "--logo-size",
        type=float,
        default=0.15,
        help="Logo max size as ratio of image dimensions. Default: 0.15.",
    )
    parser.add_argument(
        "--logo-opacity",
        type=float,
        default=0.8,
        help="Logo opacity (0.0-1.0). Default: 0.8.",
    )

    # CTA text options
    parser.add_argument(
        "--cta-text", default=None, help="CTA text to overlay on the image."
    )
    parser.add_argument(
        "--cta-position",
        default="bottom-center",
        choices=["top-left", "top-center", "top-right", "center", "bottom-left", "bottom-center", "bottom-right"],
        help="CTA text position. Default: bottom-center.",
    )
    parser.add_argument(
        "--cta-font-size",
        type=int,
        default=0,
        help="CTA font size in pixels. 0 = auto-size. Default: 0.",
    )
    parser.add_argument(
        "--cta-text-color",
        default="neutral_light",
        help="CTA text color: brand key or hex code. Default: neutral_light.",
    )
    parser.add_argument(
        "--cta-bg-color",
        default="primary",
        help="CTA background color: brand key or hex code. Default: primary.",
    )
    parser.add_argument(
        "--cta-bg-opacity",
        type=float,
        default=0.85,
        help="CTA background opacity (0.0-1.0). Default: 0.85.",
    )

    # Gradient bar
    parser.add_argument(
        "--gradient-bar",
        choices=["top", "bottom", "none"],
        default="none",
        help="Add a gradient bar for text readability. Default: none.",
    )
    parser.add_argument(
        "--gradient-height",
        type=int,
        default=0,
        help="Gradient bar height in pixels. 0 = 15%% of image height.",
    )

    args = parser.parse_args()

    # Validate input
    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)

    # Load brand
    brand = load_brand_dna(args.brand_dna)
    colors = get_brand_colors(brand)
    typography = brand.get("typography", {})

    # Load image
    print(f"Loading: {args.input}")
    img = Image.open(input_path).convert("RGB")
    print(f"Dimensions: {img.width}x{img.height}")

    def resolve_color(color_str: str) -> tuple[int, int, int]:
        """Resolve a color string to RGB tuple."""
        if color_str in colors:
            return colors[color_str]
        if color_str.startswith("#"):
            return hex_to_rgb(color_str)
        return colors.get("primary", (26, 39, 68))

    # Apply gradient bar first (under other elements)
    if args.gradient_bar != "none":
        grad_h = args.gradient_height if args.gradient_height > 0 else int(img.height * 0.15)
        img = add_gradient_bar(img, position=args.gradient_bar, height=grad_h)
        print(f"Added {args.gradient_bar} gradient bar ({grad_h}px)")

    # Apply border
    if args.border_width > 0:
        border_color = resolve_color(args.border_color)
        if args.inner_border:
            img = add_inner_border(img, args.border_width, border_color)
            print(f"Added inner border ({args.border_width}px)")
        else:
            img = add_border(img, args.border_width, border_color)
            print(f"Added outer border ({args.border_width}px, new size: {img.width}x{img.height})")

    # Apply logo
    if args.logo:
        img = add_logo(
            img,
            args.logo,
            position=args.logo_position,
            max_size_ratio=args.logo_size,
            opacity=args.logo_opacity,
        )
        print(f"Added logo at {args.logo_position}")

    # Apply CTA text
    if args.cta_text:
        text_color = resolve_color(args.cta_text_color)
        bg_color_rgb = resolve_color(args.cta_bg_color)
        bg_color_rgba = bg_color_rgb + (int(255 * args.cta_bg_opacity),)

        font_name = typography.get("heading_font", None)

        img = add_cta_text(
            img,
            text=args.cta_text,
            position=args.cta_position,
            font_name=font_name,
            font_size=args.cta_font_size,
            text_color=text_color,
            bg_color=bg_color_rgba,
        )
        print(f"Added CTA text: '{args.cta_text}' at {args.cta_position}")

    # Save output
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    save_kwargs = {}
    if output_path.suffix.lower() in [".jpg", ".jpeg"]:
        save_kwargs["quality"] = 95
        save_kwargs["optimize"] = True
    elif output_path.suffix.lower() == ".png":
        save_kwargs["optimize"] = True
    elif output_path.suffix.lower() == ".webp":
        save_kwargs["quality"] = 95

    img.save(output_path, **save_kwargs)
    print(f"\nBranded image saved: {args.output} ({img.width}x{img.height})")


if __name__ == "__main__":
    main()
