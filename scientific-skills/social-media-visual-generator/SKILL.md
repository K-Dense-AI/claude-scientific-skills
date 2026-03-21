---
name: social-media-visual-generator
description: Generate platform-optimized visual assets for social media including post images, carousels, story graphics, reel storyboards, and video content. Integrates with generate-image skill (OpenRouter/FLUX), fal.ai MCP (image/video generation), and Gamma MCP (carousel/presentation creation). All visuals are brand-consistent via brand_dna.json.
allowed-tools: [Bash, Read, Write, Edit, WebSearch, WebFetch, Agent]
license: MIT license
metadata:
    skill-author: K-Dense Inc.
---

# Social Media Visual Generator

Visual content engine that produces platform-optimized images, carousels, stories, reel storyboards, thumbnails, and video assets for social media. Every output respects brand identity via `brand_dna.json` and meets exact platform dimension requirements.

## When to Use

- Creating a single post image sized for a specific platform
- Building multi-slide carousels (educational, storytelling, listicle)
- Designing story graphics with interactive element placeholders
- Producing reel or TikTok storyboards with shot lists and timing
- Generating thumbnails with high-contrast text overlays
- Batch-resizing a single asset across multiple platforms
- Applying brand overlays (logo, CTA, color borders) to existing images

## Quick Start

```bash
# Generate a branded Instagram post image
python scripts/generate_social_image.py \
  --prompt "Minimalist flat illustration of a person meditating at sunrise" \
  --platform instagram --format post \
  --brand-dna brand_dna.json --output output/ig_post.png

# Build a 5-slide LinkedIn carousel
python scripts/generate_carousel.py \
  --slides '[{"title":"The Problem","body":"80% of startups fail in year one"},{"title":"Root Cause","body":"Poor product-market fit"},{"title":"The Fix","body":"Talk to 100 customers before building"},{"title":"Results","body":"3x higher survival rate"},{"title":"Start Today","body":"Download our free interview template"}]' \
  --platform linkedin --style minimalist \
  --brand-dna brand_dna.json --output-dir output/carousel/

# Generate a 30-second Reel storyboard
python scripts/generate_video_storyboard.py \
  --template 30s --topic "3 morning habits that changed my productivity" \
  --brand-dna brand_dna.json --output output/storyboard.json

# Resize one image for all Instagram formats
python scripts/batch_resize.py \
  --input hero.png --platforms instagram \
  --formats post,story,reel,profile --output-dir output/resized/

# Add brand overlay to an existing image
python scripts/add_brand_overlay.py \
  --input photo.png --brand-dna brand_dna.json \
  --logo logo.png --cta-text "Learn More" --output output/branded.png
```

## Visual Generation Workflows

### 1. Single Post Image

**Pipeline:** Generate base image --> Apply brand overlay --> Resize to platform

1. Construct an enhanced prompt by prepending brand style keywords and color palette from `brand_dna.json`.
2. Generate the image using one of the available tools (see Tool Selection Guide below).
3. Run `add_brand_overlay.py` to apply logo, CTA text, and color accents.
4. Run `batch_resize.py` if the image is needed on multiple platforms.

### 2. Carousel

**Pipeline:** Define slide content --> Generate slides --> Ensure visual consistency

**Option A -- Gamma MCP (recommended for speed):**
Use the Gamma MCP `generate` tool with `format: "social"` and `cardOptions.dimensions: "1x1"` (LinkedIn, Instagram) or `"4x5"` (Instagram portrait). Pass the full slide content as `inputText`.

**Option B -- Multi-image generation:**
Run `generate_carousel.py` to produce individual slide images. The script maintains a consistent color palette and layout grid across all slides.

### 3. Reel / TikTok Storyboard

**Pipeline:** Select template --> Fill topic --> Generate storyboard JSON --> Generate key frames

1. Choose a template duration (15s, 30s, 60s, 90s) from `video_storyboard_templates.md`.
2. Run `generate_video_storyboard.py` with `--template` and `--topic`.
3. The output JSON contains frame-by-frame shot descriptions, timing, text overlays, and transition cues.
4. Optionally generate key frame images from the shot descriptions.

### 4. Story

**Pipeline:** Generate vertical image --> Add interactive element placeholders

- Always use 9:16 aspect ratio (1080x1920).
- Keep critical content within the center safe zone (810x1420, offset 135px from sides, 250px from top/bottom).
- Leave space for interactive elements: poll stickers (top third), question boxes (middle), swipe-up/link areas (bottom 200px).

### 5. Thumbnail

**Pipeline:** Generate attention-grabbing base --> Add text overlay --> Ensure contrast

- Use 1280x720 for YouTube, 1080x1920 for Reels/TikTok covers.
- Apply high-contrast text with stroke or shadow for readability.
- Use expressive close-up imagery or bold graphic elements.
- Text should occupy no more than 30% of the image area.

## Gemini API Integration (Image + Video)

The `generate_with_gemini.py` script provides direct access to Google Gemini for image generation, image editing, and video storyboard frame creation.

### Setup

Set your Gemini API key in `.env`:
```
GOOGLE_API_KEY=your-api-key-here
```

### Image Generation with Gemini
```bash
python scripts/generate_with_gemini.py \
  --prompt "Modern workspace with natural lighting" \
  --type image --platform instagram --format post \
  --brand-dna brand_dna.json --output output/workspace.png
```

### Image Editing with Gemini
```bash
python scripts/generate_with_gemini.py \
  --prompt "Change the background to a gradient of brand colors" \
  --type edit --input original.png --output output/edited.png
```

### Video Storyboard Frames with Gemini
```bash
python scripts/generate_with_gemini.py \
  --prompt "Smooth camera pan across luxury skincare products on marble" \
  --type video-frames --num-frames 6 \
  --brand-dna brand_dna.json --output output/frames/
```

This generates key frames + a `storyboard.json` metadata file for video assembly.

## Tool Selection Guide

| Tool | Best For | Speed | Quality | Cost |
|------|----------|-------|---------|------|
| **Gemini API** (`generate_with_gemini.py`) | General images, edits, video frames | Fast | High | Low (free tier) |
| **generate-image skill** (OpenRouter/FLUX) | Artistic images, specific styles | Medium | Very High | Medium |
| **fal.ai MCP** (`fal-generate`) | Batch generation, video clips | Fast | High | Per-use |
| **Gamma MCP** | Carousels, presentations, social cards | Very Fast | Good | Free |

**Recommendation:** Use Gemini as the default generator (free tier, fast, good quality). Fall back to FLUX for artistic/stylistic needs, fal.ai for video clips, Gamma for carousels.

## Brand Consistency Enforcement

Every image generation prompt is enhanced with brand context before being sent to any generation tool:

1. **Load `brand_dna.json`** -- extract `primary_color`, `secondary_color`, `accent_color`, `font_family`, `style_keywords`, and `mood`.
2. **Prepend style directive** -- e.g., `"Professional minimalist style, color palette: #1A2B3C, #4D5E6F, #7A8B9C. Clean lines, generous whitespace. "` is prepended to the user's prompt.
3. **Append brand suffix** -- e.g., `" Brand mood: innovative and approachable. No text in image unless specified."` is appended.
4. **Consistent seed management** -- when generating multiple images for one campaign, use sequential seeds to maintain stylistic coherence.

## Tool Selection Guide

| Scenario | Tool | Reason |
|---|---|---|
| Single image, any style | `generate-image` skill (FLUX via OpenRouter) | Fast, high quality, good prompt adherence |
| Image with specific style model | fal.ai MCP | Access to multiple models (SDXL, FLUX variants) |
| Multi-slide carousel/deck | Gamma MCP (`format: "social"`) | Native slide layout, text + image composition |
| Video frames / animation | fal.ai MCP (video models) | Supports video generation models |
| Quick brand mockup | Gamma MCP (`format: "presentation"`) | Rapid prototyping with themes |
| Batch image editing | `batch_resize.py` + `add_brand_overlay.py` | Local Pillow-based processing |

## Platform Dimensions Quick Reference

| Platform | Format | Dimensions | Aspect Ratio |
|---|---|---|---|
| Instagram | Square Post | 1080x1080 | 1:1 |
| Instagram | Portrait Post | 1080x1350 | 4:5 |
| Instagram | Landscape Post | 1080x566 | 1.91:1 |
| Instagram | Story / Reel | 1080x1920 | 9:16 |
| Instagram | Carousel | 1080x1080 | 1:1 |
| Facebook | Post | 1200x630 | 1.91:1 |
| Facebook | Story | 1080x1920 | 9:16 |
| Facebook | Cover | 820x312 | 2.63:1 |
| LinkedIn | Post (landscape) | 1200x627 | 1.91:1 |
| LinkedIn | Post (square) | 1080x1080 | 1:1 |
| LinkedIn | Carousel PDF | 1080x1080 | 1:1 |
| LinkedIn | Cover | 1584x396 | 4:1 |
| TikTok | Video / Cover | 1080x1920 | 9:16 |
| Twitter/X | Post | 1600x900 | 16:9 |
| Twitter/X | Header | 1500x500 | 3:1 |
| YouTube | Thumbnail | 1280x720 | 16:9 |
| YouTube | Channel Art | 2560x1440 | 16:9 |
| YouTube | Short | 1080x1920 | 9:16 |

See `references/visual_specs_by_platform.md` for full specs including safe zones, text areas, and bleed margins.
