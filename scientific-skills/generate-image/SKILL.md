---
name: generate-image
description: Generate or edit images using AI models (FLUX, Gemini). Use for scientific illustrations, diagrams, schematics, infographics, concept visualizations, and artistic images. Supports image editing to modify existing images (change colors, add/remove elements, style transfer). Useful for figures, posters, and visual explanations.
---

# Generate Image

Generate and edit high-quality images using OpenRouter's image generation models including FLUX.2 Pro and Nano Banana Pro (Gemini 3 Pro).

## Quick Start

Use the `scripts/generate_image.py` script to generate or edit images:

```bash
# Generate a new image
python scripts/generate_image.py "A beautiful sunset over mountains"

# Edit an existing image
python scripts/generate_image.py "Make the sky purple" --input photo.jpg
```

This generates/edits an image and saves it as `generated_image.png` in the current directory.

## API Key Setup

**CRITICAL**: The script requires an OpenRouter API key. Before running, check if the user has configured their API key:

1. Look for a `.env` file in the project directory or parent directories
2. Check for `OPENROUTER_API_KEY=<key>` in the `.env` file
3. If not found, inform the user they need to:
   - Create a `.env` file with `OPENROUTER_API_KEY=your-api-key-here`
   - Or set the environment variable: `export OPENROUTER_API_KEY=your-api-key-here`
   - Get an API key from: https://openrouter.ai/keys

The script will automatically detect the `.env` file and provide clear error messages if the API key is missing.

## Model Selection

**Default model**: `google/gemini-3-pro-image-preview` (high quality, recommended)

**Available models for generation and editing**:
- `google/gemini-3-pro-image-preview` - High quality, supports generation + editing
- `black-forest-labs/flux.2-pro` - Fast, high quality, supports generation + editing

**Generation only**:
- `black-forest-labs/flux.2-dev` - Development version, generation only

Select based on:
- **Quality**: Use gemini-3-pro or flux.2-pro
- **Editing**: Use gemini-3-pro or flux.2-pro (both support image editing)
- **Cost**: Use flux.2-dev for generation only

## Common Usage Patterns

### Basic generation
```bash
python scripts/generate_image.py "Your prompt here"
```

### Specify model
```bash
python scripts/generate_image.py "A cat in space" --model "black-forest-labs/flux.2-pro"
```

### Custom output path
```bash
python scripts/generate_image.py "Abstract art" --output artwork.png
```

### Edit an existing image
```bash
python scripts/generate_image.py "Make the background blue" --input photo.jpg
```

### Edit with a specific model
```bash
python scripts/generate_image.py "Add sunglasses to the person" --input portrait.png --model "black-forest-labs/flux.2-pro"
```

### Edit with custom output
```bash
python scripts/generate_image.py "Remove the text from the image" --input screenshot.png --output cleaned.png
```

### Multiple images
Run the script multiple times with different prompts or output paths:
```bash
python scripts/generate_image.py "Image 1 description" --output image1.png
python scripts/generate_image.py "Image 2 description" --output image2.png
```

## Script Parameters

- `prompt` (required): Text description of the image to generate, or editing instructions
- `--input` or `-i`: Input image path for editing (enables edit mode)
- `--model` or `-m`: OpenRouter model ID (default: google/gemini-3-pro-image-preview)
- `--output` or `-o`: Output file path (default: generated_image.png)
- `--api-key`: OpenRouter API key (overrides .env file)

## Error Handling

The script provides clear error messages for:
- Missing API key (with setup instructions)
- API errors (with status codes)
- Unexpected response formats
- Missing dependencies (requests library)

If the script fails, read the error message and address the issue before retrying.

## Notes

- Images are returned as base64-encoded data URLs and automatically saved as PNG files
- The script supports both `images` and `content` response formats from different OpenRouter models
- Generation time varies by model (typically 5-30 seconds)
- For image editing, the input image is encoded as base64 and sent to the model
- Supported input image formats: PNG, JPEG, GIF, WebP
- Check OpenRouter pricing for cost information: https://openrouter.ai/models

## Image Editing Tips

- Be specific about what changes you want (e.g., "change the sky to sunset colors" vs "edit the sky")
- Reference specific elements in the image when possible
- For best results, use clear and detailed editing instructions
- Both Gemini 3 Pro and FLUX.2 Pro support image editing through OpenRouter
