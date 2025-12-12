# OpenRouter Integration for MarkItDown

## Overview

This MarkItDown skill has been configured to use **OpenRouter** instead of direct OpenAI API access. OpenRouter provides a unified API gateway to access 100+ AI models from different providers through a single, OpenAI-compatible interface.

## Why OpenRouter?

### Benefits

1. **Multiple Model Access**: Access GPT-4, Claude, Gemini, and 100+ other models through one API
2. **No Vendor Lock-in**: Switch between models without code changes
3. **Competitive Pricing**: Often better rates than going direct
4. **Simple Migration**: OpenAI-compatible API means minimal code changes
5. **Flexible Choice**: Choose the best model for each task

### Popular Models for Image Description

| Model | Provider | Use Case | Vision Support |
|-------|----------|----------|----------------|
| `anthropic/claude-sonnet-4.5` | Anthropic | **Recommended** - Best overall for scientific analysis | ✅ |
| `anthropic/claude-opus-4.5` | Anthropic | Excellent technical analysis | ✅ |
| `openai/gpt-4o` | OpenAI | Strong vision understanding | ✅ |
| `openai/gpt-4-vision` | OpenAI | GPT-4 with vision | ✅ |
| `google/gemini-pro-vision` | Google | Cost-effective option | ✅ |

See https://openrouter.ai/models for the complete list.

## Getting Started

### 1. Get an API Key

1. Visit https://openrouter.ai/keys
2. Sign up or log in
3. Create a new API key
4. Copy the key (starts with `sk-or-v1-...`)

### 2. Set Environment Variable

```bash
# Add to your environment
export OPENROUTER_API_KEY="sk-or-v1-..."

# Make it permanent
echo 'export OPENROUTER_API_KEY="sk-or-v1-..."' >> ~/.zshrc  # macOS
echo 'export OPENROUTER_API_KEY="sk-or-v1-..."' >> ~/.bashrc # Linux

# Reload shell
source ~/.zshrc  # or source ~/.bashrc
```

### 3. Use in Python

```python
from markitdown import MarkItDown
from openai import OpenAI

# Initialize OpenRouter client (OpenAI-compatible)
client = OpenAI(
    api_key="your-openrouter-api-key",  # or use env var
    base_url="https://openrouter.ai/api/v1"
)

# Create MarkItDown with AI support
md = MarkItDown(
    llm_client=client,
    llm_model="anthropic/claude-sonnet-4.5"  # Choose your model
)

# Convert with AI-enhanced descriptions
result = md.convert("presentation.pptx")
print(result.text_content)
```

## Using the Scripts

All skill scripts have been updated to use OpenRouter:

### convert_with_ai.py

```bash
# Set API key
export OPENROUTER_API_KEY="sk-or-v1-..."

# Convert with default model (advanced vision model)
python scripts/convert_with_ai.py paper.pdf output.md --prompt-type scientific

# Use GPT-4o as alternative
python scripts/convert_with_ai.py paper.pdf output.md \
  --model openai/gpt-4o \
  --prompt-type scientific

# Use Gemini Pro Vision (cost-effective)
python scripts/convert_with_ai.py slides.pptx output.md \
  --model google/gemini-pro-vision \
  --prompt-type presentation

# List available prompt types
python scripts/convert_with_ai.py --list-prompts
```

### Choosing the Right Model

```bash
# For scientific papers - use advanced vision model for technical analysis
python scripts/convert_with_ai.py research.pdf output.md \
  --model anthropic/claude-sonnet-4.5 \
  --prompt-type scientific

# For presentations - use advanced vision model
python scripts/convert_with_ai.py slides.pptx output.md \
  --model anthropic/claude-sonnet-4.5 \
  --prompt-type presentation

# For data visualizations - use advanced vision model
python scripts/convert_with_ai.py charts.pdf output.md \
  --model anthropic/claude-sonnet-4.5 \
  --prompt-type data_viz

# For medical images - use advanced vision model for detailed analysis
python scripts/convert_with_ai.py xray.jpg output.md \
  --model anthropic/claude-sonnet-4.5 \
  --prompt-type medical
```

## Code Examples

### Basic Usage

```python
from markitdown import MarkItDown
from openai import OpenAI
import os

# Initialize OpenRouter client
client = OpenAI(
    api_key=os.environ.get("OPENROUTER_API_KEY"),
    base_url="https://openrouter.ai/api/v1"
)

# Use advanced vision model for image descriptions
md = MarkItDown(
    llm_client=client,
    llm_model="anthropic/claude-sonnet-4.5"
)

result = md.convert("document.pptx")
print(result.text_content)
```

### Switching Models Dynamically

```python
from markitdown import MarkItDown
from openai import OpenAI
import os

client = OpenAI(
    api_key=os.environ["OPENROUTER_API_KEY"],
    base_url="https://openrouter.ai/api/v1"
)

# Use different models for different file types
def convert_with_best_model(filepath):
    if filepath.endswith('.pdf'):
        # Use advanced vision model for technical PDFs
        md = MarkItDown(
            llm_client=client,
            llm_model="anthropic/claude-sonnet-4.5",
            llm_prompt="Describe scientific figures with technical precision"
        )
    elif filepath.endswith('.pptx'):
        # Use advanced vision model for presentations
        md = MarkItDown(
            llm_client=client,
            llm_model="anthropic/claude-sonnet-4.5",
            llm_prompt="Describe slide content and visual elements"
        )
    else:
        # Use advanced vision model as default
        md = MarkItDown(
            llm_client=client,
            llm_model="anthropic/claude-sonnet-4.5"
        )
    
    return md.convert(filepath)

# Use it
result = convert_with_best_model("paper.pdf")
```

### Custom Prompts per Model

```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI(
    api_key="your-openrouter-api-key",
    base_url="https://openrouter.ai/api/v1"
)

# Scientific analysis with advanced vision model
scientific_prompt = """
Analyze this scientific figure. Provide:
1. Type of visualization and methodology
2. Quantitative data points and trends
3. Statistical significance
4. Technical interpretation
Be precise and use scientific terminology.
"""

md_scientific = MarkItDown(
    llm_client=client,
    llm_model="anthropic/claude-sonnet-4.5",
    llm_prompt=scientific_prompt
)

# Visual analysis with advanced vision model
visual_prompt = """
Describe this image comprehensively:
1. Main visual elements and composition
2. Colors, layout, and design
3. Text and labels
4. Overall message
"""

md_visual = MarkItDown(
    llm_client=client,
    llm_model="anthropic/claude-sonnet-4.5",
    llm_prompt=visual_prompt
)
```

## Model Comparison

### For Scientific Content

**Recommended: anthropic/claude-sonnet-4.5**
- Excellent at technical analysis
- Superior reasoning capabilities
- Best at understanding scientific figures
- Most detailed and accurate explanations
- Advanced vision capabilities

**Alternative: openai/gpt-4o**
- Good vision understanding
- Fast processing
- Good at charts and graphs

### For Presentations

**Recommended: anthropic/claude-sonnet-4.5**
- Superior vision capabilities
- Excellent at understanding slide layouts
- Fast and reliable
- Best technical comprehension

### For Cost-Effectiveness

**Recommended: google/gemini-pro-vision**
- Lower cost per request
- Good quality
- Fast processing

## Pricing Considerations

OpenRouter pricing varies by model. Check current rates at https://openrouter.ai/models

**Tips for Cost Optimization:**
1. Use advanced vision models for best quality on complex scientific content
2. Use cheaper models (Gemini) for simple images
3. Batch process similar content with the same model
4. Use appropriate prompts to get better results in fewer retries

## Troubleshooting

### API Key Issues

```bash
# Check if key is set
echo $OPENROUTER_API_KEY

# Should show: sk-or-v1-...
# If empty, set it:
export OPENROUTER_API_KEY="sk-or-v1-..."
```

### Model Not Found

If you get a "model not found" error, check:
1. Model name format: `provider/model-name`
2. Model availability: https://openrouter.ai/models
3. Vision support: Ensure model supports vision for image description

### Rate Limits

OpenRouter has rate limits. If you hit them:
1. Add delays between requests
2. Use batch processing scripts with `--workers` parameter
3. Consider upgrading your OpenRouter plan

## Migration Notes

This skill was updated from direct OpenAI API to OpenRouter. Key changes:

1. **Environment Variable**: `OPENAI_API_KEY` → `OPENROUTER_API_KEY`
2. **Client Initialization**: Added `base_url="https://openrouter.ai/api/v1"`
3. **Model Names**: `gpt-4o` → `openai/gpt-4o` (with provider prefix)
4. **Script Updates**: All scripts now use OpenRouter by default

## Resources

- **OpenRouter Website**: https://openrouter.ai
- **Get API Keys**: https://openrouter.ai/keys
- **Model List**: https://openrouter.ai/models
- **Pricing**: https://openrouter.ai/models (click on model for details)
- **Documentation**: https://openrouter.ai/docs
- **Support**: https://openrouter.ai/discord

## Example Workflow

Here's a complete workflow using OpenRouter:

```bash
# 1. Set up API key
export OPENROUTER_API_KEY="sk-or-v1-your-key-here"

# 2. Convert a scientific paper with Claude
python scripts/convert_with_ai.py \
  research_paper.pdf \
  output.md \
  --model anthropic/claude-opus-4.5 \
  --prompt-type scientific

# 3. Convert presentation with GPT-4o
python scripts/convert_with_ai.py \
  talk_slides.pptx \
  slides.md \
  --model openai/gpt-4o \
  --prompt-type presentation

# 4. Batch convert with cost-effective model
python scripts/batch_convert.py \
  images/ \
  markdown_output/ \
  --extensions .jpg .png
```

## Support

For OpenRouter-specific issues:
- Discord: https://openrouter.ai/discord
- Email: support@openrouter.ai

For MarkItDown skill issues:
- Check documentation in this skill directory
- Review examples in `assets/example_usage.md`

