# MarkItDown Installation Guide

## Prerequisites

- Python 3.10 or higher
- pip package manager
- Virtual environment (recommended)

## Basic Installation

### Install All Features (Recommended)

```bash
pip install 'markitdown[all]'
```

This installs support for all file formats and features.

### Install Specific Features

If you only need certain file formats, you can install specific dependencies:

```bash
# PDF support only
pip install 'markitdown[pdf]'

# Office documents
pip install 'markitdown[docx,pptx,xlsx]'

# Multiple formats
pip install 'markitdown[pdf,docx,pptx,xlsx,audio-transcription]'
```

### Install from Source

```bash
git clone https://github.com/microsoft/markitdown.git
cd markitdown
pip install -e 'packages/markitdown[all]'
```

## Optional Dependencies

| Feature | Installation | Use Case |
|---------|--------------|----------|
| All formats | `pip install 'markitdown[all]'` | Everything |
| PDF | `pip install 'markitdown[pdf]'` | PDF documents |
| Word | `pip install 'markitdown[docx]'` | DOCX files |
| PowerPoint | `pip install 'markitdown[pptx]'` | PPTX files |
| Excel (new) | `pip install 'markitdown[xlsx]'` | XLSX files |
| Excel (old) | `pip install 'markitdown[xls]'` | XLS files |
| Outlook | `pip install 'markitdown[outlook]'` | MSG files |
| Azure DI | `pip install 'markitdown[az-doc-intel]'` | Enhanced PDF |
| Audio | `pip install 'markitdown[audio-transcription]'` | WAV/MP3 |
| YouTube | `pip install 'markitdown[youtube-transcription]'` | YouTube videos |

## System Dependencies

### OCR Support (for scanned documents and images)

#### macOS
```bash
brew install tesseract
```

#### Ubuntu/Debian
```bash
sudo apt-get update
sudo apt-get install tesseract-ocr
```

#### Windows
Download from: https://github.com/UB-Mannheim/tesseract/wiki

### Poppler Utils (for advanced PDF operations)

#### macOS
```bash
brew install poppler
```

#### Ubuntu/Debian
```bash
sudo apt-get install poppler-utils
```

## Verification

Test your installation:

```bash
# Check version
python -c "import markitdown; print('MarkItDown installed successfully')"

# Test basic conversion
echo "Test" > test.txt
markitdown test.txt
rm test.txt
```

## Virtual Environment Setup

### Using venv

```bash
# Create virtual environment
python -m venv markitdown-env

# Activate (macOS/Linux)
source markitdown-env/bin/activate

# Activate (Windows)
markitdown-env\Scripts\activate

# Install
pip install 'markitdown[all]'
```

### Using conda

```bash
# Create environment
conda create -n markitdown python=3.12

# Activate
conda activate markitdown

# Install
pip install 'markitdown[all]'
```

### Using uv

```bash
# Create virtual environment
uv venv --python=3.12 .venv

# Activate
source .venv/bin/activate

# Install
uv pip install 'markitdown[all]'
```

## AI Enhancement Setup (Optional)

For AI-powered image descriptions using OpenRouter:

### OpenRouter API

OpenRouter provides unified access to multiple AI models (GPT-4, Claude, Gemini, etc.) through a single API.

```bash
# Install OpenAI SDK (required, already included with markitdown)
pip install openai

# Get API key from https://openrouter.ai/keys

# Set API key
export OPENROUTER_API_KEY="sk-or-v1-..."

# Add to shell profile for persistence
echo 'export OPENROUTER_API_KEY="sk-or-v1-..."' >> ~/.bashrc  # Linux
echo 'export OPENROUTER_API_KEY="sk-or-v1-..."' >> ~/.zshrc   # macOS
```

**Why OpenRouter?**
- Access to 100+ AI models through one API
- Choose between GPT-4, Claude, Gemini, and more
- Competitive pricing
- No vendor lock-in
- Simple OpenAI-compatible interface

**Popular Models for Image Description:**
- `anthropic/claude-sonnet-4.5` - **Recommended** - Best for scientific vision
- `anthropic/claude-opus-4.5` - Excellent technical analysis
- `openai/gpt-4o` - Good vision understanding
- `google/gemini-pro-vision` - Cost-effective option

See https://openrouter.ai/models for complete model list and pricing.

## Azure Document Intelligence Setup (Optional)

For enhanced PDF conversion:

1. Create Azure Document Intelligence resource in Azure Portal
2. Get endpoint and key
3. Set environment variables:

```bash
export AZURE_DOCUMENT_INTELLIGENCE_KEY="your-key"
export AZURE_DOCUMENT_INTELLIGENCE_ENDPOINT="https://your-endpoint.cognitiveservices.azure.com/"
```

## Docker Installation (Alternative)

```bash
# Clone repository
git clone https://github.com/microsoft/markitdown.git
cd markitdown

# Build image
docker build -t markitdown:latest .

# Run
docker run --rm -i markitdown:latest < input.pdf > output.md
```

## Troubleshooting

### Import Error
```
ModuleNotFoundError: No module named 'markitdown'
```

**Solution**: Ensure you're in the correct virtual environment and markitdown is installed:
```bash
pip install 'markitdown[all]'
```

### Missing Feature
```
Error: PDF conversion not supported
```

**Solution**: Install the specific feature:
```bash
pip install 'markitdown[pdf]'
```

### OCR Not Working

**Solution**: Install Tesseract OCR (see System Dependencies above)

### Permission Errors

**Solution**: Use virtual environment or install with `--user` flag:
```bash
pip install --user 'markitdown[all]'
```

## Upgrading

```bash
# Upgrade to latest version
pip install --upgrade 'markitdown[all]'

# Check version
pip show markitdown
```

## Uninstallation

```bash
pip uninstall markitdown
```

## Next Steps

After installation:
1. Read `QUICK_REFERENCE.md` for basic usage
2. See `SKILL.md` for comprehensive guide
3. Try example scripts in `scripts/` directory
4. Check `assets/example_usage.md` for practical examples

## Skill Scripts Setup

To use the skill scripts:

```bash
# Navigate to scripts directory
cd /Users/vinayak/Documents/claude-scientific-writer/.claude/skills/markitdown/scripts

# Scripts are already executable, just run them
python batch_convert.py --help
python convert_with_ai.py --help
python convert_literature.py --help
```

## Testing Installation

Create a test file to verify everything works:

```python
# test_markitdown.py
from markitdown import MarkItDown

def test_basic():
    md = MarkItDown()
    # Create a simple test file
    with open("test.txt", "w") as f:
        f.write("Hello MarkItDown!")
    
    # Convert it
    result = md.convert("test.txt")
    print("✓ Basic conversion works")
    print(result.text_content)
    
    # Cleanup
    import os
    os.remove("test.txt")

if __name__ == "__main__":
    test_basic()
```

Run it:
```bash
python test_markitdown.py
```

## Getting Help

- **Documentation**: See `SKILL.md` and `README.md`
- **GitHub Issues**: https://github.com/microsoft/markitdown/issues
- **Examples**: `assets/example_usage.md`
- **API Reference**: `references/api_reference.md`

