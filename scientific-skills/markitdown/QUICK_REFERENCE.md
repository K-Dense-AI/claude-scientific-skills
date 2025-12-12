# MarkItDown Quick Reference

## Installation

```bash
# All features
pip install 'markitdown[all]'

# Specific formats
pip install 'markitdown[pdf,docx,pptx,xlsx]'
```

## Basic Usage

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("file.pdf")
print(result.text_content)
```

## Command Line

```bash
# Simple conversion
markitdown input.pdf > output.md
markitdown input.pdf -o output.md

# With plugins
markitdown --use-plugins file.pdf -o output.md
```

## Common Tasks

### Convert PDF
```python
md = MarkItDown()
result = md.convert("paper.pdf")
```

### Convert with AI
```python
from openai import OpenAI

# Use OpenRouter for multiple model access
client = OpenAI(
    api_key="your-openrouter-api-key",
    base_url="https://openrouter.ai/api/v1"
)

md = MarkItDown(
    llm_client=client,
    llm_model="anthropic/claude-sonnet-4.5"  # recommended for vision
)
result = md.convert("slides.pptx")
```

### Batch Convert
```bash
python scripts/batch_convert.py input/ output/ --extensions .pdf .docx
```

### Literature Conversion
```bash
python scripts/convert_literature.py papers/ markdown/ --create-index
```

## Supported Formats

| Format | Extension | Notes |
|--------|-----------|-------|
| PDF | `.pdf` | Full text + OCR |
| Word | `.docx` | Tables, formatting |
| PowerPoint | `.pptx` | Slides + notes |
| Excel | `.xlsx`, `.xls` | Tables |
| Images | `.jpg`, `.png`, `.gif`, `.webp` | EXIF + OCR |
| Audio | `.wav`, `.mp3` | Transcription |
| HTML | `.html`, `.htm` | Clean conversion |
| Data | `.csv`, `.json`, `.xml` | Structured |
| Archives | `.zip` | Iterates contents |
| E-books | `.epub` | Full text |
| YouTube | URLs | Transcripts |

## Optional Dependencies

```bash
[all]                  # All features
[pdf]                  # PDF support
[docx]                 # Word documents
[pptx]                 # PowerPoint
[xlsx]                 # Excel
[xls]                  # Old Excel
[outlook]              # Outlook messages
[az-doc-intel]         # Azure Document Intelligence
[audio-transcription]  # Audio files
[youtube-transcription] # YouTube videos
```

## AI-Enhanced Conversion

### Scientific Papers
```python
from openai import OpenAI

# Initialize OpenRouter client
client = OpenAI(
    api_key="your-openrouter-api-key",
    base_url="https://openrouter.ai/api/v1"
)

md = MarkItDown(
    llm_client=client,
    llm_model="anthropic/claude-sonnet-4.5",  # recommended for scientific vision
    llm_prompt="Describe scientific figures with technical precision"
)
result = md.convert("paper.pdf")
```

### Custom Prompts
```python
prompt = """
Analyze this data visualization. Describe:
- Type of chart/graph
- Key trends and patterns
- Notable data points
"""

md = MarkItDown(
    llm_client=client,
    llm_model="anthropic/claude-sonnet-4.5",
    llm_prompt=prompt
)
```

### Available Models via OpenRouter
- `anthropic/claude-sonnet-4.5` - **Recommended for scientific vision**
- `anthropic/claude-opus-4.5` - Advanced vision model
- `openai/gpt-4o` - GPT-4 Omni (vision)
- `openai/gpt-4-vision` - GPT-4 Vision
- `google/gemini-pro-vision` - Gemini Pro Vision

See https://openrouter.ai/models for full list

## Azure Document Intelligence

```python
md = MarkItDown(docintel_endpoint="https://YOUR-ENDPOINT.cognitiveservices.azure.com/")
result = md.convert("complex_layout.pdf")
```

## Batch Processing

### Python
```python
from markitdown import MarkItDown
from pathlib import Path

md = MarkItDown()

for file in Path("input/").glob("*.pdf"):
    result = md.convert(str(file))
    output = Path("output") / f"{file.stem}.md"
    output.write_text(result.text_content)
```

### Script
```bash
# Parallel conversion
python scripts/batch_convert.py input/ output/ --workers 8

# Recursive
python scripts/batch_convert.py input/ output/ -r
```

## Error Handling

```python
try:
    result = md.convert("file.pdf")
except FileNotFoundError:
    print("File not found")
except Exception as e:
    print(f"Error: {e}")
```

## Streaming

```python
with open("large_file.pdf", "rb") as f:
    result = md.convert_stream(f, file_extension=".pdf")
```

## Common Prompts

### Scientific
```
Analyze this scientific figure. Describe:
- Type of visualization
- Key data points and trends
- Axes, labels, and legends
- Scientific significance
```

### Medical
```
Describe this medical image. Include:
- Type of imaging (X-ray, MRI, CT, etc.)
- Anatomical structures visible
- Notable findings
- Clinical relevance
```

### Data Visualization
```
Analyze this data visualization:
- Chart type
- Variables and axes
- Data ranges
- Key patterns and outliers
```

## Performance Tips

1. **Reuse instance**: Create once, use many times
2. **Parallel processing**: Use ThreadPoolExecutor for multiple files
3. **Stream large files**: Use `convert_stream()` for big files
4. **Choose right format**: Install only needed dependencies

## Environment Variables

```bash
# OpenRouter for AI-enhanced conversions
export OPENROUTER_API_KEY="sk-or-v1-..."

# Azure Document Intelligence (optional)
export AZURE_DOCUMENT_INTELLIGENCE_KEY="key..."
export AZURE_DOCUMENT_INTELLIGENCE_ENDPOINT="https://..."
```

## Scripts Quick Reference

### batch_convert.py
```bash
python scripts/batch_convert.py INPUT OUTPUT [OPTIONS]

Options:
  --extensions .pdf .docx    File types to convert
  --recursive, -r            Search subdirectories
  --workers 4                Parallel workers
  --verbose, -v              Detailed output
  --plugins, -p              Enable plugins
```

### convert_with_ai.py
```bash
python scripts/convert_with_ai.py INPUT OUTPUT [OPTIONS]

Options:
  --api-key KEY              OpenRouter API key
  --model MODEL              Model name (default: anthropic/claude-sonnet-4.5)
  --prompt-type TYPE         Preset prompt (scientific, medical, etc.)
  --custom-prompt TEXT       Custom prompt
  --list-prompts             Show available prompts
```

### convert_literature.py
```bash
python scripts/convert_literature.py INPUT OUTPUT [OPTIONS]

Options:
  --organize-by-year, -y     Organize by year
  --create-index, -i         Create index file
  --recursive, -r            Search subdirectories
```

## Troubleshooting

### Missing Dependencies
```bash
pip install 'markitdown[pdf]'  # Install PDF support
```

### Binary File Error
```python
# Wrong
with open("file.pdf", "r") as f:

# Correct
with open("file.pdf", "rb") as f:  # Binary mode
```

### OCR Not Working
```bash
# macOS
brew install tesseract

# Ubuntu
sudo apt-get install tesseract-ocr
```

## More Information

- **Full Documentation**: See `SKILL.md`
- **API Reference**: See `references/api_reference.md`
- **Format Details**: See `references/file_formats.md`
- **Examples**: See `assets/example_usage.md`
- **GitHub**: https://github.com/microsoft/markitdown

