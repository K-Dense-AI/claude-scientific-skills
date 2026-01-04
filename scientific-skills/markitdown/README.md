# MarkItDown Skill

This skill provides comprehensive support for converting various file formats to Markdown using Microsoft's MarkItDown tool.

## Overview

MarkItDown is a Python tool that converts files and office documents to Markdown format. This skill includes:

- Complete API documentation
- Format-specific conversion guides
- Utility scripts for batch processing
- AI-enhanced conversion examples
- Integration with scientific workflows

## Contents

### Main Skill File
- **SKILL.md** - Complete guide to using MarkItDown with quick start, examples, and best practices

### References
- **api_reference.md** - Detailed API documentation, class references, and method signatures
- **file_formats.md** - Format-specific details for all supported file types

### Scripts
- **batch_convert.py** - Batch convert multiple files with parallel processing
- **convert_with_ai.py** - AI-enhanced conversion with custom prompts
- **convert_literature.py** - Scientific literature conversion with metadata extraction

### Assets
- **example_usage.md** - Practical examples for common use cases

## Installation

```bash
# Install with all features
pip install 'markitdown[all]'

# Or install specific features
pip install 'markitdown[pdf,docx,pptx,xlsx]'
```

## Quick Start

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("document.pdf")
print(result.text_content)
```

## Supported Formats

- **Documents**: PDF, DOCX, PPTX, XLSX, EPUB
- **Images**: JPEG, PNG, GIF, WebP (with OCR)
- **Audio**: WAV, MP3 (with transcription)
- **Web**: HTML, YouTube URLs
- **Data**: CSV, JSON, XML
- **Archives**: ZIP files

## Key Features

### 1. AI-Enhanced Conversions
Use AI models via OpenRouter to generate detailed image descriptions:

```python
from openai import OpenAI

# OpenRouter provides access to 100+ AI models
client = OpenAI(
    api_key="your-openrouter-api-key",
    base_url="https://openrouter.ai/api/v1"
)

md = MarkItDown(
    llm_client=client,
    llm_model="anthropic/claude-sonnet-4.5"  # recommended for vision
)
result = md.convert("presentation.pptx")
```

### 2. Batch Processing
Convert multiple files efficiently:

```bash
python scripts/batch_convert.py papers/ output/ --extensions .pdf .docx
```

### 3. Scientific Literature
Convert and organize research papers:

```bash
python scripts/convert_literature.py papers/ output/ --organize-by-year --create-index
```

### 4. Azure Document Intelligence
Enhanced PDF conversion with Microsoft Document Intelligence:

```python
md = MarkItDown(docintel_endpoint="https://YOUR-ENDPOINT.cognitiveservices.azure.com/")
result = md.convert("complex_document.pdf")
```

## Use Cases

### Literature Review
Convert research papers to Markdown for easier analysis and note-taking.

### Data Extraction
Extract tables from Excel files into Markdown format.

### Presentation Processing
Convert PowerPoint slides with AI-generated descriptions.

### Document Analysis
Process documents for LLM consumption with token-efficient Markdown.

### YouTube Transcripts
Fetch and convert YouTube video transcriptions.

## Scripts Usage

### Batch Convert
```bash
# Convert all PDFs in a directory
python scripts/batch_convert.py input_dir/ output_dir/ --extensions .pdf

# Recursive with multiple formats
python scripts/batch_convert.py docs/ markdown/ --extensions .pdf .docx .pptx -r
```

### AI-Enhanced Conversion
```bash
# Convert with AI descriptions via OpenRouter
export OPENROUTER_API_KEY="sk-or-v1-..."
python scripts/convert_with_ai.py paper.pdf output.md --prompt-type scientific

# Use different models
python scripts/convert_with_ai.py image.png output.md --model anthropic/claude-sonnet-4.5

# Use custom prompt
python scripts/convert_with_ai.py image.png output.md --custom-prompt "Describe this diagram"
```

### Literature Conversion
```bash
# Convert papers with metadata extraction
python scripts/convert_literature.py papers/ markdown/ --organize-by-year --create-index
```

## Integration with Scientific Writer

This skill integrates seamlessly with the Scientific Writer CLI for:
- Converting source materials for paper writing
- Processing literature for reviews
- Extracting data from various document formats
- Preparing documents for LLM analysis

## Resources

- **MarkItDown GitHub**: https://github.com/microsoft/markitdown
- **PyPI**: https://pypi.org/project/markitdown/
- **OpenRouter**: https://openrouter.ai (AI model access)
- **OpenRouter API Keys**: https://openrouter.ai/keys
- **OpenRouter Models**: https://openrouter.ai/models
- **License**: MIT

## Requirements

- Python 3.10+
- Optional dependencies based on formats needed
- OpenRouter API key (for AI-enhanced conversions) - Get at https://openrouter.ai/keys
- Azure subscription (optional, for Document Intelligence)

## Examples

See `assets/example_usage.md` for comprehensive examples covering:
- Basic conversions
- Scientific workflows
- AI-enhanced processing
- Batch operations
- Error handling
- Integration patterns

