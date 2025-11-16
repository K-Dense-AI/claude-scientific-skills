---
name: markitdown
description: Convert various file formats (PDF, Office documents, images, audio, web content, structured data) to Markdown optimized for LLM processing. Use when converting documents to markdown, extracting text from PDFs/Office files, transcribing audio, performing OCR on images, extracting YouTube transcripts, or processing batches of files. Supports 20+ formats including DOCX, XLSX, PPTX, PDF, HTML, EPUB, CSV, JSON, images with OCR, and audio with transcription.
---

# MarkItDown

## Overview

MarkItDown is a Python utility that converts various file formats into Markdown format, optimized for use with large language models and text analysis pipelines. It preserves document structure (headings, lists, tables, hyperlinks) while producing clean, token-efficient Markdown output.

## When to Use This Skill

Use this skill when users request:
- Converting documents to Markdown format
- Extracting text from PDF, Word, PowerPoint, or Excel files
- Performing OCR on images to extract text
- Transcribing audio files to text
- Extracting YouTube video transcripts
- Processing HTML, EPUB, or web content to Markdown
- Converting structured data (CSV, JSON, XML) to readable Markdown
- Batch converting multiple files or ZIP archives
- Preparing documents for LLM analysis or RAG systems

## Core Capabilities

### 1. Document Conversion

Convert Office documents and PDFs to Markdown while preserving structure.

**Supported formats:**
- PDF files (with optional Azure Document Intelligence integration)
- Word documents (DOCX)
- PowerPoint presentations (PPTX)
- Excel spreadsheets (XLSX, XLS)

**Basic usage:**
```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("document.pdf")
print(result.text_content)
```

**Command-line:**
```bash
markitdown document.pdf -o output.md
```

See `references/document_conversion.md` for detailed documentation on document-specific features.

### 2. Media Processing

Extract text from images using OCR and transcribe audio files to text.

**Supported formats:**
- Images (JPEG, PNG, GIF, etc.) with EXIF metadata extraction
- Audio files with speech transcription (requires speech_recognition)

**Image with OCR:**
```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("image.jpg")
print(result.text_content)  # Includes EXIF metadata and OCR text
```

**Audio transcription:**
```python
result = md.convert("audio.wav")
print(result.text_content)  # Transcribed speech
```

See `references/media_processing.md` for advanced media handling options.

### 3. Web Content Extraction

Convert web-based content and e-books to Markdown.

**Supported formats:**
- HTML files and web pages
- YouTube video transcripts (via URL)
- EPUB books
- RSS feeds

**YouTube transcript:**
```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("https://youtube.com/watch?v=VIDEO_ID")
print(result.text_content)
```

See `references/web_content.md` for web extraction details.

### 4. Structured Data Handling

Convert structured data formats to readable Markdown tables.

**Supported formats:**
- CSV files
- JSON files
- XML files

**CSV to Markdown table:**
```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("data.csv")
print(result.text_content)  # Formatted as Markdown table
```

See `references/structured_data.md` for format-specific options.

### 5. Advanced Integrations

Enhance conversion quality with AI-powered features.

**Azure Document Intelligence:**
For enhanced PDF processing with better table extraction and layout analysis:
```python
from markitdown import MarkItDown

md = MarkItDown(docintel_endpoint="<endpoint>", docintel_key="<key>")
result = md.convert("complex.pdf")
```

**LLM-Powered Image Descriptions:**
Generate detailed image descriptions using GPT-4o:
```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI()
md = MarkItDown(llm_client=client, llm_model="gpt-4o")
result = md.convert("presentation.pptx")  # Images described with LLM
```

See `references/advanced_integrations.md` for integration details.

### 6. Batch Processing

Process multiple files or entire ZIP archives at once.

**ZIP file processing:**
```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("archive.zip")
print(result.text_content)  # All files converted and concatenated
```

**Batch script:**
Use the provided batch processing script for directory conversion:
```bash
python scripts/batch_convert.py /path/to/documents /path/to/output
```

See `scripts/batch_convert.py` for implementation details.

## Installation

**Full installation (all features):**
```bash
uv pip install 'markitdown[all]'
```

**Modular installation (specific features):**
```bash
uv pip install 'markitdown[pdf]'           # PDF support
uv pip install 'markitdown[docx]'          # Word support
uv pip install 'markitdown[pptx]'          # PowerPoint support
uv pip install 'markitdown[xlsx]'          # Excel support
uv pip install 'markitdown[audio]'         # Audio transcription
uv pip install 'markitdown[youtube]'       # YouTube transcripts
```

**Requirements:**
- Python 3.10 or higher

## Output Format

MarkItDown produces clean, token-efficient Markdown optimized for LLM consumption:
- Preserves headings, lists, and tables
- Maintains hyperlinks and formatting
- Includes metadata where relevant (EXIF, document properties)
- No temporary files created (streaming approach)

## Common Workflows

**Preparing documents for RAG:**
```python
from markitdown import MarkItDown

md = MarkItDown()

# Convert knowledge base documents
docs = ["manual.pdf", "guide.docx", "faq.html"]
markdown_content = []

for doc in docs:
    result = md.convert(doc)
    markdown_content.append(result.text_content)

# Now ready for embedding and indexing
```

**Document analysis pipeline:**
```bash
# Convert all PDFs in directory
for file in documents/*.pdf; do
    markitdown "$file" -o "markdown/$(basename "$file" .pdf).md"
done
```

## Plugin System

MarkItDown supports extensible plugins for custom conversion logic. Plugins are disabled by default for security:

```python
from markitdown import MarkItDown

# Enable plugins if needed
md = MarkItDown(enable_plugins=True)
```

## Resources

This skill includes comprehensive reference documentation for each capability:

- **references/document_conversion.md** - Detailed PDF, DOCX, PPTX, XLSX conversion options
- **references/media_processing.md** - Image OCR and audio transcription details
- **references/web_content.md** - HTML, YouTube, and EPUB extraction
- **references/structured_data.md** - CSV, JSON, XML conversion formats
- **references/advanced_integrations.md** - Azure Document Intelligence and LLM integration
- **scripts/batch_convert.py** - Batch processing utility for directories
