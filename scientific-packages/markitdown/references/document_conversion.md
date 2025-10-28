# Document Conversion Reference

This document provides detailed information about converting Office documents and PDFs to Markdown using MarkItDown.

## PDF Files

PDF conversion extracts text, tables, and structure from PDF documents.

### Basic PDF Conversion

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("document.pdf")
print(result.text_content)
```

### PDF with Azure Document Intelligence

For complex PDFs with tables, forms, and sophisticated layouts, use Azure Document Intelligence for enhanced extraction:

```python
from markitdown import MarkItDown

md = MarkItDown(
    docintel_endpoint="https://YOUR-ENDPOINT.cognitiveservices.azure.com/",
    docintel_key="YOUR-API-KEY"
)
result = md.convert("complex_table.pdf")
print(result.text_content)
```

**Benefits of Azure Document Intelligence:**
- Superior table extraction and reconstruction
- Better handling of multi-column layouts
- Form field recognition
- Improved text ordering in complex documents

### PDF Handling Notes

- Scanned PDFs require OCR (automatically handled if tesseract is installed)
- Password-protected PDFs are not supported
- Large PDFs may take longer to process
- Vector graphics and embedded images are extracted where possible

## Word Documents (DOCX)

Word document conversion preserves headings, paragraphs, lists, tables, and hyperlinks.

### Basic DOCX Conversion

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("document.docx")
print(result.text_content)
```

### DOCX Structure Preservation

MarkItDown preserves:
- **Headings** → Markdown headers (`#`, `##`, etc.)
- **Bold/Italic** → Markdown emphasis (`**bold**`, `*italic*`)
- **Lists** → Markdown lists (ordered and unordered)
- **Tables** → Markdown tables
- **Hyperlinks** → Markdown links `[text](url)`
- **Images** → Referenced with descriptions (can use LLM for descriptions)

### Command-Line Usage

```bash
# Basic conversion
markitdown report.docx -o report.md

# With output directory
markitdown report.docx -o output/report.md
```

### DOCX with Images

To generate descriptions for images in Word documents, use LLM integration:

```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI()
md = MarkItDown(llm_client=client, llm_model="gpt-4o")
result = md.convert("document_with_images.docx")
```

## PowerPoint Presentations (PPTX)

PowerPoint conversion extracts text from slides while preserving structure.

### Basic PPTX Conversion

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("presentation.pptx")
print(result.text_content)
```

### PPTX Structure

MarkItDown processes presentations as:
- Each slide becomes a major section
- Slide titles become headers
- Bullet points are preserved
- Tables are converted to Markdown tables
- Notes are included if present

### PPTX with Image Descriptions

Presentations often contain important visual information. Use LLM integration to describe images:

```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI()
md = MarkItDown(
    llm_client=client,
    llm_model="gpt-4o",
    llm_prompt="Describe this slide image in detail, focusing on key information"
)
result = md.convert("presentation.pptx")
```

**Custom prompts for presentations:**
- "Describe charts and graphs with their key data points"
- "Explain diagrams and their relationships"
- "Summarize visual content for accessibility"

## Excel Spreadsheets (XLSX, XLS)

Excel conversion formats spreadsheet data as Markdown tables.

### Basic XLSX Conversion

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("data.xlsx")
print(result.text_content)
```

### Multi-Sheet Workbooks

For workbooks with multiple sheets:
- Each sheet becomes a separate section
- Sheet names are used as headers
- Empty sheets are skipped
- Formulas are evaluated (values shown, not formulas)

### XLSX Conversion Details

**What's preserved:**
- Cell values (text, numbers, dates)
- Table structure (rows and columns)
- Sheet names
- Cell formatting (bold headers)

**What's not preserved:**
- Formulas (only computed values)
- Charts and graphs (use LLM integration for descriptions)
- Cell colors and conditional formatting
- Comments and notes

### Large Spreadsheets

For large spreadsheets, consider:
- Processing may be slower for files with many rows/columns
- Very wide tables may not format well in Markdown
- Consider filtering or preprocessing data if possible

### XLS (Legacy Excel) Files

Legacy `.xls` files are supported but require additional dependencies:

```bash
pip install 'markitdown[xls]'
```

Then use normally:
```python
md = MarkItDown()
result = md.convert("legacy_data.xls")
```

## Common Document Conversion Patterns

### Batch Document Processing

```python
from markitdown import MarkItDown
import os

md = MarkItDown()

# Process all documents in a directory
for filename in os.listdir("documents"):
    if filename.endswith(('.pdf', '.docx', '.pptx', '.xlsx')):
        result = md.convert(f"documents/{filename}")

        # Save to output directory
        output_name = os.path.splitext(filename)[0] + ".md"
        with open(f"markdown/{output_name}", "w") as f:
            f.write(result.text_content)
```

### Document with Mixed Content

For documents containing multiple types of content (text, tables, images):

```python
from markitdown import MarkItDown
from openai import OpenAI

# Use LLM for image descriptions + Azure for complex tables
client = OpenAI()
md = MarkItDown(
    llm_client=client,
    llm_model="gpt-4o",
    docintel_endpoint="YOUR-ENDPOINT",
    docintel_key="YOUR-KEY"
)

result = md.convert("complex_report.pdf")
```

### Error Handling

```python
from markitdown import MarkItDown

md = MarkItDown()

try:
    result = md.convert("document.pdf")
    print(result.text_content)
except Exception as e:
    print(f"Conversion failed: {e}")
    # Handle specific errors (file not found, unsupported format, etc.)
```

## Output Quality Tips

**For best results:**
1. Use Azure Document Intelligence for PDFs with complex tables
2. Enable LLM descriptions for documents with important visual content
3. Ensure source documents are well-structured (proper headings, etc.)
4. For scanned documents, ensure good scan quality for OCR accuracy
5. Test with sample documents to verify output quality

## Performance Considerations

**Conversion speed depends on:**
- Document size and complexity
- Number of images (especially with LLM descriptions)
- Use of Azure Document Intelligence
- Available system resources

**Optimization tips:**
- Disable LLM integration if image descriptions aren't needed
- Use standard extraction (not Azure) for simple documents
- Process large batches in parallel when possible
- Consider streaming for very large documents
