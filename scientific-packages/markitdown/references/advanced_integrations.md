# Advanced Integrations Reference

This document provides detailed information about advanced MarkItDown features including Azure Document Intelligence integration, LLM-powered descriptions, and plugin system.

## Azure Document Intelligence Integration

Azure Document Intelligence (formerly Form Recognizer) provides superior PDF processing with advanced table extraction and layout analysis.

### Setup

**Prerequisites:**
1. Azure subscription
2. Document Intelligence resource created in Azure
3. Endpoint URL and API key

**Create Azure Resource:**
```bash
# Using Azure CLI
az cognitiveservices account create \
  --name my-doc-intelligence \
  --resource-group my-resource-group \
  --kind FormRecognizer \
  --sku F0 \
  --location eastus
```

### Basic Usage

```python
from markitdown import MarkItDown

md = MarkItDown(
    docintel_endpoint="https://YOUR-RESOURCE.cognitiveservices.azure.com/",
    docintel_key="YOUR-API-KEY"
)

result = md.convert("complex_document.pdf")
print(result.text_content)
```

### Configuration from Environment Variables

```python
import os
from markitdown import MarkItDown

# Set environment variables
os.environ['AZURE_DOCUMENT_INTELLIGENCE_ENDPOINT'] = 'YOUR-ENDPOINT'
os.environ['AZURE_DOCUMENT_INTELLIGENCE_KEY'] = 'YOUR-KEY'

# Use without explicit credentials
md = MarkItDown(
    docintel_endpoint=os.getenv('AZURE_DOCUMENT_INTELLIGENCE_ENDPOINT'),
    docintel_key=os.getenv('AZURE_DOCUMENT_INTELLIGENCE_KEY')
)

result = md.convert("document.pdf")
```

### When to Use Azure Document Intelligence

**Use for:**
- Complex PDFs with sophisticated tables
- Multi-column layouts
- Forms and structured documents
- Scanned documents requiring OCR
- PDFs with mixed content types
- Documents with intricate formatting

**Benefits over standard extraction:**
- **Superior table extraction** - Better handling of merged cells, complex layouts
- **Layout analysis** - Understands document structure (headers, footers, columns)
- **Form fields** - Extracts key-value pairs from forms
- **Reading order** - Maintains correct text flow in complex layouts
- **OCR quality** - High-quality text extraction from scanned documents

### Comparison Example

**Standard extraction:**
```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("complex_table.pdf")
# May struggle with complex tables
```

**Azure Document Intelligence:**
```python
from markitdown import MarkItDown

md = MarkItDown(
    docintel_endpoint="YOUR-ENDPOINT",
    docintel_key="YOUR-KEY"
)
result = md.convert("complex_table.pdf")
# Better table reconstruction and layout understanding
```

### Cost Considerations

Azure Document Intelligence is a paid service:
- **Free tier**: 500 pages per month
- **Paid tiers**: Pay per page processed
- Monitor usage to control costs
- Use standard extraction for simple documents

### Error Handling

```python
from markitdown import MarkItDown

md = MarkItDown(
    docintel_endpoint="YOUR-ENDPOINT",
    docintel_key="YOUR-KEY"
)

try:
    result = md.convert("document.pdf")
    print(result.text_content)
except Exception as e:
    print(f"Document Intelligence error: {e}")
    # Common issues: authentication, quota exceeded, unsupported file
```

## LLM-Powered Image Descriptions

Generate detailed, contextual descriptions for images using large language models.

### Setup with OpenAI

```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI(api_key="YOUR-OPENAI-API-KEY")
md = MarkItDown(llm_client=client, llm_model="gpt-4o")

result = md.convert("image.jpg")
print(result.text_content)
```

### Supported Use Cases

**Images in documents:**
```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI()
md = MarkItDown(llm_client=client, llm_model="gpt-4o")

# PowerPoint with images
result = md.convert("presentation.pptx")

# Word documents with images
result = md.convert("report.docx")

# Standalone images
result = md.convert("diagram.png")
```

### Custom Prompts

Customize the LLM prompt for specific needs:

```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI()

# For diagrams
md = MarkItDown(
    llm_client=client,
    llm_model="gpt-4o",
    llm_prompt="Analyze this diagram and explain all components, connections, and relationships in detail"
)

# For charts
md = MarkItDown(
    llm_client=client,
    llm_model="gpt-4o",
    llm_prompt="Describe this chart, including the type, axes, data points, trends, and key insights"
)

# For UI screenshots
md = MarkItDown(
    llm_client=client,
    llm_model="gpt-4o",
    llm_prompt="Describe this user interface screenshot, listing all UI elements, their layout, and functionality"
)

# For scientific figures
md = MarkItDown(
    llm_client=client,
    llm_model="gpt-4o",
    llm_prompt="Describe this scientific figure in detail, including methodology, results shown, and significance"
)
```

### Model Selection

**GPT-4o (Recommended):**
- Best vision capabilities
- High-quality descriptions
- Good at understanding context
- Higher cost per image

**GPT-4o-mini:**
- Lower cost alternative
- Good for simpler images
- Faster processing
- May miss subtle details

```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI()

# High quality (more expensive)
md_quality = MarkItDown(llm_client=client, llm_model="gpt-4o")

# Budget option (less expensive)
md_budget = MarkItDown(llm_client=client, llm_model="gpt-4o-mini")
```

### Configuration from Environment

```python
import os
from markitdown import MarkItDown
from openai import OpenAI

# Set API key in environment
os.environ['OPENAI_API_KEY'] = 'YOUR-API-KEY'

client = OpenAI()  # Uses env variable
md = MarkItDown(llm_client=client, llm_model="gpt-4o")
```

### Alternative LLM Providers

**Anthropic Claude:**
```python
from markitdown import MarkItDown
from anthropic import Anthropic

# Note: Check current compatibility with MarkItDown
client = Anthropic(api_key="YOUR-API-KEY")
# May require adapter for MarkItDown compatibility
```

**Azure OpenAI:**
```python
from markitdown import MarkItDown
from openai import AzureOpenAI

client = AzureOpenAI(
    api_key="YOUR-AZURE-KEY",
    api_version="2024-02-01",
    azure_endpoint="https://YOUR-RESOURCE.openai.azure.com"
)

md = MarkItDown(llm_client=client, llm_model="gpt-4o")
```

### Cost Management

**Strategies to reduce LLM costs:**

1. **Selective processing:**
```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI()

# Only use LLM for important documents
if is_important_document(file):
    md = MarkItDown(llm_client=client, llm_model="gpt-4o")
else:
    md = MarkItDown()  # Standard processing

result = md.convert(file)
```

2. **Image filtering:**
```python
# Pre-process to identify images that need descriptions
# Only use LLM for complex/important images
```

3. **Batch processing:**
```python
# Process multiple images in batches
# Monitor costs and set limits
```

4. **Model selection:**
```python
# Use gpt-4o-mini for simple images
# Reserve gpt-4o for complex visualizations
```

### Performance Considerations

**LLM processing adds latency:**
- Each image requires an API call
- Processing time: 1-5 seconds per image
- Network dependent
- Consider parallel processing for multiple images

**Batch optimization:**
```python
from markitdown import MarkItDown
from openai import OpenAI
import concurrent.futures

client = OpenAI()
md = MarkItDown(llm_client=client, llm_model="gpt-4o")

def process_image(image_path):
    return md.convert(image_path)

# Process multiple images in parallel
images = ["img1.jpg", "img2.jpg", "img3.jpg"]
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
    results = list(executor.map(process_image, images))
```

## Combined Advanced Features

### Azure Document Intelligence + LLM Descriptions

Combine both for maximum quality:

```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI()
md = MarkItDown(
    llm_client=client,
    llm_model="gpt-4o",
    docintel_endpoint="YOUR-AZURE-ENDPOINT",
    docintel_key="YOUR-AZURE-KEY"
)

# Best possible PDF conversion with image descriptions
result = md.convert("complex_report.pdf")
```

**Use cases:**
- Research papers with figures
- Business reports with charts
- Technical documentation with diagrams
- Presentations with visual data

### Smart Document Processing Pipeline

```python
from markitdown import MarkItDown
from openai import OpenAI
import os

def smart_convert(file_path):
    """Intelligently choose processing method based on file type."""
    client = OpenAI()
    ext = os.path.splitext(file_path)[1].lower()

    # PDFs with complex tables: Use Azure
    if ext == '.pdf':
        md = MarkItDown(
            docintel_endpoint=os.getenv('AZURE_ENDPOINT'),
            docintel_key=os.getenv('AZURE_KEY')
        )

    # Documents/presentations with images: Use LLM
    elif ext in ['.pptx', '.docx']:
        md = MarkItDown(
            llm_client=client,
            llm_model="gpt-4o"
        )

    # Simple formats: Standard processing
    else:
        md = MarkItDown()

    return md.convert(file_path)

# Use it
result = smart_convert("document.pdf")
```

## Plugin System

MarkItDown supports custom plugins for extending functionality.

### Plugin Architecture

Plugins are disabled by default for security:

```python
from markitdown import MarkItDown

# Enable plugins
md = MarkItDown(enable_plugins=True)
```

### Creating Custom Plugins

**Plugin structure:**
```python
class CustomConverter:
    """Custom converter plugin for MarkItDown."""

    def can_convert(self, file_path):
        """Check if this plugin can handle the file."""
        return file_path.endswith('.custom')

    def convert(self, file_path):
        """Convert file to Markdown."""
        # Your conversion logic here
        return {
            'text_content': '# Converted Content\n\n...'
        }
```

### Plugin Registration

```python
from markitdown import MarkItDown

md = MarkItDown(enable_plugins=True)

# Register custom plugin
md.register_plugin(CustomConverter())

# Use normally
result = md.convert("file.custom")
```

### Plugin Use Cases

**Custom formats:**
- Proprietary document formats
- Specialized scientific data formats
- Legacy file formats

**Enhanced processing:**
- Custom OCR engines
- Specialized table extraction
- Domain-specific parsing

**Integration:**
- Enterprise document systems
- Custom databases
- Specialized APIs

### Plugin Security

**Important security considerations:**
- Plugins run with full system access
- Only enable for trusted plugins
- Validate plugin code before use
- Disable plugins in production unless required

## Error Handling for Advanced Features

```python
from markitdown import MarkItDown
from openai import OpenAI

def robust_convert(file_path):
    """Convert with fallback strategies."""
    try:
        # Try with all advanced features
        client = OpenAI()
        md = MarkItDown(
            llm_client=client,
            llm_model="gpt-4o",
            docintel_endpoint=os.getenv('AZURE_ENDPOINT'),
            docintel_key=os.getenv('AZURE_KEY')
        )
        return md.convert(file_path)

    except Exception as azure_error:
        print(f"Azure failed: {azure_error}")

        try:
            # Fallback: LLM only
            client = OpenAI()
            md = MarkItDown(llm_client=client, llm_model="gpt-4o")
            return md.convert(file_path)

        except Exception as llm_error:
            print(f"LLM failed: {llm_error}")

            # Final fallback: Standard processing
            md = MarkItDown()
            return md.convert(file_path)

# Use it
result = robust_convert("document.pdf")
```

## Best Practices

### Azure Document Intelligence
- Use for complex PDFs only (cost optimization)
- Monitor usage and costs
- Store credentials securely
- Handle quota limits gracefully
- Fall back to standard processing if needed

### LLM Integration
- Use appropriate models for task complexity
- Customize prompts for specific use cases
- Monitor API costs
- Implement rate limiting
- Cache results when possible
- Handle API errors gracefully

### Combined Features
- Test cost/quality tradeoffs
- Use selectively for important documents
- Implement intelligent routing
- Monitor performance and costs
- Have fallback strategies

### Security
- Store API keys securely (environment variables, secrets manager)
- Never commit credentials to code
- Disable plugins unless required
- Validate all inputs
- Use least privilege access
