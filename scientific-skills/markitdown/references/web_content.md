# Web Content Extraction Reference

This document provides detailed information about extracting content from HTML, YouTube, EPUB, and other web-based formats.

## HTML Conversion

Convert HTML files and web pages to clean Markdown format.

### Basic HTML Conversion

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("webpage.html")
print(result.text_content)
```

### HTML Processing Features

**What's preserved:**
- Headings (`<h1>` → `#`, `<h2>` → `##`, etc.)
- Paragraphs and text formatting
- Links (`<a>` → `[text](url)`)
- Lists (ordered and unordered)
- Tables → Markdown tables
- Code blocks and inline code
- Emphasis (bold, italic)

**What's removed:**
- Scripts and styles
- Navigation elements
- Advertising content
- Boilerplate markup
- HTML comments

### HTML from URLs

Convert web pages directly from URLs:

```python
from markitdown import MarkItDown
import requests

md = MarkItDown()

# Fetch and convert web page
response = requests.get("https://example.com/article")
with open("temp.html", "wb") as f:
    f.write(response.content)

result = md.convert("temp.html")
print(result.text_content)
```

### Clean Web Article Extraction

For extracting main content from web articles:

```python
from markitdown import MarkItDown
import requests
from readability import Document  # pip install readability-lxml

md = MarkItDown()

# Fetch page
url = "https://example.com/article"
response = requests.get(url)

# Extract main content
doc = Document(response.content)
html_content = doc.summary()

# Save and convert
with open("article.html", "w") as f:
    f.write(html_content)

result = md.convert("article.html")
print(result.text_content)
```

### HTML with Images

HTML files containing images can be enhanced with LLM descriptions:

```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI()
md = MarkItDown(llm_client=client, llm_model="gpt-4o")
result = md.convert("page_with_images.html")
```

## YouTube Transcripts

Extract video transcripts from YouTube videos.

### Basic YouTube Conversion

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("https://www.youtube.com/watch?v=VIDEO_ID")
print(result.text_content)
```

### YouTube Installation

```bash
pip install 'markitdown[youtube]'
```

This installs the `youtube-transcript-api` dependency.

### YouTube URL Formats

MarkItDown supports various YouTube URL formats:
- `https://www.youtube.com/watch?v=VIDEO_ID`
- `https://youtu.be/VIDEO_ID`
- `https://www.youtube.com/embed/VIDEO_ID`
- `https://m.youtube.com/watch?v=VIDEO_ID`

### YouTube Transcript Features

**What's included:**
- Full video transcript text
- Timestamps (optional, depending on availability)
- Video metadata (title, description)
- Captions in available languages

**Transcript languages:**
```python
from markitdown import MarkItDown

md = MarkItDown()

# Get transcript in specific language (if available)
# Language codes: 'en', 'es', 'fr', 'de', etc.
result = md.convert("https://youtube.com/watch?v=VIDEO_ID")
```

### YouTube Playlist Processing

Process multiple videos from a playlist:

```python
from markitdown import MarkItDown

md = MarkItDown()

video_ids = [
    "VIDEO_ID_1",
    "VIDEO_ID_2",
    "VIDEO_ID_3"
]

transcripts = []
for vid_id in video_ids:
    url = f"https://youtube.com/watch?v={vid_id}"
    result = md.convert(url)
    transcripts.append({
        'video_id': vid_id,
        'transcript': result.text_content
    })
```

### YouTube Use Cases

**Content Analysis:**
- Analyze video content without watching
- Extract key information from tutorials
- Build searchable transcript databases

**Research:**
- Process interview transcripts
- Extract lecture content
- Analyze presentation content

**Accessibility:**
- Generate text versions of video content
- Create searchable video archives

### YouTube Limitations

- Requires videos to have captions/transcripts available
- Auto-generated captions may have transcription errors
- Some videos may disable transcript access
- Rate limiting may apply for bulk processing

## EPUB Books

Convert EPUB e-books to Markdown format.

### Basic EPUB Conversion

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("book.epub")
print(result.text_content)
```

### EPUB Processing Features

**What's extracted:**
- Book text content
- Chapter structure
- Headings and formatting
- Tables of contents
- Footnotes and references

**What's preserved:**
- Heading hierarchy
- Text emphasis (bold, italic)
- Links and references
- Lists and tables

### EPUB with Images

EPUB files often contain images (covers, diagrams, illustrations):

```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI()
md = MarkItDown(llm_client=client, llm_model="gpt-4o")
result = md.convert("illustrated_book.epub")
```

### EPUB Use Cases

**Research:**
- Convert textbooks to searchable format
- Extract content for analysis
- Build digital libraries

**Content Processing:**
- Prepare books for LLM training data
- Convert to different formats
- Create summaries and extracts

**Accessibility:**
- Convert to more accessible formats
- Extract text for screen readers
- Process for text-to-speech

## RSS Feeds

Process RSS feeds to extract article content.

### Basic RSS Processing

```python
from markitdown import MarkItDown
import feedparser

md = MarkItDown()

# Parse RSS feed
feed = feedparser.parse("https://example.com/feed.xml")

# Convert each entry
for entry in feed.entries:
    # Save entry HTML
    with open("temp.html", "w") as f:
        f.write(entry.summary)

    result = md.convert("temp.html")
    print(f"## {entry.title}\n\n{result.text_content}\n\n")
```

## Combined Web Content Workflows

### Web Scraping Pipeline

```python
from markitdown import MarkItDown
import requests
from bs4 import BeautifulSoup

md = MarkItDown()

def scrape_and_convert(url):
    """Scrape webpage and convert to Markdown."""
    response = requests.get(url)
    soup = BeautifulSoup(response.content, 'html.parser')

    # Extract main content
    main_content = soup.find('article') or soup.find('main')

    if main_content:
        # Save HTML
        with open("temp.html", "w") as f:
            f.write(str(main_content))

        # Convert to Markdown
        result = md.convert("temp.html")
        return result.text_content

    return None

# Use it
markdown = scrape_and_convert("https://example.com/article")
print(markdown)
```

### YouTube Learning Content Extraction

```python
from markitdown import MarkItDown

md = MarkItDown()

# Course videos
course_videos = [
    ("https://youtube.com/watch?v=ID1", "Lesson 1: Introduction"),
    ("https://youtube.com/watch?v=ID2", "Lesson 2: Basics"),
    ("https://youtube.com/watch?v=ID3", "Lesson 3: Advanced")
]

course_content = []
for url, title in course_videos:
    result = md.convert(url)
    course_content.append(f"# {title}\n\n{result.text_content}")

# Combine into course document
full_course = "\n\n---\n\n".join(course_content)
with open("course_transcript.md", "w") as f:
    f.write(full_course)
```

### Documentation Scraping

```python
from markitdown import MarkItDown
import requests
from urllib.parse import urljoin, urlparse

md = MarkItDown()

def scrape_documentation(base_url, page_urls):
    """Scrape multiple documentation pages."""
    docs = []

    for page_url in page_urls:
        full_url = urljoin(base_url, page_url)

        # Fetch page
        response = requests.get(full_url)
        with open("temp.html", "wb") as f:
            f.write(response.content)

        # Convert
        result = md.convert("temp.html")
        docs.append({
            'url': full_url,
            'content': result.text_content
        })

    return docs

# Example usage
base = "https://docs.example.com/"
pages = ["intro.html", "getting-started.html", "api.html"]
documentation = scrape_documentation(base, pages)
```

### EPUB Library Processing

```python
from markitdown import MarkItDown
import os

md = MarkItDown()

def process_epub_library(library_path, output_path):
    """Convert all EPUB books in a directory."""
    for filename in os.listdir(library_path):
        if filename.endswith('.epub'):
            epub_path = os.path.join(library_path, filename)

            try:
                result = md.convert(epub_path)

                # Save markdown
                output_file = filename.replace('.epub', '.md')
                output_full = os.path.join(output_path, output_file)

                with open(output_full, 'w') as f:
                    f.write(result.text_content)

                print(f"Converted: {filename}")
            except Exception as e:
                print(f"Failed to convert {filename}: {e}")

# Process library
process_epub_library("books", "markdown_books")
```

## Error Handling

### HTML Conversion Errors

```python
from markitdown import MarkItDown

md = MarkItDown()

try:
    result = md.convert("webpage.html")
    print(result.text_content)
except FileNotFoundError:
    print("HTML file not found")
except Exception as e:
    print(f"Conversion error: {e}")
```

### YouTube Transcript Errors

```python
from markitdown import MarkItDown

md = MarkItDown()

try:
    result = md.convert("https://youtube.com/watch?v=VIDEO_ID")
    print(result.text_content)
except Exception as e:
    print(f"Failed to get transcript: {e}")
    # Common issues: No transcript available, video unavailable, network error
```

### EPUB Conversion Errors

```python
from markitdown import MarkItDown

md = MarkItDown()

try:
    result = md.convert("book.epub")
    print(result.text_content)
except Exception as e:
    print(f"EPUB processing error: {e}")
    # Common issues: Corrupted file, unsupported DRM, invalid format
```

## Best Practices

### HTML Processing
- Clean HTML before conversion for better results
- Use readability libraries to extract main content
- Handle different encodings appropriately
- Remove unnecessary markup

### YouTube Processing
- Check transcript availability before batch processing
- Handle API rate limits gracefully
- Store transcripts to avoid re-fetching
- Respect YouTube's terms of service

### EPUB Processing
- DRM-protected EPUBs cannot be processed
- Large EPUBs may require more memory
- Some formatting may not translate perfectly
- Test with representative samples first

### Web Scraping Ethics
- Respect robots.txt
- Add delays between requests
- Identify your scraper in User-Agent
- Cache results to minimize requests
- Follow website terms of service
