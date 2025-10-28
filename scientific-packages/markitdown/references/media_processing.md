# Media Processing Reference

This document provides detailed information about processing images and audio files with MarkItDown.

## Image Processing

MarkItDown can extract text from images using OCR and retrieve EXIF metadata.

### Basic Image Conversion

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("photo.jpg")
print(result.text_content)
```

### Image Processing Features

**What's extracted:**
1. **EXIF Metadata** - Camera settings, date, location, etc.
2. **OCR Text** - Text detected in the image (requires tesseract)
3. **Image Description** - AI-generated description (with LLM integration)

### EXIF Metadata Extraction

Images from cameras and smartphones contain EXIF metadata that's automatically extracted:

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("IMG_1234.jpg")
print(result.text_content)
```

**Example output includes:**
- Camera make and model
- Capture date and time
- GPS coordinates (if available)
- Exposure settings (ISO, shutter speed, aperture)
- Image dimensions
- Orientation

### OCR (Optical Character Recognition)

Extract text from images containing text (screenshots, scanned documents, photos of text):

**Requirements:**
- Install tesseract OCR engine:
  ```bash
  # macOS
  brew install tesseract

  # Ubuntu/Debian
  apt-get install tesseract-ocr

  # Windows
  # Download installer from https://github.com/UB-Mannheim/tesseract/wiki
  ```

**Usage:**
```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("screenshot.png")
print(result.text_content)  # Contains OCR'd text
```

**Best practices for OCR:**
- Use high-resolution images for better accuracy
- Ensure good contrast between text and background
- Straighten skewed text if possible
- Use well-lit, clear images

### LLM-Generated Image Descriptions

Generate detailed, contextual descriptions of images using GPT-4o or other vision models:

```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI()
md = MarkItDown(llm_client=client, llm_model="gpt-4o")
result = md.convert("diagram.png")
print(result.text_content)
```

**Custom prompts for specific needs:**

```python
# For diagrams
md = MarkItDown(
    llm_client=client,
    llm_model="gpt-4o",
    llm_prompt="Describe this diagram in detail, explaining all components and their relationships"
)

# For charts
md = MarkItDown(
    llm_client=client,
    llm_model="gpt-4o",
    llm_prompt="Analyze this chart and provide key data points and trends"
)

# For UI screenshots
md = MarkItDown(
    llm_client=client,
    llm_model="gpt-4o",
    llm_prompt="Describe this user interface, listing all visible elements and their layout"
)
```

### Supported Image Formats

MarkItDown supports all common image formats:
- JPEG/JPG
- PNG
- GIF
- BMP
- TIFF
- WebP
- HEIC (requires additional libraries on some platforms)

## Audio Processing

MarkItDown can transcribe audio files to text using speech recognition.

### Basic Audio Conversion

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("recording.wav")
print(result.text_content)  # Transcribed speech
```

### Audio Transcription Setup

**Installation:**
```bash
pip install 'markitdown[audio]'
```

This installs the `speech_recognition` library and dependencies.

### Supported Audio Formats

- WAV
- AIFF
- FLAC
- MP3 (requires ffmpeg or libav)
- OGG (requires ffmpeg or libav)
- Other formats supported by speech_recognition

### Audio Transcription Engines

MarkItDown uses the `speech_recognition` library, which supports multiple backends:

**Default (Google Speech Recognition):**
```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("audio.wav")
```

**Note:** Default Google Speech Recognition requires internet connection.

### Audio Quality Considerations

For best transcription accuracy:
- Use clear audio with minimal background noise
- Prefer WAV or FLAC for better quality
- Ensure speech is clear and at good volume
- Avoid multiple overlapping speakers
- Use mono audio when possible

### Audio Preprocessing Tips

For better results, consider preprocessing audio:

```python
# Example: If you have pydub installed
from pydub import AudioSegment
from pydub.effects import normalize

# Load and normalize audio
audio = AudioSegment.from_file("recording.mp3")
audio = normalize(audio)
audio.export("normalized.wav", format="wav")

# Then convert with MarkItDown
from markitdown import MarkItDown
md = MarkItDown()
result = md.convert("normalized.wav")
```

## Combined Media Workflows

### Processing Multiple Images in Batch

```python
from markitdown import MarkItDown
from openai import OpenAI
import os

client = OpenAI()
md = MarkItDown(llm_client=client, llm_model="gpt-4o")

# Process all images in directory
for filename in os.listdir("images"):
    if filename.lower().endswith(('.png', '.jpg', '.jpeg')):
        result = md.convert(f"images/{filename}")

        # Save markdown with same name
        output = filename.rsplit('.', 1)[0] + '.md'
        with open(f"output/{output}", "w") as f:
            f.write(result.text_content)
```

### Screenshot Analysis Pipeline

```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI()
md = MarkItDown(
    llm_client=client,
    llm_model="gpt-4o",
    llm_prompt="Describe this screenshot comprehensively, including UI elements, text, and layout"
)

screenshots = ["screen1.png", "screen2.png", "screen3.png"]
analysis = []

for screenshot in screenshots:
    result = md.convert(screenshot)
    analysis.append({
        'file': screenshot,
        'content': result.text_content
    })

# Now ready for further processing
```

### Document Images with OCR

For scanned documents or photos of documents:

```python
from markitdown import MarkItDown

md = MarkItDown()

# Process scanned pages
pages = ["page1.jpg", "page2.jpg", "page3.jpg"]
full_text = []

for page in pages:
    result = md.convert(page)
    full_text.append(result.text_content)

# Combine into single document
document = "\n\n---\n\n".join(full_text)
print(document)
```

### Presentation Slide Images

When you have presentation slides as images:

```python
from markitdown import MarkItDown
from openai import OpenAI

client = OpenAI()
md = MarkItDown(
    llm_client=client,
    llm_model="gpt-4o",
    llm_prompt="Describe this presentation slide, including title, bullet points, and visual elements"
)

# Process slide images
for i in range(1, 21):  # 20 slides
    result = md.convert(f"slides/slide_{i}.png")
    print(f"## Slide {i}\n\n{result.text_content}\n\n")
```

## Error Handling

### Image Processing Errors

```python
from markitdown import MarkItDown

md = MarkItDown()

try:
    result = md.convert("image.jpg")
    print(result.text_content)
except FileNotFoundError:
    print("Image file not found")
except Exception as e:
    print(f"Error processing image: {e}")
```

### Audio Processing Errors

```python
from markitdown import MarkItDown

md = MarkItDown()

try:
    result = md.convert("audio.mp3")
    print(result.text_content)
except Exception as e:
    print(f"Transcription failed: {e}")
    # Common issues: format not supported, no speech detected, network error
```

## Performance Optimization

### Image Processing

- **LLM descriptions**: Slower but more informative
- **OCR only**: Faster for text extraction
- **EXIF only**: Fastest, metadata only
- **Batch processing**: Process multiple images in parallel

### Audio Processing

- **File size**: Larger files take longer
- **Audio length**: Transcription time scales with duration
- **Format conversion**: WAV/FLAC are faster than MP3/OGG
- **Network dependency**: Default transcription requires internet

## Use Cases

### Document Digitization
Convert scanned documents or photos of documents to searchable text.

### Meeting Notes
Transcribe audio recordings of meetings to text for analysis.

### Presentation Analysis
Extract content from presentation slide images.

### Screenshot Documentation
Generate descriptions of UI screenshots for documentation.

### Image Archiving
Extract metadata and content from photo collections.

### Accessibility
Generate alt-text descriptions for images using LLM integration.

### Data Extraction
OCR text from images containing tables, forms, or structured data.
