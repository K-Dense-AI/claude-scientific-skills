# Platypus Guide - High-Level Page Layout

Platypus ("Page Layout and Typography Using Scripts") provides high-level document layout for complex, flowing documents with minimal code.

## Architecture Overview

Platypus uses a layered design:

1. **DocTemplates** - Document container with page formatting rules
2. **PageTemplates** - Specifications for different page layouts
3. **Frames** - Regions where content flows
4. **Flowables** - Content elements (paragraphs, tables, images, spacers)
5. **Canvas** - Underlying rendering engine (usually hidden)

## Quick Start

```python
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import inch

# Create document
doc = SimpleDocTemplate("output.pdf", pagesize=letter,
                       rightMargin=72, leftMargin=72,
                       topMargin=72, bottomMargin=18)

# Create story (list of flowables)
story = []
styles = getSampleStyleSheet()

# Add content
story.append(Paragraph("Title", styles['Title']))
story.append(Spacer(1, 0.2*inch))
story.append(Paragraph("Body text here", styles['BodyText']))
story.append(PageBreak())

# Build PDF
doc.build(story)
```

## Core Components

### DocTemplates

#### SimpleDocTemplate
Most common template for standard documents:

```python
doc = SimpleDocTemplate(
    filename,
    pagesize=letter,
    rightMargin=72,      # 1 inch = 72 points
    leftMargin=72,
    topMargin=72,
    bottomMargin=18,
    title=None,          # PDF metadata
    author=None,
    subject=None
)
```

#### BaseDocTemplate (Advanced)
For complex documents with multiple page layouts:

```python
from reportlab.platypus import BaseDocTemplate, PageTemplate, Frame
from reportlab.lib.pagesizes import letter

doc = BaseDocTemplate("output.pdf", pagesize=letter)

# Define frames (content regions)
frame1 = Frame(doc.leftMargin, doc.bottomMargin,
              doc.width/2-6, doc.height, id='col1')
frame2 = Frame(doc.leftMargin+doc.width/2+6, doc.bottomMargin,
              doc.width/2-6, doc.height, id='col2')

# Create page template
template = PageTemplate(id='TwoCol', frames=[frame1, frame2])
doc.addPageTemplates([template])

# Build with story
doc.build(story)
```

### Frames

Frames define regions where content flows:

```python
from reportlab.platypus import Frame

frame = Frame(
    x1, y1,              # Lower-left corner
    width, height,       # Dimensions
    leftPadding=6,       # Internal padding
    bottomPadding=6,
    rightPadding=6,
    topPadding=6,
    id=None,             # Optional identifier
    showBoundary=0       # 1 to show frame border (debugging)
)
```

### PageTemplates

Define page layouts with frames and optional functions:

```python
def header_footer(canvas, doc):
    """Called on each page for headers/footers"""
    canvas.saveState()
    canvas.setFont('Helvetica', 9)
    canvas.drawString(inch, 0.75*inch, f"Page {doc.page}")
    canvas.restoreState()

template = PageTemplate(
    id='Normal',
    frames=[frame],
    onPage=header_footer,     # Function called for each page
    onPageEnd=None,
    pagesize=letter
)
```

## Flowables

Flowables are content elements that flow through frames.

### Common Flowables

```python
from reportlab.platypus import (
    Paragraph, Spacer, PageBreak, FrameBreak,
    Image, Table, KeepTogether, CondPageBreak
)

# Spacer - vertical whitespace
Spacer(width, height)

# Page break - force new page
PageBreak()

# Frame break - move to next frame
FrameBreak()

# Conditional page break - break if less than N space remaining
CondPageBreak(height)

# Keep together - prevent splitting across pages
KeepTogether([flowable1, flowable2, ...])
```

### Paragraph Flowable
See `text_and_fonts.md` for detailed Paragraph usage.

```python
from reportlab.platypus import Paragraph
from reportlab.lib.styles import ParagraphStyle

style = ParagraphStyle(
    'CustomStyle',
    fontSize=12,
    leading=14,
    alignment=0  # 0=left, 1=center, 2=right, 4=justify
)

para = Paragraph("Text with <b>bold</b> and <i>italic</i>", style)
story.append(para)
```

### Image Flowable

```python
from reportlab.platypus import Image

# Auto-size to fit
img = Image('photo.jpg')

# Fixed size
img = Image('photo.jpg', width=2*inch, height=2*inch)

# Maintain aspect ratio with max width
img = Image('photo.jpg', width=4*inch, height=3*inch,
           kind='proportional')

story.append(img)
```

### Table Flowable
See `tables_reference.md` for detailed Table usage.

```python
from reportlab.platypus import Table

data = [['Header1', 'Header2'],
        ['Row1Col1', 'Row1Col2'],
        ['Row2Col1', 'Row2Col2']]

table = Table(data, colWidths=[2*inch, 2*inch])
story.append(table)
```

## Page Layouts

### Single Column Document

```python
doc = SimpleDocTemplate("output.pdf", pagesize=letter)
story = []
# Add flowables...
doc.build(story)
```

### Two-Column Layout

```python
from reportlab.platypus import BaseDocTemplate, PageTemplate, Frame

doc = BaseDocTemplate("output.pdf", pagesize=letter)
width, height = letter
margin = inch

# Two side-by-side frames
frame1 = Frame(margin, margin, width/2 - 1.5*margin, height - 2*margin, id='col1')
frame2 = Frame(width/2 + 0.5*margin, margin, width/2 - 1.5*margin, height - 2*margin, id='col2')

template = PageTemplate(id='TwoCol', frames=[frame1, frame2])
doc.addPageTemplates([template])

story = []
# Content flows left column first, then right column
# Add flowables...
doc.build(story)
```

### Multiple Page Templates

```python
from reportlab.platypus import NextPageTemplate

# Define templates
cover_template = PageTemplate(id='Cover', frames=[cover_frame])
body_template = PageTemplate(id='Body', frames=[body_frame])

doc.addPageTemplates([cover_template, body_template])

story = []
# Cover page content
story.append(Paragraph("Cover", title_style))
story.append(NextPageTemplate('Body'))  # Switch to body template
story.append(PageBreak())

# Body content
story.append(Paragraph("Chapter 1", heading_style))
# ...

doc.build(story)
```

## Headers and Footers

Headers and footers are added via `onPage` callback functions:

```python
def header_footer(canvas, doc):
    """Draw header and footer on each page"""
    canvas.saveState()

    # Header
    canvas.setFont('Helvetica-Bold', 12)
    canvas.drawCentredString(letter[0]/2, letter[1] - 0.5*inch,
                            "Document Title")

    # Footer
    canvas.setFont('Helvetica', 9)
    canvas.drawString(inch, 0.75*inch, "Left Footer")
    canvas.drawRightString(letter[0] - inch, 0.75*inch,
                          f"Page {doc.page}")

    canvas.restoreState()

# Apply to template
template = PageTemplate(id='Normal', frames=[frame], onPage=header_footer)
```

## Table of Contents

```python
from reportlab.platypus import TableOfContents
from reportlab.lib.styles import ParagraphStyle

# Create TOC
toc = TableOfContents()
toc.levelStyles = [
    ParagraphStyle(name='TOC1', fontSize=14, leftIndent=0),
    ParagraphStyle(name='TOC2', fontSize=12, leftIndent=20),
]

story = []
story.append(toc)
story.append(PageBreak())

# Add entries
story.append(Paragraph("Chapter 1<a name='ch1'/>", heading_style))
toc.addEntry(0, "Chapter 1", doc.page, 'ch1')

# Must call build twice for TOC to populate
doc.build(story)
```

## Document Properties

```python
from reportlab.lib.pagesizes import letter, A4
from reportlab.lib.units import inch, cm, mm

# Page sizes
letter  # US Letter (8.5" x 11")
A4      # ISO A4 (210mm x 297mm)
landscape(letter)  # Rotate to landscape

# Units
inch    # 72 points
cm      # 28.35 points
mm      # 2.835 points

# Custom page size
custom_size = (6*inch, 9*inch)
```

## Best Practices

1. **Use SimpleDocTemplate** for most documents - it handles common layouts
2. **Build story list** completely before calling `doc.build(story)`
3. **Use Spacer** for vertical spacing instead of empty Paragraphs
4. **Group related content** with KeepTogether to prevent awkward page breaks
5. **Test page breaks** early with realistic content amounts
6. **Use styles consistently** - create style once, reuse throughout document
7. **Set showBoundary=1** on Frames during development to visualize layout
8. **Headers/footers go in onPage** callback, not in story
9. **For long documents**, use BaseDocTemplate with multiple page templates
10. **Build TOC documents twice** to properly populate table of contents
