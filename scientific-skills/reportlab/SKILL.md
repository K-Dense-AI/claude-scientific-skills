---
name: reportlab
description: "PDF generation toolkit. Create invoices, reports, certificates, forms, charts, tables, barcodes, QR codes, Canvas/Platypus APIs, for professional document automation."
---

# ReportLab PDF Generation

## Overview

ReportLab is a powerful Python library for programmatic PDF generation. Create anything from simple documents to complex reports with tables, charts, images, and interactive forms.

**Two main approaches:**
- **Canvas API** (low-level): Direct drawing with coordinate-based positioning - use for precise layouts
- **Platypus** (high-level): Flowing document layout with automatic page breaks - use for multi-page documents

**Core capabilities:**
- Text with rich formatting and custom fonts
- Tables with complex styling and cell spanning
- Charts (bar, line, pie, area, scatter)
- Barcodes and QR codes (Code128, EAN, QR, etc.)
- Images with transparency
- PDF features (links, bookmarks, forms, encryption)

## Choosing the Right Approach

### Use Canvas API when:
- Creating labels, business cards, certificates
- Precise positioning is critical (x, y coordinates)
- Single-page documents or simple layouts
- Drawing graphics, shapes, and custom designs
- Adding barcodes or QR codes at specific locations

### Use Platypus when:
- Creating multi-page documents (reports, articles, books)
- Content should flow automatically across pages
- Need headers/footers that repeat on each page
- Working with paragraphs that can split across pages
- Building complex documents with table of contents

### Use Both when:
- Complex reports need both flowing content AND precise positioning
- Adding headers/footers to Platypus documents (use `onPage` callback with Canvas)
- Embedding custom graphics (Canvas) within flowing documents (Platypus)

## Quick Start Examples

### Simple Canvas Document

```python
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch

c = canvas.Canvas("output.pdf", pagesize=letter)
width, height = letter

# Draw text
c.setFont("Helvetica-Bold", 24)
c.drawString(inch, height - inch, "Hello ReportLab!")

# Draw a rectangle
c.setFillColorRGB(0.2, 0.4, 0.8)
c.rect(inch, 5*inch, 4*inch, 2*inch, fill=1)

# Save
c.showPage()
c.save()
```

### Simple Platypus Document

```python
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import inch

doc = SimpleDocTemplate("output.pdf", pagesize=letter)
story = []
styles = getSampleStyleSheet()

# Add content
story.append(Paragraph("Document Title", styles['Title']))
story.append(Spacer(1, 0.2*inch))
story.append(Paragraph("This is body text with <b>bold</b> and <i>italic</i>.", styles['BodyText']))

# Build PDF
doc.build(story)
```

## Common Tasks

### Creating Tables

Tables work with both Canvas (via Drawing) and Platypus (as Flowables):

```python
from reportlab.platypus import Table, TableStyle
from reportlab.lib import colors
from reportlab.lib.units import inch

# Define data
data = [
    ['Product', 'Q1', 'Q2', 'Q3', 'Q4'],
    ['Widget A', '100', '150', '130', '180'],
    ['Widget B', '80', '120', '110', '160'],
]

# Create table
table = Table(data, colWidths=[2*inch, 1*inch, 1*inch, 1*inch, 1*inch])

# Apply styling
style = TableStyle([
    # Header row
    ('BACKGROUND', (0, 0), (-1, 0), colors.darkblue),
    ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
    ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
    ('ALIGN', (0, 0), (-1, -1), 'CENTER'),

    # Data rows
    ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.lightgrey]),
    ('GRID', (0, 0), (-1, -1), 1, colors.black),
])

table.setStyle(style)

# Add to Platypus story
story.append(table)

# Or draw on Canvas
table.wrapOn(c, width, height)
table.drawOn(c, x, y)
```

**Detailed table reference:** See `references/tables_reference.md` for cell spanning, borders, alignment, and advanced styling.

### Creating Charts

Charts use the graphics framework and can be added to both Canvas and Platypus:

```python
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.barcharts import VerticalBarChart
from reportlab.lib import colors

# Create drawing
drawing = Drawing(400, 200)

# Create chart
chart = VerticalBarChart()
chart.x = 50
chart.y = 50
chart.width = 300
chart.height = 125

# Set data
chart.data = [[100, 150, 130, 180, 140]]
chart.categoryAxis.categoryNames = ['Q1', 'Q2', 'Q3', 'Q4', 'Q5']

# Style
chart.bars[0].fillColor = colors.blue
chart.valueAxis.valueMin = 0
chart.valueAxis.valueMax = 200

# Add to drawing
drawing.add(chart)

# Use in Platypus
story.append(drawing)

# Or render directly to PDF
from reportlab.graphics import renderPDF
renderPDF.drawToFile(drawing, 'chart.pdf', 'Chart Title')
```

**Available chart types:** Bar (vertical/horizontal), Line, Pie, Area, Scatter
**Detailed charts reference:** See `references/charts_reference.md` for all chart types, styling, legends, and customization.

### Adding Barcodes and QR Codes

```python
from reportlab.graphics.barcode import code128
from reportlab.graphics.barcode.qr import QrCodeWidget
from reportlab.graphics.shapes import Drawing
from reportlab.graphics import renderPDF

# Code128 barcode (general purpose)
barcode = code128.Code128("ABC123456789", barHeight=0.5*inch)

# On Canvas
barcode.drawOn(c, x, y)

# QR Code
qr = QrCodeWidget("https://example.com")
qr.barWidth = 2*inch
qr.barHeight = 2*inch

# Wrap in Drawing for Platypus
d = Drawing()
d.add(qr)
story.append(d)
```

**Supported formats:** Code128, Code39, EAN-13, EAN-8, UPC-A, ISBN, QR, Data Matrix, and 20+ more
**Detailed barcode reference:** See `references/barcodes_reference.md` for all formats and usage examples.

### Working with Text and Fonts

```python
from reportlab.platypus import Paragraph
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_JUSTIFY

# Create custom style
custom_style = ParagraphStyle(
    'CustomStyle',
    fontSize=12,
    leading=14,           # Line spacing
    alignment=TA_JUSTIFY,
    spaceAfter=10,
    textColor=colors.black,
)

# Paragraph with inline formatting
text = """
This paragraph has <b>bold</b>, <i>italic</i>, and <u>underlined</u> text.
You can also use <font color="blue">colors</font> and <font size="14">different sizes</font>.
Chemical formula: H<sub>2</sub>O, Einstein: E=mc<sup>2</sup>
"""

para = Paragraph(text, custom_style)
story.append(para)
```

**Using custom fonts:**

```python
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

# Register TrueType font
pdfmetrics.registerFont(TTFont('CustomFont', 'CustomFont.ttf'))

# Use in Canvas
c.setFont('CustomFont', 12)

# Use in Paragraph style
style = ParagraphStyle('Custom', fontName='CustomFont', fontSize=12)
```

**Detailed text reference:** See `references/text_and_fonts.md` for paragraph styles, font families, Asian languages, Greek letters, and formatting.

### Adding Images

```python
from reportlab.platypus import Image
from reportlab.lib.units import inch

# In Platypus
img = Image('photo.jpg', width=4*inch, height=3*inch)
story.append(img)

# Maintain aspect ratio
img = Image('photo.jpg', width=4*inch, height=3*inch, kind='proportional')

# In Canvas
c.drawImage('photo.jpg', x, y, width=4*inch, height=3*inch)

# With transparency (mask white background)
c.drawImage('logo.png', x, y, mask=[255,255,255,255,255,255])
```

### Creating Forms

```python
from reportlab.pdfgen import canvas
from reportlab.lib.colors import black, white, lightgrey

c = canvas.Canvas("form.pdf")

# Text field
c.acroForm.textfield(
    name="name",
    tooltip="Enter your name",
    x=100, y=700,
    width=200, height=20,
    borderColor=black,
    fillColor=lightgrey,
    forceBorder=True
)

# Checkbox
c.acroForm.checkbox(
    name="agree",
    x=100, y=650,
    size=20,
    buttonStyle='check',
    checked=False
)

# Dropdown
c.acroForm.choice(
    name="country",
    x=100, y=600,
    width=150, height=20,
    options=[("United States", "US"), ("Canada", "CA")],
    forceBorder=True
)

c.save()
```

**Detailed PDF features reference:** See `references/pdf_features.md` for forms, links, bookmarks, encryption, and metadata.

### Headers and Footers

For Platypus documents, use page callbacks:

```python
from reportlab.platypus import BaseDocTemplate, PageTemplate, Frame

def add_header_footer(canvas, doc):
    """Called on each page"""
    canvas.saveState()

    # Header
    canvas.setFont('Helvetica', 9)
    canvas.drawString(inch, height - 0.5*inch, "Document Title")

    # Footer
    canvas.drawRightString(width - inch, 0.5*inch, f"Page {doc.page}")

    canvas.restoreState()

# Set up document
doc = BaseDocTemplate("output.pdf")
frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id='normal')
template = PageTemplate(id='normal', frames=[frame], onPage=add_header_footer)
doc.addPageTemplates([template])

# Build with story
doc.build(story)
```

## Helper Scripts

This skill includes helper scripts for common tasks:

### Quick Document Generator

Use `scripts/quick_document.py` for rapid document creation:

```python
from scripts.quick_document import create_simple_document, create_styled_table

# Simple document from content blocks
content = [
    {'type': 'heading', 'content': 'Introduction'},
    {'type': 'paragraph', 'content': 'Your text here...'},
    {'type': 'bullet', 'content': 'Bullet point'},
]

create_simple_document("output.pdf", "My Document", content_blocks=content)

# Styled tables with presets
data = [['Header1', 'Header2'], ['Data1', 'Data2']]
table = create_styled_table(data, style_name='striped')  # 'default', 'striped', 'minimal', 'report'
```

## Template Examples

Complete working examples in `assets/`:

### Invoice Template

`assets/invoice_template.py` - Professional invoice with:
- Company and client information
- Itemized table with calculations
- Tax and totals
- Terms and notes
- Logo placement

```python
from assets.invoice_template import create_invoice

create_invoice(
    filename="invoice.pdf",
    invoice_number="INV-2024-001",
    invoice_date="January 15, 2024",
    due_date="February 15, 2024",
    company_info={'name': 'Acme Corp', 'address': '...', 'phone': '...', 'email': '...'},
    client_info={'name': 'Client Name', ...},
    items=[
        {'description': 'Service', 'quantity': 1, 'unit_price': 500.00},
        ...
    ],
    tax_rate=0.08,
    notes="Thank you for your business!",
)
```

### Report Template

`assets/report_template.py` - Multi-page business report with:
- Cover page
- Table of contents
- Multiple sections with subsections
- Charts and tables
- Headers and footers

```python
from assets.report_template import create_report

report_data = {
    'title': 'Quarterly Report',
    'subtitle': 'Q4 2023',
    'author': 'Analytics Team',
    'sections': [
        {
            'title': 'Executive Summary',
            'content': 'Report content...',
            'table_data': {...},
            'chart_data': {...}
        },
        ...
    ]
}

create_report("report.pdf", report_data)
```

## Reference Documentation

Comprehensive API references organized by feature:

- **`references/canvas_api.md`** - Low-level Canvas: drawing primitives, coordinates, transformations, state management, images, paths
- **`references/platypus_guide.md`** - High-level Platypus: document templates, frames, flowables, page layouts, TOC
- **`references/text_and_fonts.md`** - Text formatting: paragraph styles, inline markup, custom fonts, Asian languages, bullets, sequences
- **`references/tables_reference.md`** - Tables: creation, styling, cell spanning, borders, alignment, colors, gradients
- **`references/charts_reference.md`** - Charts: all chart types, data handling, axes, legends, colors, rendering
- **`references/barcodes_reference.md`** - Barcodes: Code128, QR codes, EAN, UPC, postal codes, and 20+ formats
- **`references/pdf_features.md`** - PDF features: links, bookmarks, forms, encryption, metadata, page transitions

## Best Practices

### Coordinate System (Canvas)
- Origin (0, 0) is **lower-left corner** (not top-left)
- Y-axis points **upward**
- Units are in **points** (72 points = 1 inch)
- Always specify page size explicitly

```python
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch

width, height = letter
margin = inch

# Top of page
y_top = height - margin

# Bottom of page
y_bottom = margin
```

### Choosing Page Size

```python
from reportlab.lib.pagesizes import letter, A4, landscape

# US Letter (8.5" x 11")
pagesize=letter

# ISO A4 (210mm x 297mm)
pagesize=A4

# Landscape
pagesize=landscape(letter)

# Custom
pagesize=(6*inch, 9*inch)
```

### Performance Tips

1. **Use `drawImage()` over `drawInlineImage()`** - caches images for reuse
2. **Enable compression for large files:** `canvas.Canvas("file.pdf", pageCompression=1)`
3. **Reuse styles** - create once, use throughout document
4. **Use Forms/XObjects** for repeated graphics

### Common Patterns

**Centering text on Canvas:**
```python
text = "Centered Text"
text_width = c.stringWidth(text, "Helvetica", 12)
x = (width - text_width) / 2
c.drawString(x, y, text)

# Or use built-in
c.drawCentredString(width/2, y, text)
```

**Page breaks in Platypus:**
```python
from reportlab.platypus import PageBreak

story.append(PageBreak())
```

**Keep content together (no split):**
```python
from reportlab.platypus import KeepTogether

story.append(KeepTogether([
    heading,
    paragraph1,
    paragraph2,
]))
```

**Alternate row colors:**
```python
style = TableStyle([
    ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.lightgrey]),
])
```

## Troubleshooting

**Text overlaps or disappears:**
- Check Y-coordinates - remember origin is bottom-left
- Ensure text fits within page bounds
- Verify `leading` (line spacing) is greater than `fontSize`

**Table doesn't fit on page:**
- Reduce column widths
- Decrease font size
- Use landscape orientation
- Enable table splitting with `repeatRows`

**Barcode not scanning:**
- Increase `barHeight` (try 0.5 inch minimum)
- Set `quiet=1` for quiet zones
- Test print quality (300+ DPI recommended)
- Validate data format for barcode type

**Font not found:**
- Register TrueType fonts with `pdfmetrics.registerFont()`
- Use font family name exactly as registered
- Check font file path is correct

**Images have white background:**
- Use `mask` parameter to make white transparent
- Provide RGB range to mask: `mask=[255,255,255,255,255,255]`
- Or use PNG with alpha channel

## Example Workflows

### Creating an Invoice

1. Start with invoice template from `assets/invoice_template.py`
2. Customize company info, logo path
3. Add items with descriptions, quantities, prices
4. Set tax rate if applicable
5. Add notes and payment terms
6. Generate PDF

### Creating a Report

1. Start with report template from `assets/report_template.py`
2. Define sections with titles and content
3. Add tables for data using `create_styled_table()`
4. Add charts using graphics framework
5. Build with `doc.multiBuild(story)` for TOC

### Creating a Certificate

1. Use Canvas API for precise positioning
2. Load custom fonts for elegant typography
3. Add border graphics or image background
4. Position text elements (name, date, achievement)
5. Optional: Add QR code for verification

### Creating Labels with Barcodes

1. Use Canvas with custom page size (label dimensions)
2. Calculate grid positions for multiple labels per page
3. Draw label content (text, images)
4. Add barcode at specific position
5. Use `showPage()` between labels or grids

## Installation

```bash
uv pip install reportlab

# For image support
uv pip install pillow

# For charts
uv pip install reportlab[renderPM]

# For barcode support (included in reportlab)
# QR codes require: uv pip install qrcode
```

## When to Use This Skill

This skill should be used when:
- Generating PDF documents programmatically
- Creating invoices, receipts, or billing documents
- Building reports with tables and charts
- Generating certificates, badges, or credentials
- Creating shipping labels or product labels with barcodes
- Designing forms or fillable PDFs
- Producing multi-page documents with consistent formatting
- Converting data to PDF format for archival or distribution
- Creating custom layouts that require precise positioning

This skill provides comprehensive guidance for all ReportLab capabilities, from simple documents to complex multi-page reports with charts, tables, and interactive elements.
