# Barcodes Reference

Comprehensive guide to creating barcodes and QR codes in ReportLab.

## Available Barcode Types

ReportLab supports a wide range of 1D and 2D barcode formats.

### 1D Barcodes (Linear)

- **Code128** - Compact, encodes full ASCII
- **Code39** (Standard39) - Alphanumeric, widely supported
- **Code93** (Standard93) - Compressed Code39
- **EAN-13** - European Article Number (retail)
- **EAN-8** - Short form of EAN
- **EAN-5** - 5-digit add-on (pricing)
- **UPC-A** - Universal Product Code (North America)
- **ISBN** - International Standard Book Number
- **Code11** - Telecommunications
- **Codabar** - Blood banks, FedEx, libraries
- **I2of5** (Interleaved 2 of 5) - Warehouse/distribution
- **MSI** - Inventory control
- **POSTNET** - US Postal Service
- **USPS_4State** - US Postal Service
- **FIM** (A, B, C, D) - Facing Identification Mark (mail sorting)

### 2D Barcodes

- **QR** - QR Code (widely used for URLs, contact info)
- **ECC200DataMatrix** - Data Matrix format

## Using Barcodes with Canvas

### Code128 (Recommended for General Use)

Code128 is versatile and compact - encodes full ASCII character set with mandatory checksum.

```python
from reportlab.pdfgen import canvas
from reportlab.graphics.barcode import code128
from reportlab.lib.units import inch

c = canvas.Canvas("barcode.pdf")

# Create barcode
barcode = code128.Code128("HELLO123")

# Draw on canvas
barcode.drawOn(c, 1*inch, 5*inch)

c.save()
```

### Code128 Options

```python
barcode = code128.Code128(
    value="ABC123",      # Required: data to encode
    barWidth=0.01*inch,  # Width of narrowest bar
    barHeight=0.5*inch,  # Height of bars
    quiet=1,             # Add quiet zones (margins)
    lquiet=None,         # Left quiet zone width
    rquiet=None,         # Right quiet zone width
    stop=1,              # Show stop symbol
)

# Draw with specific size
barcode.drawOn(canvas, x, y)

# Get dimensions
width = barcode.width
height = barcode.height
```

### Code39 (Standard39)

Supports: 0-9, A-Z (uppercase), space, and special chars (-.$/+%*).

```python
from reportlab.graphics.barcode import code39

barcode = code39.Standard39(
    value="HELLO",
    barWidth=0.01*inch,
    barHeight=0.5*inch,
    quiet=1,
    checksum=0,  # 0 or 1
)

barcode.drawOn(canvas, x, y)
```

### Extended Code39

Encodes full ASCII (pairs of Code39 characters).

```python
from reportlab.graphics.barcode import code39

barcode = code39.Extended39(
    value="Hello World!",  # Can include lowercase and symbols
    barWidth=0.01*inch,
    barHeight=0.5*inch,
)

barcode.drawOn(canvas, x, y)
```

### Code93

```python
from reportlab.graphics.barcode import code93

# Standard93 - uppercase, digits, some symbols
barcode = code93.Standard93(
    value="HELLO93",
    barWidth=0.01*inch,
    barHeight=0.5*inch,
)

# Extended93 - full ASCII
barcode = code93.Extended93(
    value="Hello 93!",
    barWidth=0.01*inch,
    barHeight=0.5*inch,
)

barcode.drawOn(canvas, x, y)
```

### EAN-13 (European Article Number)

13-digit barcode for retail products.

```python
from reportlab.graphics.barcode import eanbc

# Must be exactly 12 digits (13th is calculated checksum)
barcode = eanbc.Ean13BarcodeWidget(
    value="123456789012"
)

# Draw
from reportlab.graphics import renderPDF
from reportlab.graphics.shapes import Drawing

d = Drawing()
d.add(barcode)
renderPDF.draw(d, canvas, x, y)
```

### EAN-8

Short form, 8 digits.

```python
from reportlab.graphics.barcode import eanbc

# Must be exactly 7 digits (8th is calculated)
barcode = eanbc.Ean8BarcodeWidget(
    value="1234567"
)
```

### UPC-A

12-digit barcode used in North America.

```python
from reportlab.graphics.barcode import usps

# 11 digits (12th is checksum)
barcode = usps.UPCA(
    value="01234567890"
)

barcode.drawOn(canvas, x, y)
```

### ISBN (Books)

```python
from reportlab.graphics.barcode.widgets import ISBNBarcodeWidget

# 10 or 13 digit ISBN
barcode = ISBNBarcodeWidget(
    value="978-0-123456-78-9"
)

# With pricing (EAN-5 add-on)
barcode = ISBNBarcodeWidget(
    value="978-0-123456-78-9",
    price=True,
)
```

### QR Codes

Most versatile 2D barcode - can encode URLs, text, contact info, etc.

```python
from reportlab.graphics.barcode.qr import QrCodeWidget
from reportlab.graphics.shapes import Drawing
from reportlab.graphics import renderPDF

# Create QR code
qr = QrCodeWidget("https://example.com")

# Size in pixels (QR codes are square)
qr.barWidth = 100  # Width in points
qr.barHeight = 100  # Height in points

# Error correction level
# L = 7% recovery, M = 15%, Q = 25%, H = 30%
qr.qrVersion = 1  # Auto-size (1-40, or None for auto)
qr.errorLevel = 'M'  # L, M, Q, H

# Draw
d = Drawing()
d.add(qr)
renderPDF.draw(d, canvas, x, y)
```

### QR Code - More Options

```python
# URL QR Code
qr = QrCodeWidget("https://example.com")

# Contact information (vCard)
vcard_data = """BEGIN:VCARD
VERSION:3.0
FN:John Doe
TEL:+1-555-1234
EMAIL:john@example.com
END:VCARD"""
qr = QrCodeWidget(vcard_data)

# WiFi credentials
wifi_data = "WIFI:T:WPA;S:NetworkName;P:Password;;"
qr = QrCodeWidget(wifi_data)

# Plain text
qr = QrCodeWidget("Any text here")
```

### Data Matrix (ECC200)

Compact 2D barcode for small items.

```python
from reportlab.graphics.barcode.datamatrix import DataMatrixWidget

barcode = DataMatrixWidget(
    value="DATA123"
)

d = Drawing()
d.add(barcode)
renderPDF.draw(d, canvas, x, y)
```

### Postal Barcodes

```python
from reportlab.graphics.barcode import usps

# POSTNET (older format)
barcode = usps.POSTNET(
    value="55555-1234",  # ZIP or ZIP+4
)

# USPS 4-State (newer)
barcode = usps.USPS_4State(
    value="12345678901234567890",  # 20-digit routing code
    routing="12345678901"
)

barcode.drawOn(canvas, x, y)
```

### FIM (Facing Identification Mark)

Used for mail sorting.

```python
from reportlab.graphics.barcode import usps

# FIM-A, FIM-B, FIM-C, or FIM-D
barcode = usps.FIM(
    value="A"  # A, B, C, or D
)

barcode.drawOn(canvas, x, y)
```

## Using Barcodes with Platypus

For flowing documents, wrap barcodes in Flowables.

### Simple Approach - Drawing Flowable

```python
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.barcode.qr import QrCodeWidget
from reportlab.lib.units import inch

# Create drawing
d = Drawing(2*inch, 2*inch)

# Create barcode
qr = QrCodeWidget("https://example.com")
qr.barWidth = 2*inch
qr.barHeight = 2*inch
qr.x = 0
qr.y = 0

d.add(qr)

# Add to story
story.append(d)
```

### Custom Flowable Wrapper

```python
from reportlab.platypus import Flowable
from reportlab.graphics.barcode import code128
from reportlab.lib.units import inch

class BarcodeFlowable(Flowable):
    def __init__(self, code, barcode_type='code128', width=2*inch, height=0.5*inch):
        Flowable.__init__(self)
        self.code = code
        self.barcode_type = barcode_type
        self.width_val = width
        self.height_val = height

        # Create barcode
        if barcode_type == 'code128':
            self.barcode = code128.Code128(code, barWidth=width/100, barHeight=height)
        # Add other types as needed

    def draw(self):
        self.barcode.drawOn(self.canv, 0, 0)

    def wrap(self, availWidth, availHeight):
        return (self.barcode.width, self.barcode.height)

# Use in story
story.append(BarcodeFlowable("PRODUCT123"))
```

## Complete Examples

### Product Label with Barcode

```python
from reportlab.pdfgen import canvas
from reportlab.graphics.barcode import code128
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch

def create_product_label(filename, product_code, product_name):
    c = canvas.Canvas(filename, pagesize=(4*inch, 2*inch))

    # Product name
    c.setFont("Helvetica-Bold", 14)
    c.drawCentredString(2*inch, 1.5*inch, product_name)

    # Barcode
    barcode = code128.Code128(product_code)
    barcode_width = barcode.width
    barcode_height = barcode.height

    # Center barcode
    x = (4*inch - barcode_width) / 2
    y = 0.5*inch

    barcode.drawOn(c, x, y)

    # Code text
    c.setFont("Courier", 10)
    c.drawCentredString(2*inch, 0.3*inch, product_code)

    c.save()

create_product_label("label.pdf", "ABC123456789", "Premium Widget")
```

### QR Code Contact Card

```python
from reportlab.pdfgen import canvas
from reportlab.graphics.barcode.qr import QrCodeWidget
from reportlab.graphics.shapes import Drawing
from reportlab.graphics import renderPDF
from reportlab.lib.units import inch

def create_contact_card(filename, name, phone, email):
    c = canvas.Canvas(filename, pagesize=(3.5*inch, 2*inch))

    # Contact info
    c.setFont("Helvetica-Bold", 12)
    c.drawString(0.5*inch, 1.5*inch, name)
    c.setFont("Helvetica", 10)
    c.drawString(0.5*inch, 1.3*inch, phone)
    c.drawString(0.5*inch, 1.1*inch, email)

    # Create vCard data
    vcard = f"""BEGIN:VCARD
VERSION:3.0
FN:{name}
TEL:{phone}
EMAIL:{email}
END:VCARD"""

    # QR code
    qr = QrCodeWidget(vcard)
    qr.barWidth = 1.5*inch
    qr.barHeight = 1.5*inch

    d = Drawing()
    d.add(qr)

    renderPDF.draw(d, c, 1.8*inch, 0.2*inch)

    c.save()

create_contact_card("contact.pdf", "John Doe", "+1-555-1234", "john@example.com")
```

### Shipping Label with Multiple Barcodes

```python
from reportlab.pdfgen import canvas
from reportlab.graphics.barcode import code128
from reportlab.lib.units import inch

def create_shipping_label(filename, tracking_code, zip_code):
    c = canvas.Canvas(filename, pagesize=(6*inch, 4*inch))

    # Title
    c.setFont("Helvetica-Bold", 16)
    c.drawString(0.5*inch, 3.5*inch, "SHIPPING LABEL")

    # Tracking barcode
    c.setFont("Helvetica", 10)
    c.drawString(0.5*inch, 2.8*inch, "Tracking Number:")

    tracking_barcode = code128.Code128(tracking_code, barHeight=0.5*inch)
    tracking_barcode.drawOn(c, 0.5*inch, 2*inch)

    c.setFont("Courier", 9)
    c.drawString(0.5*inch, 1.8*inch, tracking_code)

    # Additional info can be added

    c.save()

create_shipping_label("shipping.pdf", "1Z999AA10123456784", "12345")
```

## Barcode Selection Guide

**Choose Code128 when:**
- General purpose encoding
- Need to encode numbers and letters
- Want compact size
- Widely supported

**Choose Code39 when:**
- Older systems require it
- Don't need lowercase letters
- Want maximum compatibility

**Choose QR Code when:**
- Need to encode URLs
- Want mobile device scanning
- Need high data capacity
- Want error correction

**Choose EAN/UPC when:**
- Retail product identification
- Need industry-standard format
- Global distribution

**Choose Data Matrix when:**
- Very limited space
- Small items (PCB, electronics)
- Need 2D compact format

## Best Practices

1. **Test scanning** early with actual barcode scanners/readers
2. **Add quiet zones** (white space) around barcodes - set `quiet=1`
3. **Choose appropriate height** - taller barcodes are easier to scan
4. **Include human-readable text** below barcode for manual entry
5. **Use Code128** as default for general purpose - it's compact and versatile
6. **For URLs, use QR codes** - much easier for mobile users
7. **Check barcode standards** for your industry (retail uses EAN/UPC)
8. **Test print quality** - low DPI can make barcodes unscannable
9. **Validate data** before encoding - wrong check digits cause issues
10. **Consider error correction** for QR codes - use 'M' or 'H' for important data
