# Canvas API Reference

The Canvas API provides low-level, precise control over PDF generation using coordinate-based drawing.

## Coordinate System

- Origin (0, 0) is at the **lower-left corner** (not top-left like web graphics)
- X-axis points right, Y-axis points upward
- Units are in points (72 points = 1 inch)
- Default page size is A4; explicitly specify page size for consistency

## Basic Setup

```python
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter, A4
from reportlab.lib.units import inch

# Create canvas
c = canvas.Canvas("output.pdf", pagesize=letter)

# Get page dimensions
width, height = letter

# Draw content
c.drawString(100, 100, "Hello World")

# Finish page and save
c.showPage()  # Complete current page
c.save()      # Write PDF to disk
```

## Text Drawing

### Basic String Methods
```python
# Basic text placement
c.drawString(x, y, text)           # Left-aligned at x, y
c.drawRightString(x, y, text)      # Right-aligned at x, y
c.drawCentredString(x, y, text)    # Center-aligned at x, y

# Font control
c.setFont(fontname, size)          # e.g., "Helvetica", 12
c.setFillColor(color)              # Text color
```

### Text Objects (Advanced)
For complex text operations with multiple lines and precise control:

```python
t = c.beginText(x, y)
t.setFont("Times-Roman", 14)
t.textLine("First line")
t.textLine("Second line")
t.setTextOrigin(x, y)  # Reset position
c.drawText(t)
```

## Drawing Primitives

### Lines
```python
c.line(x1, y1, x2, y2)                    # Single line
c.lines([(x1,y1,x2,y2), (x3,y3,x4,y4)])  # Multiple lines
c.grid(xlist, ylist)                      # Grid from coordinate lists
```

### Shapes
```python
c.rect(x, y, width, height, stroke=1, fill=0)
c.roundRect(x, y, width, height, radius, stroke=1, fill=0)
c.circle(x_ctr, y_ctr, r, stroke=1, fill=0)
c.ellipse(x1, y1, x2, y2, stroke=1, fill=0)
c.wedge(x, y, radius, startAng, extent, stroke=1, fill=0)
```

### Bezier Curves
```python
c.bezier(x1, y1, x2, y2, x3, y3, x4, y4)
```

## Path Objects

For complex shapes, use path objects:

```python
p = c.beginPath()
p.moveTo(x, y)              # Move without drawing
p.lineTo(x, y)              # Draw line to point
p.curveTo(x1, y1, x2, y2, x3, y3)  # Bezier curve
p.arc(x1, y1, x2, y2, startAng, extent)
p.arcTo(x1, y1, x2, y2, startAng, extent)
p.close()                   # Close path to start point

# Draw the path
c.drawPath(p, stroke=1, fill=0)
```

## Colors

### RGB (Screen Display)
```python
from reportlab.lib.colors import red, blue, Color

c.setFillColorRGB(r, g, b)      # r, g, b are 0-1
c.setStrokeColorRGB(r, g, b)
c.setFillColor(red)             # Named colors
c.setStrokeColor(blue)

# Custom with alpha transparency
c.setFillColor(Color(0.5, 0, 0, alpha=0.5))
```

### CMYK (Professional Printing)
```python
from reportlab.lib.colors import CMYKColor, PCMYKColor

c.setFillColorCMYK(c, m, y, k)  # 0-1 range
c.setStrokeColorCMYK(c, m, y, k)

# Integer percentages (0-100)
c.setFillColor(PCMYKColor(100, 50, 0, 0))
```

## Line Styling

```python
c.setLineWidth(width)           # Thickness in points
c.setLineCap(mode)              # 0=butt, 1=round, 2=square
c.setLineJoin(mode)             # 0=miter, 1=round, 2=bevel
c.setDash(array, phase)         # e.g., [3, 3] for dotted line
```

## Coordinate Transformations

**IMPORTANT:** Transformations are incremental and cumulative.

```python
# Translation (move origin)
c.translate(dx, dy)

# Rotation (in degrees, counterclockwise)
c.rotate(theta)

# Scaling
c.scale(xscale, yscale)

# Skewing
c.skew(alpha, beta)
```

### State Management
```python
# Save current graphics state
c.saveState()

# ... apply transformations and draw ...

# Restore previous state
c.restoreState()
```

**Note:** State cannot be preserved across `showPage()` calls.

## Images

```python
from reportlab.lib.utils import ImageReader

# Preferred method (with caching)
c.drawImage(image_source, x, y, width=None, height=None,
            mask=None, preserveAspectRatio=False)

# image_source can be:
# - Filename string
# - PIL Image object
# - ImageReader object

# For transparency, specify RGB mask range
c.drawImage("logo.png", 100, 500, mask=[255, 255, 255, 255, 255, 255])

# Inline (inefficient, no caching)
c.drawInlineImage(image_source, x, y, width=None, height=None)
```

## Page Management

```python
# Complete current page
c.showPage()

# Set page size for next page
c.setPageSize(size)  # e.g., letter, A4

# Page compression (smaller files, slower generation)
c = canvas.Canvas("output.pdf", pageCompression=1)
```

## Common Patterns

### Margins and Layout
```python
from reportlab.lib.units import inch
from reportlab.lib.pagesizes import letter

width, height = letter
margin = inch

# Draw within margins
content_width = width - 2*margin
content_height = height - 2*margin

# Text at top margin
c.drawString(margin, height - margin, "Header")

# Text at bottom margin
c.drawString(margin, margin, "Footer")
```

### Headers and Footers
```python
def draw_header_footer(c, width, height):
    c.saveState()
    c.setFont("Helvetica", 9)
    c.drawString(inch, height - 0.5*inch, "Company Name")
    c.drawRightString(width - inch, 0.5*inch, f"Page {c.getPageNumber()}")
    c.restoreState()

# Call on each page
draw_header_footer(c, width, height)
c.showPage()
```

## Best Practices

1. **Always specify page size** - Different platforms have different defaults
2. **Use variables for measurements** - `margin = inch` instead of hardcoded values
3. **Match saveState/restoreState** - Always balance these calls
4. **Apply transformations externally** for engineering drawings to prevent line width scaling
5. **Use drawImage over drawInlineImage** for better performance with repeated images
6. **Draw from bottom-up** - Remember Y-axis points upward
