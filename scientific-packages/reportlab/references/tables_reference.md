# Tables Reference

Comprehensive guide to creating and styling tables in ReportLab.

## Basic Table Creation

```python
from reportlab.platypus import Table, TableStyle
from reportlab.lib import colors

# Simple data (list of lists or tuples)
data = [
    ['Header 1', 'Header 2', 'Header 3'],
    ['Row 1, Col 1', 'Row 1, Col 2', 'Row 1, Col 3'],
    ['Row 2, Col 1', 'Row 2, Col 2', 'Row 2, Col 3'],
]

# Create table
table = Table(data)

# Add to story
story.append(table)
```

## Table Constructor

```python
table = Table(
    data,                    # Required: list of lists/tuples
    colWidths=None,          # List of column widths or single value
    rowHeights=None,         # List of row heights or single value
    style=None,              # TableStyle object
    splitByRow=1,            # Split across pages by rows (not columns)
    repeatRows=0,            # Number of header rows to repeat
    repeatCols=0,            # Number of header columns to repeat
    rowSplitRange=None,      # Tuple (start, end) of splittable rows
    spaceBefore=None,        # Space before table
    spaceAfter=None,         # Space after table
    cornerRadii=None,        # [TL, TR, BL, BR] for rounded corners
)
```

### Column Widths

```python
from reportlab.lib.units import inch

# Equal widths
table = Table(data, colWidths=2*inch)

# Different widths per column
table = Table(data, colWidths=[1.5*inch, 2*inch, 1*inch])

# Auto-calculate widths (default)
table = Table(data)

# Percentage-based (of available width)
table = Table(data, colWidths=[None, None, None])  # Equal auto-sizing
```

## Cell Content Types

### Text and Newlines

```python
# Newlines work in cells
data = [
    ['Line 1\nLine 2', 'Single line'],
    ['Another\nmulti-line\ncell', 'Text'],
]
```

### Paragraph Objects

```python
from reportlab.platypus import Paragraph
from reportlab.lib.styles import getSampleStyleSheet

styles = getSampleStyleSheet()

data = [
    [Paragraph("Formatted <b>bold</b> text", styles['Normal']),
     Paragraph("More <i>italic</i> text", styles['Normal'])],
]

table = Table(data)
```

### Images

```python
from reportlab.platypus import Image

data = [
    ['Description', Image('logo.png', width=1*inch, height=1*inch)],
    ['Product', Image('product.jpg', width=2*inch, height=1.5*inch)],
]

table = Table(data)
```

### Nested Tables

```python
# Create inner table
inner_data = [['A', 'B'], ['C', 'D']]
inner_table = Table(inner_data)

# Use in outer table
outer_data = [
    ['Label', inner_table],
    ['Other', 'Content'],
]

outer_table = Table(outer_data)
```

## TableStyle

Styles are applied using command lists:

```python
from reportlab.platypus import TableStyle
from reportlab.lib import colors

style = TableStyle([
    # Command format: ('COMMAND', (startcol, startrow), (endcol, endrow), *args)
    ('GRID', (0, 0), (-1, -1), 1, colors.black),  # Grid over all cells
    ('BACKGROUND', (0, 0), (-1, 0), colors.grey),  # Header background
    ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),  # Header text color
])

table = Table(data)
table.setStyle(style)
```

### Cell Coordinate System

- Columns and rows are 0-indexed: `(col, row)`
- Negative indices count from end: `-1` is last column/row
- `(0, 0)` is top-left cell
- `(-1, -1)` is bottom-right cell

```python
# Examples:
(0, 0), (2, 0)      # First three cells of header row
(0, 1), (-1, -1)    # All cells except header
(0, 0), (-1, -1)    # Entire table
```

## Styling Commands

### Text Formatting

```python
style = TableStyle([
    # Font name
    ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),

    # Font size
    ('FONTSIZE', (0, 0), (-1, -1), 10),

    # Text color
    ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
    ('TEXTCOLOR', (0, 1), (-1, -1), colors.black),

    # Combined font command
    ('FONT', (0, 0), (-1, 0), 'Helvetica-Bold', 12),  # name, size
])
```

### Alignment

```python
style = TableStyle([
    # Horizontal alignment: LEFT, CENTER, RIGHT, DECIMAL
    ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
    ('ALIGN', (0, 1), (0, -1), 'LEFT'),      # First column left
    ('ALIGN', (1, 1), (-1, -1), 'RIGHT'),    # Other columns right

    # Vertical alignment: TOP, MIDDLE, BOTTOM
    ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
    ('VALIGN', (0, 0), (-1, 0), 'BOTTOM'),   # Header bottom-aligned
])
```

### Cell Padding

```python
style = TableStyle([
    # Individual padding
    ('LEFTPADDING', (0, 0), (-1, -1), 12),
    ('RIGHTPADDING', (0, 0), (-1, -1), 12),
    ('TOPPADDING', (0, 0), (-1, -1), 6),
    ('BOTTOMPADDING', (0, 0), (-1, -1), 6),

    # Or set all at once by setting each
])
```

### Background Colors

```python
style = TableStyle([
    # Solid background
    ('BACKGROUND', (0, 0), (-1, 0), colors.blue),
    ('BACKGROUND', (0, 1), (-1, -1), colors.lightgrey),

    # Alternating row colors
    ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.lightblue]),

    # Alternating column colors
    ('COLBACKGROUNDS', (0, 0), (-1, -1), [colors.white, colors.lightgrey]),
])
```

### Gradient Backgrounds

```python
from reportlab.lib.colors import Color

style = TableStyle([
    # Vertical gradient (top to bottom)
    ('BACKGROUND', (0, 0), (-1, 0), colors.blue),
    ('VERTICALGRADIENT', (0, 0), (-1, 0),
     [colors.blue, colors.lightblue]),

    # Horizontal gradient (left to right)
    ('HORIZONTALGRADIENT', (0, 1), (-1, 1),
     [colors.red, colors.yellow]),
])
```

### Lines and Borders

```python
style = TableStyle([
    # Complete grid
    ('GRID', (0, 0), (-1, -1), 1, colors.black),

    # Box/outline only
    ('BOX', (0, 0), (-1, -1), 2, colors.black),
    ('OUTLINE', (0, 0), (-1, -1), 2, colors.black),  # Same as BOX

    # Inner grid only
    ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.grey),

    # Directional lines
    ('LINEABOVE', (0, 0), (-1, 0), 2, colors.black),   # Header border
    ('LINEBELOW', (0, 0), (-1, 0), 1, colors.black),   # Header bottom
    ('LINEBEFORE', (0, 0), (0, -1), 1, colors.black),  # Left border
    ('LINEAFTER', (-1, 0), (-1, -1), 1, colors.black), # Right border

    # Thickness and color
    ('LINEABOVE', (0, 1), (-1, 1), 0.5, colors.grey),  # Thin grey line
])
```

### Cell Spanning

```python
data = [
    ['Spanning Header', '', ''],           # Span will merge these
    ['A', 'B', 'C'],
    ['D', 'E', 'F'],
]

style = TableStyle([
    # Span 3 columns in first row
    ('SPAN', (0, 0), (2, 0)),

    # Center the spanning cell
    ('ALIGN', (0, 0), (2, 0), 'CENTER'),
])

table = Table(data)
table.setStyle(style)
```

**Important:** Cells that are spanned over must contain empty strings `''`.

### Advanced Spanning Examples

```python
# Span multiple rows and columns
data = [
    ['A', 'B', 'B', 'C'],
    ['A', 'D', 'E', 'F'],
    ['A', 'G', 'H', 'I'],
]

style = TableStyle([
    # Span rows in column 0
    ('SPAN', (0, 0), (0, 2)),  # Merge A cells vertically

    # Span columns in row 0
    ('SPAN', (1, 0), (2, 0)),  # Merge B cells horizontally

    ('GRID', (0, 0), (-1, -1), 1, colors.black),
])
```

## Special Commands

### Rounded Corners

```python
table = Table(data, cornerRadii=[5, 5, 5, 5])  # [TL, TR, BL, BR]

# Or in style
style = TableStyle([
    ('ROUNDEDCORNERS', [10, 10, 0, 0]),  # Rounded top corners only
])
```

### No Split

Prevent table from splitting at specific locations:

```python
style = TableStyle([
    # Don't split between rows 0 and 2
    ('NOSPLIT', (0, 0), (-1, 2)),
])
```

### Split-Specific Styling

Apply styles only to first or last part when table splits:

```python
style = TableStyle([
    # Style for first part after split
    ('LINEBELOW', (0, 'splitfirst'), (-1, 'splitfirst'), 2, colors.red),

    # Style for last part after split
    ('LINEABOVE', (0, 'splitlast'), (-1, 'splitlast'), 2, colors.blue),
])
```

## Repeating Headers

```python
# Repeat first row on each page
table = Table(data, repeatRows=1)

# Repeat first 2 rows
table = Table(data, repeatRows=2)
```

## Complete Examples

### Styled Report Table

```python
from reportlab.platypus import Table, TableStyle
from reportlab.lib import colors
from reportlab.lib.units import inch

data = [
    ['Product', 'Quantity', 'Unit Price', 'Total'],
    ['Widget A', '10', '$5.00', '$50.00'],
    ['Widget B', '5', '$12.00', '$60.00'],
    ['Widget C', '20', '$3.00', '$60.00'],
    ['', '', 'Subtotal:', '$170.00'],
]

table = Table(data, colWidths=[2.5*inch, 1*inch, 1*inch, 1*inch])

style = TableStyle([
    # Header row
    ('BACKGROUND', (0, 0), (-1, 0), colors.darkblue),
    ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
    ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
    ('FONTSIZE', (0, 0), (-1, 0), 12),
    ('ALIGN', (0, 0), (-1, 0), 'CENTER'),
    ('BOTTOMPADDING', (0, 0), (-1, 0), 12),

    # Data rows
    ('BACKGROUND', (0, 1), (-1, -2), colors.beige),
    ('GRID', (0, 0), (-1, -2), 0.5, colors.grey),
    ('ALIGN', (1, 1), (-1, -1), 'RIGHT'),
    ('ALIGN', (0, 1), (0, -1), 'LEFT'),

    # Total row
    ('BACKGROUND', (0, -1), (-1, -1), colors.lightgrey),
    ('LINEABOVE', (0, -1), (-1, -1), 2, colors.black),
    ('FONTNAME', (2, -1), (-1, -1), 'Helvetica-Bold'),
])

table.setStyle(style)
```

### Alternating Row Colors

```python
data = [
    ['Name', 'Age', 'City'],
    ['Alice', '30', 'New York'],
    ['Bob', '25', 'Boston'],
    ['Charlie', '35', 'Chicago'],
    ['Diana', '28', 'Denver'],
]

table = Table(data, colWidths=[2*inch, 1*inch, 1.5*inch])

style = TableStyle([
    # Header
    ('BACKGROUND', (0, 0), (-1, 0), colors.darkslategray),
    ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
    ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),

    # Alternating rows (zebra striping)
    ('ROWBACKGROUNDS', (0, 1), (-1, -1),
     [colors.white, colors.lightgrey]),

    # Borders
    ('BOX', (0, 0), (-1, -1), 2, colors.black),
    ('LINEBELOW', (0, 0), (-1, 0), 2, colors.black),

    # Padding
    ('LEFTPADDING', (0, 0), (-1, -1), 12),
    ('RIGHTPADDING', (0, 0), (-1, -1), 12),
    ('TOPPADDING', (0, 0), (-1, -1), 6),
    ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
])

table.setStyle(style)
```

## Best Practices

1. **Set colWidths explicitly** for consistent layout
2. **Use repeatRows** for multi-page tables with headers
3. **Apply padding** for better readability (especially LEFTPADDING and RIGHTPADDING)
4. **Use ROWBACKGROUNDS** for alternating colors instead of styling each row
5. **Put empty strings** in cells that will be spanned
6. **Test page breaks** early with realistic data amounts
7. **Use Paragraph objects** in cells for complex formatted text
8. **Set VALIGN to MIDDLE** for better appearance with varying row heights
9. **Keep tables simple** - complex nested tables are hard to maintain
10. **Use consistent styling** - define once, apply to all tables
