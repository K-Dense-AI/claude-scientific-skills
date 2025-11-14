# Text and Fonts Reference

Comprehensive guide to text formatting, paragraph styles, and font handling in ReportLab.

## Text Encoding

**IMPORTANT:** All text input should be UTF-8 encoded or Python Unicode objects (since ReportLab 2.0).

```python
# Correct - UTF-8 strings
text = "Hello 世界 مرحبا"
para = Paragraph(text, style)

# For legacy data, convert first
import codecs
decoded_text = codecs.decode(legacy_bytes, 'latin-1')
```

## Paragraph Styles

### Creating Styles

```python
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_RIGHT, TA_JUSTIFY
from reportlab.lib.colors import black, blue, red
from reportlab.lib.units import inch

# Get default styles
styles = getSampleStyleSheet()
normal = styles['Normal']
heading = styles['Heading1']

# Create custom style
custom_style = ParagraphStyle(
    'CustomStyle',
    parent=normal,           # Inherit from another style

    # Font properties
    fontName='Helvetica',
    fontSize=12,
    leading=14,              # Line spacing (should be > fontSize)

    # Indentation (in points)
    leftIndent=0,
    rightIndent=0,
    firstLineIndent=0,       # Positive = indent, negative = outdent

    # Spacing
    spaceBefore=0,
    spaceAfter=0,

    # Alignment
    alignment=TA_LEFT,       # TA_LEFT, TA_CENTER, TA_RIGHT, TA_JUSTIFY

    # Colors
    textColor=black,
    backColor=None,          # Background color

    # Borders
    borderWidth=0,
    borderColor=None,
    borderPadding=0,
    borderRadius=None,

    # Bullets
    bulletFontName='Helvetica',
    bulletFontSize=12,
    bulletIndent=0,
    bulletText=None,         # Text for bullets (e.g., '•')

    # Advanced
    wordWrap=None,           # 'CJK' for Asian languages
    allowWidows=1,           # Allow widow lines
    allowOrphans=0,          # Prevent orphan lines
    endDots=None,            # Trailing dots for TOC entries
    splitLongWords=1,
    hyphenationLang=None,    # 'en_US', etc. (requires pyphen)
)

# Add to stylesheet
styles.add(custom_style)
```

### Built-in Styles

```python
styles = getSampleStyleSheet()

# Common styles
styles['Normal']         # Body text
styles['BodyText']       # Similar to Normal
styles['Heading1']       # Top-level heading
styles['Heading2']       # Second-level heading
styles['Heading3']       # Third-level heading
styles['Title']          # Document title
styles['Bullet']         # Bulleted list items
styles['Definition']     # Definition text
styles['Code']           # Code samples
```

## Paragraph Formatting

### Basic Paragraph

```python
from reportlab.platypus import Paragraph

para = Paragraph("This is a paragraph.", style)
story.append(para)
```

### Inline Formatting Tags

```python
text = """
<b>Bold text</b>
<i>Italic text</i>
<u>Underlined text</u>
<strike>Strikethrough text</strike>
<strong>Strong (bold) text</strong>
"""

para = Paragraph(text, normal_style)
```

### Font Control

```python
text = """
<font face="Courier" size="14" color="blue">
Custom font, size, and color
</font>

<font color="#FF0000">Hex color codes work too</font>
"""

para = Paragraph(text, normal_style)
```

### Superscripts and Subscripts

```python
text = """
H<sub>2</sub>O is water.
E=mc<super>2</super> or E=mc<sup>2</sup>
X<sub><i>i</i></sub> for subscripted variables
"""

para = Paragraph(text, normal_style)
```

### Greek Letters

```python
text = """
<greek>alpha</greek>, <greek>beta</greek>, <greek>gamma</greek>
<greek>epsilon</greek>, <greek>pi</greek>, <greek>omega</greek>
"""

para = Paragraph(text, normal_style)
```

### Links

```python
# External link
text = '<link href="https://example.com" color="blue">Click here</link>'

# Internal link (to bookmark)
text = '<link href="#section1" color="blue">Go to Section 1</link>'

# Anchor for internal links
text = '<a name="section1"/>Section 1 Heading'

para = Paragraph(text, normal_style)
```

### Inline Images

```python
text = """
Here is an inline image: <img src="icon.png" width="12" height="12" valign="middle"/>
"""

para = Paragraph(text, normal_style)
```

### Line Breaks

```python
text = """
First line<br/>
Second line<br/>
Third line
"""

para = Paragraph(text, normal_style)
```

## Font Handling

### Standard Fonts

ReportLab includes 14 standard PDF fonts (no embedding needed):

```python
# Helvetica family
'Helvetica'
'Helvetica-Bold'
'Helvetica-Oblique'
'Helvetica-BoldOblique'

# Times family
'Times-Roman'
'Times-Bold'
'Times-Italic'
'Times-BoldItalic'

# Courier family
'Courier'
'Courier-Bold'
'Courier-Oblique'
'Courier-BoldOblique'

# Symbol and Dingbats
'Symbol'
'ZapfDingbats'
```

### TrueType Fonts

```python
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

# Register single font
pdfmetrics.registerFont(TTFont('CustomFont', 'CustomFont.ttf'))

# Use in Canvas
canvas.setFont('CustomFont', 12)

# Use in Paragraph style
style = ParagraphStyle('Custom', fontName='CustomFont', fontSize=12)
```

### Font Families

Register related fonts as a family for bold/italic support:

```python
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.fonts import addMapping

# Register fonts
pdfmetrics.registerFont(TTFont('Vera', 'Vera.ttf'))
pdfmetrics.registerFont(TTFont('VeraBd', 'VeraBd.ttf'))
pdfmetrics.registerFont(TTFont('VeraIt', 'VeraIt.ttf'))
pdfmetrics.registerFont(TTFont('VeraBI', 'VeraBI.ttf'))

# Map family (normal, bold, italic, bold-italic)
addMapping('Vera', 0, 0, 'Vera')       # normal
addMapping('Vera', 1, 0, 'VeraBd')     # bold
addMapping('Vera', 0, 1, 'VeraIt')     # italic
addMapping('Vera', 1, 1, 'VeraBI')     # bold-italic

# Now <b> and <i> tags work with this family
style = ParagraphStyle('VeraStyle', fontName='Vera', fontSize=12)
para = Paragraph("Normal <b>Bold</b> <i>Italic</i> <b><i>Both</i></b>", style)
```

### Font Search Paths

```python
from reportlab.pdfbase.ttfonts import TTFSearchPath

# Add custom font directory
TTFSearchPath.append('/path/to/fonts/')

# Now fonts in this directory can be found by name
pdfmetrics.registerFont(TTFont('MyFont', 'MyFont.ttf'))
```

### Asian Language Support

#### Using Adobe Language Packs (no embedding)

```python
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.cidfonts import UnicodeCIDFont

# Register CID fonts
pdfmetrics.registerFont(UnicodeCIDFont('HeiseiMin-W3'))    # Japanese
pdfmetrics.registerFont(UnicodeCIDFont('STSong-Light'))    # Chinese (Simplified)
pdfmetrics.registerFont(UnicodeCIDFont('MSung-Light'))     # Chinese (Traditional)
pdfmetrics.registerFont(UnicodeCIDFont('HYSMyeongJo-Medium'))  # Korean

# Use in styles
style = ParagraphStyle('Japanese', fontName='HeiseiMin-W3', fontSize=12)
para = Paragraph("日本語テキスト", style)
```

#### Using TrueType Fonts with Asian Characters

```python
# Register TrueType font with full Unicode support
pdfmetrics.registerFont(TTFont('SimSun', 'simsun.ttc'))

style = ParagraphStyle('Chinese', fontName='SimSun', fontSize=12, wordWrap='CJK')
para = Paragraph("中文文本", style)
```

Note: Set `wordWrap='CJK'` for proper line breaking in Asian languages.

## Numbering and Sequences

Auto-numbering using `<seq>` tags:

```python
# Simple numbering
text = "<seq id='chapter'/> Introduction"  # Outputs: 1 Introduction
text = "<seq id='chapter'/> Methods"       # Outputs: 2 Methods

# Reset counter
text = "<seq id='figure' reset='yes'/>"

# Formatting templates
text = "Figure <seq template='%(chapter)s-%(figure+)s' id='figure'/>"
# Outputs: Figure 1-1, Figure 1-2, etc.

# Multi-level numbering
text = "Section <seq template='%(chapter)s.%(section+)s' id='section'/>"
```

## Bullets and Lists

### Using Bullet Style

```python
bullet_style = ParagraphStyle(
    'Bullet',
    parent=normal_style,
    leftIndent=20,
    bulletIndent=10,
    bulletText='•',          # Unicode bullet
    bulletFontName='Helvetica',
)

story.append(Paragraph("First item", bullet_style))
story.append(Paragraph("Second item", bullet_style))
story.append(Paragraph("Third item", bullet_style))
```

### Custom Bullet Characters

```python
# Different bullet styles
bulletText='•'     # Filled circle
bulletText='◦'     # Open circle
bulletText='▪'     # Square
bulletText='▸'     # Triangle
bulletText='→'     # Arrow
bulletText='1.'    # Numbers
bulletText='a)'    # Letters
```

## Text Measurement

```python
from reportlab.pdfbase.pdfmetrics import stringWidth

# Measure string width
width = stringWidth("Hello World", "Helvetica", 12)

# Check if text fits in available width
max_width = 200
if stringWidth(text, font_name, font_size) > max_width:
    # Text is too wide
    pass
```

## Best Practices

1. **Always use UTF-8** for text input
2. **Set leading > fontSize** for readability (typically fontSize + 2)
3. **Register font families** for proper bold/italic support
4. **Escape HTML** if displaying user content: use `<` for < and `>` for >
5. **Use getSampleStyleSheet()** as a starting point, don't create all styles from scratch
6. **Test Asian fonts** early if supporting multi-language content
7. **Set wordWrap='CJK'** for Chinese/Japanese/Korean text
8. **Use stringWidth()** to check if text fits before rendering
9. **Define styles once** at document start, reuse throughout
10. **Enable hyphenation** for justified text: `hyphenationLang='en_US'` (requires pyphen package)
