# PDF Features Reference

Advanced PDF capabilities: links, bookmarks, forms, encryption, and metadata.

## Document Metadata

Set PDF document properties viewable in PDF readers.

```python
from reportlab.pdfgen import canvas

c = canvas.Canvas("output.pdf")

# Set metadata
c.setAuthor("John Doe")
c.setTitle("Annual Report 2024")
c.setSubject("Financial Analysis")
c.setKeywords("finance, annual, report, 2024")
c.setCreator("MyApp v1.0")

# ... draw content ...

c.save()
```

With Platypus:

```python
from reportlab.platypus import SimpleDocTemplate

doc = SimpleDocTemplate(
    "output.pdf",
    title="Annual Report 2024",
    author="John Doe",
    subject="Financial Analysis",
)

doc.build(story)
```

## Bookmarks and Destinations

Create internal navigation structure.

### Simple Bookmarks

```python
from reportlab.pdfgen import canvas

c = canvas.Canvas("output.pdf")

# Create bookmark for current page
c.bookmarkPage("intro")  # Internal key
c.addOutlineEntry("Introduction", "intro", level=0)

c.showPage()

# Another bookmark
c.bookmarkPage("chapter1")
c.addOutlineEntry("Chapter 1", "chapter1", level=0)

# Sub-sections
c.bookmarkPage("section1_1")
c.addOutlineEntry("Section 1.1", "section1_1", level=1)  # Nested

c.save()
```

### Bookmark Levels

```python
# Create hierarchical outline
c.bookmarkPage("ch1")
c.addOutlineEntry("Chapter 1", "ch1", level=0)

c.bookmarkPage("ch1_s1")
c.addOutlineEntry("Section 1.1", "ch1_s1", level=1)

c.bookmarkPage("ch1_s1_1")
c.addOutlineEntry("Subsection 1.1.1", "ch1_s1_1", level=2)

c.bookmarkPage("ch2")
c.addOutlineEntry("Chapter 2", "ch2", level=0)
```

### Destination Fit Modes

Control how the page displays when navigating:

```python
# bookmarkPage with fit mode
c.bookmarkPage(
    key="chapter1",
    fit="Fit"  # Fit entire page in window
)

# Or use bookmarkHorizontalAbsolute
c.bookmarkHorizontalAbsolute(key="section", top=500)

# Available fit modes:
# "Fit" - Fit whole page
# "FitH" - Fit horizontally
# "FitV" - Fit vertically
# "FitR" - Fit rectangle
# "XYZ" - Specific position and zoom
```

## Hyperlinks

### External Links

```python
from reportlab.pdfgen import canvas
from reportlab.lib.units import inch

c = canvas.Canvas("output.pdf")

# Draw link rectangle
c.linkURL(
    "https://www.example.com",
    rect=(1*inch, 5*inch, 3*inch, 5.5*inch),  # (x1, y1, x2, y2)
    relative=0,  # 0 for absolute positioning
    thickness=1,
    color=(0, 0, 1),  # Blue
    dashArray=None
)

# Draw text over link area
c.setFillColorRGB(0, 0, 1)  # Blue text
c.drawString(1*inch, 5.2*inch, "Click here to visit example.com")

c.save()
```

### Internal Links

Link to bookmarked locations within the document:

```python
# Create destination
c.bookmarkPage("target_section")

# Later, create link to that destination
c.linkRect(
    "Link Text",
    "target_section",  # Bookmark key
    rect=(1*inch, 3*inch, 2*inch, 3.2*inch),
    relative=0
)
```

### Links in Paragraphs

For Platypus documents:

```python
from reportlab.platypus import Paragraph

# External link
text = '<link href="https://example.com" color="blue">Visit our website</link>'
para = Paragraph(text, style)

# Internal link (to anchor)
text = '<link href="#section1" color="blue">Go to Section 1</link>'
para1 = Paragraph(text, style)

# Create anchor
text = '<a name="section1"/>Section 1 Heading'
para2 = Paragraph(text, heading_style)

story.append(para1)
story.append(para2)
```

## Interactive Forms

Create fillable PDF forms.

### Text Fields

```python
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfform
from reportlab.lib.colors import black, white

c = canvas.Canvas("form.pdf")

# Create text field
c.acroForm.textfield(
    name="name",
    tooltip="Enter your name",
    x=100,
    y=700,
    width=200,
    height=20,
    borderColor=black,
    fillColor=white,
    textColor=black,
    forceBorder=True,
    fontSize=12,
    maxlen=100,  # Maximum character length
)

# Label
c.drawString(100, 725, "Name:")

c.save()
```

### Checkboxes

```python
# Create checkbox
c.acroForm.checkbox(
    name="agree",
    tooltip="I agree to terms",
    x=100,
    y=650,
    size=20,
    buttonStyle='check',  # 'check', 'circle', 'cross', 'diamond', 'square', 'star'
    borderColor=black,
    fillColor=white,
    textColor=black,
    forceBorder=True,
    checked=False,  # Initial state
)

c.drawString(130, 655, "I agree to the terms and conditions")
```

### Radio Buttons

```python
# Radio button group - only one can be selected
c.acroForm.radio(
    name="payment",  # Same name for group
    tooltip="Credit Card",
    value="credit",  # Value when selected
    x=100,
    y=600,
    size=15,
    selected=False,
)
c.drawString(125, 603, "Credit Card")

c.acroForm.radio(
    name="payment",  # Same name
    tooltip="PayPal",
    value="paypal",
    x=100,
    y=580,
    size=15,
    selected=False,
)
c.drawString(125, 583, "PayPal")
```

### List Boxes

```python
# Listbox with multiple options
c.acroForm.listbox(
    name="country",
    tooltip="Select your country",
    value="US",  # Default selected
    x=100,
    y=500,
    width=150,
    height=80,
    borderColor=black,
    fillColor=white,
    textColor=black,
    forceBorder=True,
    options=[
        ("United States", "US"),
        ("Canada", "CA"),
        ("Mexico", "MX"),
        ("Other", "OTHER"),
    ],  # List of (label, value) tuples
    multiple=False,  # Allow multiple selections
)
```

### Choice (Dropdown)

```python
# Dropdown menu
c.acroForm.choice(
    name="state",
    tooltip="Select state",
    value="CA",
    x=100,
    y=450,
    width=150,
    height=20,
    borderColor=black,
    fillColor=white,
    textColor=black,
    forceBorder=True,
    options=[
        ("California", "CA"),
        ("New York", "NY"),
        ("Texas", "TX"),
    ],
)
```

### Complete Form Example

```python
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter
from reportlab.lib.colors import black, white, lightgrey
from reportlab.lib.units import inch

def create_registration_form(filename):
    c = canvas.Canvas(filename, pagesize=letter)
    c.setFont("Helvetica-Bold", 16)
    c.drawString(inch, 10*inch, "Registration Form")

    y = 9*inch
    c.setFont("Helvetica", 12)

    # Name field
    c.drawString(inch, y, "Full Name:")
    c.acroForm.textfield(
        name="fullname",
        x=2*inch, y=y-5, width=4*inch, height=20,
        borderColor=black, fillColor=lightgrey, forceBorder=True
    )

    # Email field
    y -= 0.5*inch
    c.drawString(inch, y, "Email:")
    c.acroForm.textfield(
        name="email",
        x=2*inch, y=y-5, width=4*inch, height=20,
        borderColor=black, fillColor=lightgrey, forceBorder=True
    )

    # Age dropdown
    y -= 0.5*inch
    c.drawString(inch, y, "Age Group:")
    c.acroForm.choice(
        name="age_group",
        x=2*inch, y=y-5, width=2*inch, height=20,
        borderColor=black, fillColor=lightgrey, forceBorder=True,
        options=[("18-25", "18-25"), ("26-35", "26-35"),
                ("36-50", "36-50"), ("51+", "51+")]
    )

    # Newsletter checkbox
    y -= 0.5*inch
    c.acroForm.checkbox(
        name="newsletter",
        x=inch, y=y-5, size=15,
        buttonStyle='check', borderColor=black, forceBorder=True
    )
    c.drawString(inch + 25, y, "Subscribe to newsletter")

    c.save()

create_registration_form("registration.pdf")
```

## Encryption and Security

Protect PDFs with passwords and permissions.

### Basic Encryption

```python
from reportlab.pdfgen import canvas

c = canvas.Canvas("secure.pdf")

# Encrypt with user password
c.encrypt(
    userPassword="user123",    # Password to open
    ownerPassword="owner456",  # Password to change permissions
    canPrint=1,                # Allow printing
    canModify=0,               # Disallow modifications
    canCopy=1,                 # Allow text copying
    canAnnotate=0,             # Disallow annotations
    strength=128,              # 40 or 128 bit encryption
)

# ... draw content ...

c.save()
```

### Permission Settings

```python
c.encrypt(
    userPassword="user123",
    ownerPassword="owner456",
    canPrint=1,        # 1 = allow, 0 = deny
    canModify=0,       # Prevent content modification
    canCopy=1,         # Allow text/graphics copying
    canAnnotate=0,     # Prevent comments/annotations
    strength=128,      # Use 128-bit encryption
)
```

### Advanced Encryption

```python
from reportlab.lib.pdfencrypt import StandardEncryption

# Create encryption object
encrypt = StandardEncryption(
    userPassword="user123",
    ownerPassword="owner456",
    canPrint=1,
    canModify=0,
    canCopy=1,
    canAnnotate=1,
    strength=128,
)

# Use with canvas
c = canvas.Canvas("secure.pdf")
c._doc.encrypt = encrypt

# ... draw content ...

c.save()
```

### Platypus with Encryption

```python
from reportlab.platypus import SimpleDocTemplate

doc = SimpleDocTemplate("secure.pdf")

# Set encryption
doc.encrypt = True
doc.canPrint = 1
doc.canModify = 0

# Or use encrypt() method
doc.encrypt = encrypt_object

doc.build(story)
```

## Page Transitions

Add visual effects for presentations.

```python
from reportlab.pdfgen import canvas

c = canvas.Canvas("presentation.pdf")

# Set transition for current page
c.setPageTransition(
    effectname="Wipe",  # Transition effect
    duration=1,         # Duration in seconds
    direction=0         # Direction (effect-specific)
)

# Available effects:
# "Split", "Blinds", "Box", "Wipe", "Dissolve",
# "Glitter", "R" (Replace), "Fly", "Push", "Cover",
# "Uncover", "Fade"

# Direction values (effect-dependent):
# 0, 90, 180, 270 for most directional effects

# Example: Slide with fade transition
c.setFont("Helvetica-Bold", 24)
c.drawString(100, 400, "Slide 1")
c.setPageTransition("Fade", 0.5)
c.showPage()

c.drawString(100, 400, "Slide 2")
c.setPageTransition("Wipe", 1, 90)
c.showPage()

c.save()
```

## PDF/A Compliance

Create archival-quality PDFs.

```python
from reportlab.pdfgen import canvas

c = canvas.Canvas("pdfa.pdf")

# Enable PDF/A-1b compliance
c.setPageCompression(0)  # PDF/A requires uncompressed
# Note: Full PDF/A requires additional XMP metadata
# This is simplified - full compliance needs more setup

# ... draw content ...

c.save()
```

## Compression

Control file size vs generation speed.

```python
# Enable page compression
c = canvas.Canvas("output.pdf", pageCompression=1)

# Compression reduces file size but slows generation
# 0 = no compression (faster, larger files)
# 1 = compression (slower, smaller files)
```

## Forms and XObjects

Reusable graphics elements.

```python
from reportlab.pdfgen import canvas

c = canvas.Canvas("output.pdf")

# Begin form (reusable object)
c.beginForm("logo")
c.setFillColorRGB(0, 0, 1)
c.rect(0, 0, 100, 50, fill=1)
c.setFillColorRGB(1, 1, 1)
c.drawString(10, 20, "LOGO")
c.endForm()

# Use form multiple times
c.doForm("logo")  # At current position
c.translate(200, 0)
c.doForm("logo")  # At translated position
c.translate(200, 0)
c.doForm("logo")

c.save()

# Benefits: Smaller file size, faster rendering
```

## Best Practices

1. **Always set metadata** for professional documents
2. **Use bookmarks** for documents > 10 pages
3. **Make links visually distinct** (blue, underlined)
4. **Test forms** in multiple PDF readers (behavior varies)
5. **Use strong encryption (128-bit)** for sensitive data
6. **Set both user and owner passwords** for full security
7. **Enable printing** unless specifically restricted
8. **Test page transitions** - some readers don't support all effects
9. **Use meaningful bookmark titles** for navigation
10. **Consider PDF/A** for long-term archival needs
11. **Validate form field names** - must be unique and valid identifiers
12. **Add tooltips** to form fields for better UX
