# Academic Presentation Styling Guide

This document provides styling guidelines for creating professional academic presentations.

## Color Schemes

### Professional Academic Palette
```css
:root {
  /* Primary colors */
  --primary-dark: #1a365d;      /* Navy blue - headers */
  --primary-medium: #2c5282;    /* Medium blue - accents */
  --primary-light: #4299e1;     /* Light blue - highlights */
  
  /* Secondary colors */
  --secondary-dark: #2d3748;    /* Dark gray - text */
  --secondary-medium: #4a5568;  /* Medium gray - subtext */
  --secondary-light: #cbd5e0;   /* Light gray - borders */
  
  /* Background colors */
  --bg-white: #ffffff;
  --bg-light: #f7fafc;
  --bg-accent: #ebf8ff;
  
  /* Accent colors for data */
  --accent-1: #48bb78;          /* Green */
  --accent-2: #ed8936;          /* Orange */
  --accent-3: #9f7aea;          /* Purple */
  --accent-4: #f56565;          /* Red */
}
```

### Alternative: Minimalist Palette
```css
:root {
  --primary: #000000;
  --secondary: #424242;
  --accent: #0066cc;
  --bg: #ffffff;
  --bg-light: #f5f5f5;
}
```

## Typography

### Font Families
```css
:root {
  /* REQUIRED: Use Arial for all text */
  --font-family: 'Arial', sans-serif;
}

/* Apply to all elements */
body, h1, h2, h3, h4, h5, h6, p, li, span, div {
  font-family: var(--font-family);
}
```

### Font Sizes (16:9 at 960x540px)
```css
:root {
  --text-xs: 14px;
  --text-sm: 16px;
  --text-base: 20px;
  --text-lg: 24px;
  --text-xl: 32px;
  --text-2xl: 40px;
  --text-3xl: 48px;
}
```

### Typography Classes
```css
h1 {
  font-size: var(--text-3xl);
  font-weight: 700;
  color: var(--primary-dark);
  line-height: 1.2;
  margin-bottom: 20px;
}

h2 {
  font-size: var(--text-2xl);
  font-weight: 600;
  color: var(--primary-medium);
  line-height: 1.3;
  margin-bottom: 16px;
}

h3 {
  font-size: var(--text-xl);
  font-weight: 600;
  color: var(--primary-medium);
  line-height: 1.4;
}

p, li {
  font-size: var(--text-base);
  line-height: 1.6;
  color: var(--secondary-dark);
}

.small-text {
  font-size: var(--text-sm);
  color: var(--secondary-medium);
}

.citation {
  font-size: var(--text-xs);
  color: var(--secondary-medium);
  font-style: italic;
}
```

## Layout Patterns

### Title Slide
```html
<div class="title-slide">
  <h1>Presentation Title</h1>
  <p class="subtitle">Subtitle or Topic Area</p>
  <div class="author-info">
    <p>Presenter Name</p>
    <p class="small-text">Institution/Lab</p>
    <p class="small-text">Date</p>
  </div>
</div>
```

### Content Slide Layouts

#### Single Column
```html
<div class="content-slide">
  <h2>Slide Title</h2>
  <div class="content">
    <!-- Main content -->
  </div>
  <div class="citation">Source: Author et al. (2024)</div>
</div>
```

#### Two Column
```html
<div class="content-slide">
  <h2>Slide Title</h2>
  <div class="row">
    <div class="col">
      <!-- Left column -->
    </div>
    <div class="col">
      <!-- Right column -->
    </div>
  </div>
  <div class="citation">Source: Author et al. (2024)</div>
</div>
```

#### Figure + Text
```html
<div class="content-slide">
  <h2>Slide Title</h2>
  <div class="row">
    <div class="col" style="flex: 1.5">
      <img src="figure.png" class="figure">
    </div>
    <div class="col" style="flex: 1">
      <ul>
        <li>Key point 1</li>
        <li>Key point 2</li>
        <li>Key point 3</li>
      </ul>
    </div>
  </div>
  <div class="citation">Source: Author et al. (2024), Figure 2</div>
</div>
```

## Spacing and Margins

```css
:root {
  --spacing-xs: 8px;
  --spacing-sm: 12px;
  --spacing-md: 20px;
  --spacing-lg: 32px;
  --spacing-xl: 48px;
  
  --slide-padding: 40px;
  --content-margin: 20px;
}

body {
  padding: var(--slide-padding);
}

.content {
  margin-top: var(--content-margin);
}
```

## Citation Display

### In-slide Citation (Reference Numbers)
```html
<!-- In-text citation with reference numbers -->
<p>Recent studies demonstrate enzyme cascade efficiency [1,2,5].</p>

<!-- Top-right corner for slide-level citation -->
<div class="citation" style="position: absolute; top: 10px; right: 40px;">
  [1,3]
</div>
```

### Figure Captions
```html
<!-- Figure with original caption and citation -->
<div class="figure-container">
  <img src="figure.png" class="figure">
  <p class="figure-caption">
    <strong>Figure 2B.</strong> Original caption text from the paper describing 
    the experimental results and methodology. (Figure 2B from [3])
  </p>
</div>
```

```css
.figure-caption {
  font-size: var(--text-sm);
  line-height: 1.4;
  color: var(--secondary-dark);
  margin-top: 8px;
  text-align: left;
}

.figure-caption strong {
  font-weight: 600;
}
```

### Schematic Diagrams
```html
<!-- Section summary schematic -->
<div class="schematic-container">
  <svg class="schematic" viewBox="0 0 800 400">
    <!-- SVG content for workflow/mechanism diagram -->
  </svg>
  <p class="schematic-caption">
    Schematic summary of enzymatic cascade mechanism integrating findings 
    from multiple studies [1,2,5].
  </p>
</div>
```

```css
.schematic-container {
  margin: 20px 0;
}

.schematic {
  width: 100%;
  max-height: 350px;
}

.schematic-caption {
  font-size: var(--text-sm);
  line-height: 1.4;
  color: var(--secondary-dark);
  margin-top: 12px;
  font-style: italic;
}
```

### References Slide
```css
.references-slide {
  font-size: var(--text-sm);
  line-height: 1.5;
}

.references-slide h2 {
  margin-bottom: 24px;
}

.reference-item {
  margin-bottom: 12px;
  padding-left: 24px;
  text-indent: -24px;
}

.reference-number {
  font-weight: 600;
  color: var(--primary-medium);
}
```

Example:
```html
<div class="references-slide">
  <h2>References</h2>
  <div class="reference-item">
    <span class="reference-number">1.</span> Kim J, et al. (2024). 
    Enzymatic cascade for rare sugar biosynthesis. Nature 123:456-789. 
    DOI: 10.1038/xxxxx
  </div>
  <div class="reference-item">
    <span class="reference-number">2.</span> Lee S, et al. (2023). 
    Optimization of multi-enzyme systems. Cell 456:123-456. 
    DOI: 10.1016/xxxxx
  </div>
</div>
```

## Special Elements

### Callout Boxes
```css
.callout {
  background: var(--bg-accent);
  border-left: 4px solid var(--primary-medium);
  padding: 16px 20px;
  margin: 20px 0;
}

.callout-warning {
  background: #fff5f5;
  border-left-color: var(--accent-4);
}
```

### Code Blocks
```css
.code {
  background: var(--bg-light);
  padding: 16px;
  font-family: var(--font-mono);
  font-size: var(--text-sm);
  border-radius: 4px;
  overflow-x: auto;
}
```

### Data Tables
```css
table {
  width: 100%;
  border-collapse: collapse;
  font-size: var(--text-sm);
}

th {
  background: var(--primary-dark);
  color: white;
  padding: 12px;
  text-align: left;
}

td {
  padding: 10px 12px;
  border-bottom: 1px solid var(--secondary-light);
}

tr:hover {
  background: var(--bg-light);
}
```

## Best Practices

### Presentation Structure
- Maximum 30 slides total
- Include section summary schematics
- Use original figure captions with citations
- Number references in order of first appearance

### Content Density
- Maximum 6-8 bullet points per slide
- Maximum 40-50 words per slide (excluding title and captions)
- Use images/figures to reduce text
- Include complete figure captions from original papers

### Visual Hierarchy
1. Title: Largest, bold
2. Section headers: Medium, semi-bold
3. Body text: Base size, regular weight
4. Citations: Smallest, italic

### Color Usage
- Use primary colors for structure (headers, borders)
- Use accent colors sparingly for emphasis
- Maintain high contrast for readability
- Avoid pure black (#000000), use dark gray instead

### Figure Guidelines
- Minimum 300 DPI for raster images
- Vector graphics (SVG) preferred when possible
- Include scale bars and labels
- Caption below figure with citation
