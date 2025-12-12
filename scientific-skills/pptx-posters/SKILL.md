---
name: latex-posters
description: "Create professional research posters in LaTeX using beamerposter, tikzposter, or baposter. Support for conference presentations, academic posters, and scientific communication. Includes layout design, color schemes, multi-column formats, figure integration, and poster-specific best practices for visual communication."
allowed-tools: [Read, Write, Edit, Bash]
---

# LaTeX Research Posters

## Overview

Research posters are a critical medium for scientific communication at conferences, symposia, and academic events. This skill provides comprehensive guidance for creating professional, visually appealing research posters using LaTeX packages. Generate publication-quality posters with proper layout, typography, color schemes, and visual hierarchy.

## When to Use This Skill

This skill should be used when:
- Creating research posters for conferences, symposia, or poster sessions
- Designing academic posters for university events or thesis defenses
- Preparing visual summaries of research for public engagement
- Converting scientific papers into poster format
- Creating template posters for research groups or departments
- Designing posters that comply with specific conference size requirements (A0, A1, 36×48", etc.)
- Building posters with complex multi-column layouts
- Integrating figures, tables, equations, and citations in poster format

## Visual Enhancement with Scientific Schematics

**⚠️ MANDATORY: Every research poster MUST include at least 2-3 AI-generated figures using the scientific-schematics skill.**

This is not optional. Posters are primarily visual media - text-heavy posters fail to communicate effectively. Before finalizing any poster:
1. Generate at minimum TWO schematics or diagrams
2. Target 3-4 figures for comprehensive posters (methodology flowchart, key results visualization, conceptual framework)
3. Figures should occupy 40-50% of poster area

**How to generate figures:**
- Use the **scientific-schematics** skill to generate AI-powered publication-quality diagrams
- Simply describe your desired diagram in natural language
- Nano Banana Pro will automatically generate, review, and refine the schematic

**How to generate schematics:**
```bash
python scripts/generate_schematic.py "your diagram description" -o figures/output.png
```

The AI will automatically:
- Create publication-quality images with proper formatting
- Review and refine through multiple iterations
- Ensure accessibility (colorblind-friendly, high contrast)
- Save outputs in the figures/ directory

**When to add schematics:**
- Research methodology flowcharts for poster content
- Conceptual framework diagrams
- Experimental design visualizations
- Data analysis pipeline diagrams
- System architecture diagrams
- Biological pathway illustrations
- Any complex concept that benefits from visualization

For detailed guidance on creating schematics, refer to the scientific-schematics skill documentation.

---

## Core Capabilities

### 1. LaTeX Poster Packages

Support for three major LaTeX poster packages, each with distinct advantages. For detailed comparison and package-specific guidance, refer to `references/latex_poster_packages.md`.

**beamerposter**:
- Extension of the Beamer presentation class
- Familiar syntax for Beamer users
- Excellent theme support and customization
- Best for: Traditional academic posters, institutional branding

**tikzposter**:
- Modern, flexible design with TikZ integration
- Built-in color themes and layout templates
- Extensive customization through TikZ commands
- Best for: Colorful, modern designs, custom graphics

**baposter**:
- Box-based layout system
- Automatic spacing and positioning
- Professional-looking default styles
- Best for: Multi-column layouts, consistent spacing

### 2. Poster Layout and Structure

Create effective poster layouts following visual communication principles. For comprehensive layout guidance, refer to `references/poster_layout_design.md`.

**Common Poster Sections**:
- **Header/Title**: Title, authors, affiliations, logos
- **Introduction/Background**: Research context and motivation
- **Methods/Approach**: Methodology and experimental design
- **Results**: Key findings with figures and data visualizations
- **Conclusions**: Main takeaways and implications
- **References**: Key citations (typically abbreviated)
- **Acknowledgments**: Funding, collaborators, institutions

**Layout Strategies**:
- **Column-based layouts**: 2-column, 3-column, or 4-column grids
- **Block-based layouts**: Flexible arrangement of content blocks
- **Z-pattern flow**: Guide readers through content logically
- **Visual hierarchy**: Use size, color, and spacing to emphasize key points

### 3. Design Principles for Research Posters

Apply evidence-based design principles for maximum impact. For detailed design guidance, refer to `references/poster_design_principles.md`.

**Typography**:
- Title: 72-120pt for visibility from distance
- Section headers: 48-72pt
- Body text: 24-36pt minimum for readability from 4-6 feet
- Use sans-serif fonts (Arial, Helvetica, Calibri) for clarity
- Limit to 2-3 font families maximum

**Color and Contrast**:
- Use high-contrast color schemes for readability
- Institutional color palettes for branding
- Color-blind friendly palettes (avoid red-green combinations)
- White space is active space—don't overcrowd

**Visual Elements**:
- High-resolution figures (300 DPI minimum for print)
- Large, clear labels on all figures
- Consistent figure styling throughout
- Strategic use of icons and graphics
- Balance text with visual content (40-50% visual recommended)

**Content Guidelines**:
- **Less is more**: 300-800 words total recommended
- Bullet points over paragraphs for scannability
- Clear, concise messaging
- Self-explanatory figures with minimal text explanation
- QR codes for supplementary materials or online resources

### 4. Standard Poster Sizes

Support for international and conference-specific poster dimensions:

**International Standards**:
- A0 (841 × 1189 mm / 33.1 × 46.8 inches) - Most common European standard
- A1 (594 × 841 mm / 23.4 × 33.1 inches) - Smaller format
- A2 (420 × 594 mm / 16.5 × 23.4 inches) - Compact posters

**North American Standards**:
- 36 × 48 inches (914 × 1219 mm) - Common US conference size
- 42 × 56 inches (1067 × 1422 mm) - Large format
- 48 × 72 inches (1219 × 1829 mm) - Extra large

**Orientation**:
- Portrait (vertical) - Most common, traditional
- Landscape (horizontal) - Better for wide content, timelines

### 5. Package-Specific Templates

Provide ready-to-use templates for each major package. Templates available in `assets/` directory.

**beamerposter Templates**:
- `beamerposter_classic.tex` - Traditional academic style
- `beamerposter_modern.tex` - Clean, minimal design
- `beamerposter_colorful.tex` - Vibrant theme with blocks

**tikzposter Templates**:
- `tikzposter_default.tex` - Standard tikzposter layout
- `tikzposter_rays.tex` - Modern design with ray theme
- `tikzposter_wave.tex` - Professional wave-style theme

**baposter Templates**:
- `baposter_portrait.tex` - Classic portrait layout
- `baposter_landscape.tex` - Landscape multi-column
- `baposter_minimal.tex` - Minimalist design

### 6. Figure and Image Integration

Optimize visual content for poster presentations:

**Best Practices**:
- Use vector graphics (PDF, SVG) when possible for scalability
- Raster images: minimum 300 DPI at final print size
- Consistent image styling (borders, captions, sizes)
- Group related figures together
- Use subfigures for comparisons

**LaTeX Figure Commands**:
```latex
% Include graphics package
\usepackage{graphicx}

% Simple figure
\includegraphics[width=0.8\linewidth]{figure.pdf}

% Figure with caption in tikzposter
\block{Results}{
  \begin{tikzfigure}
    \includegraphics[width=0.9\linewidth]{results.png}
  \end{tikzfigure}
}

% Multiple subfigures
\usepackage{subcaption}
\begin{figure}
  \begin{subfigure}{0.48\linewidth}
    \includegraphics[width=\linewidth]{fig1.pdf}
    \caption{Condition A}
  \end{subfigure}
  \begin{subfigure}{0.48\linewidth}
    \includegraphics[width=\linewidth]{fig2.pdf}
    \caption{Condition B}
  \end{subfigure}
\end{figure}
```

### 7. Color Schemes and Themes

Provide professional color palettes for various contexts:

**Academic Institution Colors**:
- Match university or department branding
- Use official color codes (RGB, CMYK, or LaTeX color definitions)

**Scientific Color Palettes** (color-blind friendly):
- Viridis: Professional gradient from purple to yellow
- ColorBrewer: Research-tested palettes for data visualization
- IBM Color Blind Safe: Accessible corporate palette

**Package-Specific Theme Selection**:

**beamerposter**:
```latex
\usetheme{Berlin}
\usecolortheme{beaver}
```

**tikzposter**:
```latex
\usetheme{Rays}
\usecolorstyle{Denmark}
```

**baposter**:
```latex
\begin{poster}{
  background=plain,
  bgColorOne=white,
  headerColorOne=blue!70,
  textborder=rounded
}
```

### 8. Typography and Text Formatting

Ensure readability and visual appeal:

**Font Selection**:
```latex
% Sans-serif fonts recommended for posters
\usepackage{helvet}      % Helvetica
\usepackage{avant}       % Avant Garde
\usepackage{sfmath}      % Sans-serif math fonts

% Set default to sans-serif
\renewcommand{\familydefault}{\sfdefault}
```

**Text Sizing**:
```latex
% Adjust text sizes for visibility
\setbeamerfont{title}{size=\VeryHuge}
\setbeamerfont{author}{size=\Large}
\setbeamerfont{institute}{size=\normalsize}
```

**Emphasis and Highlighting**:
- Use bold for key terms: `\textbf{important}`
- Color highlights sparingly: `\textcolor{blue}{highlight}`
- Boxes for critical information
- Avoid italics (harder to read from distance)

### 9. QR Codes and Interactive Elements

Enhance poster interactivity for modern conferences:

**QR Code Integration**:
```latex
\usepackage{qrcode}

% Link to paper, code repository, or supplementary materials
\qrcode[height=2cm]{https://github.com/username/project}

% QR code with caption
\begin{center}
  \qrcode[height=3cm]{https://doi.org/10.1234/paper}\\
  \small Scan for full paper
\end{center}
```

**Digital Enhancements**:
- Link to GitHub repositories for code
- Link to video presentations or demos
- Link to interactive web visualizations
- Link to supplementary data or appendices

### 10. Compilation and Output

Generate high-quality PDF output for printing or digital display:

**Compilation Commands**:
```bash
# Basic compilation
pdflatex poster.tex

# With bibliography
pdflatex poster.tex
bibtex poster
pdflatex poster.tex
pdflatex poster.tex

# For beamer-based posters
lualatex poster.tex  # Better font support
xelatex poster.tex   # Unicode and modern fonts
```

**Ensuring Full Page Coverage**:

Posters should use the entire page without excessive margins. Configure packages correctly:

**beamerposter - Full Page Setup**:
```latex
\documentclass[final,t]{beamer}
\usepackage[size=a0,scale=1.4,orientation=portrait]{beamerposter}

% Remove default beamer margins
\setbeamersize{text margin left=0mm, text margin right=0mm}

% Use geometry for precise control
\usepackage[margin=10mm]{geometry}  % 10mm margins all around

% Remove navigation symbols
\setbeamertemplate{navigation symbols}{}

% Remove footline and headline if not needed
\setbeamertemplate{footline}{}
\setbeamertemplate{headline}{}
```

**tikzposter - Full Page Setup**:
```latex
\documentclass[
  25pt,                      % Font scaling
  a0paper,                   % Paper size
  portrait,                  % Orientation
  margin=10mm,               % Outer margins (minimal)
  innermargin=15mm,          % Space inside blocks
  blockverticalspace=15mm,   % Space between blocks
  colspace=15mm,             % Space between columns
  subcolspace=8mm            % Space between subcolumns
]{tikzposter}

% This ensures content fills the page
```

**baposter - Full Page Setup**:
```latex
\documentclass[a0paper,portrait,fontscale=0.285]{baposter}

\begin{poster}{
  grid=false,
  columns=3,
  colspacing=1.5em,          % Space between columns
  eyecatcher=true,
  background=plain,
  bgColorOne=white,
  borderColor=blue!50,
  headerheight=0.12\textheight,  % 12% for header
  textborder=roundedleft,
  headerborder=closed,
  boxheaderheight=2em        % Consistent box header heights
}
% Content here
\end{poster}
```

**Common Issues and Fixes**:

**Problem**: Large white margins around poster
```latex
% Fix for beamerposter
\setbeamersize{text margin left=5mm, text margin right=5mm}

% Fix for tikzposter
\documentclass[..., margin=5mm, innermargin=10mm]{tikzposter}

% Fix for baposter - adjust in document class
\documentclass[a0paper, margin=5mm]{baposter}
```

**Problem**: Content doesn't fill vertical space
```latex
% Use \vfill between sections to distribute space
\block{Introduction}{...}
\vfill
\block{Methods}{...}
\vfill
\block{Results}{...}

% Or manually adjust block spacing
\vspace{1cm}  % Add space between specific blocks
```

**Problem**: Poster extends beyond page boundaries
```latex
% Check total width calculation
% For 3 columns with spacing:
% Total = 3×columnwidth + 2×colspace + 2×margins
% Ensure this equals \paperwidth

% Debug by adding visible page boundary
\usepackage{eso-pic}
\AddToShipoutPictureBG{
  \AtPageLowerLeft{
    \put(0,0){\framebox(\LenToUnit{\paperwidth},\LenToUnit{\paperheight}){}}
  }
}
```

**Print Preparation**:
- Generate PDF/X-1a for professional printing
- Embed all fonts
- Convert colors to CMYK if required
- Check resolution of all images (minimum 300 DPI)
- Add bleed area if required by printer (usually 3-5mm)
- Verify page size matches requirements exactly

**Digital Display**:
- RGB color space for screen display
- Optimize file size for email/web
- Test readability on different screens

### 11. PDF Review and Quality Control

**CRITICAL**: Always review the generated PDF before printing or presenting. Use this systematic checklist:

**Step 1: Page Size Verification**
```bash
# Check PDF dimensions (should match poster size exactly)
pdfinfo poster.pdf | grep "Page size"

# Expected outputs:
# A0: 2384 x 3370 points (841 x 1189 mm)
# 36x48": 2592 x 3456 points
# A1: 1684 x 2384 points (594 x 841 mm)
```

**Step 2: Visual Inspection Checklist**

Open PDF at 100% zoom and check:

**Layout and Spacing**:
- [ ] Content fills entire page (no large white margins)
- [ ] Consistent spacing between columns
- [ ] Consistent spacing between blocks/sections
- [ ] All elements aligned properly (use ruler tool)
- [ ] No overlapping text or figures
- [ ] White space evenly distributed

**Typography**:
- [ ] Title clearly visible and large (72pt+)
- [ ] Section headers readable (48-72pt)
- [ ] Body text readable at 100% zoom (24-36pt minimum)
- [ ] No text cutoff or running off edges
- [ ] Consistent font usage throughout
- [ ] All special characters render correctly (symbols, Greek letters)

**Visual Elements**:
- [ ] All figures display correctly
- [ ] No pixelated or blurry images
- [ ] Figure captions present and readable
- [ ] Colors render as expected (not washed out or too dark)
- [ ] Logos display clearly
- [ ] QR codes visible and scannable

**Content Completeness**:
- [ ] Title and authors complete
- [ ] All sections present (Intro, Methods, Results, Conclusions)
- [ ] References included
- [ ] Contact information visible
- [ ] Acknowledgments (if applicable)
- [ ] No placeholder text remaining (Lorem ipsum, TODO, etc.)

**Technical Quality**:
- [ ] No LaTeX compilation warnings in important areas
- [ ] All citations resolved (no [?] marks)
- [ ] All cross-references working
- [ ] Page boundaries correct (no content cut off)

**Step 3: Reduced-Scale Print Test**

**Essential Pre-Printing Test**:
```bash
# Create reduced-size test print (25% of final size)
# This simulates viewing full poster from ~8-10 feet

# For A0 poster, print on A4 paper (24.7% scale)
# For 36x48" poster, print on letter paper (~25% scale)
```

**Print Test Checklist**:
- [ ] Title readable from 6 feet away
- [ ] Section headers readable from 4 feet away
- [ ] Body text readable from 2 feet away
- [ ] Figures clear and understandable
- [ ] Colors printed accurately
- [ ] No obvious design flaws

**Step 4: Digital Quality Checks**

**Font Embedding Verification**:
```bash
# Check that all fonts are embedded (required for printing)
pdffonts poster.pdf

# All fonts should show "yes" in "emb" column
# If any show "no", recompile with:
pdflatex -dEmbedAllFonts=true poster.tex
```

**Image Resolution Check**:
```bash
# Extract image information
pdfimages -list poster.pdf

# Check that all images are at least 300 DPI
# Formula: DPI = pixels / (inches in poster)
# For A0 width (33.1"): 300 DPI = 9930 pixels minimum
```

**File Size Optimization**:
```bash
# For email/web, compress if needed (>50MB)
gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 \
   -dPDFSETTINGS=/printer -dNOPAUSE -dQUIET -dBATCH \
   -sOutputFile=poster_compressed.pdf poster.pdf

# For printing, keep original (no compression)
```

**Step 5: Accessibility Check**

**Color Contrast Verification**:
- [ ] Text-background contrast ratio ≥ 4.5:1 (WCAG AA)
- [ ] Important elements contrast ratio ≥ 7:1 (WCAG AAA)
- Test online: https://webaim.org/resources/contrastchecker/

**Color Blindness Simulation**:
- [ ] View PDF through color blindness simulator
- [ ] Information not lost with red-green simulation
- [ ] Use Coblis (color-blindness.com) or similar tool

**Step 6: Content Proofreading**

**Systematic Review**:
- [ ] Spell-check all text
- [ ] Verify all author names and affiliations
- [ ] Check all numbers and statistics for accuracy
- [ ] Confirm all citations are correct
- [ ] Review figure labels and captions
- [ ] Check for typos in headers and titles

**Peer Review**:
- [ ] Ask colleague to review poster
- [ ] 30-second test: Can they identify main message?
- [ ] 5-minute review: Do they understand conclusions?
- [ ] Note any confusing elements

**Step 7: Technical Validation**

**LaTeX Compilation Log Review**:
```bash
# Check for warnings in .log file
grep -i "warning\|error\|overfull\|underfull" poster.log

# Common issues to fix:
# - Overfull hbox: Text extending beyond margins
# - Underfull hbox: Excessive spacing
# - Missing references: Citations not resolved
# - Missing figures: Image files not found
```

**Fix Common Warnings**:
```latex
% Overfull hbox (text too wide)
\usepackage{microtype}  % Better spacing
\sloppy  % Allow slightly looser spacing
\hyphenation{long-word}  % Manual hyphenation

% Missing fonts
\usepackage[T1]{fontenc}  % Better font encoding

% Image not found
% Ensure paths are correct and files exist
\graphicspath{{./figures/}{./images/}}
```

**Step 8: Final Pre-Print Checklist**

**Before Sending to Printer**:
- [ ] PDF size exactly matches requirements (check with pdfinfo)
- [ ] All fonts embedded (check with pdffonts)
- [ ] Color mode correct (RGB for screen, CMYK for print if required)
- [ ] Bleed area added if required (usually 3-5mm)
- [ ] Crop marks visible if required
- [ ] Test print completed and reviewed
- [ ] File naming clear: [LastName]_[Conference]_Poster.pdf
- [ ] Backup copy saved

**Printing Specifications to Confirm**:
- [ ] Paper type (matte vs. glossy)
- [ ] Printing method (inkjet, large format, fabric)
- [ ] Color profile (provided to printer if required)
- [ ] Delivery deadline and shipping address
- [ ] Tube or flat packaging preference

**Digital Presentation Checklist**:
- [ ] PDF size optimized (<10MB for email)
- [ ] Tested on multiple PDF viewers (Adobe, Preview, etc.)
- [ ] Displays correctly on different screens
- [ ] QR codes tested and functional
- [ ] Alternative formats prepared (PNG for social media)

**Review Script** (Available in `scripts/review_poster.sh`):
```bash
#!/bin/bash
# Automated poster PDF review script

echo "Poster PDF Quality Check"
echo "======================="

# Check file exists
if [ ! -f "$1" ]; then
    echo "Error: File not found"
    exit 1
fi

echo "File: $1"
echo ""

# Check page size
echo "1. Page Dimensions:"
pdfinfo "$1" | grep "Page size"
echo ""

# Check fonts
echo "2. Font Embedding:"
pdffonts "$1" | head -20
echo ""

# Check file size
echo "3. File Size:"
ls -lh "$1" | awk '{print $5}'
echo ""

# Count pages (should be 1 for poster)
echo "4. Page Count:"
pdfinfo "$1" | grep "Pages"
echo ""

echo "Manual checks required:"
echo "- Visual inspection at 100% zoom"
echo "- Reduced-scale print test (25%)"
echo "- Color contrast verification"
echo "- Proofreading for typos"
```

**Common PDF Issues and Solutions**:

| Issue | Cause | Solution |
|-------|-------|----------|
| Large white margins | Incorrect margin settings | Reduce margin in documentclass |
| Content cut off | Exceeds page boundaries | Check total width/height calculations |
| Blurry images | Low resolution (<300 DPI) | Replace with higher resolution images |
| Missing fonts | Fonts not embedded | Compile with -dEmbedAllFonts=true |
| Wrong page size | Incorrect paper size setting | Verify documentclass paper size |
| Colors look wrong | RGB vs CMYK mismatch | Convert color space for print |
| File too large (>50MB) | Uncompressed images | Optimize images or compress PDF |
| QR codes don't work | Too small or low resolution | Minimum 2×2cm, high contrast |

### 11. Common Poster Content Patterns

Effective content organization for different research types:

**Experimental Research Poster**:
1. Title and authors
2. Introduction: Problem and hypothesis
3. Methods: Experimental design (with diagram)
4. Results: Key findings (2-4 main figures)
5. Conclusions: Main takeaways (3-5 bullet points)
6. Future work (optional)
7. References and acknowledgments

**Computational/Modeling Poster**:
1. Title and authors
2. Motivation: Problem statement
3. Approach: Algorithm or model (with flowchart)
4. Implementation: Technical details
5. Results: Performance metrics and comparisons
6. Applications: Use cases
7. Code availability (QR code to GitHub)
8. References

**Review/Survey Poster**:
1. Title and authors
2. Scope: Topic overview
3. Methods: Literature search strategy
4. Key findings: Main themes (organized by category)
5. Trends: Visualizations of publication patterns
6. Gaps: Identified research needs
7. Conclusions: Summary and implications
8. References

### 12. Accessibility and Inclusive Design

Design posters that are accessible to diverse audiences:

**Color Blindness Considerations**:
- Avoid red-green combinations (most common color blindness)
- Use patterns or shapes in addition to color
- Test with color-blindness simulators
- Provide high contrast (WCAG AA standard: 4.5:1 minimum)

**Visual Impairment Accommodations**:
- Large, clear fonts (minimum 24pt body text)
- High contrast text and background
- Clear visual hierarchy
- Avoid complex textures or patterns in backgrounds

**Language and Content**:
- Clear, concise language
- Define acronyms and jargon
- International audience considerations
- Consider multilingual QR code options for global conferences

### 13. Poster Presentation Best Practices

Guidance beyond LaTeX for effective poster sessions:

**Content Strategy**:
- Tell a story, don't just list facts
- Focus on 1-3 main messages
- Use visual abstract or graphical summary
- Leave room for conversation (don't over-explain)

**Physical Presentation Tips**:
- Bring printed handouts or business cards with QR code
- Prepare 30-second, 2-minute, and 5-minute verbal summaries
- Stand to the side, not blocking the poster
- Engage viewers with open-ended questions

**Digital Backups**:
- Save poster as PDF on mobile device
- Prepare digital version for email sharing
- Create social media-friendly image version
- Have backup printed copy or digital display option

## Workflow for Poster Creation

### Stage 1: Planning and Content Development

1. **Determine poster requirements**:
   - Conference size specifications (A0, 36×48", etc.)
   - Orientation (portrait vs. landscape)
   - Submission deadlines and format requirements

2. **Develop content outline**:
   - Identify 1-3 core messages
   - Select key figures (typically 3-6 main visuals)
   - Draft concise text for each section (bullet points preferred)
   - Aim for 300-800 words total

3. **Choose LaTeX package**:
   - beamerposter: If familiar with Beamer, need institutional themes
   - tikzposter: For modern, colorful designs with flexibility
   - baposter: For structured, professional multi-column layouts

### Stage 2: Design and Layout

1. **Select or create template**:
   - Start with provided templates in `assets/`
   - Customize color scheme to match branding
   - Configure page size and orientation

2. **Design layout structure**:
   - Plan column structure (2, 3, or 4 columns)
   - Map content flow (typically left-to-right, top-to-bottom)
   - Allocate space for title (10-15%), content (70-80%), footer (5-10%)

3. **Set typography**:
   - Configure font sizes for different hierarchy levels
   - Ensure minimum 24pt body text
   - Test readability from 4-6 feet distance

### Stage 3: Content Integration

1. **Create poster header**:
   - Title (concise, descriptive, 10-15 words)
   - Authors and affiliations
   - Institution logos (high-resolution)
   - Conference logo if required

2. **Populate content sections**:
   - Keep text minimal and scannable
   - Use bullet points, not paragraphs
   - Write in active voice
   - Integrate figures with clear captions

3. **Add visual elements**:
   - High-resolution figures (300 DPI minimum)
   - Consistent styling across all figures
   - Color-coded elements for emphasis
   - QR codes for supplementary materials

4. **Include references**:
   - Cite key papers only (5-10 references typical)
   - Use abbreviated citation style
   - Consider QR code to full bibliography

### Stage 4: Refinement and Testing

1. **Review and iterate**:
   - Check for typos and errors
   - Verify all figures are high resolution
   - Ensure consistent formatting
   - Confirm color scheme works well together

2. **Test readability**:
   - Print at 25% scale and read from 2-3 feet (simulates poster from 8-12 feet)
   - Check color on different monitors
   - Verify QR codes function correctly
   - Ask colleague to review

3. **Optimize for printing**:
   - Embed all fonts in PDF
   - Verify image resolution
   - Check PDF size requirements
   - Include bleed area if required

### Stage 5: Compilation and Delivery

1. **Compile final PDF**:
   ```bash
   pdflatex poster.tex
   # Or for better font support:
   lualatex poster.tex
   ```

2. **Verify output quality**:
   - Check all elements are visible and correctly positioned
   - Zoom to 100% and inspect figure quality
   - Verify colors match expectations
   - Confirm PDF opens correctly on different viewers

3. **Prepare for printing**:
   - Export as PDF/X-1a if required
   - Save backup copies
   - Get test print on regular paper first
   - Order professional printing 2-3 days before deadline

4. **Create supplementary materials**:
   - Save PNG/JPG version for social media
   - Create handout version (8.5×11" summary)
   - Prepare digital version for email sharing

## Integration with Other Skills

This skill works effectively with:
- **Scientific Writing**: For developing poster content from papers
- **Figure Creation**: For generating high-quality visualizations
- **Literature Review**: For contextualizing research
- **Data Analysis**: For creating result figures and charts

## Common Pitfalls to Avoid

**Design Mistakes**:
- ❌ Too much text (over 1000 words)
- ❌ Font sizes too small (under 24pt body text)
- ❌ Low-contrast color combinations
- ❌ Cluttered layout with no white space
- ❌ Inconsistent styling across sections
- ❌ Poor quality or pixelated images

**Content Mistakes**:
- ❌ No clear narrative or message
- ❌ Too many research questions or objectives
- ❌ Overuse of jargon without definitions
- ❌ Results without context or interpretation
- ❌ Missing author contact information

**Technical Mistakes**:
- ❌ Wrong poster dimensions for conference requirements
- ❌ RGB colors sent to CMYK printer (color shift)
- ❌ Fonts not embedded in PDF
- ❌ File size too large for submission portal
- ❌ QR codes too small or not tested

**Best Practices**:
- ✅ Follow conference size specifications exactly
- ✅ Test print at reduced scale before final printing
- ✅ Use high-contrast, accessible color schemes
- ✅ Keep text minimal and highly scannable
- ✅ Include clear contact information and QR codes
- ✅ Balance text and visuals (40-50% visual content)
- ✅ Proofread carefully (errors are magnified on posters!)

## Package Installation

Ensure required LaTeX packages are installed:

```bash
# For TeX Live (Linux/Mac)
tlmgr install beamerposter tikzposter baposter

# For MiKTeX (Windows)
# Packages typically auto-install on first use

# Additional recommended packages
tlmgr install qrcode graphics xcolor tcolorbox subcaption
```

## Scripts and Automation

Helper scripts available in `scripts/` directory:

- `compile_poster.sh`: Automated compilation with error handling
- `generate_template.py`: Interactive template generator
- `resize_images.py`: Batch image optimization for posters
- `poster_checklist.py`: Pre-submission validation tool

## References

Comprehensive reference files for detailed guidance:

- `references/latex_poster_packages.md`: Detailed comparison of beamerposter, tikzposter, and baposter with examples
- `references/poster_layout_design.md`: Layout principles, grid systems, and visual flow
- `references/poster_design_principles.md`: Typography, color theory, visual hierarchy, and accessibility
- `references/poster_content_guide.md`: Content organization, writing style, and section-specific guidance

## Templates

Ready-to-use poster templates in `assets/` directory:

- beamerposter templates (classic, modern, colorful)
- tikzposter templates (default, rays, wave, envelope)
- baposter templates (portrait, landscape, minimal)
- Example posters from various scientific disciplines
- Color scheme definitions and institutional templates

Load these templates and customize for your specific research and conference requirements.

