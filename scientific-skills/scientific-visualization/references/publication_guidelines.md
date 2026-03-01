# Publication-Ready Figure Guidelines

## Core Principles

Scientific figures must be clear, accurate, and accessible. Publication-ready figures follow these fundamental principles:

1. **Clarity**: Information should be immediately understandable
2. **Accuracy**: Data representation must be truthful and unmanipulated
3. **Accessibility**: Figures should be interpretable by all readers, including those with visual impairments
4. **Professional**: Clean, polished appearance suitable for peer-reviewed journals

## Resolution and File Format

### Resolution Requirements
- **Raster images (photos, microscopy)**: 300-600 DPI at final print size
- **Line art and graphs**: 600-1200 DPI (or vector format)
- **Combined figures**: 300-600 DPI

### Publisher DPI Requirements Comparison

| Image Type | Nature | ACS | Elsevier |
|---|---|---|---|
| **Line art** (graphs, diagrams) | 1200 DPI | 1200 DPI | 1000 DPI |
| **Halftone** (photographs, micrographs) | 300 DPI | 300 DPI | 300 DPI |
| **Combination** (line art + halftone) | 600 DPI | 600 DPI | 500 DPI |

> **Note**: These are *minimum* requirements. When in doubt, export at the highest DPI listed for your image type. Vector formats (PDF, EPS) bypass DPI constraints entirely for line art elements.

### File Formats

#### Vector Format Priority: PDF > EPS > SVG
- **PDF** (preferred): Universal compatibility, embeds fonts reliably, accepted by all major publishers
- **EPS**: Legacy standard, widely supported; some journals still prefer EPS over PDF
- **SVG**: Excellent for web/supplementary; limited support for print submissions (not accepted by Nature, ACS)
- Infinitely scalable without quality loss
- Smaller file sizes for line art
- Best for: plots, diagrams, schematics

#### Raster Format Priority: TIFF > PNG >> JPEG (avoid)
- **TIFF** (preferred): Lossless compression (LZW), universally accepted for print, supports CMYK
- **PNG**: Lossless, good for web and supplementary materials; some journals do not accept PNG for print
- Use for: photographs, microscopy, images with continuous tone

> **WARNING -- JPEG/JPG**: Lossy compression introduces irreversible artifacts, especially around sharp edges, text, and thin lines. JPEG is **unacceptable** for line art, graphs, and any figure containing text. If you must use JPEG (e.g., very large photographic images where TIFF file size is prohibitive), use **maximum quality (95-100%)** and document the reason. Preferred alternatives:
> - Convert to **TIFF** (LZW compression) for equivalent file size with no quality loss
> - Use **PNG** for web/digital submissions
> - Export line art elements as **PDF/EPS** vectors and composite photographics separately

### Size Specifications
- **Single column**: 85-90 mm (3.35-3.54 inches) width
- **1.5 column**: 114-120 mm (4.49-4.72 inches) width
- **Double column**: 174-180 mm (6.85-7.08 inches) width
- **Maximum height**: Usually 230-240 mm (9-9.5 inches)

### Publisher Figure Size Comparison

| Dimension | Nature | ACS | Elsevier |
|---|---|---|---|
| **Single column width** | 89 mm (3.5 in) | 84 mm (3.3 in) | 90 mm (3.5 in) |
| **1.5 column width** | 120 mm (4.7 in) | 110 mm (4.3 in) | 140 mm (5.5 in) |
| **Double column width** | 183 mm (7.2 in) | 175 mm (6.9 in) | 190 mm (7.5 in) |
| **Maximum height** | 247 mm (9.7 in) | 234 mm (9.2 in) | 240 mm (9.4 in) |
| **Minimum width** | 30 mm (1.2 in) | 42 mm (1.65 in) | 30 mm (1.2 in) |

> **Practical tip**: Design at the single-column width (89 mm) for maximum compatibility. Scale up only when detail is lost at the smaller size.

## Typography

### Font Guidelines
- **Font family**: Sans-serif fonts (Arial, Helvetica, Calibri) for most journals
  - Some journals prefer specific fonts (check guidelines)
  - Consistency across all figures in manuscript

- **Font sizes at final print size**:
  - Axis labels: 7-9 pt minimum
  - Tick labels: 6-8 pt minimum
  - Legends: 6-8 pt
  - Panel labels (A, B, C): 8-12 pt, bold
  - Title: Generally avoided in multi-panel figures

- **Font weight**: Regular weight for most text; bold for panel labels only

### Text Embedding Rules

When exporting figures as PDF or EPS, font embedding is critical for ensuring text renders identically across systems. Incorrect font embedding is one of the most common reasons for figure rejection.

- **`pdf.fonttype: 42` (TrueType)**: Always set `matplotlib.rcParams['pdf.fonttype'] = 42` to embed fonts as TrueType outlines. This ensures text is searchable, selectable, and renders correctly on all systems.
  ```python
  import matplotlib
  matplotlib.rcParams['pdf.fonttype'] = 42  # TrueType
  matplotlib.rcParams['ps.fonttype'] = 42   # TrueType for EPS
  ```
- **Type 3 fonts are prohibited**: Nature and ACS journals explicitly reject figures containing Type 3 (bitmap) fonts. Type 3 fonts render poorly at different zoom levels and are not searchable. Matplotlib defaults to Type 3 unless overridden with the setting above.
- **Do NOT convert text to outlines (paths)**: While outlining text avoids font embedding issues, it makes text non-searchable and non-selectable, which violates accessibility requirements of most publishers. Nature, ACS, and Elsevier all require embedded (not outlined) text.
- **Font subsetting**: Most export tools subset fonts automatically. Ensure the full character set used in your figure is included in the embedded subset.

### Text Best Practices
- Use sentence case for axis labels ("Time (hours)" not "TIME (HOURS)")
- Include units in parentheses
- Avoid abbreviations unless space-constrained (define in caption)
- No text smaller than 5-6 pt at final size

## Color Usage

### Color Selection Principles
1. **Colorblind-friendly**: ~8% of males have color vision deficiency
   - Avoid red/green combinations
   - Use blue/orange, blue/yellow, or add texture/pattern
   - Test with colorblindness simulators

2. **Purposeful color**: Color should convey meaning, not just aesthetics
   - Use color to distinguish categories or highlight key data
   - Maintain consistency across figures (same treatment = same color)

3. **Print considerations**:
   - Colors may appear different in print vs. screen
   - Use CMYK color space for print, RGB for digital
   - Ensure sufficient contrast (especially for grayscale conversion)

### Recommended Color Palettes
- **Qualitative (categories)**: ColorBrewer, Okabe-Ito palette
- **Sequential (low to high)**: Viridis, Cividis, Blues, Oranges
- **Diverging (negative to positive)**: RdBu, PuOr, BrBG (ensure colorblind-safe)

### Grayscale Compatibility
- All figures should be interpretable in grayscale
- Use different line styles (solid, dashed, dotted) and markers
- Add patterns/hatching to bars and areas

## Layout and Composition

### Multi-Panel Figures
- **Panel labels**: Use bold uppercase letters (A, B, C) in top-left corner
- **Spacing**: Adequate white space between panels
- **Alignment**: Align panels along edges or axes where possible
- **Sizing**: Related panels should have consistent sizes
- **Arrangement**: Logical flow (left-to-right, top-to-bottom)

### Plot Elements

#### Axes
- **Axis lines**: 0.5-1 pt thickness
- **Tick marks**: Point inward or outward consistently
- **Tick frequency**: Enough to read values, not cluttered (typically 4-7 major ticks)
- **Axis labels**: Required on all plots; state units
- **Axis ranges**: Start from zero for bar charts (unless scientifically inappropriate)

#### Lines and Markers
- **Line width**: 1-2 pt for data lines; 0.5-1 pt for reference lines
- **Marker size**: 3-6 pt, larger than line width
- **Marker types**: Differentiate when multiple series (circles, squares, triangles)
- **Error bars**: 0.5-1 pt width; include caps if appropriate

#### Legends
- **Position**: Inside plot area if space permits, outside otherwise
- **Frame**: Optional; if used, thin line (0.5 pt)
- **Order**: Match order of data appearance (top to bottom or left to right)
- **Content**: Concise descriptions; full details in caption

### White Space and Margins
- Remove unnecessary white space around plots
- Maintain consistent margins
- `tight_layout()` or `constrained_layout=True` in matplotlib

## Data Representation Best Practices

### Statistical Rigor
- **Error bars**: Always show uncertainty (SD, SEM, CI) and state which in caption
- **Sample size**: Indicate n in figure or caption
- **Significance**: Mark statistical significance clearly (*, **, ***)
- **Replicates**: Show individual data points when possible, not just summary statistics

### Appropriate Chart Types
- **Bar plots**: Comparing discrete categories; always start y-axis at zero
- **Line plots**: Time series or continuous relationships
- **Scatter plots**: Correlation between variables; add regression line if appropriate
- **Box plots**: Distribution comparisons; show outliers
- **Heatmaps**: Matrix data, correlations, expression patterns
- **Violin plots**: Distribution shape comparison (better than box plots for bimodal data)

### Avoiding Distortion
- **No 3D effects**: Distorts perception of values
- **No unnecessary decorations**: No gradients, shadows, or chart junk
- **Consistent scales**: Use same scale for comparable panels
- **No truncated axes**: Unless clearly indicated and scientifically justified
- **Linear vs. log scales**: Choose appropriate scale; always label clearly

## Accessibility

### Colorblind Considerations
- Test with online simulators (e.g., Coblis, Color Oracle)
- Use patterns/textures in addition to color
- Provide alternative representations in supplementary materials if needed

### Visual Impairment
- High contrast between elements
- Thick enough lines (minimum 0.5 pt)
- Clear, uncluttered layouts

### Data Availability
- Include data tables in supplementary materials
- Provide source data files for graphs
- Consider interactive figures for online supplementary materials

## Ten Simple Rules for Better Figures

Adapted from Rougier, Droettboom & Bourne (2014), *PLOS Computational Biology* 10(9): e1003833. These rules provide a practical framework for creating effective scientific figures.

### Rule 1: Know Your Audience
- Adapt complexity and terminology to the target readership (specialist journal vs. general audience)
- **Implementation**: Before designing, ask: "Will a non-specialist in my subfield understand this figure without reading the caption?"

### Rule 2: Identify Your Message
- Each figure should convey one clear, specific message
- **Implementation**: Write a single sentence describing the figure's takeaway *before* creating it. If you need more than one sentence, consider splitting into multiple panels or figures.

### Rule 3: Adapt the Figure to the Support Medium
- Design differently for print (static, high-res, CMYK) vs. screen (interactive, RGB) vs. presentation (large fonts, high contrast)
- **Implementation**: Set DPI, color space, and font sizes based on the final medium. A figure designed for a journal may need redesigning for a poster or slide.

### Rule 4: Captions Are Not Optional
- The caption should make the figure self-contained; readers should understand it without reading the full text
- **Implementation**: Include: what is shown, how to read it (axes, colors, symbols), key statistical details (n, p-values, error bar type), and the main conclusion.

### Rule 5: Do Not Trust the Defaults
- Default settings in most plotting libraries (matplotlib, R/ggplot2) are not optimized for publication
- **Implementation**: Create a custom style sheet or `rcParams` configuration that enforces publication standards (font sizes, line widths, DPI, color palette). Apply it at the start of every project.

### Rule 6: Use Color Effectively
- Color should encode data, not decorate; must be accessible to colorblind readers
- **Implementation**: Use perceptually uniform colormaps (viridis, cividis). Test every figure with a colorblindness simulator before submission.

### Rule 7: Do Not Mislead the Reader
- Avoid visual distortions: truncated axes, dual y-axes with mismatched scales, area-encoding errors, 3D effects
- **Implementation**: Start bar charts at zero. Use consistent scales across comparable panels. Never use 3D perspective for 2D data.

### Rule 8: Avoid "Chartjunk"
- Remove every element that does not contribute to understanding: background colors, grid lines, 3D effects, unnecessary borders
- **Implementation**: Apply the "ink-to-data ratio" principle -- maximize the proportion of ink used to display data. Use `ax.spines` and `ax.tick_params` to strip unnecessary decoration.

### Rule 9: Message Trumps Beauty
- Prioritize clarity of communication over aesthetic appeal; a clear figure is always better than a beautiful but confusing one
- **Implementation**: After creating a figure, show it to a colleague without context. If they cannot state the main message within 10 seconds, simplify.

### Rule 10: Get the Right Tool
- Choose the right software for the task; no single tool is best for everything
- **Implementation**: Use matplotlib/seaborn for statistical plots, Inkscape/Illustrator for post-processing and annotations, ImageJ/FIJI for microscopy, and BioRender for biological schematics. Avoid doing everything in one tool when a combination yields better results.

## Common Mistakes to Avoid

1. **Font too small**: Text unreadable at final print size
2. **Low resolution**: Pixelated or blurry images
3. **Chart junk**: Unnecessary grid lines, 3D effects, decorations
4. **Poor color choices**: Red/green combinations, low contrast
5. **Missing elements**: No axis labels, no units, no error bars
6. **Inconsistent styling**: Different fonts/sizes within figure or between figures
7. **Data distortion**: Truncated axes, inappropriate scales, 3D effects
8. **JPEG compression**: Artifacts around text and lines
9. **Too much information**: Cramming too many data series into one plot
10. **Inaccessible legends**: Legends outside the figure boundary after export

## Figure Checklist

Before submission, verify:

- [ ] Resolution meets journal requirements (300+ DPI for raster)
- [ ] File format is acceptable (vector for plots, TIFF/PNG for images)
- [ ] Figure dimensions match journal specifications
- [ ] All text is readable at final size (minimum 6-7 pt)
- [ ] Fonts are consistent and embedded (for PDF/EPS)
- [ ] Colors are colorblind-friendly
- [ ] Figure is interpretable in grayscale
- [ ] All axes are labeled with units
- [ ] Error bars or uncertainty indicators are present
- [ ] Statistical significance is marked if applicable
- [ ] Panel labels are present and consistent (A, B, C)
- [ ] Legend is clear and complete
- [ ] No chart junk or unnecessary elements
- [ ] File naming follows journal conventions
- [ ] Figure caption is comprehensive
- [ ] Source data is available

## Automated Validation Checklist

Use this checklist programmatically before saving/exporting any publication figure. These items can be verified automatically in code (e.g., via a `validate_figure()` function).

### Resolution and Dimensions
- [ ] **DPI check**: Raster output DPI >= 300 (halftone) / >= 600 (combination) / >= 1200 (line art)
- [ ] **Figure width**: Within target journal's column width (single: 85-90 mm, double: 174-190 mm)
- [ ] **Figure height**: Does not exceed maximum page height (230-247 mm depending on journal)
- [ ] **Aspect ratio**: Reasonable aspect ratio (not excessively wide or tall)

### File Format and Size
- [ ] **Format**: Vector (PDF/EPS) for plots; TIFF/PNG for photographic content
- [ ] **No JPEG**: Output is not JPEG/JPG (or if JPEG is unavoidable, quality >= 95%)
- [ ] **File size**: Under journal limit (typically 10-50 MB per figure; check specific journal)
- [ ] **Color space**: CMYK for print, RGB for digital-only (verify `rcParams` or export settings)

### Font and Text
- [ ] **Font type embedded**: `pdf.fonttype = 42` (TrueType); no Type 3 fonts present
- [ ] **Minimum font size**: All text >= 6 pt at final print dimensions
- [ ] **Maximum font size**: No text excessively large (panel labels typically <= 12 pt)
- [ ] **Font consistency**: Single font family used throughout the figure
- [ ] **Text is not outlined**: Fonts are embedded as text, not converted to paths

### Visual Elements
- [ ] **Line width minimum**: All lines >= 0.5 pt at final size
- [ ] **Colorblind safe**: Palette passes deuteranopia/protanopia simulation check
- [ ] **Grayscale legible**: Figure remains interpretable when converted to grayscale
- [ ] **No transparency issues**: If transparent elements are used, they render correctly when flattened

### Data Integrity
- [ ] **Axis labels present**: All axes have labels with units
- [ ] **Error bars / uncertainty**: Statistical uncertainty is shown where applicable
- [ ] **Panel labels**: Multi-panel figures have A, B, C labels
- [ ] **Legend present**: All data series are identified in the legend or caption
- [ ] **No clipping**: Data points, labels, and legends are not clipped by figure boundaries

### Implementation Example (matplotlib)
```python
def validate_figure(fig, target_dpi=300, max_width_mm=183, max_height_mm=247):
    """Validate a matplotlib figure against publication standards."""
    issues = []

    # Check DPI
    if fig.dpi < target_dpi:
        issues.append(f"DPI ({fig.dpi}) below target ({target_dpi})")

    # Check dimensions
    w_in, h_in = fig.get_size_inches()
    w_mm, h_mm = w_in * 25.4, h_in * 25.4
    if w_mm > max_width_mm:
        issues.append(f"Width ({w_mm:.0f} mm) exceeds max ({max_width_mm} mm)")
    if h_mm > max_height_mm:
        issues.append(f"Height ({h_mm:.0f} mm) exceeds max ({max_height_mm} mm)")

    # Check font sizes
    for ax in fig.get_axes():
        for text in ax.get_xticklabels() + ax.get_yticklabels():
            if text.get_fontsize() < 6:
                issues.append(f"Tick label font size ({text.get_fontsize()} pt) below 6 pt")
                break
        if ax.xaxis.label.get_fontsize() < 6:
            issues.append("X-axis label font size below 6 pt")
        if ax.yaxis.label.get_fontsize() < 6:
            issues.append("Y-axis label font size below 6 pt")

    # Check rcParams for font embedding
    import matplotlib
    if matplotlib.rcParams.get('pdf.fonttype', 3) != 42:
        issues.append("pdf.fonttype is not 42 (TrueType) -- Type 3 fonts will be generated")
    if matplotlib.rcParams.get('ps.fonttype', 3) != 42:
        issues.append("ps.fonttype is not 42 (TrueType) -- Type 3 fonts in EPS output")

    if issues:
        print("VALIDATION WARNINGS:")
        for i, issue in enumerate(issues, 1):
            print(f"  {i}. {issue}")
    else:
        print("All checks passed.")

    return len(issues) == 0
```

## Journal-Specific Considerations

Always consult the specific journal's author guidelines. Common variations include:

- **Nature journals**: RGB, 300 DPI minimum, specific size requirements
- **Science**: EPS or high-res TIFF, specific font requirements
- **Cell Press**: PDF or EPS preferred, Arial or Helvetica fonts
- **PLOS**: TIFF or EPS, specific color space requirements
- **ACS journals**: Application files (AI, EPS) or high-res TIFF

See `journal_requirements.md` for detailed specifications from major publishers.
