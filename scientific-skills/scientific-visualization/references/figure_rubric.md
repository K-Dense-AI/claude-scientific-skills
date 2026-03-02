# Figure Quality Rubric

## Overview

A quantitative scoring rubric for evaluating scientific figures before submission. Based on the ScientificFigures framework (Rokhsareh et al.) and adapted for publication, poster, and presentation contexts. Each figure is scored on six criteria (0-2 points each) for a maximum of 12 points.

## Six Evaluation Criteria

### 1. Scale & Resolution (0-2 points)

Does the figure have appropriate dimensions and resolution for the target medium?

| Score | Criterion |
|-------|-----------|
| **0** | Wrong dimensions, low DPI (<150), pixelated at print size |
| **1** | Acceptable dimensions but minor issues (slightly over/under-sized, DPI adequate but not optimal) |
| **2** | Exact journal dimensions, DPI exceeds requirements, vector format for line art |

**What to check:**
- Figure width matches journal column specification (single/double)
- DPI meets minimum: 300 (photos), 600 (combination), 1000+ (line art)
- Exported at final size (not rescaled in LaTeX/Word)
- Vector format used for plots/diagrams

### 2. Units & Labels (0-2 points)

Are all axes, legends, and annotations complete and correctly formatted?

| Score | Criterion |
|-------|-----------|
| **0** | Missing axis labels, no units, unlabeled data series |
| **1** | Labels present but incomplete (missing units, inconsistent formatting, abbreviations undefined) |
| **2** | All axes labeled with units in parentheses, consistent formatting, all abbreviations defined in caption |

**What to check:**
- Every axis has a label with unit: "Time (hours)", "Concentration (mM)"
- Sentence case for labels (not ALL CAPS or Title Case)
- Panel labels present (A, B, C) for multi-panel figures
- Legend entries match data series
- Scale bars present for microscopy images
- Statistical annotations defined (*, **, ***)

### 3. Color Usage (0-2 points)

Are colors accessible, meaningful, and consistent?

| Score | Criterion |
|-------|-----------|
| **0** | Red-green combinations, rainbow/jet colormap, color as only encoding |
| **1** | Acceptable palette but not optimized (some accessibility issues, inconsistent across panels) |
| **2** | Colorblind-safe palette (Okabe-Ito/Paul Tol), redundant encoding (markers + line styles), consistent across all figures in manuscript |

**What to check:**
- No red-green combinations
- Perceptually uniform colormap for continuous data (viridis, cividis)
- Colorblind-safe palette for categories (Okabe-Ito, Paul Tol)
- Same condition = same color across all figures
- Redundant encoding (different markers, line styles, or patterns)
- Tested with colorblind simulator

### 4. Emphasis & Visual Hierarchy (0-2 points)

Does the figure direct the viewer's attention to the key message?

| Score | Criterion |
|-------|-----------|
| **0** | No clear focal point, all elements equally weighted, viewer cannot identify the main message |
| **1** | Some visual hierarchy but key finding not immediately obvious |
| **2** | Clear focal point, main message identifiable within 5 seconds, supporting elements subordinated |

**What to check:**
- Can a colleague identify the main message in <10 seconds?
- Key data emphasized (color, size, position, annotation)
- Supporting elements de-emphasized (lighter colors, smaller size)
- Visual flow guides the eye from context to result
- Title/annotation highlights the key finding (if appropriate for figure type)

### 5. Ink:Content Ratio (0-2 points)

Is every visual element necessary? (Based on Tufte's data-ink ratio principle)

| Score | Criterion |
|-------|-----------|
| **0** | Excessive chart junk: gridlines, 3D effects, shadows, background colors, decorative borders |
| **1** | Some unnecessary elements remain (minor gridlines, redundant borders) |
| **2** | Every pixel serves a purpose: data, labels, or essential structure only |

**What to check:**
- Top and right spines removed (unless essential)
- No background fill or gradient
- No 3D effects on 2D data
- No drop shadows or decorative elements
- Gridlines only if they aid value reading
- Legend frame removed or minimal
- No duplicate information (same data shown in both color and text)

### 6. Accessibility (0-2 points)

Can the figure be understood by all readers, including those with visual impairments?

| Score | Criterion |
|-------|-----------|
| **0** | Fails in grayscale, unreadable at print size, relies solely on color |
| **1** | Partially accessible (works in grayscale OR has redundant encoding, but not both) |
| **2** | Fully accessible: grayscale-compatible, colorblind-safe, text readable at final size (>=6 pt), high contrast |

**What to check:**
- Figure interpretable in grayscale
- Passes colorblind simulation (deuteranopia + protanopia)
- All text >= 6 pt at final print size
- Sufficient contrast between adjacent elements
- Line width >= 0.5 pt
- Patterns/hatching used in addition to color for bars/areas

## Scoring Thresholds

| Total Score | Rating | Action |
|-------------|--------|--------|
| **11-12** | Excellent | Ready for submission |
| **9-10** | Good | Minor revisions recommended |
| **7-8** | Acceptable | Revisions needed before submission |
| **5-6** | Poor | Significant revisions required |
| **0-4** | Unacceptable | Redesign the figure |

### Context-Specific Minimum Scores

| Context | Minimum Score | Rationale |
|---------|--------------|-----------|
| Journal (Nature, Science, Cell) | 10/12 | Highest standards; peer review scrutiny |
| Journal (general) | 9/12 | Professional publication standards |
| Conference paper | 9/12 | Peer-reviewed proceedings |
| Thesis/dissertation | 8/12 | Formal academic document |
| Grant proposal | 8/12 | Competitive review process |
| Preprint (arXiv, bioRxiv) | 7/12 | Public but not formally peer-reviewed |
| Poster | 7/12 | Viewed at distance; simpler requirements |
| Presentation slides | 6/12 | Complemented by verbal explanation |
| Internal report | 6/12 | Limited audience |

## Automated Rubric Scorer

```python
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


def score_figure(fig, target_journal='nature', verbose=True):
    """
    Score a matplotlib figure against the 6-criteria rubric.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to evaluate.
    target_journal : str
        Target journal for size checks.
    verbose : bool
        Print detailed scoring breakdown.

    Returns
    -------
    dict
        Scores per criterion and total.
    """
    scores = {}

    # --- 1. Scale & Resolution ---
    score_scale = 2
    w_in, h_in = fig.get_size_inches()
    dpi = fig.dpi

    journal_widths = {
        'nature': 3.5, 'science': 2.17, 'acs': 3.25,
        'cell': 3.35, 'elsevier': 3.54, 'plos': 3.27,
    }
    target_w = journal_widths.get(target_journal, 3.5)

    if dpi < 150:
        score_scale = 0
    elif dpi < 300:
        score_scale = 1
    if w_in > target_w * 2.2:
        score_scale = min(score_scale, 1)

    scores['scale_resolution'] = score_scale

    # --- 2. Units & Labels ---
    score_labels = 2
    for ax in fig.get_axes():
        xlabel = ax.get_xlabel()
        ylabel = ax.get_ylabel()
        if not xlabel.strip() or not ylabel.strip():
            score_labels = 0
            break
        if '(' not in xlabel or '(' not in ylabel:
            score_labels = min(score_labels, 1)

    scores['units_labels'] = score_labels

    # --- 3. Color Usage ---
    score_color = 2
    bad_cmaps = {'jet', 'rainbow', 'hsv', 'gist_rainbow', 'gist_ncar'}
    for ax in fig.get_axes():
        for img in ax.get_images():
            if img.get_cmap().name in bad_cmaps:
                score_color = 0
                break

    # Check for red-green combinations
    for ax in fig.get_axes():
        colors_rgb = []
        for line in ax.get_lines():
            try:
                c = mpl.colors.to_rgba(line.get_color())
                colors_rgb.append(c)
            except (ValueError, TypeError):
                pass
        has_red = any(c[0] > 0.6 and c[1] < 0.4 and c[2] < 0.4 for c in colors_rgb)
        has_green = any(c[1] > 0.6 and c[0] < 0.4 and c[2] < 0.4 for c in colors_rgb)
        if has_red and has_green:
            score_color = min(score_color, 0)

    scores['color_usage'] = score_color

    # --- 4. Emphasis ---
    # Heuristic: check if figure has a title or annotation highlighting key finding
    score_emphasis = 1  # Default: partial (cannot fully automate)
    for ax in fig.get_axes():
        if ax.get_title() or ax.texts:
            score_emphasis = 2
            break

    scores['emphasis'] = score_emphasis

    # --- 5. Ink:Content Ratio ---
    score_ink = 2
    for ax in fig.get_axes():
        if ax.spines['top'].get_visible() and ax.spines['right'].get_visible():
            score_ink = min(score_ink, 1)
        if ax.get_facecolor() != (1.0, 1.0, 1.0, 1.0):
            score_ink = min(score_ink, 1)
        if ax.grid:
            # Check if grid is enabled via lines
            if any(line.get_visible() for line in ax.get_xgridlines()):
                score_ink = min(score_ink, 1)

    scores['ink_content_ratio'] = score_ink

    # --- 6. Accessibility ---
    score_access = 2
    for ax in fig.get_axes():
        for text in ax.get_xticklabels() + ax.get_yticklabels():
            if text.get_fontsize() < 5:
                score_access = 0
                break
            elif text.get_fontsize() < 6:
                score_access = min(score_access, 1)

    # Check pdf.fonttype setting
    if mpl.rcParams.get('pdf.fonttype', 3) != 42:
        score_access = min(score_access, 1)

    scores['accessibility'] = score_access

    # --- Total ---
    total = sum(scores.values())
    scores['total'] = total
    scores['max'] = 12

    if verbose:
        print(f"\n{'=' * 50}")
        print(f"FIGURE QUALITY RUBRIC SCORE")
        print(f"{'=' * 50}")
        criteria_names = {
            'scale_resolution': '1. Scale & Resolution',
            'units_labels': '2. Units & Labels',
            'color_usage': '3. Color Usage',
            'emphasis': '4. Emphasis',
            'ink_content_ratio': '5. Ink:Content Ratio',
            'accessibility': '6. Accessibility',
        }
        for key, name in criteria_names.items():
            score = scores[key]
            bar = '*' * score + '.' * (2 - score)
            print(f"  {name:30s} [{bar}] {score}/2")
        print(f"{'─' * 50}")
        print(f"  {'TOTAL':30s}       {total}/12")

        # Rating
        if total >= 11:
            rating = 'EXCELLENT - Ready for submission'
        elif total >= 9:
            rating = 'GOOD - Minor revisions recommended'
        elif total >= 7:
            rating = 'ACCEPTABLE - Revisions needed'
        elif total >= 5:
            rating = 'POOR - Significant revisions required'
        else:
            rating = 'UNACCEPTABLE - Redesign the figure'
        print(f"  Rating: {rating}")
        print(f"{'=' * 50}\n")

    return scores
```

## Manual Evaluation Worksheet

Use this worksheet when reviewing figures that cannot be scored automatically (e.g., images, external graphics).

```
Figure: ________________    Target: ________________

1. Scale & Resolution    [ ] 0  [ ] 1  [ ] 2    Notes: ________
2. Units & Labels        [ ] 0  [ ] 1  [ ] 2    Notes: ________
3. Color Usage           [ ] 0  [ ] 1  [ ] 2    Notes: ________
4. Emphasis              [ ] 0  [ ] 1  [ ] 2    Notes: ________
5. Ink:Content Ratio     [ ] 0  [ ] 1  [ ] 2    Notes: ________
6. Accessibility         [ ] 0  [ ] 1  [ ] 2    Notes: ________
                                          TOTAL: ___ / 12

Rating:  [ ] Excellent (11-12)  [ ] Good (9-10)  [ ] Acceptable (7-8)
         [ ] Poor (5-6)  [ ] Unacceptable (0-4)

Action required: ________________________________________
```

## Resources

- Rokhsareh et al., ScientificFigures: https://github.com/nrokh/ScientificFigures
- Rougier et al. (2014), "Ten Simple Rules for Better Figures", PLOS Computational Biology
- Tufte, E. (2001), *The Visual Display of Quantitative Information*
- See `publication_guidelines.md` for the full Ten Simple Rules breakdown
- See `figure_antipatterns.md` for common mistakes and automated detection
