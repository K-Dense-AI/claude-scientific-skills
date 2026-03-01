# Scientific Figure Antipatterns

## Overview

This reference catalogs common antipatterns in scientific figure design, organized into five categories. Each antipattern includes a description of the problem, why it matters, the recommended solution, and (where applicable) code to detect or fix it automatically. Avoiding these mistakes improves clarity, accessibility, and the likelihood of acceptance during peer review.

## Category 1: Resolution and Format Mistakes

### 1.1 Using JPEG for Graphs and Line Art

**Problem:** JPEG uses lossy compression that introduces visible artifacts around sharp edges, text, and thin lines. These artifacts degrade readability, especially after journal typesetting rescales the image.

**Solution:** Use vector formats (PDF, EPS, SVG) for all graphs, plots, and diagrams. Use PNG or TIFF (lossless) for raster images like photographs or microscopy.

```python
def save_figure_correct(fig, basename, dpi=300):
    """Save figure in appropriate formats, never JPEG."""
    # Vector format for graphs (scalable, no artifacts)
    fig.savefig(f'{basename}.pdf', dpi=dpi, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    # Lossless raster for backup / web use
    fig.savefig(f'{basename}.png', dpi=dpi, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    # Never: fig.savefig(f'{basename}.jpg', ...)
```

### 1.2 Submitting Screenshots Instead of Exported Figures

**Problem:** Screen captures have fixed resolution (typically 72-96 DPI) far below publication requirements (300-1200 DPI). They also include window chrome, menu bars, and compression artifacts.

**Solution:** Always export directly from the plotting library. Set DPI explicitly before saving.

```python
import matplotlib.pyplot as plt

# Correct: export from code
fig, ax = plt.subplots(figsize=(3.5, 2.5))
ax.plot([1, 2, 3], [4, 5, 6])
fig.savefig('figure.pdf', dpi=300, bbox_inches='tight')

# Never: take a screenshot of the matplotlib window
```

### 1.3 Resizing Figures in LaTeX or Word Instead of at Source

**Problem:** Scaling a raster image down in LaTeX (`\includegraphics[width=...]`) or Word does not increase its effective resolution. Scaling up makes it blurry. Either way, fonts become inconsistent with the rest of the manuscript.

**Solution:** Set the figure to its final dimensions in the plotting code. Font sizes and line widths should be correct at the intended print size.

```python
# Set figure size to match journal column width BEFORE plotting
# Nature single column = 89 mm = 3.5 inches
fig, ax = plt.subplots(figsize=(3.5, 2.5))

# Fonts are set for this exact size
ax.set_xlabel('Time (hours)', fontsize=8)
ax.set_ylabel('Concentration (mM)', fontsize=8)
ax.tick_params(labelsize=7)

fig.savefig('figure.pdf', dpi=300, bbox_inches='tight')
```

```latex
% In LaTeX: include at natural size, no scaling
\includegraphics{figure.pdf}

% Avoid this:
% \includegraphics[width=0.5\textwidth]{figure.pdf}  % Rescales everything
```

### 1.4 Insufficient Resolution for Raster Elements

**Problem:** Raster images (microscopy, photographs) saved below 300 DPI appear pixelated in print. Line art below 600 DPI shows jagged edges.

**Solution:** Check and enforce minimum resolution before submission.

```python
from PIL import Image

def check_image_resolution(filepath, min_dpi=300):
    """Check if an image meets minimum resolution requirements."""
    img = Image.open(filepath)
    dpi = img.info.get('dpi', (72, 72))

    width_px, height_px = img.size
    width_in = width_px / dpi[0]
    height_in = height_px / dpi[1]

    print(f"Image: {filepath}")
    print(f"  Pixel dimensions: {width_px} x {height_px}")
    print(f"  DPI: {dpi[0]} x {dpi[1]}")
    print(f"  Print size at current DPI: {width_in:.2f} x {height_in:.2f} inches")

    if dpi[0] < min_dpi or dpi[1] < min_dpi:
        print(f"  WARNING: DPI below minimum ({min_dpi}). "
              f"Increase resolution or reduce print size.")
        return False
    print(f"  OK: Meets minimum DPI requirement.")
    return True
```

## Category 2: Color Mistakes

### 2.1 Using Rainbow (Jet) Colormap

**Problem:** The rainbow/jet colormap is not perceptually uniform: equal data differences produce unequal visual differences. It creates artificial boundaries, obscures real gradients, and fails completely in grayscale. Despite being a historical default, it is scientifically misleading.

**Solution:** Use perceptually uniform colormaps: viridis, cividis, plasma, inferno, or magma.

```python
import matplotlib.pyplot as plt
import numpy as np

data = np.random.randn(20, 20)

fig, axes = plt.subplots(1, 2, figsize=(7, 3))

# Wrong: jet/rainbow
axes[0].imshow(data, cmap='jet')
axes[0].set_title('Jet (misleading)', fontsize=8)

# Correct: viridis
axes[1].imshow(data, cmap='viridis')
axes[1].set_title('Viridis (perceptually uniform)', fontsize=8)

fig.tight_layout()
```

### 2.2 Red-Green Color Combinations

**Problem:** Approximately 8% of males and 0.5% of females have red-green color vision deficiency (deuteranopia or protanopia). Red-green combinations become indistinguishable for these readers, making the figure uninterpretable.

**Solution:** Replace red-green with blue-orange or use the Okabe-Ito palette. Always test with a colorblind simulator.

```python
# Wrong
colors_bad = ['#FF0000', '#00FF00']  # Red and green

# Correct: Okabe-Ito blue and orange
colors_good = ['#0072B2', '#E69F00']  # Blue and orange

# Correct: Add redundant encoding
# Different markers + different line styles + different colors
plt.plot(x, y1, color='#0072B2', marker='o', linestyle='-', label='Control')
plt.plot(x, y2, color='#E69F00', marker='s', linestyle='--', label='Treatment')
```

### 2.3 Too Many Colors (More Than 6 Categorical Colors)

**Problem:** Humans can reliably distinguish only about 6-8 colors simultaneously. Beyond this limit, colors become confusable and the legend becomes a decoding exercise rather than an aid.

**Solution:** Limit categorical colors to 6 or fewer. For more categories, use faceting (small multiples), grouping, or a sequential colormap with labeled annotations.

```python
import matplotlib.pyplot as plt

# Wrong: 12 categories each with a different color
# categories = [f'Group {i}' for i in range(12)]
# plt.plot(..., color=colors[i]) for each

# Correct: Facet into subplots
fig, axes = plt.subplots(2, 3, figsize=(7, 4), sharex=True, sharey=True)
groups = ['A', 'B', 'C', 'D', 'E', 'F']
for ax, group in zip(axes.flat, groups):
    # Plot each group in its own panel
    ax.plot(x, data[group], color='#0072B2')
    ax.set_title(group, fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
fig.tight_layout()
```

### 2.4 Using Color as the Only Encoding

**Problem:** If color is the sole channel for distinguishing data series, the figure fails for colorblind readers and in grayscale printouts.

**Solution:** Always provide a redundant visual channel: line style, marker shape, or pattern/hatching.

```python
styles = [
    {'color': '#0072B2', 'marker': 'o', 'linestyle': '-'},
    {'color': '#E69F00', 'marker': 's', 'linestyle': '--'},
    {'color': '#009E73', 'marker': '^', 'linestyle': '-.'},
]

for i, (label, style) in enumerate(zip(labels, styles)):
    ax.plot(x, y[i], label=label, markersize=4, **style)
```

### 2.5 Inconsistent Color Mapping Across Figures

**Problem:** Using blue for "Control" in Figure 1 and red for "Control" in Figure 3 forces readers to re-decode every figure and increases the risk of misinterpretation.

**Solution:** Define a color dictionary once and reuse throughout the manuscript.

```python
# Define once, use everywhere
CONDITION_COLORS = {
    'Control': '#0072B2',
    'Treatment A': '#E69F00',
    'Treatment B': '#009E73',
    'Mutant': '#D55E00',
}

# In every figure
ax.bar(conditions, values, color=[CONDITION_COLORS[c] for c in conditions])
```

## Category 3: Chart Type Mistakes

### 3.1 Pie Charts for Quantitative Comparison

**Problem:** Humans are poor at comparing areas and angles. Pie charts make it nearly impossible to judge small differences between slices, and completely fail when there are more than 4-5 categories.

**Solution:** Use horizontal bar charts for composition data, or stacked bar charts for comparing compositions across groups.

```python
import matplotlib.pyplot as plt

categories = ['Protein', 'Lipid', 'Carbohydrate', 'Nucleic acid', 'Other']
values = [42, 25, 18, 10, 5]

fig, axes = plt.subplots(1, 2, figsize=(7, 2.5))

# Wrong: Pie chart
axes[0].pie(values, labels=categories, autopct='%1.0f%%')
axes[0].set_title('Pie chart (hard to compare)', fontsize=8)

# Correct: Horizontal bar chart
axes[1].barh(categories, values, color='#0072B2', height=0.6)
axes[1].set_xlabel('Composition (%)', fontsize=8)
axes[1].set_title('Bar chart (easy to compare)', fontsize=8)
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)
axes[1].invert_yaxis()

fig.tight_layout()
```

### 3.2 3D Effects on 2D Data

**Problem:** 3D bar charts, 3D pie charts, and 3D surface plots of 2D data distort perceived values through perspective foreshortening. Bars at the back appear smaller; bars at the front occlude bars behind them.

**Solution:** Always use 2D representations for 2D data. Reserve 3D only for genuine 3D datasets where the third dimension carries information (e.g., molecular surfaces, topography).

```python
# Wrong: 3D bar chart
# from mpl_toolkits.mplot3d import Axes3D
# ax = fig.add_subplot(111, projection='3d')
# ax.bar3d(...)

# Correct: 2D grouped bar chart
x = np.arange(len(categories))
width = 0.35
ax.bar(x - width/2, values_a, width, label='Condition A', color='#0072B2')
ax.bar(x + width/2, values_b, width, label='Condition B', color='#E69F00')
```

### 3.3 Dual Y-Axes (Two Different Scales on Left and Right)

**Problem:** Dual y-axis charts are inherently misleading because the relationship between the two scales is arbitrary. By choosing different scale ranges, any two datasets can be made to appear correlated or anti-correlated. Readers also frequently misread which data series belongs to which axis.

**Solution:** Use two separate panels (subplots) with independent y-axes, placed side by side or stacked.

```python
import matplotlib.pyplot as plt

# Wrong: Dual y-axes
# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()
# ax1.plot(x, y1, color='blue')
# ax2.plot(x, y2, color='red')

# Correct: Two panels
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3.5, 4), sharex=True)

ax1.plot(x, y1, color='#0072B2')
ax1.set_ylabel('Temperature (K)', fontsize=8)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

ax2.plot(x, y2, color='#E69F00')
ax2.set_ylabel('Pressure (kPa)', fontsize=8)
ax2.set_xlabel('Time (min)', fontsize=8)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

fig.tight_layout()
```

### 3.4 Truncated Y-Axis on Bar Charts

**Problem:** Starting the y-axis above zero on a bar chart exaggerates differences between groups. A bar that appears twice as tall may represent only a 5% difference. This is one of the most common sources of visual deception in scientific and popular graphics.

**Solution:** Always start bar chart y-axes at zero. If the effect is real, it will be visible. If differences are too small to see, consider a dot plot or use an inset to zoom in.

```python
import matplotlib.pyplot as plt

values = [98, 100, 103, 101]
categories = ['A', 'B', 'C', 'D']

fig, axes = plt.subplots(1, 2, figsize=(7, 2.5))

# Wrong: Truncated y-axis
axes[0].bar(categories, values, color='#0072B2')
axes[0].set_ylim(95, 105)  # Exaggerates small differences
axes[0].set_title('Truncated (misleading)', fontsize=8)

# Correct: Y-axis starts at zero
axes[1].bar(categories, values, color='#0072B2')
axes[1].set_ylim(0, 110)
axes[1].set_title('Full range (honest)', fontsize=8)

for ax in axes:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

fig.tight_layout()
```

### 3.5 Dynamite Plots (Bar Chart + Error Bar for Continuous Data)

**Problem:** Bar charts with error bars ("dynamite plots") hide the underlying data distribution. They obscure whether the distribution is normal, bimodal, or skewed, and can mask outliers or unequal sample sizes.

**Solution:** Use box plots, violin plots, or strip/swarm plots that show the actual data distribution.

```python
import matplotlib.pyplot as plt
import numpy as np

np.random.seed(42)
data_a = np.random.normal(100, 15, 30)
data_b = np.concatenate([np.random.normal(80, 5, 15),
                          np.random.normal(120, 5, 15)])  # Bimodal!

fig, axes = plt.subplots(1, 2, figsize=(7, 3))

# Wrong: Bar + error bar hides bimodality
for ax_idx, (d, label) in enumerate([(data_a, 'Unimodal'), (data_b, 'Bimodal')]):
    axes[0].bar(ax_idx, np.mean(d), yerr=np.std(d), capsize=5,
                color=['#0072B2', '#E69F00'][ax_idx], alpha=0.7)
axes[0].set_xticks([0, 1])
axes[0].set_xticklabels(['A', 'B'])
axes[0].set_title('Bar plot (hides distribution)', fontsize=8)

# Correct: Strip plot shows bimodality
for ax_idx, (d, label) in enumerate([(data_a, 'A'), (data_b, 'B')]):
    x = np.random.normal(ax_idx, 0.05, len(d))
    axes[1].scatter(x, d, s=12, alpha=0.5,
                    color=['#0072B2', '#E69F00'][ax_idx])
    axes[1].plot([ax_idx - 0.15, ax_idx + 0.15],
                 [np.mean(d), np.mean(d)], 'k-', linewidth=1.5)
axes[1].set_xticks([0, 1])
axes[1].set_xticklabels(['A', 'B'])
axes[1].set_title('Strip plot (reveals distribution)', fontsize=8)

for ax in axes:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

fig.tight_layout()
```

## Category 4: Layout Mistakes

### 4.1 Using Software Default Styling

**Problem:** Default matplotlib, Excel, or R settings produce figures with gray backgrounds, excessive grid lines, thick borders, and garish colors that look unprofessional and often violate journal requirements.

**Solution:** Apply publication-quality settings before creating any plots.

```python
import matplotlib.pyplot as plt
import matplotlib as mpl

def set_publication_style():
    """Apply publication-ready defaults."""
    mpl.rcParams.update({
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'font.size': 8,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica'],
        'axes.labelsize': 9,
        'axes.titlesize': 9,
        'axes.linewidth': 0.5,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.major.size': 3,
        'ytick.major.size': 3,
        'legend.fontsize': 7,
        'legend.frameon': False,
        'lines.linewidth': 1.5,
        'lines.markersize': 4,
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'axes.grid': False,
        'image.cmap': 'viridis',
    })
    # Set colorblind-friendly color cycle
    okabe_ito = ['#E69F00', '#56B4E9', '#009E73', '#F0E442',
                 '#0072B2', '#D55E00', '#CC79A7']
    mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=okabe_ito)

# Call at the start of every script
set_publication_style()
```

### 4.2 Missing or Inadequate Axis Labels

**Problem:** Axes without labels, or labels without units, force readers to guess what is being plotted. This is the single most common feedback item from reviewers and editors.

**Solution:** Every axis must have a label that includes the variable name and its unit in parentheses.

```python
# Wrong
ax.set_xlabel('x')
ax.set_ylabel('y')

# Wrong: no units
ax.set_xlabel('Temperature')
ax.set_ylabel('Rate')

# Correct
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Reaction rate (mol L$^{-1}$ s$^{-1}$)')
```

### 4.3 Chart Junk (Unnecessary Visual Elements)

**Problem:** Gridlines, background colors, 3D effects, drop shadows, gradients, borders, and decorative elements add visual clutter without conveying data. Every non-data pixel reduces the data-ink ratio (Tufte).

**Solution:** Remove all non-essential elements. Keep only data, axes, labels, and legends.

```python
# Remove chart junk
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(False)                          # Remove gridlines
ax.set_facecolor('white')              # White background
ax.legend(frameon=False)               # No legend box
# No: shadows, gradients, 3D effects, background images
```

### 4.4 Legends That Obscure Data

**Problem:** Large legends placed inside the plot area cover data points. Legends with colored boxes instead of line samples are harder to decode.

**Solution:** Place legends where they do not overlap data. Use line/marker samples that match the actual plot elements.

```python
# Option 1: Inside the plot in an empty region
ax.legend(frameon=False, loc='upper left')

# Option 2: Outside the plot
ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc='upper left',
          borderaxespad=0)

# Option 3: Direct labeling (best for 2-3 series)
ax.text(x[-1] + 0.2, y1[-1], 'Control', fontsize=7, color='#0072B2',
        va='center')
ax.text(x[-1] + 0.2, y2[-1], 'Treatment', fontsize=7, color='#E69F00',
        va='center')
```

### 4.5 Inconsistent Panel Sizing and Spacing

**Problem:** Panels of different sizes, inconsistent spacing, or misaligned axes make multi-panel figures look unprofessional and are harder to compare visually.

**Solution:** Use GridSpec or constrained_layout for precise control.

```python
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(7, 4))
gs = gridspec.GridSpec(2, 3, figure=fig,
                       hspace=0.4, wspace=0.4,
                       left=0.08, right=0.95,
                       top=0.95, bottom=0.1)

# Consistent panel sizes
ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[0, 2])
ax_d = fig.add_subplot(gs[1, :2])  # Span two columns intentionally
ax_e = fig.add_subplot(gs[1, 2])
```

## Category 5: Data Representation Mistakes

### 5.1 Showing Only Mean Without Variability

**Problem:** A bar or point showing only the mean value provides no information about variability, reliability, or distribution shape. Readers cannot assess whether differences are meaningful.

**Solution:** Always show variability. State in the caption whether error bars represent SD, SEM, or confidence intervals.

```python
# Wrong: mean only
ax.bar(groups, means)

# Correct: mean + SEM + individual data points
ax.bar(groups, means, yerr=sems, capsize=3, color='#0072B2',
       alpha=0.4, edgecolor='#0072B2')
for i, (g, d) in enumerate(zip(groups, raw_data)):
    x_jitter = np.random.normal(i, 0.05, len(d))
    ax.scatter(x_jitter, d, s=10, color='#0072B2', alpha=0.6, zorder=3)
```

### 5.2 Not Reporting Sample Size (n)

**Problem:** Without n, readers cannot assess statistical power or whether error bars are meaningful. An SEM from n=3 and n=300 look very different in interpretation even if the bar height is similar.

**Solution:** Report n in the figure (annotations), caption, or both.

```python
# Add n to bar labels
for i, (group, n) in enumerate(zip(groups, sample_sizes)):
    ax.text(i, -5, f'n={n}', ha='center', fontsize=6, color='gray')

# Or in the caption:
# "Values represent mean +/- SEM. Sample sizes: Control (n=15),
#  Treatment A (n=12), Treatment B (n=18)."
```

### 5.3 Hiding Individual Data Points Behind Summary Statistics

**Problem:** Summary statistics (bar + error bar) can hide important patterns: outliers, bimodal distributions, unequal group sizes, or individual trends in longitudinal data.

**Solution:** Overlay individual data points on all summary plots when n < 100.

```python
import matplotlib.pyplot as plt
import numpy as np

def plot_with_individual_points(ax, groups, data_list, colors=None):
    """Plot summary statistics with overlaid individual data points."""
    if colors is None:
        colors = ['#0072B2'] * len(groups)

    for i, (group, data, color) in enumerate(zip(groups, data_list, colors)):
        # Summary: box or bar
        bp = ax.boxplot([data], positions=[i], widths=0.4,
                        patch_artist=True, showfliers=False,
                        boxprops=dict(facecolor='lightgray', edgecolor='black'),
                        medianprops=dict(color='black', linewidth=1.5),
                        whiskerprops=dict(linewidth=0.8),
                        capprops=dict(linewidth=0.8))

        # Individual points with jitter
        x_jitter = np.random.normal(i, 0.06, len(data))
        ax.scatter(x_jitter, data, s=10, alpha=0.5, color=color, zorder=3)

    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups)
```

### 5.4 SEM vs. SD Confusion

**Problem:** SEM (Standard Error of the Mean) and SD (Standard Deviation) answer different questions. SEM reflects precision of the mean estimate; SD reflects spread of the data. Using SEM makes error bars look deceptively small, especially with large n. Authors sometimes choose SEM specifically because bars look better.

**Solution:** Choose the appropriate measure and state it explicitly.

- **SD**: Use when showing data spread/variability
- **SEM**: Use when showing precision of the mean
- **95% CI**: Often the best choice for inferential comparisons

```python
import numpy as np

def compute_error_measures(data):
    """Compute SD, SEM, and 95% CI for a dataset."""
    n = len(data)
    mean = np.mean(data)
    sd = np.std(data, ddof=1)
    sem = sd / np.sqrt(n)
    ci_95 = 1.96 * sem  # Approximate for large n

    return {
        'mean': mean,
        'sd': sd,
        'sem': sem,
        'ci_95': ci_95,
        'n': n,
    }
```

### 5.5 Misleading Connecting Lines

**Problem:** Connecting data points with lines implies continuity or a functional relationship. Lines between discrete, unrelated categories (e.g., different treatments) suggest a trend that does not exist.

**Solution:** Use lines only for time series, continuous variables, or when interpolation between points is meaningful. Use bars or dot plots for discrete categories.

```python
# Wrong: lines between discrete categories
# ax.plot(['Control', 'Drug A', 'Drug B'], [10, 15, 13], '-o')

# Correct: dot plot for discrete categories
ax.scatter(['Control', 'Drug A', 'Drug B'], [10, 15, 13],
           s=40, color='#0072B2', zorder=3)
ax.errorbar(['Control', 'Drug A', 'Drug B'], [10, 15, 13],
            yerr=[1, 2, 1.5], fmt='none', color='black', capsize=3)
```

## Automated Validation Functions

### Comprehensive Figure Checker

```python
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import warnings


def validate_figure(fig, target_journal='general', verbose=True):
    """
    Validate a matplotlib figure against common antipatterns.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to validate.
    target_journal : str
        Target journal for size/format checks. Options:
        'nature', 'science', 'acs', 'general'.
    verbose : bool
        If True, print detailed diagnostics.

    Returns
    -------
    dict
        Dictionary with 'errors' (must fix) and 'warnings' (should fix).
    """
    errors = []
    warn = []

    # Journal size limits (width in inches)
    size_limits = {
        'nature': {'single': 3.5, 'double': 7.2, 'max_height': 9.7},
        'science': {'single': 2.17, 'double': 6.89, 'max_height': 9.17},
        'acs': {'single': 3.25, 'double': 7.0, 'max_height': 9.5},
        'general': {'single': 3.5, 'double': 7.5, 'max_height': 10.0},
    }
    limits = size_limits.get(target_journal, size_limits['general'])

    # --- Check figure dimensions ---
    fig_w, fig_h = fig.get_size_inches()
    if fig_w > limits['double']:
        errors.append(
            f"Figure width ({fig_w:.1f} in) exceeds maximum "
            f"({limits['double']} in) for {target_journal}."
        )
    if fig_h > limits['max_height']:
        errors.append(
            f"Figure height ({fig_h:.1f} in) exceeds maximum "
            f"({limits['max_height']} in) for {target_journal}."
        )

    for ax in fig.get_axes():
        # --- Check axis labels ---
        xlabel = ax.get_xlabel()
        ylabel = ax.get_ylabel()
        if not xlabel and ax.get_xticks().size > 0:
            errors.append(f"Axis missing x-label: {ax.get_title() or 'unnamed'}")
        if not ylabel and ax.get_yticks().size > 0:
            errors.append(f"Axis missing y-label: {ax.get_title() or 'unnamed'}")

        # Check for units in labels
        for label_text, axis_name in [(xlabel, 'x'), (ylabel, 'y')]:
            if label_text and '(' not in label_text and label_text.strip():
                warn.append(
                    f"{axis_name}-label '{label_text}' may be missing units. "
                    f"Use format: 'Variable name (unit)'."
                )

        # --- Check font sizes ---
        for text_obj in ax.get_xticklabels() + ax.get_yticklabels():
            size = text_obj.get_fontsize()
            if size < 5:
                errors.append(
                    f"Tick label font size ({size} pt) below minimum readable "
                    f"size (5 pt)."
                )
            elif size < 6:
                warn.append(
                    f"Tick label font size ({size} pt) may be too small for "
                    f"some journals (recommend 6-8 pt)."
                )

        # --- Check for colorbar on images ---
        for img in ax.get_images():
            has_colorbar = any(
                isinstance(child, mpl.colorbar.Colorbar)
                for child in fig.get_children()
            )
            # Heuristic: check if any axes look like a colorbar
            if not has_colorbar:
                warn.append(
                    "Image/heatmap detected but no colorbar found. "
                    "Add a colorbar with labeled ticks."
                )

        # --- Check for excessive colors ---
        all_colors = set()
        for line in ax.get_lines():
            all_colors.add(line.get_color())
        for coll in ax.collections:
            fc = coll.get_facecolor()
            if hasattr(fc, '__len__') and len(fc) > 0:
                for c in fc:
                    all_colors.add(tuple(c) if hasattr(c, '__len__') else c)
        if len(all_colors) > 8:
            warn.append(
                f"Found {len(all_colors)} distinct colors in one axes. "
                f"Consider reducing to 6-8 or using faceting."
            )

    # --- Check colormap usage ---
    for ax in fig.get_axes():
        for img in ax.get_images():
            cmap_name = img.get_cmap().name
            bad_cmaps = ['jet', 'rainbow', 'hsv', 'gist_rainbow',
                         'gist_ncar', 'nipy_spectral']
            if cmap_name in bad_cmaps:
                errors.append(
                    f"Non-perceptually-uniform colormap '{cmap_name}' detected. "
                    f"Use viridis, cividis, plasma, or inferno instead."
                )

    # --- Report ---
    if verbose:
        if errors:
            print("ERRORS (must fix):")
            for e in errors:
                print(f"  [X] {e}")
        if warn:
            print("WARNINGS (should fix):")
            for w in warn:
                print(f"  [!] {w}")
        if not errors and not warn:
            print("All checks passed.")

    return {'errors': errors, 'warnings': warn}


def check_red_green(fig):
    """
    Check if a figure uses red-green color combinations.

    Returns list of warnings about potential colorblind issues.
    """
    issues = []

    def is_red(color_rgb):
        r, g, b = color_rgb[:3]
        return r > 0.6 and g < 0.4 and b < 0.4

    def is_green(color_rgb):
        r, g, b = color_rgb[:3]
        return g > 0.6 and r < 0.4 and b < 0.4

    for ax in fig.get_axes():
        colors_in_plot = []
        for line in ax.get_lines():
            c = mpl.colors.to_rgba(line.get_color())
            colors_in_plot.append(c)
        for coll in ax.collections:
            fc = coll.get_facecolor()
            if hasattr(fc, '__len__'):
                for c in fc:
                    colors_in_plot.append(tuple(c))

        has_red = any(is_red(c) for c in colors_in_plot)
        has_green = any(is_green(c) for c in colors_in_plot)

        if has_red and has_green:
            issues.append(
                f"Red-green color combination detected in "
                f"'{ax.get_title() or 'axes'}'. "
                f"This is problematic for ~8% of male readers. "
                f"Replace with blue-orange or Okabe-Ito palette."
            )

    return issues
```

### Usage Example

```python
# After creating your figure:
fig, ax = plt.subplots(figsize=(3.5, 2.5))
ax.plot(x, y1, color='red', label='A')
ax.plot(x, y2, color='green', label='B')
ax.set_xlabel('Time')  # Missing units!

# Run validation
results = validate_figure(fig, target_journal='nature')
# Output:
#   ERRORS (must fix):
#     [X] Non-perceptually-uniform colormap... (if heatmap present)
#   WARNINGS (should fix):
#     [!] x-label 'Time' may be missing units.

# Check colorblind safety
issues = check_red_green(fig)
# Output: Red-green color combination detected...
```

## Quick Reference: Antipattern Checklist

Before submitting any figure, verify:

- [ ] No JPEG artifacts (saved as PDF/PNG/TIFF, not JPG)
- [ ] Not a screenshot (exported programmatically at correct DPI)
- [ ] Sized correctly at source (not rescaled in LaTeX/Word)
- [ ] Resolution meets journal requirements (300+ DPI for raster)
- [ ] No rainbow/jet colormap (use viridis, cividis, etc.)
- [ ] No red-green color pairs (use blue-orange, Okabe-Ito)
- [ ] 6 or fewer categorical colors (facet if more needed)
- [ ] Redundant encoding beyond color (markers, line styles)
- [ ] Consistent colors across all figures in manuscript
- [ ] No pie charts (use bar charts instead)
- [ ] No 3D effects on 2D data
- [ ] No dual y-axes (use separate panels)
- [ ] Bar chart y-axes start at zero
- [ ] No dynamite plots (show distributions instead)
- [ ] Publication-quality styling applied (not software defaults)
- [ ] All axes labeled with units
- [ ] No chart junk (no gridlines, shadows, gradients)
- [ ] Legend does not obscure data
- [ ] Consistent panel sizes and spacing
- [ ] Error bars shown with variability measure stated
- [ ] Sample sizes (n) reported
- [ ] Individual data points visible where appropriate
- [ ] SD vs. SEM chosen appropriately and stated in caption

## Resources

- **Weissgerber et al. (2015)**: "Beyond Bar and Line Graphs" -- PLOS Biology article on data presentation
- **Rougier et al. (2014)**: "Ten Simple Rules for Better Figures" -- PLOS Computational Biology
- **Tufte, Edward (2001)**: *The Visual Display of Quantitative Information* -- Foundational text on data-ink ratio
- **Nature Methods Points of View**: Column series on visualization best practices
- **Crameri et al. (2020)**: "The misuse of colour in science communication" -- Nature Communications
- See `color_palettes.md` for colorblind-friendly palette details
- See `publication_guidelines.md` for general figure design principles
- See `journal_requirements.md` for publisher-specific technical specifications
