# Figure Reproducibility Workflow

## Overview

Every figure in a scientific manuscript should be reproducible from source data and code. This reference describes a workflow for scripting all figures so they can be regenerated with a single command. Based on best practices from Wookai/paper-tips-and-tricks and reproducible research templates.

## Directory Structure

Organize figure source files alongside the manuscript:

```
project/
  manuscript/
    main.tex
    references.bib
  figures/
    config.py            # Shared style configuration
    figure1_kinetics.py  # One script per figure
    figure2_comparison.py
    figure3_heatmap.py
    figure_S1_controls.py
    Makefile             # Regenerate all figures
  data/
    raw/                 # Unprocessed data files
    processed/           # Cleaned/aggregated data
  output/
    pdf/                 # Generated figure PDFs
    png/                 # Generated figure PNGs
    tiff/                # Generated figure TIFFs
```

## Shared Configuration File

Create a single `config.py` that enforces consistent styling across all figures:

```python
"""
figures/config.py
Shared configuration for all manuscript figures.

Usage:
    from config import setup_figure, save_figure, COLORS
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path

# --- Project paths ---
PROJECT_ROOT = Path(__file__).parent.parent
FIGURE_DIR = Path(__file__).parent
OUTPUT_DIR = PROJECT_ROOT / 'output'
DATA_DIR = PROJECT_ROOT / 'data'

# Ensure output directories exist
for subdir in ['pdf', 'png', 'tiff']:
    (OUTPUT_DIR / subdir).mkdir(parents=True, exist_ok=True)

# --- Color palette (Okabe-Ito) ---
COLORS = {
    'control': '#0072B2',       # Blue
    'treatment_a': '#E69F00',   # Orange
    'treatment_b': '#009E73',   # Bluish green
    'treatment_c': '#D55E00',   # Vermillion
    'highlight': '#CC79A7',     # Reddish purple
    'neutral': '#999999',       # Gray
}
COLOR_LIST = list(COLORS.values())

# --- Marker and line style cycles ---
MARKERS = ['o', 's', '^', 'v', 'D', 'P']
LINESTYLES = ['-', '--', '-.', ':', (0, (3, 1, 1, 1))]

# --- Journal configuration ---
JOURNAL = 'nature'  # Change per manuscript

JOURNAL_WIDTHS = {  # inches
    'nature': {'single': 3.5, 'double': 7.2},
    'science': {'single': 2.17, 'double': 6.89},
    'acs': {'single': 3.25, 'double': 7.0},
    'cell': {'single': 3.35, 'double': 7.01},
    'elsevier': {'single': 3.54, 'double': 7.48},
}


def setup_figure(width='single', height_ratio=0.75, journal=None):
    """
    Create a figure with correct dimensions for the target journal.

    Parameters
    ----------
    width : str or float
        'single', 'double', or width in inches.
    height_ratio : float
        Height as fraction of width (default 0.75).
    journal : str, optional
        Override the global JOURNAL setting.

    Returns
    -------
    fig, ax : matplotlib Figure and Axes
    """
    j = journal or JOURNAL
    widths = JOURNAL_WIDTHS.get(j, JOURNAL_WIDTHS['nature'])

    if isinstance(width, str):
        w = widths.get(width, widths['single'])
    else:
        w = float(width)

    h = w * height_ratio

    _apply_style()
    fig, ax = plt.subplots(figsize=(w, h))
    return fig, ax


def setup_multipanel(nrows, ncols, width='double', height_ratio=0.5,
                     journal=None, **gridspec_kw):
    """
    Create a multi-panel figure with correct dimensions.

    Returns
    -------
    fig, axes : matplotlib Figure and array of Axes
    """
    j = journal or JOURNAL
    widths = JOURNAL_WIDTHS.get(j, JOURNAL_WIDTHS['nature'])

    if isinstance(width, str):
        w = widths.get(width, widths['double'])
    else:
        w = float(width)

    h = w * height_ratio

    _apply_style()
    fig, axes = plt.subplots(nrows, ncols, figsize=(w, h),
                              gridspec_kw=gridspec_kw)
    return fig, axes


def save_figure(fig, name, formats=None, dpi=300):
    """
    Save figure in multiple formats to the output directory.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    name : str
        Base filename (without extension).
    formats : list of str
        File formats to save. Default: ['pdf', 'png'].
    dpi : int
        DPI for raster formats.
    """
    if formats is None:
        formats = ['pdf', 'png']

    for fmt in formats:
        outpath = OUTPUT_DIR / fmt / f'{name}.{fmt}'
        save_kwargs = {
            'dpi': dpi,
            'bbox_inches': 'tight',
            'facecolor': 'white',
            'edgecolor': 'none',
        }
        if fmt == 'tiff':
            save_kwargs['pil_kwargs'] = {'compression': 'tiff_lzw'}

        fig.savefig(outpath, **save_kwargs)
        print(f"  Saved: {outpath}")

    plt.close(fig)


def add_panel_labels(axes, labels=None, fontsize=10, x=-0.15, y=1.05):
    """
    Add panel labels (A, B, C, ...) to a list of axes.

    Parameters
    ----------
    axes : list of matplotlib.axes.Axes
    labels : list of str, optional
        Custom labels. Default: uppercase A, B, C, ...
    """
    if labels is None:
        labels = [chr(65 + i) for i in range(len(axes))]

    for ax, label in zip(axes, labels):
        ax.text(x, y, label, transform=ax.transAxes,
                fontsize=fontsize, fontweight='bold', va='top')


def _apply_style():
    """Apply publication-quality matplotlib defaults."""
    mpl.rcParams.update({
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'font.size': 8,
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
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
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })
    mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=COLOR_LIST)
```

## Per-Figure Script Template

Each figure gets its own script that imports from `config.py`:

```python
#!/usr/bin/env python3
"""
figures/figure1_kinetics.py
Figure 1: Michaelis-Menten kinetics of wild-type and mutant enzymes.
"""

import numpy as np
import pandas as pd
from config import setup_figure, save_figure, COLORS, DATA_DIR


def main():
    # --- Load data ---
    df = pd.read_csv(DATA_DIR / 'processed' / 'kinetics.csv')

    # --- Create figure ---
    fig, ax = setup_figure(width='single', height_ratio=0.8)

    # --- Plot ---
    for enzyme, color in [('WT', COLORS['control']),
                           ('MutA', COLORS['treatment_a']),
                           ('MutB', COLORS['treatment_b'])]:
        subset = df[df['enzyme'] == enzyme]
        ax.errorbar(subset['substrate_mM'], subset['rate_mean'],
                    yerr=subset['rate_sem'], fmt='o-', color=color,
                    markersize=4, capsize=2, label=enzyme)

    # --- Labels ---
    ax.set_xlabel('Substrate concentration (mM)')
    ax.set_ylabel('Reaction rate (nmol min$^{-1}$ mg$^{-1}$)')
    ax.legend()

    # --- Save ---
    save_figure(fig, 'figure1_kinetics', formats=['pdf', 'png', 'tiff'])


if __name__ == '__main__':
    main()
```

## Makefile for Batch Regeneration

```makefile
# figures/Makefile
# Regenerate all manuscript figures from source data.

PYTHON = python
SCRIPTS = $(wildcard figure*.py)
TARGETS = $(patsubst %.py,../output/pdf/%.pdf,$(SCRIPTS))

.PHONY: all clean

all: $(TARGETS)
	@echo "All figures generated."

../output/pdf/%.pdf: %.py config.py
	$(PYTHON) $<

clean:
	rm -rf ../output/pdf/*.pdf ../output/png/*.png ../output/tiff/*.tiff
	@echo "Cleaned all generated figures."
```

**Usage:**
```bash
cd figures/
make        # Regenerate all figures
make clean  # Remove all generated figures
```

## Version Control Best Practices

### What to commit
- `figures/config.py` and all `figure*.py` scripts
- `data/processed/` (cleaned data used by scripts)
- `figures/Makefile`

### What to .gitignore
```
# Generated outputs (can be regenerated from scripts)
output/pdf/
output/png/
output/tiff/
```

### Commit message convention
```
fig: update Figure 2 colors to Okabe-Ito palette
fig: add error bars to Figure 3
fig: regenerate all figures for Nature resubmission
```

## LaTeX Integration

Reference generated figures from the output directory:

```latex
\begin{figure}[htbp]
  \centering
  \includegraphics{output/pdf/figure1_kinetics.pdf}
  \caption{Michaelis-Menten kinetics of wild-type and mutant enzymes.
  Data represent mean $\pm$ SEM ($n = 3$ independent experiments).}
  \label{fig:kinetics}
\end{figure}
```

**Key principle:** Never rescale figures in LaTeX. The `config.py` ensures figures are exported at exactly the correct dimensions for the target journal.

## Workflow Summary

1. **Setup**: Create `figures/config.py` with journal and palette settings
2. **Script**: Write one `figure_N_name.py` per figure, importing from config
3. **Generate**: Run `make` to regenerate all figures
4. **Validate**: Use `figure_validator.py` and `figure_rubric.md` to check quality
5. **Iterate**: Modify script, re-run, figures update automatically
6. **Submit**: Generated outputs are always consistent and publication-ready

## Resources

- Wookai, "Paper Tips and Tricks": https://github.com/Wookai/paper-tips-and-tricks
- Mouret, "Matplotlib for Papers": https://github.com/jbmouret/matplotlib_for_papers
- Chure, "Reproducible Research Template": https://github.com/gchure/reproducible_research
- See `config.py` template above for ready-to-use project configuration
