# Scientific Table Formatting Guide

## Overview

Well-formatted tables are essential for presenting quantitative data in scientific publications. The booktabs three-line style (top rule, mid rule, bottom rule) is the standard for professional typesetting, replacing the cluttered grid lines of default table styles. This reference covers LaTeX table creation, pandas-to-LaTeX conversion, number formatting with siunitx, and journal-specific requirements.

## Booktabs Three-Line Tables

### The Booktabs Philosophy

The booktabs package eliminates vertical rules and excessive horizontal rules in favor of three key lines:
- **Top rule** (`\toprule`): Heavy line at the top of the table
- **Mid rule** (`\midrule`): Medium line separating header from data
- **Bottom rule** (`\bottomrule`): Heavy line at the bottom of the table

Additional `\cmidrule` lines can separate subgroups within the header or data.

### Basic Table

```latex
\usepackage{booktabs}

\begin{table}[htbp]
  \centering
  \caption{Kinetic parameters for enzymatic reactions at 25\textdegree{}C.}
  \label{tab:kinetics}
  \begin{tabular}{lSS}
    \toprule
    Enzyme & {$K_\mathrm{m}$ (mM)} & {$k_\mathrm{cat}$ (s$^{-1}$)} \\
    \midrule
    Wild type    & 0.45  & 120  \\
    Mutant A     & 1.23  & 85   \\
    Mutant B     & 0.31  & 210  \\
    Mutant C     & 2.10  & 42   \\
    \bottomrule
  \end{tabular}
\end{table}
```

### Grouped Header Table

```latex
\usepackage{booktabs}

\begin{table}[htbp]
  \centering
  \caption{Comparison of catalyst performance under different conditions.}
  \label{tab:catalyst}
  \begin{tabular}{l SS SS}
    \toprule
    & \multicolumn{2}{c}{Condition A} & \multicolumn{2}{c}{Condition B} \\
    \cmidrule(lr){2-3} \cmidrule(lr){4-5}
    Catalyst & {Yield (\%)} & {$T_{1/2}$ (h)} & {Yield (\%)} & {$T_{1/2}$ (h)} \\
    \midrule
    Pd/C       & 92.3  & 1.2  & 88.1  & 1.8  \\
    Pt/Al$_2$O$_3$ & 78.5  & 2.4  & 81.2  & 2.1  \\
    Rh/SiO$_2$ & 95.1  & 0.8  & 93.7  & 1.0  \\
    \bottomrule
  \end{tabular}
\end{table}
```

### Complex Table with Subgroups and Notes

```latex
\usepackage{booktabs}
\usepackage{siunitx}
\usepackage{threeparttable}

\begin{table}[htbp]
  \centering
  \begin{threeparttable}
    \caption{Growth parameters of bacterial strains under aerobic and anaerobic conditions.}
    \label{tab:growth}
    \begin{tabular}{l l S[table-format=1.2] S[table-format=2.1] S[table-format=1.2]}
      \toprule
      Strain & Condition & {$\mu_\mathrm{max}$ (h$^{-1}$)} & {$Y_{X/S}$ (\%)} & {Lag (h)} \\
      \midrule
      \multirow{2}{*}{\textit{E.~coli} K-12}
        & Aerobic   & 0.87 & 42.3 & 1.20 \\
        & Anaerobic & 0.31 & 18.7 & 3.45 \\
      \addlinespace
      \multirow{2}{*}{\textit{E.~coli} BL21}
        & Aerobic   & 0.92 & 45.1 & 0.95 \\
        & Anaerobic & 0.28 & 16.2 & 4.10 \\
      \addlinespace
      \multirow{2}{*}{\textit{B.~subtilis} 168\tnote{a}}
        & Aerobic   & 0.65 & 38.9 & 1.80 \\
        & Anaerobic & 0.12 & 8.30 & 6.20 \\
      \bottomrule
    \end{tabular}
    \begin{tablenotes}
      \small
      \item[a] Sporulation-deficient mutant.
    \end{tablenotes}
  \end{threeparttable}
\end{table}
```

## pandas to LaTeX Conversion

### Basic Conversion Function

```python
import pandas as pd
import numpy as np

def dataframe_to_booktabs(df, caption, label, filename=None,
                          sig_figs=3, column_format=None):
    """
    Convert a pandas DataFrame to a booktabs-style LaTeX table.

    Parameters
    ----------
    df : pd.DataFrame
        Data to convert.
    caption : str
        Table caption text.
    label : str
        LaTeX label for cross-referencing.
    filename : str, optional
        If provided, write LaTeX to this file.
    sig_figs : int
        Number of significant figures for float columns.
    column_format : str, optional
        Custom LaTeX column format string. If None, auto-generated.

    Returns
    -------
    str
        LaTeX table code.
    """
    # Round float columns to significant figures
    df_formatted = df.copy()
    for col in df_formatted.select_dtypes(include=[np.floating]).columns:
        df_formatted[col] = df_formatted[col].apply(
            lambda x: f'{x:.{sig_figs}g}' if pd.notna(x) else ''
        )

    # Generate column format if not provided
    if column_format is None:
        formats = ['l']  # index column
        for col in df.columns:
            if pd.api.types.is_numeric_dtype(df[col]):
                formats.append('S')
            else:
                formats.append('l')
        column_format = ' '.join(formats)

    # Build LaTeX string
    latex = df_formatted.to_latex(
        index=True,
        escape=False,
        column_format=column_format,
        caption=caption,
        label=label,
        position='htbp',
    )

    # Replace default rules with booktabs rules
    latex = latex.replace('\\hline', '')
    lines = latex.split('\n')
    result = []
    header_found = False
    for line in lines:
        result.append(line)
        if '\\centering' in line or ('\\begin{tabular}' in line and not header_found):
            continue
        if not header_found and '\\\\' in line:
            result.append('    \\midrule')
            header_found = True

    latex = '\n'.join(result)

    if filename:
        with open(filename, 'w') as f:
            f.write(latex)

    return latex
```

### Conversion with siunitx Formatting

```python
def dataframe_to_siunitx_table(df, caption, label,
                                 uncertainty_cols=None,
                                 si_units=None):
    """
    Convert DataFrame to LaTeX table with siunitx number formatting.

    Parameters
    ----------
    df : pd.DataFrame
        Data to convert.
    caption : str
        Table caption.
    label : str
        LaTeX label.
    uncertainty_cols : dict, optional
        Mapping of value column to uncertainty column, e.g.,
        {'yield': 'yield_err'} produces "92.3 +/- 1.2" formatting.
    si_units : dict, optional
        Mapping of column name to SI unit string, e.g.,
        {'temperature': 'K', 'pressure': 'kPa'}.

    Returns
    -------
    str
        LaTeX table code with siunitx formatting.
    """
    if uncertainty_cols is None:
        uncertainty_cols = {}
    if si_units is None:
        si_units = {}

    # Build header with units
    header_parts = []
    for col in df.columns:
        if col in uncertainty_cols.values():
            continue  # Skip uncertainty columns (merged with value)
        if col in si_units:
            header_parts.append(f'{{{col} (\\si{{{si_units[col]}}})}}}')
        else:
            header_parts.append(f'{{{col}}}')

    # Build column format
    col_formats = []
    for col in df.columns:
        if col in uncertainty_cols.values():
            continue
        if col in uncertainty_cols:
            col_formats.append('S@{\\,\\(\\pm\\)\\,}S')
        elif pd.api.types.is_numeric_dtype(df[col]):
            col_formats.append('S')
        else:
            col_formats.append('l')

    # Build rows
    rows = []
    for _, row in df.iterrows():
        cells = []
        for col in df.columns:
            if col in uncertainty_cols.values():
                continue
            if col in uncertainty_cols:
                val = row[col]
                unc = row[uncertainty_cols[col]]
                cells.append(f'{val} & {unc}')
            else:
                cells.append(str(row[col]))
        rows.append(' & '.join(cells) + ' \\\\')

    # Assemble table
    col_format_str = ' '.join(col_formats)
    header_str = ' & '.join(header_parts) + ' \\\\'

    latex = f"""\\begin{{table}}[htbp]
  \\centering
  \\caption{{{caption}}}
  \\label{{{label}}}
  \\begin{{tabular}}{{{col_format_str}}}
    \\toprule
    {header_str}
    \\midrule
"""
    for r in rows:
        latex += f'    {r}\n'
    latex += """    \\bottomrule
  \\end{tabular}
\\end{table}"""

    return latex


# Example usage
data = pd.DataFrame({
    'Catalyst': ['Pd/C', 'Pt/Al2O3', 'Rh/SiO2'],
    'Yield': [92.3, 78.5, 95.1],
    'Yield_err': [1.2, 2.4, 0.8],
    'TOF': [1250, 890, 1580],
})

latex_code = dataframe_to_siunitx_table(
    data,
    caption='Catalytic performance with uncertainties.',
    label='tab:catalysis',
    uncertainty_cols={'Yield': 'Yield_err'},
    si_units={'TOF': '\\per\\hour'},
)
print(latex_code)
```

## siunitx Package Usage

### Essential Preamble

```latex
\usepackage{siunitx}

% Common configuration
\sisetup{
  separate-uncertainty = true,    % 1.23 +/- 0.05 instead of 1.23(5)
  table-align-uncertainty = true, % Align uncertainties in columns
  detect-weight = true,           % Detect bold in headers
  detect-family = true,           % Detect font family
  group-separator = {,},          % Thousands separator
  group-minimum-digits = 4,       % Apply separator for 4+ digit numbers
}
```

### Number Formatting Examples

```latex
% Basic numbers
\num{12345}          % 12,345
\num{0.00034}        % 0.000 34
\num{1.23e4}         % 1.23 x 10^4
\num{6.022e23}       % 6.022 x 10^23

% Numbers with uncertainty
\num{1.23(4)}        % 1.23 +/- 0.04 (compact notation)
\num{1.23 +- 0.04}   % 1.23 +/- 0.04 (explicit notation)

% Units
\SI{273.15}{\kelvin}          % 273.15 K
\SI{101.325}{\kilo\pascal}    % 101.325 kPa
\SI{3.0e8}{\metre\per\second} % 3.0 x 10^8 m/s

% Ranges
\SIrange{20}{30}{\celsius}    % 20 to 30 C
\numrange{1.2}{3.4}           % 1.2 to 3.4
```

### S Column for Aligned Numbers

```latex
% S column aligns numbers on the decimal point
\begin{tabular}{l S[table-format=3.2] S[table-format=1.2e1]}
  \toprule
  Sample & {Mass (\si{\gram})} & {$k$ (\si{\per\second})} \\
  \midrule
  A & 12.34  & 3.45e-3 \\
  B & 123.40 & 1.20e-2 \\
  C & 1.23   & 9.87e-4 \\
  \bottomrule
\end{tabular}
```

**table-format specifier**: `3.2` means up to 3 digits before decimal, 2 after. Use `1.2e1` for scientific notation (1 digit, 2 decimals, 1-digit exponent).

## Alignment and Number Formatting Rules

### Alignment Conventions

| Data Type | Alignment | Rationale |
|-----------|-----------|-----------|
| Text / labels | Left (`l`) | Natural reading direction |
| Integers | Right (`r`) or `S` column | Aligns ones, tens, hundreds places |
| Decimals | Decimal-aligned (`S`) | Aligns magnitudes visually |
| Uncertainties | Aligned with value (`S@{\,}S`) | Paired with corresponding value |
| Mixed text/number | Left (`l`) | Treat as text when mixed |

### Significant Figures

Follow these conventions for reporting numerical results:

1. **Match precision to uncertainty**: Report values to the same decimal place as the uncertainty
   - Correct: 12.34 +/- 0.05
   - Wrong: 12.3421 +/- 0.05

2. **Consistent decimal places within a column**: All values in a column should have the same number of decimal places
   - Correct: 1.20, 3.45, 0.98
   - Wrong: 1.2, 3.45, 0.978

3. **Trailing zeros are meaningful**: 1.20 implies precision to hundredths; do not strip to 1.2

4. **Scientific notation for very large/small numbers**: Use when values span more than 3 orders of magnitude
   - Correct: 3.45 x 10^-3
   - Wrong: 0.00345

### Uncertainty Notation

```latex
% Preferred styles (choose one and be consistent):

% Style 1: Plus-minus notation (most common in chemistry/biology)
$K_\mathrm{m} = \SI{0.45 +- 0.03}{\milli\molar}$

% Style 2: Parenthetical notation (common in physics)
$K_\mathrm{m} = \SI{0.45(3)}{\milli\molar}$

% Style 3: Separate column (for tables)
% Value column | Uncertainty column
% 0.45         | 0.03
```

## Journal-Specific Table Requirements

### Nature Portfolio
- Tables placed after references in manuscript
- No vertical rules; booktabs style required
- Short informative title above the table
- Footnotes below using superscript lowercase letters (a, b, c)
- Units in column headers, not repeated in cells
- Maximum width: 183 mm (double column)

### Science (AAAS)
- Tables embedded in manuscript text
- Booktabs style strongly preferred
- Title as first row of table or above
- Footnotes use symbols: *, dagger, double-dagger, section
- Concise: data better suited to figures should be plotted

### ACS Journals
- Booktabs style with three-line format
- Title above table; double-spaced
- Footnotes with superscript italic lowercase letters
- Chemical formulas in table cells should use proper formatting
- Numbers aligned on decimal points

### Elsevier
- Tables on separate pages at end of manuscript
- Three-line format (booktabs)
- Footnotes: superscript lowercase letters
- Avoid overly wide tables (prefer splitting into two)
- Units in column headers in parentheses

### PLOS
- Embedded in manuscript
- No color in tables (data only)
- Three-line style preferred
- Footnotes with superscript lowercase letters
- Submit as editable text (not images)

### Quick Comparison

| Feature | Nature | Science | ACS | Elsevier | PLOS |
|---------|--------|---------|-----|----------|------|
| Placement | After refs | In text | In text | End of MS | In text |
| Vertical rules | No | No | No | No | No |
| Footnote style | a, b, c | *, dagger | ^a, ^b | a, b, c | a, b, c |
| Max width | 183 mm | 17.5 cm | 7 in | 190 mm | 17.3 cm |
| Booktabs | Required | Preferred | Required | Preferred | Preferred |

## Table Caption Writing

### Structure

A good table caption has three components:

1. **Title**: Brief descriptive statement (sentence case, ends with period)
2. **Method context** (if needed): Experimental conditions or data source
3. **Abbreviation definitions**: Define any abbreviations used in the table

### Examples

**Good captions:**
```
Table 1. Kinetic parameters for wild-type and mutant enzymes at pH 7.4 and 25 C.
Values represent mean +/- SD from three independent experiments (n = 3).
K_m, Michaelis constant; k_cat, catalytic rate constant; ND, not determined.
```

```
Table 2. Elemental analysis of synthesized nanoparticles by ICP-OES.
All values in weight percent. Detection limit: 0.01 wt%.
```

**Bad captions:**
```
Table 1. Results.                          % Too vague
Table 2. This table shows the data         % Starts with "This table shows"
         from our experiments.
Table 3. Comprehensive Overview of All     % Too long; title case
         Parameters Measured During the
         Complete Experimental Campaign
```

### Caption Checklist
- [ ] Sentence case (not Title Case)
- [ ] Ends with a period
- [ ] Describes content, not just "Results" or "Data"
- [ ] Defines all abbreviations
- [ ] States statistical measure (mean +/- SD, median, etc.)
- [ ] States sample size (n)
- [ ] Does not repeat information obvious from column headers

## Common Table Mistakes

### 1. Using Grid Lines Instead of Booktabs

**Wrong:**
```latex
\begin{tabular}{|l|c|c|}
  \hline
  Sample & Value 1 & Value 2 \\
  \hline
  A & 1.23 & 4.56 \\
  \hline
  B & 7.89 & 0.12 \\
  \hline
\end{tabular}
```

**Correct:**
```latex
\begin{tabular}{lSS}
  \toprule
  Sample & {Value 1} & {Value 2} \\
  \midrule
  A & 1.23 & 4.56 \\
  B & 7.89 & 0.12 \\
  \bottomrule
\end{tabular}
```

### 2. Misaligned Decimal Points

**Wrong:** Numbers aligned left or center, making magnitude comparison difficult.

**Correct:** Use `S` columns from siunitx for automatic decimal alignment.

### 3. Repeating Units in Every Cell

**Wrong:**
```
Mass: 12.3 g, 45.6 g, 78.9 g
```

**Correct:**
```
Mass (g): 12.3, 45.6, 78.9  % Unit in header only
```

### 4. Inconsistent Significant Figures

**Wrong:** Mixing precision levels within a column (1.2, 3.456, 78).

**Correct:** All values in a column use the same number of decimal places.

### 5. Overly Wide Tables

**Wrong:** Table exceeds page width, gets clipped or scaled down illegibly.

**Solutions:**
- Abbreviate column headers (define in footnotes)
- Split into two tables
- Rotate to landscape (`\begin{sidewaystable}`)
- Use `\resizebox{\textwidth}{!}{...}` as a last resort (can distort font sizes)

### 6. Missing Uncertainty Information

**Wrong:** Reporting only mean values without any indication of variability.

**Correct:** Include SD, SEM, or confidence intervals; state n and which measure is used.

### 7. Table as Screenshot or Image

**Wrong:** Pasting a table as a raster image from Excel or Word.

**Correct:** Always submit tables as editable text (LaTeX, Word table, or CSV). Image tables cannot be typeset, searched, or reformatted by publishers.

### 8. Color-Dependent Tables

**Wrong:** Using color shading as the only way to convey grouping or significance.

**Correct:** Use `\addlinespace`, `\cmidrule`, or footnote markers instead. Color may not survive print or be accessible to colorblind readers.

### 9. Excessive Precision from Software Output

**Wrong:** Copying raw floating-point output (e.g., 3.141592653589793).

**Correct:** Round to meaningful precision based on measurement uncertainty. Use `\num{3.14}` or format in pandas before export.

### 10. No Cross-Reference in Text

**Wrong:** "The results are shown below" with no `\ref{}`.

**Correct:** "The kinetic parameters are summarized in Table~\ref{tab:kinetics}."

## Resources

- **booktabs package documentation**: CTAN booktabs
- **siunitx package documentation**: CTAN siunitx
- **threeparttable package**: For proper table footnotes
- **pandas to_latex()**: https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_latex.html
- **Small Guide to Making Nice Tables** (Markus Pueschel): Classic reference for LaTeX table design
- See `journal_requirements.md` for publisher-specific figure and table specifications
