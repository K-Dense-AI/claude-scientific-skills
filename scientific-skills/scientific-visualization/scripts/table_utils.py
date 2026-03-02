#!/usr/bin/env python3
"""
Table Utilities for Publication-Ready LaTeX Tables

This module provides utilities to convert pandas DataFrames into
publication-quality LaTeX tables with booktabs formatting, uncertainty
notation, and best-value highlighting.
"""

import re
import math
from typing import List, Optional, Union

import pandas as pd


def df_to_booktabs(
    df: pd.DataFrame,
    caption: str,
    label: str,
    column_format: Optional[str] = None,
    escape: bool = True,
    index: bool = False,
) -> str:
    """
    Convert a pandas DataFrame to a booktabs-style LaTeX table.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to convert
    caption : str
        Table caption text
    label : str
        LaTeX label for cross-referencing (e.g., 'tab:results')
    column_format : str, optional
        LaTeX column alignment string (e.g., 'lccc'). If None, auto-generated
        with left-aligned first column and centered remaining columns.
    escape : bool, default True
        If True, escape special LaTeX characters in cell values
    index : bool, default False
        If True, include the DataFrame index as the first column

    Returns
    -------
    str
        LaTeX table string with booktabs formatting

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({'Method': ['A', 'B'], 'RMSE': [0.123, 0.098]})
    >>> print(df_to_booktabs(df, caption='Results', label='tab:results'))
    """
    if column_format is None:
        n_cols = len(df.columns) + (1 if index else 0)
        column_format = "l" + "c" * (n_cols - 1)

    # Build header row
    headers = []
    if index:
        headers.append(df.index.name or "")
    headers.extend(df.columns.tolist())

    if escape:
        headers = [_escape_latex(str(h)) for h in headers]

    header_line = " & ".join(headers) + r" \\"

    # Build data rows
    rows = []
    for idx, row in df.iterrows():
        cells = []
        if index:
            cells.append(str(idx))
        cells.extend([str(v) for v in row])
        if escape:
            cells = [_escape_latex(c) for c in cells]
        rows.append(" & ".join(cells) + r" \\")

    data_lines = "\n        ".join(rows)

    table = rf"""\begin{{table}}[htbp]
    \centering
    \caption{{{caption}}}
    \label{{{label}}}
    \begin{{tabular}}{{{column_format}}}
        \toprule
        {header_line}
        \midrule
        {data_lines}
        \bottomrule
    \end{{tabular}}
\end{{table}}"""

    return table


def format_uncertainty(
    value: float,
    uncertainty: float,
    sig_figs: int = 3,
) -> str:
    """
    Format a value with its uncertainty using standard notation.

    The uncertainty is rounded to 1 significant figure, then the value
    is rounded to match the decimal place of the uncertainty.

    Parameters
    ----------
    value : float
        The measured value
    uncertainty : float
        The uncertainty (standard deviation, standard error, etc.)
    sig_figs : int, default 3
        Number of significant figures for the value when uncertainty is zero

    Returns
    -------
    str
        Formatted string in the form 'value +/- uncertainty' or
        'value \\pm uncertainty' suitable for LaTeX

    Examples
    --------
    >>> format_uncertainty(3.14159, 0.025)
    '3.14 \\\\pm 0.03'
    >>> format_uncertainty(1234, 56)
    '1230 \\\\pm 60'
    >>> format_uncertainty(0.00321, 0.00004)
    '0.00321 \\\\pm 0.00004'
    """
    if uncertainty == 0:
        return _format_sig_figs(value, sig_figs)

    uncertainty = abs(uncertainty)

    # Determine the order of magnitude of the uncertainty
    unc_order = math.floor(math.log10(uncertainty))

    # Round uncertainty to 1 significant figure
    unc_rounded = round(uncertainty, -unc_order)

    # Determine number of decimal places
    if unc_order >= 0:
        decimal_places = 0
    else:
        decimal_places = -unc_order

    # Round value to match
    val_rounded = round(value, decimal_places)

    # Format strings
    if decimal_places > 0:
        val_str = f"{val_rounded:.{decimal_places}f}"
        unc_str = f"{unc_rounded:.{decimal_places}f}"
    else:
        val_str = f"{int(val_rounded)}"
        unc_str = f"{int(unc_rounded)}"

    return rf"{val_str} \pm {unc_str}"


def highlight_best(
    df: pd.DataFrame,
    column: str,
    criterion: str = "max",
) -> pd.DataFrame:
    """
    Highlight the best value in a column by wrapping it in LaTeX bold.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to process
    column : str
        Column name to evaluate
    criterion : str, default 'max'
        Criterion for best value: 'max' or 'min'

    Returns
    -------
    pd.DataFrame
        A copy of the DataFrame with the best value in the specified
        column wrapped in '\\textbf{...}'

    Examples
    --------
    >>> df = pd.DataFrame({'Method': ['A', 'B'], 'Accuracy': [0.95, 0.98]})
    >>> result = highlight_best(df, 'Accuracy', criterion='max')
    >>> print(result.loc[1, 'Accuracy'])
    \\textbf{0.98}
    """
    if column not in df.columns:
        raise ValueError(f"Column '{column}' not found in DataFrame")

    if criterion not in ("max", "min"):
        raise ValueError(f"Criterion must be 'max' or 'min', got '{criterion}'")

    df_out = df.copy()

    # Convert column to numeric for comparison
    numeric_vals = pd.to_numeric(df_out[column], errors="coerce")

    if numeric_vals.isna().all():
        return df_out

    if criterion == "max":
        best_idx = numeric_vals.idxmax()
    else:
        best_idx = numeric_vals.idxmin()

    best_val = df_out.loc[best_idx, column]
    df_out.loc[best_idx, column] = rf"\textbf{{{best_val}}}"

    return df_out


def validate_table(df: pd.DataFrame) -> dict:
    """
    Validate a DataFrame for publication-readiness.

    Checks include:
    - Consistent significant figures within each numeric column
    - Detection of unit information in column headers (e.g., 'Mass (kg)')
    - Missing values
    - Duplicate rows

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to validate

    Returns
    -------
    dict
        Validation report with keys:
        - 'passed' : bool, True if all checks pass
        - 'warnings' : list of str, warning messages
        - 'errors' : list of str, error messages
        - 'column_info' : dict, per-column information

    Examples
    --------
    >>> df = pd.DataFrame({'Mass (kg)': [1.23, 4.5], 'Temp': [300, 310]})
    >>> report = validate_table(df)
    >>> print(report['passed'])
    """
    warnings: List[str] = []
    errors: List[str] = []
    column_info: dict = {}

    # Check for missing values
    missing = df.isnull().sum()
    for col, count in missing.items():
        if count > 0:
            warnings.append(f"Column '{col}' has {count} missing value(s)")

    # Check for duplicate rows
    n_dupes = df.duplicated().sum()
    if n_dupes > 0:
        warnings.append(f"Found {n_dupes} duplicate row(s)")

    # Per-column analysis
    unit_pattern = re.compile(r"\(([^)]+)\)$")

    for col in df.columns:
        info: dict = {"has_unit": False, "unit": None, "sig_figs_consistent": True}

        # Check for unit in header
        match = unit_pattern.search(str(col).strip())
        if match:
            info["has_unit"] = True
            info["unit"] = match.group(1)
        else:
            # Only warn for numeric columns without units
            if pd.api.types.is_numeric_dtype(df[col]):
                warnings.append(
                    f"Column '{col}' appears numeric but has no unit in header"
                )

        # Check significant figures consistency for numeric columns
        if pd.api.types.is_numeric_dtype(df[col]):
            sig_figs_list = []
            for val in df[col].dropna():
                sf = _count_sig_figs(val)
                if sf is not None:
                    sig_figs_list.append(sf)

            if sig_figs_list and len(set(sig_figs_list)) > 1:
                info["sig_figs_consistent"] = False
                warnings.append(
                    f"Column '{col}' has inconsistent significant figures: "
                    f"{sorted(set(sig_figs_list))}"
                )

            info["sig_figs"] = sorted(set(sig_figs_list)) if sig_figs_list else []

        column_info[col] = info

    passed = len(errors) == 0

    report = {
        "passed": passed,
        "warnings": warnings,
        "errors": errors,
        "column_info": column_info,
    }

    # Print summary
    status = "PASSED" if passed else "FAILED"
    print(f"\nTable Validation: {status}")
    print(f"  {len(errors)} error(s), {len(warnings)} warning(s)")
    for e in errors:
        print(f"  [ERROR] {e}")
    for w in warnings:
        print(f"  [WARN]  {w}")

    return report


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _escape_latex(text: str) -> str:
    """Escape special LaTeX characters in a string."""
    special_chars = {
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    for char, replacement in special_chars.items():
        text = text.replace(char, replacement)
    return text


def _format_sig_figs(value: float, sig_figs: int) -> str:
    """Format a float to a given number of significant figures."""
    if value == 0:
        return f"0.{'0' * (sig_figs - 1)}"
    magnitude = math.floor(math.log10(abs(value)))
    decimal_places = sig_figs - 1 - magnitude
    if decimal_places > 0:
        return f"{value:.{decimal_places}f}"
    rounded = round(value, -int(magnitude - sig_figs + 1))
    return f"{int(rounded)}" if decimal_places <= 0 else f"{rounded}"


def _count_sig_figs(value: float) -> Optional[int]:
    """
    Estimate the number of significant figures from a float value.

    Returns None if the value cannot be analysed (e.g., zero).
    """
    if value == 0:
        return None

    # Convert to string, strip leading/trailing zeros carefully
    text = f"{value:g}"  # compact representation
    text = text.lstrip("-")  # ignore sign

    if "e" in text or "E" in text:
        mantissa = text.split("e")[0].split("E")[0]
        mantissa = mantissa.replace(".", "")
        return len(mantissa)

    # Remove leading zeros and decimal point
    stripped = text.replace(".", "")
    stripped = stripped.lstrip("0")
    return len(stripped) if stripped else None


if __name__ == "__main__":
    # Example usage
    df = pd.DataFrame(
        {
            "Method": ["Random Forest", "XGBoost", "Neural Net"],
            "RMSE": [0.123, 0.098, 0.105],
            "R^2": [0.89, 0.94, 0.92],
            "Time (s)": [1.2, 3.4, 15.7],
        }
    )

    print("=== DataFrame to Booktabs ===")
    print(df_to_booktabs(df, caption="Model comparison", label="tab:models"))

    print("\n=== Uncertainty Formatting ===")
    print(format_uncertainty(3.14159, 0.025))
    print(format_uncertainty(1234, 56))
    print(format_uncertainty(0.00321, 0.00004))

    print("\n=== Highlight Best ===")
    highlighted = highlight_best(df, "RMSE", criterion="min")
    print(highlighted)

    print("\n=== Validate Table ===")
    validate_table(df)
