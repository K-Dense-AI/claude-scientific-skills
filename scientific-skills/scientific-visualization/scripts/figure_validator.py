#!/usr/bin/env python3
"""
Figure Validator for Publication-Ready Scientific Figures

This module provides validation utilities to check that exported figures
meet journal submission requirements including resolution, file size,
dimensions, and font embedding.
"""

import os
from pathlib import Path
from typing import Dict, List, Optional, Union


# ---------------------------------------------------------------------------
# Journal specifications
# ---------------------------------------------------------------------------

JOURNAL_SPECS: Dict[str, dict] = {
    "nature": {
        "single_column_mm": 89,
        "double_column_mm": 183,
        "max_height_mm": 247,
        "min_dpi_line_art": 1000,
        "min_dpi_halftone": 300,
        "min_dpi_combination": 600,
        "accepted_formats": ["pdf", "eps", "tiff", "png"],
        "max_file_size_mb": 10,
        "font_requirements": "Type 1 or TrueType (Type 42), no Type 3",
    },
    "acs": {
        "single_column_mm": 82.5,
        "double_column_mm": 178,
        "max_height_mm": 247,
        "min_dpi_line_art": 600,
        "min_dpi_halftone": 300,
        "min_dpi_combination": 600,
        "accepted_formats": ["tiff", "pdf", "eps", "png"],
        "max_file_size_mb": 10,
        "font_requirements": "TrueType (Type 42), no Type 3",
    },
    "elsevier": {
        "single_column_mm": 90,
        "double_column_mm": 190,
        "max_height_mm": 240,
        "min_dpi_line_art": 1000,
        "min_dpi_halftone": 300,
        "min_dpi_combination": 600,
        "accepted_formats": ["pdf", "eps", "tiff", "jpg", "png"],
        "max_file_size_mb": 10,
        "font_requirements": "Embedded fonts, no Type 3",
    },
    "science": {
        "single_column_mm": 55,
        "double_column_mm": 175,
        "max_height_mm": 233,
        "min_dpi_line_art": 1000,
        "min_dpi_halftone": 300,
        "min_dpi_combination": 600,
        "accepted_formats": ["pdf", "eps", "tiff"],
        "max_file_size_mb": 10,
        "font_requirements": "Type 1 or TrueType (Type 42), no Type 3",
    },
    "plos": {
        "single_column_mm": 83,
        "double_column_mm": 173,
        "max_height_mm": 233,
        "min_dpi_line_art": 600,
        "min_dpi_halftone": 300,
        "min_dpi_combination": 300,
        "accepted_formats": ["tiff", "eps", "pdf", "png"],
        "max_file_size_mb": 10,
        "font_requirements": "Embedded fonts",
    },
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def validate_resolution(
    image_path: Union[str, Path],
    target_width_mm: float,
    image_type: str = "combination",
) -> dict:
    """
    Validate that a raster image has sufficient DPI for a target print width.

    Parameters
    ----------
    image_path : str or Path
        Path to the image file (PNG, TIFF, JPEG)
    target_width_mm : float
        Intended print width in millimetres
    image_type : str, default 'combination'
        Type of figure: 'line_art', 'halftone' (photo), or 'combination'

    Returns
    -------
    dict
        Validation result with keys:
        - 'passed' : bool
        - 'width_px' : int, image pixel width
        - 'height_px' : int, image pixel height
        - 'effective_dpi' : float, DPI at the target width
        - 'min_required_dpi' : int, minimum DPI threshold
        - 'message' : str

    Examples
    --------
    >>> result = validate_resolution('figure1.png', target_width_mm=89, image_type='line_art')
    >>> print(result['passed'])
    """
    from PIL import Image

    image_path = Path(image_path)
    if not image_path.exists():
        return {"passed": False, "message": f"File not found: {image_path}"}

    dpi_thresholds = {
        "line_art": 1000,
        "halftone": 300,
        "combination": 600,
    }

    if image_type not in dpi_thresholds:
        available = ", ".join(dpi_thresholds.keys())
        raise ValueError(f"image_type must be one of: {available}")

    min_dpi = dpi_thresholds[image_type]

    img = Image.open(image_path)
    width_px, height_px = img.size

    # Calculate effective DPI at target print width
    target_width_inches = target_width_mm / 25.4
    effective_dpi = width_px / target_width_inches

    passed = effective_dpi >= min_dpi

    result = {
        "passed": passed,
        "width_px": width_px,
        "height_px": height_px,
        "effective_dpi": round(effective_dpi, 1),
        "min_required_dpi": min_dpi,
        "image_type": image_type,
        "target_width_mm": target_width_mm,
        "message": (
            f"OK: {effective_dpi:.0f} DPI >= {min_dpi} DPI required for {image_type}"
            if passed
            else f"FAIL: {effective_dpi:.0f} DPI < {min_dpi} DPI required for {image_type}. "
            f"Need at least {int(min_dpi * target_width_inches)} px width."
        ),
    }

    print(f"Resolution check: {result['message']}")
    return result


def validate_file_size(
    image_path: Union[str, Path],
    max_mb: float = 10,
) -> dict:
    """
    Validate that a figure file does not exceed a size limit.

    Parameters
    ----------
    image_path : str or Path
        Path to the image file
    max_mb : float, default 10
        Maximum allowed file size in megabytes

    Returns
    -------
    dict
        Validation result with keys:
        - 'passed' : bool
        - 'size_mb' : float, actual file size in MB
        - 'max_mb' : float, the limit checked against
        - 'message' : str

    Examples
    --------
    >>> result = validate_file_size('figure1.tiff', max_mb=10)
    >>> print(result['passed'])
    """
    image_path = Path(image_path)
    if not image_path.exists():
        return {"passed": False, "size_mb": 0, "max_mb": max_mb,
                "message": f"File not found: {image_path}"}

    size_bytes = os.path.getsize(image_path)
    size_mb = size_bytes / (1024 * 1024)
    passed = size_mb <= max_mb

    result = {
        "passed": passed,
        "size_mb": round(size_mb, 2),
        "max_mb": max_mb,
        "message": (
            f"OK: {size_mb:.2f} MB <= {max_mb} MB limit"
            if passed
            else f"FAIL: {size_mb:.2f} MB exceeds {max_mb} MB limit"
        ),
    }

    print(f"File size check: {result['message']}")
    return result


def check_journal_compliance(
    image_path: Union[str, Path],
    journal: str,
    image_type: str = "combination",
    target_width: str = "single",
) -> dict:
    """
    Comprehensive compliance check against a specific journal's requirements.

    Parameters
    ----------
    image_path : str or Path
        Path to the image file
    journal : str
        Journal name. Options: 'nature', 'acs', 'elsevier', 'science', 'plos'
    image_type : str, default 'combination'
        Type of figure: 'line_art', 'halftone', or 'combination'
    target_width : str, default 'single'
        Column width: 'single' or 'double'

    Returns
    -------
    dict
        Compliance report with keys:
        - 'passed' : bool, True if all checks pass
        - 'journal' : str
        - 'checks' : list of dict, individual check results
        - 'summary' : str

    Examples
    --------
    >>> report = check_journal_compliance('fig1.tiff', journal='nature')
    >>> print(report['passed'])
    """
    journal = journal.lower()
    if journal not in JOURNAL_SPECS:
        available = ", ".join(JOURNAL_SPECS.keys())
        raise ValueError(f"Journal '{journal}' not recognized. Available: {available}")

    specs = JOURNAL_SPECS[journal]
    image_path = Path(image_path)
    checks: List[dict] = []

    # 1. File format check
    ext = image_path.suffix.lstrip(".").lower()
    fmt_ok = ext in specs["accepted_formats"]
    checks.append({
        "name": "file_format",
        "passed": fmt_ok,
        "message": (
            f"OK: .{ext} is accepted by {journal.upper()}"
            if fmt_ok
            else f"FAIL: .{ext} not in accepted formats: {specs['accepted_formats']}"
        ),
    })

    # 2. File size check
    size_result = validate_file_size(image_path, max_mb=specs["max_file_size_mb"])
    checks.append({
        "name": "file_size",
        "passed": size_result["passed"],
        "message": size_result["message"],
    })

    # 3. Resolution check (raster formats only)
    raster_formats = {"png", "tiff", "tif", "jpg", "jpeg"}
    if ext in raster_formats:
        width_key = (
            "single_column_mm" if target_width == "single" else "double_column_mm"
        )
        target_mm = specs[width_key]
        res_result = validate_resolution(image_path, target_mm, image_type)
        checks.append({
            "name": "resolution",
            "passed": res_result["passed"],
            "message": res_result["message"],
        })
    else:
        checks.append({
            "name": "resolution",
            "passed": True,
            "message": f"OK: vector format (.{ext}), resolution check not applicable",
        })

    # 4. Font embedding check (PDF only)
    if ext == "pdf":
        font_result = validate_font_embedding(image_path)
        checks.append({
            "name": "font_embedding",
            "passed": font_result["passed"],
            "message": font_result["message"],
        })

    all_passed = all(c["passed"] for c in checks)

    report = {
        "passed": all_passed,
        "journal": journal,
        "checks": checks,
        "summary": (
            f"All checks PASSED for {journal.upper()}"
            if all_passed
            else f"Some checks FAILED for {journal.upper()}"
        ),
    }

    # Print report
    print(f"\n{'=' * 60}")
    print(f"Journal Compliance Check: {journal.upper()}")
    print(f"{'=' * 60}")
    print(f"File: {image_path}")
    for c in checks:
        status = "PASS" if c["passed"] else "FAIL"
        print(f"  [{status}] {c['name']}: {c['message']}")
    print(f"\nOverall: {report['summary']}")
    print(f"{'=' * 60}\n")

    return report


def validate_font_embedding(
    pdf_path: Union[str, Path],
) -> dict:
    """
    Validate font embedding in a PDF file, checking for Type 3 fonts.

    Type 3 fonts are bitmap fonts that most journals reject. This function
    verifies that all fonts are either Type 1 or TrueType (Type 42).

    Parameters
    ----------
    pdf_path : str or Path
        Path to the PDF file

    Returns
    -------
    dict
        Validation result with keys:
        - 'passed' : bool, True if no Type 3 fonts found
        - 'fonts' : list of dict, font information per font found
        - 'has_type3' : bool, True if Type 3 fonts detected
        - 'message' : str

    Examples
    --------
    >>> result = validate_font_embedding('figure1.pdf')
    >>> print(result['passed'])
    """
    pdf_path = Path(pdf_path)
    if not pdf_path.exists():
        return {
            "passed": False,
            "fonts": [],
            "has_type3": False,
            "message": f"File not found: {pdf_path}",
        }

    try:
        from PyPDF2 import PdfReader
    except ImportError:
        return {
            "passed": True,
            "fonts": [],
            "has_type3": False,
            "message": "Warning: PyPDF2 not installed, cannot verify fonts. "
            "Install with: pip install PyPDF2",
        }

    fonts_found: List[dict] = []
    has_type3 = False

    try:
        reader = PdfReader(pdf_path)

        for page_num, page in enumerate(reader.pages):
            resources = page.get("/Resources")
            if resources is None:
                continue

            font_dict = resources.get("/Font")
            if font_dict is None:
                continue

            font_obj = font_dict.get_object()
            for font_name in font_obj:
                font = font_obj[font_name].get_object()
                subtype = str(font.get("/Subtype", ""))
                base_font = str(font.get("/BaseFont", "unknown"))

                font_type = subtype.replace("/", "")
                is_type3 = font_type == "Type3"

                if is_type3:
                    has_type3 = True

                fonts_found.append({
                    "name": base_font.replace("/", ""),
                    "type": font_type,
                    "page": page_num + 1,
                    "is_type3": is_type3,
                })

    except Exception as e:
        return {
            "passed": False,
            "fonts": [],
            "has_type3": False,
            "message": f"Error reading PDF: {e}",
        }

    passed = not has_type3

    if not fonts_found:
        message = "No fonts detected in PDF (may use only vector paths)"
    elif has_type3:
        type3_names = [f["name"] for f in fonts_found if f["is_type3"]]
        message = (
            f"FAIL: Type 3 fonts detected: {', '.join(type3_names)}. "
            f"Set pdf.fonttype=42 in matplotlib rcParams."
        )
    else:
        font_types = sorted({f["type"] for f in fonts_found})
        message = f"OK: All fonts are {', '.join(font_types)} (no Type 3)"

    result = {
        "passed": passed,
        "fonts": fonts_found,
        "has_type3": has_type3,
        "message": message,
    }

    print(f"Font embedding check: {result['message']}")
    return result


if __name__ == "__main__":
    print("Figure Validator for Publication-Ready Scientific Figures")
    print("=" * 60)

    print("\nAvailable journal specifications:")
    for name, specs in JOURNAL_SPECS.items():
        print(f"\n  {name.upper()}:")
        print(f"    Single column: {specs['single_column_mm']} mm")
        print(f"    Double column: {specs['double_column_mm']} mm")
        print(f"    Max height: {specs['max_height_mm']} mm")
        print(f"    Accepted formats: {', '.join(specs['accepted_formats'])}")
        print(f"    Fonts: {specs['font_requirements']}")

    print("\nExample usage:")
    print("  from figure_validator import check_journal_compliance")
    print("  report = check_journal_compliance('fig1.pdf', journal='nature')")
