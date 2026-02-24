#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from pathlib import Path
import pandas as pd


def pick_metric(path: Path, col: str):
    if not path.exists():
        return "N/A"
    df = pd.read_csv(path)
    if col not in df.columns or df.empty:
        return "N/A"
    return f"{df[col].iloc[0]}"


def main():
    root = Path(__file__).resolve().parents[3]
    ex = root / "scientific-skills" / "qtl-analysis" / "examples"
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)

    lambda_after = pick_metric(
        ex / "genomic-control" / "output" / "lambda_summary.csv", "value"
    )
    qtl_count = "N/A"
    mtrait = ex / "multi-trait-gwas" / "output" / "gwas_results.csv"
    if mtrait.exists():
        df = pd.read_csv(mtrait)
        if "P_MultiTrait" in df.columns:
            qtl_count = int((df["P_MultiTrait"] < 0.05 / len(df)).sum())

    html = f"""
<html><head><title>QTL Analysis Report</title></head><body>
<h1>QTL Analysis Report</h1>
<p>Generated from committed example outputs.</p>
<h2>Key Metrics</h2>
<ul>
  <li>Multi-trait significant hits (Bonferroni): {qtl_count}</li>
  <li>Genomic control lambda (after): {lambda_after}</li>
</ul>
<h2>Included Modules</h2>
<ul>
  <li>Bayesian GP</li>
  <li>Advanced GWAS (multi-trait, GxE, covariates, thresholds, genomic control, rare variants)</li>
  <li>Kinship & relatedness (pedigree NRM, genomic NRM, IBS)</li>
</ul>
</body></html>
"""

    md = f"""# QTL Analysis Report

Generated from example outputs.

## Key Metrics
- Multi-trait significant hits (Bonferroni): {qtl_count}
- Genomic control lambda (after): {lambda_after}

## Included Modules
- Bayesian GP
- Advanced GWAS
- Kinship and relatedness
"""

    (out / "analysis_report.html").write_text(html.strip() + "\n", encoding="utf-8")
    (out / "analysis_report.md").write_text(md, encoding="utf-8")
    print("Saved analysis_report.html and analysis_report.md")


if __name__ == "__main__":
    main()
