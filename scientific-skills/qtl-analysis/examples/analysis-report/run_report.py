#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def pick_metric(path: Path, col: str):
    if not path.exists():
        return "N/A"
    df = pd.read_csv(path)
    if col not in df.columns or df.empty:
        return "N/A"
    return f"{df[col].iloc[0]}"


def main():
    root = Path(__file__).resolve().parents[4]
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

    metrics = pd.DataFrame(
        {
            "metric": ["multi_trait_hits", "lambda_after_gc"],
            "value": [
                float(qtl_count) if qtl_count != "N/A" else 0.0,
                float(lambda_after) if lambda_after != "N/A" else 0.0,
            ],
        }
    )
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(9.2, 4.2), gridspec_kw={"width_ratios": [1.6, 1]}
    )
    vals = metrics["value"].astype(float).tolist()
    ax1.bar(metrics["metric"], [max(v, 1.0) for v in vals], color="#e6e6e6")
    ax1.bar(metrics["metric"], vals, color=["#1f77b4", "#2ca02c"])
    ax1.set_title("QTL Report Key Metrics")
    ax1.set_ylabel("Value")
    ax1.grid(axis="y", alpha=0.2)
    for i, v in enumerate(vals):
        ax1.text(i, v + 0.03, f"{v:.3g}", ha="center", fontsize=9)

    ax2.axis("off")
    ax2.text(
        0.02,
        0.95,
        f"Bonferroni hits: {qtl_count}\nLambda after GC: {lambda_after}\nData source: committed example outputs",
        va="top",
        fontsize=10,
        bbox={
            "facecolor": "#f8f8f8",
            "edgecolor": "#cccccc",
            "boxstyle": "round,pad=0.4",
        },
    )

    fig.tight_layout()
    fig.savefig(out / "analysis_report_metrics.png", dpi=150)
    plt.close(fig)

    (out / "analysis_report.html").write_text(html.strip() + "\n", encoding="utf-8")
    (out / "analysis_report.md").write_text(md, encoding="utf-8")
    print(
        "Saved analysis_report.html, analysis_report.md, and analysis_report_metrics.png"
    )


if __name__ == "__main__":
    main()
