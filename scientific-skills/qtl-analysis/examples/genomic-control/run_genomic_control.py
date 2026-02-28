#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

"""Genomic control example with lambda adjustment and QQ plots."""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import chi2


def simulate_pvalues(n_markers: int = 5000, inflation: float = 1.35, seed: int = 42):
    rng = np.random.default_rng(seed)
    null_chi = rng.chisquare(df=1, size=n_markers)
    inflated_chi = null_chi * inflation
    pvals = 1 - chi2.cdf(inflated_chi, df=1)
    return pvals, inflated_chi


def lambda_gc(chi_stats: np.ndarray) -> float:
    return float(np.median(chi_stats) / 0.456)


def adjust_with_gc(chi_stats: np.ndarray, lam: float):
    adj_chi = chi_stats / lam
    adj_p = 1 - chi2.cdf(adj_chi, df=1)
    return adj_p, adj_chi


def qq_points(pvals: np.ndarray):
    n = len(pvals)
    expected = -np.log10(np.linspace(1 / (n + 1), n / (n + 1), n))
    observed = -np.log10(np.sort(np.clip(pvals, 1e-300, 1)))
    return expected, observed


def plot_qq(before_p: np.ndarray, after_p: np.ndarray, out_path: Path):
    exp_b, obs_b = qq_points(before_p)
    exp_a, obs_a = qq_points(after_p)
    max_x = max(np.max(exp_b), np.max(exp_a))
    max_y = max(np.max(obs_b), np.max(obs_a))

    plt.figure(figsize=(8, 6))
    plt.scatter(exp_b, obs_b, s=10, alpha=0.45, label="Before GC", color="#e74c3c")
    plt.scatter(exp_a, obs_a, s=10, alpha=0.45, label="After GC", color="#2ecc71")
    plt.plot([0, max_x], [0, max_x], "k--", linewidth=1.5, label="Expected")
    plt.xlim(0, max_x * 1.05)
    plt.ylim(0, max(max_x, max_y) * 1.05)
    plt.xlabel("Expected -log10(p)")
    plt.ylabel("Observed -log10(p)")
    plt.title("Genomic Control Adjustment")
    plt.legend()
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()


def main():
    out_dir = Path(__file__).parent / "output"
    out_dir.mkdir(exist_ok=True)

    before_p, chi_stats = simulate_pvalues()
    lam = lambda_gc(chi_stats)
    after_p, adj_chi = adjust_with_gc(chi_stats, lam)
    lam_after = lambda_gc(adj_chi)

    plot_qq(before_p, after_p, out_dir / "qq_before_after_gc.png")

    pd.DataFrame(
        {
            "marker": np.arange(len(before_p)),
            "p_before": before_p,
            "p_after": after_p,
            "chi_before": chi_stats,
            "chi_after": adj_chi,
        }
    ).to_csv(out_dir / "genomic_control_results.csv", index=False)

    pd.DataFrame(
        {
            "metric": ["lambda_before", "lambda_after"],
            "value": [lam, lam_after],
        }
    ).to_csv(out_dir / "lambda_summary.csv", index=False)

    print(f"lambda before: {lam:.3f}")
    print(f"lambda after:  {lam_after:.3f}")
    print("Saved: output/qq_before_after_gc.png")
    print("Saved: output/genomic_control_results.csv")
    print("Saved: output/lambda_summary.csv")


if __name__ == "__main__":
    main()
