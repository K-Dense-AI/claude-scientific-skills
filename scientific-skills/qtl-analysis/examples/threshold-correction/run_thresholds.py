#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under Apache License 2.0

"""
Example: Multiple Testing Correction Comparison

Compares Bonferroni, Benjamini-Hochberg FDR, and permutation-based thresholds.
Shows tradeoffs between false positive control and statistical power.

Auto-installs: pandas, numpy, matplotlib, scipy
"""

import sys
import subprocess
from pathlib import Path


def install_packages():
    packages = ["pandas", "numpy", "matplotlib", "scipy"]
    for pkg in packages:
        try:
            __import__(pkg.replace("-", "_"))
        except ImportError:
            print(f"Installing {pkg}...")
            subprocess.check_call(
                [sys.executable, "-m", "pip", "install", "--user", "-q", pkg]
            )


def simulate_gwas_pvalues(n_markers=10000, n_qtl=50, h2=0.5, seed=42):
    import numpy as np
    from scipy.stats import beta, uniform

    rng = np.random.default_rng(seed)

    # Null markers (no association) - uniform p-values
    null_p = rng.uniform(0, 1, n_markers - n_qtl)

    # True QTLs (some association) - biased toward small p-values
    # Beta distribution skewed toward 0
    qtl_p = beta.rvs(0.3, 10, size=n_qtl, random_state=rng)

    # Combine and shuffle
    all_p = np.concatenate([null_p, qtl_p])
    rng.shuffle(all_p)

    # Track which are QTLs
    is_qtl = np.zeros(n_markers, dtype=bool)
    is_qtl[:n_qtl] = True
    rng.shuffle(is_qtl)

    return all_p, is_qtl


def bonferroni_correction(p_values, alpha=0.05):
    """Bonferroni: p < alpha/m"""
    n = len(p_values)
    threshold = alpha / n
    rejected = p_values < threshold
    return rejected, threshold


def benjamini_hochberg(p_values, fdr=0.05):
    """Benjamini-Hochberg FDR control"""
    import numpy as np

    n = len(p_values)
    sorted_idx = np.argsort(p_values)
    sorted_p = p_values[sorted_idx]

    # Find largest k such that p(k) <= (k/n) * fdr
    thresholds = np.arange(1, n + 1) / n * fdr
    rejected_sorted = sorted_p <= thresholds

    # Find largest k where condition holds
    if np.any(rejected_sorted):
        k = np.where(rejected_sorted)[0][-1]
        threshold = sorted_p[k]
        rejected = p_values <= threshold
    else:
        threshold = 0
        rejected = np.zeros(n, dtype=bool)

    return rejected, threshold


def permutation_threshold(p_values, n_perm=100, alpha=0.05, seed=42):
    """Permutation-based threshold"""
    import numpy as np
    from scipy.stats import percentileofscore

    rng = np.random.default_rng(seed)
    n = len(p_values)

    # Simulate null distribution by permuting
    min_p_perm = []
    for _ in range(n_perm):
        perm_p = rng.uniform(0, 1, n)
        min_p_perm.append(np.min(perm_p))

    # Threshold is alpha percentile of null min p-values
    threshold = np.percentile(min_p_perm, alpha * 100)
    rejected = p_values < threshold

    return rejected, threshold


def plot_comparison(p_values, bonf_rej, bh_rej, perm_rej, is_qtl, output_dir):
    import matplotlib.pyplot as plt
    import numpy as np

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    n = len(p_values)
    sorted_idx = np.argsort(p_values)

    # Plot 1: P-value distribution
    ax = axes[0, 0]
    ax.hist(
        p_values[~is_qtl], bins=50, alpha=0.7, label="Null markers", color="#3498db"
    )
    ax.hist(p_values[is_qtl], bins=20, alpha=0.7, label="True QTLs", color="#e74c3c")
    ax.set_xlabel("P-value", fontsize=11)
    ax.set_ylabel("Count", fontsize=11)
    ax.set_title("P-value Distribution", fontsize=12, fontweight="bold")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: QQ plot
    ax = axes[0, 1]
    from scipy.stats import uniform

    theoretical = -np.log10(np.linspace(1 / n, 1, n))
    observed = -np.log10(np.sort(p_values))
    ax.scatter(theoretical, observed, alpha=0.3, s=10)
    ax.plot([0, max(theoretical)], [0, max(theoretical)], "r--", label="Expected")
    ax.set_xlabel("Expected -log10(p)", fontsize=11)
    ax.set_ylabel("Observed -log10(p)", fontsize=11)
    ax.set_title("QQ Plot", fontsize=12, fontweight="bold")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 3: Manhattan-style with thresholds
    ax = axes[1, 0]
    x = np.arange(n)
    colors = ["#3498db" if not q else "#e74c3c" for q in is_qtl]
    ax.scatter(x, -np.log10(p_values + 1e-300), c=colors, s=8, alpha=0.5)

    # Add threshold lines
    bonf_thresh = -np.log10(0.05 / n)
    bh_thresh = -np.log10(max(p_values[bh_rej]) if np.any(bh_rej) else 1 / n)
    perm_thresh = -np.log10(max(p_values[perm_rej]) if np.any(perm_rej) else 1)

    ax.axhline(
        y=bonf_thresh,
        color="red",
        linestyle="--",
        linewidth=2,
        label=f"Bonferroni: {bonf_thresh:.2f}",
    )
    ax.axhline(
        y=bh_thresh,
        color="green",
        linestyle="--",
        linewidth=2,
        label=f"BH FDR: {bh_thresh:.2f}",
    )
    ax.axhline(
        y=perm_thresh,
        color="purple",
        linestyle="--",
        linewidth=2,
        label=f"Permutation: {perm_thresh:.2f}",
    )

    ax.set_xlabel("Marker Rank", fontsize=11)
    ax.set_ylabel("-log10(p-value)", fontsize=11)
    ax.set_title(
        "Multiple Testing Correction Thresholds", fontsize=12, fontweight="bold"
    )
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 4: Performance comparison
    ax = axes[1, 1]

    methods = ["Bonferroni", "BH FDR", "Permutation"]
    true_positives = [
        np.sum(bonf_rej & is_qtl),
        np.sum(bh_rej & is_qtl),
        np.sum(perm_rej & is_qtl),
    ]
    false_positives = [
        np.sum(bonf_rej & ~is_qtl),
        np.sum(bh_rej & ~is_qtl),
        np.sum(perm_rej & ~is_qtl),
    ]

    x = np.arange(len(methods))
    width = 0.35

    ax.bar(
        x - width / 2, true_positives, width, label="True Positives", color="#2ecc71"
    )
    ax.bar(
        x + width / 2, false_positives, width, label="False Positives", color="#e74c3c"
    )

    ax.set_ylabel("Count", fontsize=11)
    ax.set_title("Detection Performance", fontsize=12, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels(methods)
    ax.legend()
    ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    plt.savefig(f"{output_dir}/threshold_comparison.png", dpi=150, bbox_inches="tight")
    print(f"Saved: {output_dir}/threshold_comparison.png")


def save_results(p_values, bonf_rej, bh_rej, perm_rej, is_qtl, output_dir):
    import pandas as pd
    import numpy as np

    results = pd.DataFrame(
        {
            "Marker": range(len(p_values)),
            "P_Value": p_values,
            "Is_QTL": is_qtl,
            "Bonferroni_Significant": bonf_rej,
            "BH_FDR_Significant": bh_rej,
            "Permutation_Significant": perm_rej,
        }
    )

    results.to_csv(f"{output_dir}/threshold_results.csv", index=False)

    # Summary statistics
    summary = pd.DataFrame(
        {
            "Method": ["Bonferroni", "BH FDR", "Permutation"],
            "True_Positives": [
                np.sum(bonf_rej & is_qtl),
                np.sum(bh_rej & is_qtl),
                np.sum(perm_rej & is_qtl),
            ],
            "False_Positives": [
                np.sum(bonf_rej & ~is_qtl),
                np.sum(bh_rej & ~is_qtl),
                np.sum(perm_rej & ~is_qtl),
            ],
            "True_Negatives": [
                np.sum(~bonf_rej & ~is_qtl),
                np.sum(~bh_rej & ~is_qtl),
                np.sum(~perm_rej & ~is_qtl),
            ],
            "False_Negatives": [
                np.sum(~bonf_rej & is_qtl),
                np.sum(~bh_rej & is_qtl),
                np.sum(~perm_rej & is_qtl),
            ],
            "Power": [
                np.sum(bonf_rej & is_qtl) / np.sum(is_qtl),
                np.sum(bh_rej & is_qtl) / np.sum(is_qtl),
                np.sum(perm_rej & is_qtl) / np.sum(is_qtl),
            ],
            "FDR": [
                np.sum(bonf_rej & ~is_qtl) / max(np.sum(bonf_rej), 1),
                np.sum(bh_rej & ~is_qtl) / max(np.sum(bh_rej), 1),
                np.sum(perm_rej & ~is_qtl) / max(np.sum(perm_rej), 1),
            ],
        }
    )

    summary.to_csv(f"{output_dir}/threshold_summary.csv", index=False)
    print(f"Saved: {output_dir}/threshold_summary.csv")

    return summary


def main():
    import numpy as np

    print("=" * 70)
    print("Example: Multiple Testing Correction Comparison")
    print("=" * 70)

    print("\n[1/3] Installing dependencies...")
    install_packages()

    output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(exist_ok=True)

    print("\n[2/3] Simulating GWAS p-values...")
    p_values, is_qtl = simulate_gwas_pvalues()
    print(f" Total markers: {len(p_values)}")
    print(f" True QTLs: {np.sum(is_qtl)}")

    print("\n[3/3] Applying multiple testing corrections...")

    # Bonferroni
    bonf_rej, bonf_thresh = bonferroni_correction(p_values)
    print(f"\n Bonferroni threshold: {bonf_thresh:.2e}")
    print(f"   Significant: {np.sum(bonf_rej)}")
    print(f"   True positives: {np.sum(bonf_rej & is_qtl)}")
    print(f"   False positives: {np.sum(bonf_rej & ~is_qtl)}")

    # Benjamini-Hochberg
    bh_rej, bh_thresh = benjamini_hochberg(p_values)
    print(f"\n BH FDR threshold: {bh_thresh:.2e}")
    print(f"   Significant: {np.sum(bh_rej)}")
    print(f"   True positives: {np.sum(bh_rej & is_qtl)}")
    print(f"   False positives: {np.sum(bh_rej & ~is_qtl)}")

    # Permutation
    perm_rej, perm_thresh = permutation_threshold(p_values)
    print(f"\n Permutation threshold: {perm_thresh:.2e}")
    print(f"   Significant: {np.sum(perm_rej)}")
    print(f"   True positives: {np.sum(perm_rej & is_qtl)}")
    print(f"   False positives: {np.sum(perm_rej & ~is_qtl)}")

    print("\n[4/4] Saving results and creating plots...")
    summary = save_results(
        p_values, bonf_rej, bh_rej, perm_rej, is_qtl, str(output_dir)
    )
    plot_comparison(p_values, bonf_rej, bh_rej, perm_rej, is_qtl, str(output_dir))

    print("\n" + "=" * 70)
    print("Threshold Correction Example Complete!")
    print("=" * 70)
    print("\nSummary:")
    print(summary.to_string(index=False))
    print("\nKey Finding: BH FDR offers better power than Bonferroni")
    print("while controlling false discovery rate")


if __name__ == "__main__":
    main()
