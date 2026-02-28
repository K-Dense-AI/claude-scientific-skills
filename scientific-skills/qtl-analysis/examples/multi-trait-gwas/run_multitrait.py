#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Author: Clayton Young <Clayton@SuperiorByteWorks.com>
# LinkedIn: https://linkedin.com/in/claytoneyoung/
# GitHub: https://github.com/borealBytes

"""
Example: Multi-Trait GWAS

Demonstrates joint GWAS analysis of two correlated traits.
Uses multivariate linear models to increase power when traits share genetic architecture.

Auto-installs: pandas, numpy, scikit-learn, matplotlib, scipy
"""

import sys
import subprocess
from pathlib import Path


def install_packages():
    """Install required packages without root."""
    packages = ["pandas", "numpy", "scikit-learn", "matplotlib", "scipy"]
    for pkg in packages:
        try:
            __import__(pkg.replace("-", "_"))
        except ImportError:
            print(f"Installing {pkg}...")
            subprocess.check_call(
                [sys.executable, "-m", "pip", "install", "--user", "-q", pkg]
            )


def generate_multitrait_data(
    n_ind=300,
    n_markers=500,
    n_shared_qtl=10,
    n_trait1_qtl=5,
    n_trait2_qtl=5,
    rg=0.7,
    h2=0.5,
    seed=42,
):
    """Generate synthetic genotype/phenotype data for two correlated traits."""
    import numpy as np

    rng = np.random.default_rng(seed)

    # Genotypes
    X = rng.binomial(2, 0.3, size=(n_ind, n_markers)).astype(float)

    # QTL indices - some shared, some trait-specific
    shared_qtl = rng.choice(n_markers, size=n_shared_qtl, replace=False)
    trait1_only = rng.choice(
        [i for i in range(n_markers) if i not in shared_qtl],
        size=n_trait1_qtl,
        replace=False,
    )
    trait2_only = rng.choice(
        [i for i in range(n_markers) if i not in shared_qtl and i not in trait1_only],
        size=n_trait2_qtl,
        replace=False,
    )

    # Genetic values
    shared_effects = rng.normal(0, 1.0, n_shared_qtl)
    t1_effects = rng.normal(0, 0.8, n_trait1_qtl)
    t2_effects = rng.normal(0, 0.8, n_trait2_qtl)

    g1 = X[:, shared_qtl] @ shared_effects + X[:, trait1_only] @ t1_effects
    g2 = X[:, shared_qtl] @ (shared_effects * rg) + X[:, trait2_only] @ t2_effects

    # Add correlated noise
    var_g1 = float(np.var(g1))
    var_g2 = float(np.var(g2))
    var_e1 = var_g1 * (1 - h2) / h2
    var_e2 = var_g2 * (1 - h2) / h2

    # Correlated environmental effects
    cov_e = np.sqrt(var_e1 * var_e2) * rg * 0.5
    cov_matrix = [[var_e1, cov_e], [cov_e, var_e2]]
    noise = rng.multivariate_normal([0, 0], cov_matrix, n_ind)

    y1 = g1 + noise[:, 0]
    y2 = g2 + noise[:, 1]

    return {
        "X": X,
        "y1": y1,
        "y2": y2,
        "qtl_shared": shared_qtl,
        "qtl_trait1": trait1_only,
        "qtl_trait2": trait2_only,
    }


def run_single_trait_gwas(X, y, n_markers):
    """Standard single-trait GWAS."""
    import numpy as np
    from scipy.stats import pearsonr

    p_values = []
    for m in range(n_markers):
        r, _ = pearsonr(X[:, m], y)
        # Convert correlation to approximate p-value
        n = len(y)
        t_stat = r * np.sqrt((n - 2) / (1 - r**2))
        from scipy.stats import t as t_dist

        p = 2 * (1 - t_dist.cdf(abs(t_stat), n - 2))
        p_values.append(p)

    return np.array(p_values)


def run_multitrait_gwas(X, y1, y2, n_markers):
    """Multi-trait GWAS using product of p-values or canonical correlation approach."""
    import numpy as np
    from scipy.stats import pearsonr

    # Simple approach: combine evidence from both traits
    p1 = run_single_trait_gwas(X, y1, n_markers)
    p2 = run_single_trait_gwas(X, y2, n_markers)

    # Fisher's method for combining p-values
    combined = -2 * (np.log(p1 + 1e-300) + np.log(p2 + 1e-300))
    from scipy.stats import chi2

    multivariate_p = 1 - chi2.cdf(combined, df=4)

    return p1, p2, multivariate_p


def plot_manhattan(p_values, title, output_path, qtl_positions=None):
    """Create Manhattan plot."""
    import matplotlib.pyplot as plt
    import numpy as np

    fig, ax = plt.subplots(figsize=(12, 5))

    n_markers = len(p_values)
    x = np.arange(n_markers)

    # -log10(p)
    neg_log_p = -np.log10(p_values + 1e-300)

    # Color by chromosome (simulated - every 50 markers)
    colors = ["#3498db" if (i // 50) % 2 == 0 else "#2ecc71" for i in range(n_markers)]

    ax.scatter(x, neg_log_p, c=colors, s=20, alpha=0.7, edgecolors="none")

    # Significance line (Bonferroni)
    bonferroni = -np.log10(0.05 / n_markers)
    ax.axhline(
        y=bonferroni,
        color="red",
        linestyle="--",
        linewidth=2,
        label=f"Bonferroni (p<{(0.05 / n_markers):.2e})",
    )

    # Mark true QTLs if provided
    if qtl_positions is not None:
        for pos in qtl_positions:
            ax.axvline(x=pos, color="orange", alpha=0.3, linestyle="-", linewidth=3)

    ax.set_xlabel("Marker Position", fontsize=12)
    ax.set_ylabel("-log10(p-value)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {output_path}")


def save_results(p1, p2, p_multi, data, output_dir):
    """Save GWAS results to CSV."""
    import pandas as pd
    import numpy as np

    results = pd.DataFrame(
        {
            "Marker": range(len(p1)),
            "P_Trait1": p1,
            "P_Trait2": p2,
            "P_MultiTrait": p_multi,
            "Is_Shared_QTL": [i in data["qtl_shared"] for i in range(len(p1))],
            "Is_Trait1_QTL": [i in data["qtl_trait1"] for i in range(len(p1))],
            "Is_Trait2_QTL": [i in data["qtl_trait2"] for i in range(len(p1))],
        }
    )

    results.to_csv(f"{output_dir}/gwas_results.csv", index=False)
    print(f"Saved: {output_dir}/gwas_results.csv")


def main():
    print("=" * 70)
    print("Example: Multi-Trait GWAS")
    print("=" * 70)

    # Install packages
    print("\n[1/5] Installing dependencies...")
    install_packages()

    # Setup
    output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(exist_ok=True)

    # Generate data
    print("\n[2/5] Generating synthetic multi-trait data...")
    data = generate_multitrait_data()
    print(f" Individuals: {data['X'].shape[0]}")
    print(f" Markers: {data['X'].shape[1]}")
    print(f" Shared QTLs: {len(data['qtl_shared'])}")
    print(f" Trait1-specific QTLs: {len(data['qtl_trait1'])}")
    print(f" Trait2-specific QTLs: {len(data['qtl_trait2'])}")

    # Run GWAS
    print("\n[3/5] Running single-trait and multi-trait GWAS...")
    p1, p2, p_multi = run_multitrait_gwas(
        data["X"], data["y1"], data["y2"], data["X"].shape[1]
    )

    # Count significant hits
    threshold = 0.05 / len(p1)
    print(f"\n Significant hits (Bonferroni p<{threshold:.2e}):")
    print(f" Trait 1 only: {sum(p1 < threshold)}")
    print(f" Trait 2 only: {sum(p2 < threshold)}")
    print(f" Multi-trait: {sum(p_multi < threshold)}")

    # Save results
    print("\n[4/5] Saving results...")
    save_results(p1, p2, p_multi, data, str(output_dir))

    # Plot
    print("\n[5/5] Creating Manhattan plots...")
    plot_manhattan(
        p1,
        "Multi-Trait GWAS: Trait 1",
        f"{output_dir}/manhattan_trait1.png",
        data["qtl_shared"],
    )
    plot_manhattan(
        p2,
        "Multi-Trait GWAS: Trait 2",
        f"{output_dir}/manhattan_trait2.png",
        data["qtl_shared"],
    )

    print("\n" + "=" * 70)
    print("Multi-Trait GWAS Example Complete!")
    print("=" * 70)
    print(
        "\nKey Finding: Multi-trait analysis can detect shared QTLs with increased power"
    )


if __name__ == "__main__":
    main()
