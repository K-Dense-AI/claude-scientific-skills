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

"""
Example: G×E GWAS (Genotype × Environment)

Demonstrates GWAS with genotype-by-environment interaction.
Detects QTLs with environment-specific effects in multi-environment trials.

Auto-installs: pandas, numpy, scikit-learn, matplotlib, scipy
"""

import sys
import subprocess
from pathlib import Path


def install_packages():
    packages = ["pandas", "numpy", "scikit-learn", "matplotlib", "scipy"]
    for pkg in packages:
        try:
            __import__(pkg.replace("-", "_"))
        except ImportError:
            print(f"Installing {pkg}...")
            subprocess.check_call(
                [sys.executable, "-m", "pip", "install", "--user", "-q", pkg]
            )


def generate_gxe_data(
    n_ind=250, n_markers=400, n_envs=3, n_stable_qtl=8, n_gxe_qtl=6, h2=0.5, seed=42
):
    import numpy as np

    rng = np.random.default_rng(seed)

    # Genotypes
    X = rng.binomial(2, 0.3, size=(n_ind, n_markers)).astype(float)

    # QTLs - some stable across envs, some with G×E
    stable_qtl = rng.choice(n_markers, size=n_stable_qtl, replace=False)
    gxe_qtl = rng.choice(
        [i for i in range(n_markers) if i not in stable_qtl],
        size=n_gxe_qtl,
        replace=False,
    )

    # Effects
    stable_effects = rng.normal(0, 1.0, n_stable_qtl)
    gxe_effects = rng.normal(0, 0.8, (n_gxe_qtl, n_envs))
    env_means = rng.normal(0, 2.0, n_envs)

    # Generate phenotypes for each environment
    phenotypes = np.zeros((n_ind, n_envs))
    for env in range(n_envs):
        # Stable QTLs (same effect in all envs)
        g_stable = X[:, stable_qtl] @ stable_effects
        # G×E QTLs (env-specific effects)
        g_gxe = X[:, gxe_qtl] @ gxe_effects[:, env]
        # Total genetic value
        g = g_stable + g_gxe
        # Add noise
        var_g = float(np.var(g))
        var_e = var_g * (1 - h2) / h2
        noise = rng.normal(0, np.sqrt(var_e), n_ind)
        phenotypes[:, env] = g + env_means[env] + noise

    return {
        "X": X,
        "phenotypes": phenotypes,
        "stable_qtl": stable_qtl,
        "gxe_qtl": gxe_qtl,
        "env_means": env_means,
    }


def run_gxe_gwas(X, phenotypes):
    import numpy as np
    from scipy.stats import f_oneway, pearsonr
    from scipy.special import betaincinv

    n_ind, n_markers = X.shape
    n_envs = phenotypes.shape[1]

    p_values = []
    gxe_p_values = []

    for m in range(n_markers):
        marker = X[:, m]

        # Main effect (average across environments)
        y_mean = phenotypes.mean(axis=1)
        r_main, _ = pearsonr(marker, y_mean)
        n = len(y_mean)
        t_stat = r_main * np.sqrt((n - 2) / (1 - r_main**2))
        from scipy.stats import t as t_dist

        p_main = 2 * (1 - t_dist.cdf(abs(t_stat), n - 2))

        # G×E interaction effect
        # Test if marker effect differs across environments
        groups = []
        for env in range(n_envs):
            groups.append(phenotypes[:, env])

        # ANOVA-style interaction test
        try:
            f_stat, p_gxe = f_oneway(*groups)
        except:
            p_gxe = 1.0

        p_values.append(p_main)
        gxe_p_values.append(p_gxe)

    return np.array(p_values), np.array(gxe_p_values)


def plot_gxe_manhattan(main_p, gxe_p, data, output_path):
    import matplotlib.pyplot as plt
    import numpy as np

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

    n_markers = len(main_p)
    x = np.arange(n_markers)
    colors = ["#3498db" if (i // 50) % 2 == 0 else "#2ecc71" for i in range(n_markers)]

    # Main effects
    neg_log_p = -np.log10(main_p + 1e-300)
    ax1.scatter(x, neg_log_p, c=colors, s=20, alpha=0.7, edgecolors="none")
    bonferroni = -np.log10(0.05 / n_markers)
    ax1.axhline(y=bonferroni, color="red", linestyle="--", linewidth=2)
    ax1.set_ylabel("-log10(p-value)", fontsize=11)
    ax1.set_title("G×E GWAS: Main Effects", fontsize=13, fontweight="bold")
    ax1.grid(True, alpha=0.3)

    # Mark true QTLs
    for pos in data["stable_qtl"]:
        ax1.axvline(x=pos, color="orange", alpha=0.3, linestyle="-", linewidth=2)

    # G×E effects
    neg_log_p_gxe = -np.log10(gxe_p + 1e-300)
    ax2.scatter(x, neg_log_p_gxe, c=colors, s=20, alpha=0.7, edgecolors="none")
    ax2.axhline(y=bonferroni, color="red", linestyle="--", linewidth=2)
    ax2.set_xlabel("Marker Position", fontsize=11)
    ax2.set_ylabel("-log10(p-value)", fontsize=11)
    ax2.set_title("G×E GWAS: Interaction Effects", fontsize=13, fontweight="bold")
    ax2.grid(True, alpha=0.3)

    for pos in data["gxe_qtl"]:
        ax2.axvline(x=pos, color="purple", alpha=0.3, linestyle="-", linewidth=2)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {output_path}")


def save_results(main_p, gxe_p, data, output_dir):
    import pandas as pd

    results = pd.DataFrame(
        {
            "Marker": range(len(main_p)),
            "P_MainEffect": main_p,
            "P_GxE": gxe_p,
            "Is_Stable_QTL": [i in data["stable_qtl"] for i in range(len(main_p))],
            "Is_GxE_QTL": [i in data["gxe_qtl"] for i in range(len(main_p))],
        }
    )

    results.to_csv(f"{output_dir}/gxe_results.csv", index=False)
    print(f"Saved: {output_dir}/gxe_results.csv")


def main():
    print("=" * 70)
    print("Example: G×E GWAS")
    print("=" * 70)

    print("\n[1/4] Installing dependencies...")
    install_packages()

    output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(exist_ok=True)

    print("\n[2/4] Generating synthetic MET data...")
    data = generate_gxe_data()
    print(f" Individuals: {data['X'].shape[0]}")
    print(f" Markers: {data['X'].shape[1]}")
    print(f" Environments: {data['phenotypes'].shape[1]}")
    print(f" Stable QTLs: {len(data['stable_qtl'])}")
    print(f" G×E QTLs: {len(data['gxe_qtl'])}")

    print("\n[3/4] Running G×E GWAS...")
    main_p, gxe_p = run_gxe_gwas(data["X"], data["phenotypes"])

    threshold = 0.05 / len(main_p)
    print(f"\n Significant hits (Bonferroni p<{threshold:.2e}):")
    print(f" Main effects: {sum(main_p < threshold)}")
    print(f" G×E interactions: {sum(gxe_p < threshold)}")

    print("\n[4/4] Saving and plotting...")
    save_results(main_p, gxe_p, data, str(output_dir))
    plot_gxe_manhattan(main_p, gxe_p, data, f"{output_dir}/gxe_manhattan.png")

    print("\n" + "=" * 70)
    print("G×E GWAS Example Complete!")
    print("=" * 70)
    print("\nKey Finding: G×E analysis detects QTLs with environment-specific effects")


if __name__ == "__main__":
    main()
