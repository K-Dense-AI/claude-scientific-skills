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

"""
Example: Covariate-Adjusted GWAS

Demonstrates how to include covariates (population structure, sex, batch) in GWAS
to avoid confounding and reduce false positives.

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


def generate_structured_data(n_ind=300, n_markers=500, n_qtl=15, h2=0.5, seed=42):
    import numpy as np
    from sklearn.decomposition import PCA

    rng = np.random.default_rng(seed)

    # Create 3 subpopulations with different allele frequencies
    n_per_pop = n_ind // 3
    subpops = np.repeat([0, 1, 2], n_per_pop)
    if len(subpops) < n_ind:
        subpops = np.concatenate([subpops, np.zeros(n_ind - len(subpops))])
    rng.shuffle(subpops)

    # Generate genotypes with population structure
    X = np.zeros((n_ind, n_markers))
    base_freqs = rng.uniform(0.1, 0.9, n_markers)

    for i in range(n_ind):
        pop = subpops[i]
        # Shift frequencies by population
        pop_shift = [-0.15, 0, 0.15][int(pop)]
        freqs = np.clip(base_freqs + pop_shift, 0.05, 0.95)
        X[i] = rng.binomial(2, freqs)

    # Calculate PCs for population structure
    pca = PCA(n_components=3)
    pcs = pca.fit_transform(X)

    # Generate covariates
    sex = rng.binomial(1, 0.5, n_ind)  # 0=female, 1=male
    batch = rng.choice([0, 1, 2], n_ind)  # 3 batches

    # Generate phenotype with structure and covariate effects
    qtl_idx = rng.choice(n_markers, size=n_qtl, replace=False)
    qtl_effects = rng.normal(0, 1.0, n_qtl)

    # Genetic value (includes pop structure because QTLs are correlated with structure)
    g = X[:, qtl_idx] @ qtl_effects

    # Add population structure effect (confounder)
    pop_effects = np.array([5, 0, -5])  # Strong population effect
    pop_contrib = pop_effects[subpops.astype(int)]

    # Add covariate effects
    sex_contrib = sex * 2.0  # Sex effect
    batch_contrib = np.array([-1, 0, 1])[batch]  # Batch effect

    # Add noise
    var_g = float(np.var(g))
    var_e = var_g * (1 - h2) / h2
    noise = rng.normal(0, np.sqrt(var_e), n_ind)

    # Total phenotype (includes confounding population structure)
    phenotype = g + pop_contrib + sex_contrib + batch_contrib + noise

    return {
        "X": X,
        "phenotype": phenotype,
        "sex": sex,
        "batch": batch,
        "pcs": pcs,
        "subpops": subpops,
        "qtl_idx": qtl_idx,
    }


def run_gwas_simple(X, y, n_markers):
    import numpy as np
    from scipy.stats import pearsonr, t as t_dist

    p_values = []
    for m in range(n_markers):
        r, _ = pearsonr(X[:, m], y)
        n = len(y)
        if abs(r) < 1:
            t_stat = r * np.sqrt((n - 2) / (1 - r**2))
            p = 2 * (1 - t_dist.cdf(abs(t_stat), n - 2))
        else:
            p = 1e-300
        p_values.append(p)
    return np.array(p_values)


def run_gwas_covariates(X, y, covariates, n_markers):
    import numpy as np
    from scipy.stats import pearsonr, t as t_dist
    from sklearn.linear_model import LinearRegression

    # Residualize phenotype and each marker for covariates
    reg_y = LinearRegression().fit(covariates, y)
    y_resid = y - reg_y.predict(covariates)

    p_values = []
    for m in range(n_markers):
        reg_m = LinearRegression().fit(covariates, X[:, m])
        m_resid = X[:, m] - reg_m.predict(covariates)

        # Correlation of residuals
        r, _ = pearsonr(m_resid, y_resid)
        n = len(y_resid)
        if abs(r) < 1:
            t_stat = r * np.sqrt((n - 2) / (1 - r**2))
            p = 2 * (1 - t_dist.cdf(abs(t_stat), n - 2))
        else:
            p = 1e-300
        p_values.append(p)

    return np.array(p_values)


def plot_comparison(p_simple, p_cov, data, output_dir):
    import matplotlib.pyplot as plt
    import numpy as np

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8))

    n_markers = len(p_simple)
    x = np.arange(n_markers)
    colors = ["#3498db" if (i // 50) % 2 == 0 else "#2ecc71" for i in range(n_markers)]
    bonferroni = -np.log10(0.05 / n_markers)

    # Plot 1: Simple GWAS (confounded)
    neg_log_p = -np.log10(p_simple + 1e-300)
    ax1.scatter(x, neg_log_p, c=colors, s=15, alpha=0.6, edgecolors="none")
    ax1.axhline(
        y=bonferroni,
        color="red",
        linestyle="--",
        linewidth=2,
        label=f"Bonferroni (p={0.05 / n_markers:.2e})",
    )
    ax1.set_ylabel("-log10(p-value)", fontsize=11)
    ax1.set_title(
        "GWAS WITHOUT Covariates: Confounded by Population Structure",
        fontsize=12,
        fontweight="bold",
        color="darkred",
    )
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, max(neg_log_p) * 1.1)

    # Mark true QTLs
    for pos in data["qtl_idx"]:
        ax1.axvline(x=pos, color="orange", alpha=0.3, linestyle="-", linewidth=2)

    # Plot 2: Covariate-adjusted GWAS
    neg_log_p_cov = -np.log10(p_cov + 1e-300)
    ax2.scatter(x, neg_log_p_cov, c=colors, s=15, alpha=0.6, edgecolors="none")
    ax2.axhline(
        y=bonferroni,
        color="red",
        linestyle="--",
        linewidth=2,
        label=f"Bonferroni (p={0.05 / n_markers:.2e})",
    )
    ax2.set_xlabel("Marker Position", fontsize=11)
    ax2.set_ylabel("-log10(p-value)", fontsize=11)
    ax2.set_title(
        "GWAS WITH Covariates (PCs + Sex + Batch): Corrected",
        fontsize=12,
        fontweight="bold",
        color="darkgreen",
    )
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0, max(neg_log_p) * 1.1)

    for pos in data["qtl_idx"]:
        ax2.axvline(x=pos, color="orange", alpha=0.3, linestyle="-", linewidth=2)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/covariate_comparison.png", dpi=150, bbox_inches="tight")
    print(f"Saved: {output_dir}/covariate_comparison.png")


def save_results(p_simple, p_cov, data, output_dir):
    import pandas as pd
    import numpy as np

    results = pd.DataFrame(
        {
            "Marker": range(len(p_simple)),
            "P_No_Covariates": p_simple,
            "P_With_Covariates": p_cov,
            "LogP_Difference": -np.log10(p_simple + 1e-300) - -np.log10(p_cov + 1e-300),
            "Is_QTL": [i in data["qtl_idx"] for i in range(len(p_simple))],
        }
    )

    results.to_csv(f"{output_dir}/covariate_results.csv", index=False)
    print(f"Saved: {output_dir}/covariate_results.csv")


def main():
    import numpy as np

    print("=" * 70)
    print("Example: Covariate-Adjusted GWAS")
    print("=" * 70)

    print("\n[1/5] Installing dependencies...")
    install_packages()

    output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(exist_ok=True)

    print("\n[2/5] Generating structured data with population stratification...")
    data = generate_structured_data()
    print(f" Individuals: {data['X'].shape[0]}")
    print(f" Markers: {data['X'].shape[1]}")
    print(f" QTLs: {len(data['qtl_idx'])}")
    print(f" Subpopulations: {len(np.unique(data['subpops']))}")
    print(f" Sex balance: {np.sum(data['sex'] == 1)}M / {np.sum(data['sex'] == 0)}F")

    print("\n[3/5] Running simple GWAS (no covariates)...")
    p_simple = run_gwas_simple(data["X"], data["phenotype"], data["X"].shape[1])

    print("\n[4/5] Running covariate-adjusted GWAS...")
    covariates = np.column_stack([data["pcs"], data["sex"], data["batch"]])
    p_cov = run_gwas_covariates(
        data["X"], data["phenotype"], covariates, data["X"].shape[1]
    )

    # Compare
    threshold = 0.05 / len(p_simple)
    print(f"\n Significant hits (p < {threshold:.2e}):")
    print(
        f"  Without covariates: {np.sum(p_simple < threshold)} (includes false positives)"
    )
    print(f"  With covariates: {np.sum(p_cov < threshold)} (better calibrated)")
    print(
        f"  True QTLs detected: {np.sum(p_cov[data['qtl_idx']] < threshold)}/{len(data['qtl_idx'])}"
    )

    print("\n[5/5] Saving and plotting...")
    save_results(p_simple, p_cov, data, str(output_dir))
    plot_comparison(p_simple, p_cov, data, str(output_dir))

    print("\n" + "=" * 70)
    print("Covariate GWAS Example Complete!")
    print("=" * 70)
    print("\nKey Finding: Covariate adjustment reduces false positives from")
    print("population structure while maintaining power to detect true QTLs")


if __name__ == "__main__":
    main()
