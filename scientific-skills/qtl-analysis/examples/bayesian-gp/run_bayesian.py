#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
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
Example: Bayesian Genomic Prediction

Demonstrates Bayesian genomic prediction methods: BayesA, BayesB, BayesCpi, and GBLUP.
Uses simplified Python implementations (no external R/BGLR dependency).

Equivalent to QTLmax: "Bayesian GP"
https://open.qtlmax.com/guide/index.php/2025/07/12/bayesian-gp/

Auto-installs: pandas, numpy, scikit-learn, matplotlib, scipy
"""

from __future__ import annotations

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


def generate_synthetic_data(
    n_individuals: int = 300,
    n_markers: int = 500,
    n_qtl: int = 30,
    heritability: float = 0.5,
    seed: int = 42,
) -> dict:
    """Generate synthetic genotype/phenotype data for Bayesian GP."""
    import numpy as np

    rng = np.random.default_rng(seed)

    # Generate genotypes (0/1/2, coded as additive)
    genotypes = rng.binomial(2, 0.3, size=(n_individuals, n_markers)).astype(float)

    # Select QTL indices
    qtl_idx = rng.choice(n_markers, size=n_qtl, replace=False)

    # True breeding values = sum of marker effects
    qtl_effects = rng.normal(0, 1.0, n_qtl)
    tbv = genotypes[:, qtl_idx] @ qtl_effects

    # Phenotypes = TBV + noise
    var_tbv = float(np.var(tbv))
    var_noise = var_tbv * (1 - heritability) / heritability
    phenotypes = tbv + rng.normal(0, np.sqrt(var_noise), n_individuals)

    # Build GRM for GBLUP (VanRaden method 1)
    p = genotypes.mean(axis=0) / 2.0
    denom = 2.0 * float(np.sum(p * (1 - p)))
    Z = (genotypes - 2.0 * p) / np.sqrt(denom)
    GRM = (Z @ Z.T) / n_markers

    return {
        "genotypes": genotypes,
        "phenotypes": phenotypes,
        "tbv": tbv,
        "GRM": GRM,
        "qtl_idx": qtl_idx,
    }


def bayes_a(X_train, y_train, X_test, n_iter: int = 100, burn_in: int = 20):
    """BayesA: marker-specific variance (t prior on effects)."""
    import numpy as np
    from sklearn.linear_model import Ridge

    # Simplified: Use Ridge with very small alpha as proxy for BayesA
    # True BayesA uses t-distribution priors with marker-specific variance
    model = Ridge(alpha=0.01)
    model.fit(X_train, y_train)
    return model.predict(X_test)


def bayes_b(X_train, y_train, X_test, n_iter: int = 100, pi: float = 0.9):
    """BayesB: spike-slab prior (some markers have zero effect)."""
    import numpy as np
    from sklearn.linear_model import LassoCV

    # Simplified: Use LassoCV as proxy for BayesB (sparsity)
    # True BayesB uses spike-slab mixture
    model = LassoCV(cv=3, random_state=42, max_iter=1000)
    model.fit(X_train, y_train)
    return model.predict(X_test)


def bayes_cpi(X_train, y_train, X_test, n_iter: int = 100, pi_init: float = 0.5):
    """BayesCpi: mixture with unknown proportion."""
    import numpy as np
    from sklearn.linear_model import ElasticNetCV

    # Simplified: Use ElasticNet as proxy for BayesCpi
    model = ElasticNetCV(cv=3, random_state=42, l1_ratio=0.5, max_iter=1000)
    model.fit(X_train, y_train)
    return model.predict(X_test)


def gblup_ridge(y_train, G_train, G_test_train, alpha: float = 1.0):
    """GBLUP via kernel ridge regression on GRM."""
    import numpy as np
    from sklearn.linear_model import Ridge

    n = y_train.shape[0]
    mu = np.mean(y_train)

    # GBLUP solution: alpha = (G + lambda*I)^-1 * (y - mu)
    A = G_train + alpha * np.eye(n)
    alpha_vec = np.linalg.solve(A, y_train - mu)

    # Predict: y_hat = mu + G_test_train @ alpha_vec
    y_pred = mu + G_test_train @ alpha_vec
    return y_pred


def run_cv_bayesian(data: dict, n_splits: int = 5) -> dict:
    """Run 5-fold CV for all Bayesian methods."""
    import numpy as np
    import pandas as pd
    from sklearn.model_selection import KFold
    from scipy.stats import pearsonr

    X = data["genotypes"]
    y = data["phenotypes"]
    tbv = data["tbv"]
    G = data["GRM"]

    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

    results = {
        method: {"correlation": [], "rmse": [], "r2": []}
        for method in ["BayesA", "BayesB", "BayesCpi", "GBLUP"]
    }

    for fold, (train_idx, test_idx) in enumerate(kf.split(X), 1):
        X_tr, X_te = X[train_idx], X[test_idx]
        y_tr, y_te = y[train_idx], y[test_idx]
        tbv_te = tbv[test_idx]

        # GBLUP needs GRM
        G_tr = G[np.ix_(train_idx, train_idx)]
        G_te_tr = G[np.ix_(test_idx, train_idx)]

        # BayesA
        pred = bayes_a(X_tr, y_tr, X_te)
        r, _ = pearsonr(pred, tbv_te)
        rmse = np.sqrt(np.mean((pred - tbv_te) ** 2))
        results["BayesA"]["correlation"].append(r)
        results["BayesA"]["rmse"].append(rmse)
        results["BayesA"]["r2"].append(r**2)

        # BayesB
        pred = bayes_b(X_tr, y_tr, X_te)
        r, _ = pearsonr(pred, tbv_te)
        rmse = np.sqrt(np.mean((pred - tbv_te) ** 2))
        results["BayesB"]["correlation"].append(r)
        results["BayesB"]["rmse"].append(rmse)
        results["BayesB"]["r2"].append(r**2)

        # BayesCpi
        pred = bayes_cpi(X_tr, y_tr, X_te)
        r, _ = pearsonr(pred, tbv_te)
        rmse = np.sqrt(np.mean((pred - tbv_te) ** 2))
        results["BayesCpi"]["correlation"].append(r)
        results["BayesCpi"]["rmse"].append(rmse)
        results["BayesCpi"]["r2"].append(r**2)

        # GBLUP
        pred = gblup_ridge(y_tr, G_tr, G_te_tr)
        r, _ = pearsonr(pred, tbv_te)
        rmse = np.sqrt(np.mean((pred - tbv_te) ** 2))
        results["GBLUP"]["correlation"].append(r)
        results["GBLUP"]["rmse"].append(rmse)
        results["GBLUP"]["r2"].append(r**2)

    # Summarize
    summary = {}
    for method in results:
        summary[method] = {
            "correlation_mean": np.mean(results[method]["correlation"]),
            "correlation_std": np.std(results[method]["correlation"]),
            "rmse_mean": np.mean(results[method]["rmse"]),
            "r2_mean": np.mean(results[method]["r2"]),
        }

    return summary, results


def create_comparison_plot(summary: dict, output_dir: str):
    """Create bar chart comparing methods."""
    import matplotlib.pyplot as plt
    import numpy as np

    methods = list(summary.keys())
    correlations = [summary[m]["correlation_mean"] for m in methods]
    errors = [summary[m]["correlation_std"] for m in methods]

    fig, ax = plt.subplots(figsize=(10, 6))

    colors = ["#3498db", "#e74c3c", "#2ecc71", "#9b59b6"]
    bars = ax.bar(
        methods, correlations, yerr=errors, capsize=5, color=colors, edgecolor="black"
    )

    ax.set_ylabel("Prediction Accuracy (Correlation)", fontsize=12)
    ax.set_xlabel("Method", fontsize=12)
    ax.set_title(
        "Bayesian Genomic Prediction: Method Comparison", fontsize=14, fontweight="bold"
    )
    ax.set_ylim(0, 1)
    ax.grid(axis="y", alpha=0.3)

    for bar, corr in zip(bars, correlations):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.03,
            f"{corr:.3f}",
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
        )

    plt.tight_layout()
    plt.savefig(f"{output_dir}/bayesian_comparison.png", dpi=150, bbox_inches="tight")
    print(f"Saved: {output_dir}/bayesian_comparison.png")


def save_results(summary: dict, fold_results: dict, output_dir: str):
    """Save results to CSV."""
    import pandas as pd

    # Summary
    summary_df = pd.DataFrame(summary).T
    summary_df.columns = ["Correlation_Mean", "Correlation_SD", "RMSE_Mean", "R2_Mean"]
    summary_df.to_csv(f"{output_dir}/method_summary.csv")

    # Fold-level
    rows = []
    for method in fold_results:
        for fold in range(len(fold_results[method]["correlation"])):
            rows.append(
                {
                    "Method": method,
                    "Fold": fold + 1,
                    "Correlation": fold_results[method]["correlation"][fold],
                    "RMSE": fold_results[method]["rmse"][fold],
                    "R2": fold_results[method]["r2"][fold],
                }
            )
    fold_df = pd.DataFrame(rows)
    fold_df.to_csv(f"{output_dir}/fold_results.csv", index=False)

    print(f"Saved: {output_dir}/method_summary.csv")
    print(f"Saved: {output_dir}/fold_results.csv")


def main():
    print("=" * 70)
    print("Example: Bayesian Genomic Prediction")
    print("=" * 70)

    # Install packages
    print("\n[1/5] Installing dependencies...")
    install_packages()

    # Setup
    output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(exist_ok=True)

    # Generate data
    print("\n[2/5] Generating synthetic data...")
    data = generate_synthetic_data(
        n_individuals=300, n_markers=500, n_qtl=30, heritability=0.5
    )
    print(f"  Individuals: {data['genotypes'].shape[0]}")
    print(f"  Markers: {data['genotypes'].shape[1]}")
    print(f"  QTL: {len(data['qtl_idx'])}")

    # Run CV
    print("\n[3/5] Running 5-fold cross-validation...")
    summary, fold_results = run_cv_bayesian(data, n_splits=5)

    print("\n  Results:")
    for method in summary:
        print(
            f"    {method:10s}: r = {summary[method]['correlation_mean']:.3f} (±{summary[method]['correlation_std']:.3f})"
        )

    # Save results
    print("\n[4/5] Saving results...")
    save_results(summary, fold_results, str(output_dir))

    # Plot
    print("\n[5/5] Creating visualization...")
    create_comparison_plot(summary, str(output_dir))

    print("\n" + "=" * 70)
    print("Bayesian GP Example Complete!")
    print("=" * 70)
    print("\nOutputs:")
    print(f"  - {output_dir}/bayesian_comparison.png")
    print(f"  - {output_dir}/method_summary.csv")
    print(f"  - {output_dir}/fold_results.csv")


if __name__ == "__main__":
    main()
