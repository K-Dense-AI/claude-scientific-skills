#!/usr/bin/env python3
# Copyright 2026 Clayton Young
# Licensed under Apache License 2.0

"""
Example: Elastic Net Cross-Validation for SNP Selection

Demonstrates elastic net regularization for marker selection in genomic prediction.
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


def generate_data(n_ind=300, n_markers=500, n_qtl=20, h2=0.5, seed=42):
    import numpy as np

    rng = np.random.default_rng(seed)
    X = rng.binomial(2, 0.3, size=(n_ind, n_markers)).astype(float)
    qtl_idx = rng.choice(n_markers, size=n_qtl, replace=False)
    qtl_effects = rng.normal(0, 1.0, n_qtl)
    tbv = X[:, qtl_idx] @ qtl_effects
    var_tbv = float(np.var(tbv))
    var_e = var_tbv * (1 - h2) / h2
    y = tbv + rng.normal(0, np.sqrt(var_e), n_ind)
    return {"X": X, "y": y, "tbv": tbv, "qtl_idx": qtl_idx}


def run_elastic_net_cv(data, alphas=None, n_splits=5):
    import numpy as np
    from sklearn.linear_model import ElasticNetCV, RidgeCV
    from sklearn.model_selection import KFold
    from scipy.stats import pearsonr

    X, y, tbv = data["X"], data["y"], data["tbv"]

    if alphas is None:
        alphas = np.logspace(-3, 1, 20)

    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

    # ElasticNet with CV for alpha selection
    model = ElasticNetCV(
        l1_ratio=0.5, alphas=alphas, cv=n_splits, random_state=42, max_iter=5000
    )

    results = {"alphas": [], "correlations": [], "n_selected": []}

    for train_idx, test_idx in kf.split(X):
        X_tr, X_te = X[train_idx], X[test_idx]
        y_tr, y_te = y[train_idx], y[test_idx]
        tbv_te = tbv[test_idx]

        model.fit(X_tr, y_tr)
        pred = model.predict(X_te)

        r, _ = pearsonr(pred, tbv_te)
        n_sel = np.sum(np.abs(model.coef_) > 1e-6)

        results["alphas"].append(model.alpha_)
        results["correlations"].append(r)
        results["n_selected"].append(n_sel)

    # Fit final model on all data
    final_model = ElasticNetCV(
        l1_ratio=0.5, alphas=alphas, cv=5, random_state=42, max_iter=5000
    )
    final_model.fit(X, y)
    selected = np.where(np.abs(final_model.coef_) > 1e-6)[0]

    return results, final_model, selected


def plot_lambda_optimization(cv_results, output_dir):
    import matplotlib.pyplot as plt
    import numpy as np

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Alpha selection
    ax1.plot(range(len(cv_results["alphas"])), cv_results["alphas"], "bo-")
    ax1.set_xlabel("Fold")
    ax1.set_ylabel("Selected Alpha")
    ax1.set_title("Lambda (Alpha) Selection per Fold")
    ax1.set_yscale("log")
    ax1.grid(True, alpha=0.3)

    # Accuracy
    ax2.bar(
        range(len(cv_results["correlations"])),
        cv_results["correlations"],
        color="#3498db",
        edgecolor="black",
    )
    ax2.set_xlabel("Fold")
    ax2.set_ylabel("Correlation")
    ax2.set_title("Prediction Accuracy per Fold")
    ax2.set_ylim(0, 1)
    ax2.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/lambda_optimization.png", dpi=150, bbox_inches="tight")
    print(f"Saved: {output_dir}/lambda_optimization.png")


def plot_selected_snaps(selected, data, output_dir):
    import matplotlib.pyplot as plt
    import numpy as np

    qtl_idx = data["qtl_idx"]

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.hist(selected, bins=50, alpha=0.7, color="#2ecc71", edgecolor="black")
    ax.axvline(
        qtl_idx.min(), color="red", linestyle="--", linewidth=2, label=f"True QTL range"
    )
    ax.axvline(qtl_idx.max(), color="red", linestyle="--", linewidth=2)
    ax.set_xlabel("Marker Index")
    ax.set_ylabel("Frequency")
    ax.set_title(f"Selected Markers: {len(selected)} total")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/selected_markers.png", dpi=150, bbox_inches="tight")
    print(f"Saved: {output_dir}/selected_markers.png")


def save_results(cv_results, selected, data, output_dir):
    import pandas as pd
    import numpy as np

    # CV summary
    cv_df = pd.DataFrame(
        {
            "Fold": range(1, len(cv_results["correlations"]) + 1),
            "Alpha": cv_results["alphas"],
            "Correlation": cv_results["correlations"],
            "N_Selected": cv_results["n_selected"],
        }
    )
    cv_df.to_csv(f"{output_dir}/cv_results.csv", index=False)

    # Selected SNPs
    snp_df = pd.DataFrame(
        {
            "Marker_Index": selected,
            "Effect": [
                data["X"][:, i].dot(data["y"]) / len(data["y"]) for i in selected
            ],
        }
    )
    snp_df.to_csv(f"{output_dir}/selected_snps.csv", index=False)

    print(f"Saved: {output_dir}/cv_results.csv")
    print(f"Saved: {output_dir}/selected_snps.csv")


def main():
    print("=" * 70)
    print("Example: Elastic Net Cross-Validation for SNP Selection")
    print("=" * 70)

    print("\n[1/5] Installing dependencies...")
    install_packages()

    output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(exist_ok=True)

    print("\n[2/5] Generating synthetic data...")
    import numpy as np
    data = generate_data()
    print(
        f"  Individuals: {data['X'].shape[0]}, Markers: {data['X'].shape[1]}, QTL: {len(data['qtl_idx'])}"
    )

    print("\n[3/5] Running Elastic Net CV...")
    cv_results, model, selected = run_elastic_net_cv(data)
    print(f"  Mean correlation: {np.mean(cv_results['correlations']):.3f}")
    print(f"  Selected markers: {len(selected)}")

    print("\n[4/5] Saving results...")
    save_results(cv_results, selected, data, str(output_dir))

    print("\n[5/5] Creating visualizations...")
    plot_lambda_optimization(cv_results, str(output_dir))
    plot_selected_snaps(selected, data, str(output_dir))

    print("\n" + "=" * 70)
    print("Elastic Net CV Complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
