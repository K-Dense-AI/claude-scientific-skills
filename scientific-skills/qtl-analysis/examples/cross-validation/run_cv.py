#!/usr/bin/env python3
# Copyright 2026 Clayton Young
# Licensed under Apache License 2.0

"""
Example: Cross-Validation for Genomic Prediction

Demonstrates multiple CV strategies to avoid information leakage.
Promoted from research/genomic_cv_python.py

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


def simulate_data(
    n_ind=400,
    n_markers=800,
    n_qtl=40,
    n_subpops=3,
    n_families=60,
    n_years=6,
    h2=0.6,
    seed=42,
):
    import numpy as np

    rng = np.random.default_rng(seed)

    family_subpops = np.repeat(np.arange(n_subpops), n_families // n_subpops)
    family_subpops = family_subpops[:n_families]
    rng.shuffle(family_subpops)

    base = n_ind // n_families
    family_sizes = np.full(n_families, base, dtype=int)
    family_sizes[: (n_ind - base * n_families)] += 1
    rng.shuffle(family_sizes)

    family_ids = np.concatenate(
        [np.full(family_sizes[f], f, dtype=int) for f in range(n_families)]
    )
    rng.shuffle(family_ids)
    subpops = family_subpops[family_ids]

    base_freqs = rng.uniform(0.1, 0.9, n_markers)
    subpop_freqs = []
    for _ in range(n_subpops):
        deviation = rng.normal(0, 0.10, n_markers)
        subpop_freqs.append(np.clip(base_freqs + deviation, 0.05, 0.95))

    M = np.zeros((n_ind, n_markers), dtype=float)
    for f in range(n_families):
        sp = int(family_subpops[f])
        freqs = subpop_freqs[sp]
        family_freqs = np.clip(freqs + rng.normal(0, 0.03, n_markers), 0.05, 0.95)
        idx = np.where(family_ids == f)[0]
        M[idx, :] = rng.binomial(2, family_freqs, size=(idx.size, n_markers))

    M_centered = M - M.mean(axis=0)
    sd = M_centered.std(axis=0, ddof=1)
    sd[sd == 0] = 1.0
    X = M_centered / sd

    qtl_idx = rng.choice(n_markers, size=n_qtl, replace=False)
    qtl_effects = rng.normal(0, 1, n_qtl)
    g = X[:, qtl_idx] @ qtl_effects

    var_g = float(np.var(g))
    var_e = var_g * (1 - h2) / h2

    n_per_year = n_ind // n_years
    time = np.repeat(np.arange(n_years), n_per_year)
    if len(time) < n_ind:
        time = np.concatenate([time, np.zeros(n_ind - len(time), dtype=int)])
    time = time[:n_ind]
    rng.shuffle(time)
    year_effect = rng.normal(0, 0.5, n_years)

    e = rng.normal(0, np.sqrt(var_e), n_ind)
    y = g + year_effect[time] + e

    p = M.mean(axis=0) / 2.0
    denom = 2.0 * float(np.sum(p * (1 - p)))
    W = (M - 2.0 * p) / np.sqrt(denom)
    G = (W @ W.T) / n_markers

    return {
        "X": X,
        "y": y,
        "subpops": subpops,
        "families": family_ids,
        "time": time,
        "G": G,
        "qtl_idx": qtl_idx,
    }


def rrblup_ridge(X_train, y_train, X_test, alpha=1.0):
    from sklearn.linear_model import Ridge

    model = Ridge(alpha=alpha)
    model.fit(X_train, y_train)
    return model.predict(X_test)


def run_all_cv_strategies(data):
    import numpy as np
    from sklearn.model_selection import KFold, StratifiedKFold, GroupKFold
    from scipy.stats import pearsonr
    from collections import deque

    X, y, subpops, families, time, G = (
        data["X"],
        data["y"],
        data["subpops"],
        data["families"],
        data["time"],
        data["G"],
    )

    def get_metrics(y_true, y_pred):
        r, _ = pearsonr(y_true, y_pred)
        rmse = np.sqrt(np.mean((y_true - y_pred) ** 2))
        return {"r": r, "rmse": rmse}

    results = {}

    # 1. Standard K-fold
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    results["Standard K-fold"] = [
        get_metrics(y[te], rrblup_ridge(X[tr], y[tr], X[te]))["r"]
        for tr, te in kf.split(X)
    ]

    # 2. Stratified by subpopulation
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    results["Stratified K-fold"] = [
        get_metrics(y[te], rrblup_ridge(X[tr], y[tr], X[te]))["r"]
        for tr, te in skf.split(X, subpops)
    ]

    # 3. GroupKFold by family
    gkf = GroupKFold(n_splits=5)
    results["GroupKFold (family)"] = [
        get_metrics(y[te], rrblup_ridge(X[tr], y[tr], X[te]))["r"]
        for tr, te in gkf.split(X, y, groups=families)
    ]

    # 4. Forward validation (time)
    uniq = np.sort(np.unique(time))
    fwd_preds = []
    for i in range(len(uniq) - 1):
        tr = np.where(time <= uniq[i])[0]
        te = np.where(time == uniq[i + 1])[0]
        if len(tr) > 0 and len(te) > 0:
            pred = rrblup_ridge(X[tr], y[tr], X[te])
            fwd_preds.append(get_metrics(y[te], pred)["r"])
    results["Forward validation"] = fwd_preds

    # 5. GBLUP
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    gblup_accs = []
    for tr, te in kf.split(np.arange(len(y))):
        G_tr = G[np.ix_(tr, tr)]
        G_te_tr = G[np.ix_(te, tr)]
        n = len(tr)
        mu = np.mean(y[tr])
        A = G_tr + 1.0 * np.eye(n)
        alpha_vec = np.linalg.solve(A, y[tr] - mu)
        pred = mu + G_te_tr @ alpha_vec
        gblup_accs.append(get_metrics(y[te], pred)["r"])
    results["GBLUP"] = gblup_accs

    return results


def plot_cv_comparison(results, output_dir):
    import matplotlib.pyplot as plt
    import numpy as np

    methods = list(results.keys())
    means = [np.mean(results[m]) for m in methods]
    stds = [np.std(results[m]) for m in methods]

    fig, ax = plt.subplots(figsize=(12, 6))
    colors = ["#3498db", "#e74c3c", "#2ecc71", "#f39c12", "#9b59b6"]
    bars = ax.bar(methods, means, yerr=stds, capsize=5, color=colors, edgecolor="black")

    ax.set_ylabel("Prediction Accuracy (Correlation)", fontsize=12)
    ax.set_xlabel("CV Strategy", fontsize=12)
    ax.set_title("Cross-Validation Strategy Comparison", fontsize=14, fontweight="bold")
    ax.set_ylim(0, 1)
    ax.grid(axis="y", alpha=0.3)

    for bar, mean in zip(bars, means):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.04,
            f"{mean:.3f}",
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )

    plt.xticks(rotation=15, ha="right")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/cv_comparison.png", dpi=150, bbox_inches="tight")
    print(f"Saved: {output_dir}/cv_comparison.png")


def save_results(results, output_dir):
    import pandas as pd
    import numpy as np

    rows = []
    for method, accs in results.items():
        for fold, acc in enumerate(accs, 1):
            rows.append({"Method": method, "Fold": fold, "Correlation": acc})

    df = pd.DataFrame(rows)
    df.to_csv(f"{output_dir}/cv_results.csv", index=False)

    summary = pd.DataFrame(
        {
            "Method": list(results.keys()),
            "Mean_Correlation": [np.mean(results[m]) for m in results],
            "Std_Correlation": [np.std(results[m]) for m in results],
        }
    )
    summary.to_csv(f"{output_dir}/cv_summary.csv", index=False)

    print(f"Saved: {output_dir}/cv_results.csv")
    print(f"Saved: {output_dir}/cv_summary.csv")


def main():
    print("=" * 70)
    print("Example: Cross-Validation for Genomic Prediction")
    print("=" * 70)

    print("\n[1/5] Installing dependencies...")
    install_packages()

    output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(exist_ok=True)

    print("\n[2/5] Generating synthetic data with structure...")
    import numpy as np
    data = simulate_data()
    print(f"  Individuals: {data['X'].shape[0]}, Markers: {data['X'].shape[1]}")
    print(
        f"  Families: {len(np.unique(data['families']))}, Subpops: {len(np.unique(data['subpops']))}"
    )

    print("\n[3/5] Running all CV strategies...")
    results = run_all_cv_strategies(data)

    print("\n  Results:")
    import numpy as np

    for method in results:
        print(
            f"    {method:25s}: r = {np.mean(results[method]):.3f} (±{np.std(results[method]):.3f})"
        )

    print("\n[4/5] Saving results...")
    save_results(results, str(output_dir))

    print("\n[5/5] Creating visualization...")
    plot_cv_comparison(results, str(output_dir))

    print("\n" + "=" * 70)
    print("Cross-Validation Example Complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
