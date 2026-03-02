#!/usr/bin/env python3
# Copyright 2026 Clayton Young
# Licensed under Apache License 2.0

"""
Example: G×E (Genotype × Environment) Genomic Prediction

Demonstrates multi-environment genomic prediction with reaction norm visualization.
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


def generate_met_data(n_ind=200, n_envs=3, n_markers=300, n_qtl=25, h2=0.5, seed=42):
    import numpy as np

    rng = np.random.default_rng(seed)

    X = rng.binomial(2, 0.3, size=(n_ind, n_markers)).astype(float)

    qtl_idx = rng.choice(n_markers, size=n_qtl, replace=False)
    main_effects = rng.normal(0, 1.0, n_qtl)
    gxe_effects = rng.normal(0, 0.5, (n_qtl, n_envs))

    env_means = rng.normal(0, 1.0, n_envs)

    phenotypes = np.zeros((n_ind, n_envs))
    for env in range(n_envs):
        tbv = X[:, qtl_idx] @ (main_effects + gxe_effects[:, env])
        var_tbv = float(np.var(tbv))
        var_e = var_tbv * (1 - h2) / h2
        phenotypes[:, env] = tbv + env_means[env] + rng.normal(0, np.sqrt(var_e), n_ind)

    p = X.mean(axis=0) / 2.0
    denom = 2.0 * float(np.sum(p * (1 - p)))
    Z = (X - 2.0 * p) / np.sqrt(denom)
    G = (Z @ Z.T) / n_markers

    return {
        "X": X,
        "phenotypes": phenotypes,
        "G": G,
        "qtl_idx": qtl_idx,
        "env_means": env_means,
    }


def run_single_env_model(data):
    import numpy as np
    from sklearn.linear_model import Ridge
    from sklearn.model_selection import KFold
    from scipy.stats import pearsonr

    X, Y, G = data["X"], data["phenotypes"], data["G"]
    n_ind, n_envs = Y.shape
    kf = KFold(n_splits=5, shuffle=True, random_state=42)

    results = {env: [] for env in range(n_envs)}

    for fold, (tr_idx, te_idx) in enumerate(kf.split(np.arange(n_ind))):
        for env in range(n_envs):
            y_tr = Y[tr_idx, env]
            y_te = Y[te_idx, env]

            model = Ridge(alpha=1.0)
            model.fit(X[tr_idx], y_tr)
            pred = model.predict(X[te_idx])

            r, _ = pearsonr(pred, y_te)
            results[env].append(r)

    return results


def run_multienv_model(data):
    import numpy as np
    from sklearn.linear_model import Ridge
    from sklearn.model_selection import KFold
    from scipy.stats import pearsonr

    X, Y, G = data["X"], data["phenotypes"], data["G"]
    n_ind, n_envs = Y.shape
    kf = KFold(n_splits=5, shuffle=True, random_state=42)

    results = {env: [] for env in range(n_envs)}

    for fold, (tr_idx, te_idx) in enumerate(kf.split(np.arange(n_ind))):
        y_tr = Y[tr_idx]
        y_te = Y[te_idx]

        y_mean = y_tr.mean()
        y_centered = y_tr - y_mean

        G_tr = G[np.ix_(tr_idx, tr_idx)]
        G_te_tr = G[np.ix_(te_idx, tr_idx)]

        n = len(tr_idx)
        A = G_tr + 1.0 * np.eye(n)
        alpha_vec = np.linalg.solve(A, y_centered)

        for env in range(n_envs):
            pred = y_mean + G_te_tr @ alpha_vec[:, env]
            r, _ = pearsonr(pred, y_te[:, env])
            results[env].append(r)

    return results


def plot_reaction_norms(data, output_dir):
    import matplotlib.pyplot as plt
    import numpy as np

    X, phenotypes = data["X"], data["phenotypes"]
    n_ind, n_envs = phenotypes.shape
    env_names = [f"Env {i + 1}" for i in range(n_envs)]

    sample_idx = np.random.choice(n_ind, min(30, n_ind), replace=False)

    fig, ax = plt.subplots(figsize=(10, 6))

    for idx in sample_idx:
        ax.plot(range(n_envs), phenotypes[idx], "o-", alpha=0.4, linewidth=1)

    env_means = phenotypes.mean(axis=0)
    ax.plot(
        range(n_envs), env_means, "ko-", linewidth=3, markersize=8, label="Env Mean"
    )

    ax.set_xticks(range(n_envs))
    ax.set_xticklabels(env_names)
    ax.set_xlabel("Environment")
    ax.set_ylabel("Phenotype")
    ax.set_title("Reaction Norms: G×E Interaction")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/reaction_norms.png", dpi=150, bbox_inches="tight")
    print(f"Saved: {output_dir}/reaction_norms.png")


def plot_accuracy_comparison(single_env, multi_env, output_dir):
    import matplotlib.pyplot as plt
    import numpy as np

    n_envs = len(single_env)
    envs = [f"Env {i + 1}" for i in range(n_envs)]

    single_means = [np.mean(single_env[i]) for i in range(n_envs)]
    multi_means = [np.mean(multi_env[i]) for i in range(n_envs)]

    x = np.arange(n_envs)
    width = 0.35

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(
        x - width / 2,
        single_means,
        width,
        label="Single-Env",
        color="#3498db",
        edgecolor="black",
    )
    ax.bar(
        x + width / 2,
        multi_means,
        width,
        label="Multi-Env",
        color="#e74c3c",
        edgecolor="black",
    )

    ax.set_xticks(x)
    ax.set_xticklabels(envs)
    ax.set_ylabel("Prediction Accuracy (r)")
    ax.set_title("Single-Env vs Multi-Env Model Accuracy")
    ax.legend()
    ax.set_ylim(0, 1)
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/gxe_accuracy.png", dpi=150, bbox_inches="tight")
    print(f"Saved: {output_dir}/gxe_accuracy.png")


def save_results(single_env, multi_env, output_dir):
    import pandas as pd
    import numpy as np

    rows = []
    for env in range(len(single_env)):
        for fold, (s, m) in enumerate(zip(single_env[env], multi_env[env]), 1):
            rows.append(
                {
                    "Environment": f"Env_{env + 1}",
                    "Fold": fold,
                    "Single_Env": s,
                    "Multi_Env": m,
                }
            )

    df = pd.DataFrame(rows)
    df.to_csv(f"{output_dir}/gxe_results.csv", index=False)

    summary = pd.DataFrame(
        {
            "Environment": [f"Env_{i + 1}" for i in range(len(single_env))],
            "Single_Env_Mean": [np.mean(single_env[i]) for i in range(len(single_env))],
            "Multi_Env_Mean": [np.mean(multi_env[i]) for i in range(len(multi_env))],
        }
    )
    summary.to_csv(f"{output_dir}/gxe_summary.csv", index=False)

    print(f"Saved: {output_dir}/gxe_results.csv")
    print(f"Saved: {output_dir}/gxe_summary.csv")


def main():
    print("=" * 70)
    print("Example: G×E (Genotype × Environment) Genomic Prediction")
    print("=" * 70)

    print("\n[1/5] Installing dependencies...")
    install_packages()

    output_dir = Path(__file__).parent / "output"
    output_dir.mkdir(exist_ok=True)

    print("\n[2/5] Generating multi-environment data...")
    data = generate_met_data()
    print(f"  Individuals: {data['X'].shape[0]}, Markers: {data['X'].shape[1]}")
    print(f"  Environments: {data['phenotypes'].shape[1]}")

    print("\n[3/5] Running single-environment models...")
    single_env = run_single_env_model(data)

    print("\n[4/5] Running multi-environment model...")
    multi_env = run_multienv_model(data)

    import numpy as np

    print("\n  Results:")
    for env in range(data["phenotypes"].shape[1]):
        print(
            f"    Env {env + 1}: Single={np.mean(single_env[env]):.3f}, Multi={np.mean(multi_env[env]):.3f}"
        )

    print("\n[5/5] Saving and plotting...")
    save_results(single_env, multi_env, str(output_dir))
    plot_reaction_norms(data, str(output_dir))
    plot_accuracy_comparison(single_env, multi_env, str(output_dir))

    print("\n" + "=" * 70)
    print("G×E Prediction Example Complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
