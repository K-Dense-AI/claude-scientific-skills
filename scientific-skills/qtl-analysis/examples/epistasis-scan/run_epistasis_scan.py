#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)
    rng = np.random.default_rng(314)

    n_samples = 300
    n_snps = 18
    snp_names = [f"SNP{i + 1:03d}" for i in range(n_snps)]
    g = rng.integers(0, 3, size=(n_samples, n_snps))

    y = (
        0.35 * g[:, 2]
        + 0.45 * g[:, 7]
        + 0.70 * (g[:, 2] * g[:, 7])
        + rng.normal(0, 1.0, n_samples)
    )

    rows = []
    for i, j in combinations(range(n_snps), 2):
        x = g[:, i] * g[:, j]
        r = np.corrcoef(x, y)[0, 1]
        score = abs(r) * np.sqrt(n_samples)
        p_proxy = float(np.exp(-score))
        rows.append(
            {
                "snp_a": snp_names[i],
                "snp_b": snp_names[j],
                "interaction_score": round(float(score), 5),
                "p_value_proxy": p_proxy,
            }
        )

    res = pd.DataFrame(rows).sort_values("p_value_proxy").reset_index(drop=True)
    res.to_csv(out / "epistasis_pair_scan.csv", index=False)
    res.head(20).to_csv(out / "epistasis_top_pairs.csv", index=False)

    top = res.head(20).iloc[::-1]
    labels = [f"{a} x {b}" for a, b in zip(top["snp_a"], top["snp_b"])]
    plt.figure(figsize=(8.5, 6))
    plt.barh(labels, -np.log10(top["p_value_proxy"] + 1e-12), color="#1f77b4")
    plt.xlabel("-log10(p-value proxy)")
    plt.title("Top Epistatic SNP Pair Signals")
    plt.tight_layout()
    plt.savefig(out / "epistasis_top_pairs_plot.png", dpi=160)
    plt.close()

    matrix = np.full((n_snps, n_snps), np.nan)
    for _, r in res.iterrows():
        a = int(r["snp_a"].replace("SNP", "")) - 1
        b = int(r["snp_b"].replace("SNP", "")) - 1
        matrix[a, b] = -np.log10(float(r["p_value_proxy"]) + 1e-12)
        matrix[b, a] = matrix[a, b]
    np.fill_diagonal(matrix, 0)
    plt.figure(figsize=(6.8, 5.8))
    plt.imshow(matrix, cmap="magma")
    plt.colorbar(label="-log10(p-value proxy)")
    plt.title("Genome-wide Epistasis Signal Matrix")
    plt.xlabel("SNP index")
    plt.ylabel("SNP index")
    plt.tight_layout()
    plt.savefig(out / "epistasis_signal_matrix.png", dpi=160)
    plt.close()

    print("Saved epistasis scan CSVs and plots")


if __name__ == "__main__":
    main()
