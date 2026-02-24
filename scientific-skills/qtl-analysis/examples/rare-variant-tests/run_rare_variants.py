#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

"""Rare variant burden and SKAT-like set test example."""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import chi2


def simulate_data(
    n_samples: int = 600, n_genes: int = 25, variants_per_gene: int = 12, seed: int = 42
):
    rng = np.random.default_rng(seed)
    n_variants = n_genes * variants_per_gene
    mafs = rng.uniform(0.001, 0.04, size=n_variants)
    geno = np.vstack([rng.binomial(2, maf, size=n_samples) for maf in mafs]).T
    genes = np.repeat([f"GENE_{i + 1:02d}" for i in range(n_genes)], variants_per_gene)

    causal = set(rng.choice(np.arange(n_genes), size=5, replace=False).tolist())
    beta = np.zeros(n_variants)
    for g in causal:
        idx = np.where(genes == f"GENE_{g + 1:02d}")[0]
        beta[idx] = rng.normal(0.8, 0.2, size=len(idx))

    signal = geno @ beta
    noise = rng.normal(0, np.std(signal) * 1.25 + 1e-6, size=n_samples)
    phenotype = signal + noise

    return geno, phenotype, genes, mafs, causal


def burden_test(geno_set: np.ndarray, y: np.ndarray, weights: np.ndarray):
    score = geno_set @ weights
    if np.std(score) < 1e-12:
        return 1.0
    r = np.corrcoef(score, y)[0, 1]
    if np.isnan(r):
        return 1.0
    z2 = (len(y) - 2) * (r * r) / max(1 - r * r, 1e-12)
    return float(1 - chi2.cdf(z2, df=1))


def skat_like_test(geno_set: np.ndarray, y: np.ndarray, weights: np.ndarray):
    y0 = y - np.mean(y)
    xw = geno_set * weights
    q = float(np.sum((xw.T @ y0) ** 2))
    if q <= 0:
        return 1.0
    lam = np.linalg.svd(xw, compute_uv=False) ** 2
    lam = lam[lam > 1e-12]
    if lam.size == 0:
        return 1.0
    df = max(1, int(round((np.sum(lam) ** 2) / np.sum(lam**2))))
    scale = np.sum(lam) / df
    stat = q / max(scale, 1e-12)
    return float(1 - chi2.cdf(stat, df=df))


def run_tests(geno, y, genes, mafs):
    rows = []
    for gene in np.unique(genes):
        idx = np.where(genes == gene)[0]
        g = geno[:, idx]
        w = 1.0 / np.sqrt(np.clip(mafs[idx], 1e-4, None))
        p_burden = burden_test(g, y, w)
        p_skat = skat_like_test(g, y, w)
        rows.append(
            {
                "gene": gene,
                "variants": int(len(idx)),
                "burden_p": p_burden,
                "skat_like_p": p_skat,
                "min_p": min(p_burden, p_skat),
            }
        )
    return pd.DataFrame(rows).sort_values("min_p").reset_index(drop=True)


def plot_results(df: pd.DataFrame, out_path: Path):
    x = np.arange(len(df))
    plt.figure(figsize=(10, 5))
    plt.scatter(
        x,
        -np.log10(np.clip(df["burden_p"].to_numpy(), 1e-300, 1)),
        label="Burden",
        alpha=0.7,
    )
    plt.scatter(
        x,
        -np.log10(np.clip(df["skat_like_p"].to_numpy(), 1e-300, 1)),
        label="SKAT-like",
        alpha=0.7,
    )
    plt.xlabel("Gene rank")
    plt.ylabel("-log10(p)")
    plt.title("Rare Variant Set Tests")
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()


def main():
    out_dir = Path(__file__).parent / "output"
    out_dir.mkdir(exist_ok=True)

    geno, y, genes, mafs, causal = simulate_data()
    res = run_tests(geno, y, genes, mafs)

    n_genes = len(np.unique(genes))
    threshold = 0.05 / n_genes
    res["burden_sig"] = res["burden_p"] < threshold
    res["skat_sig"] = res["skat_like_p"] < threshold
    res["is_causal"] = (
        res["gene"].str.extract(r"(\d+)").astype(int).iloc[:, 0].sub(1).isin(causal)
    )

    plot_results(res, out_dir / "rare_variant_method_comparison.png")
    res.to_csv(out_dir / "rare_variant_gene_results.csv", index=False)

    summary = pd.DataFrame(
        {
            "metric": [
                "bonferroni_threshold",
                "burden_significant_genes",
                "skat_significant_genes",
                "causal_genes_total",
            ],
            "value": [
                threshold,
                int(res["burden_sig"].sum()),
                int(res["skat_sig"].sum()),
                int(len(causal)),
            ],
        }
    )
    summary.to_csv(out_dir / "rare_variant_summary.csv", index=False)

    print(f"Causal genes: {sorted([f'GENE_{i + 1:02d}' for i in causal])}")
    print(f"Bonferroni threshold: {threshold:.3e}")
    print("Saved: output/rare_variant_method_comparison.png")
    print("Saved: output/rare_variant_gene_results.csv")
    print("Saved: output/rare_variant_summary.csv")


if __name__ == "__main__":
    main()
