#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)
    rng = np.random.default_rng(512)

    n = 220
    samples = [f"S{i + 1:03d}" for i in range(n)]
    chroms = [f"chr{i}" for i in range(1, 11)]

    calls = pd.DataFrame(
        {
            "sample": rng.choice(samples, size=900, replace=True),
            "chrom": rng.choice(chroms, size=900, replace=True),
            "start": rng.integers(1_000_000, 40_000_000, size=900),
            "length_bp": rng.integers(50_000, 1_500_000, size=900),
            "cn_state": rng.choice([0, 1, 3, 4], size=900, p=[0.08, 0.27, 0.56, 0.09]),
        }
    )
    calls["end"] = calls["start"] + calls["length_bp"]
    calls.to_csv(out / "cnv_calls.csv", index=False)

    burden = calls.groupby(["sample", "chrom"]).size().reset_index(name="cnv_burden")
    pheno = pd.DataFrame({"sample": samples, "trait": rng.normal(0, 1, n)})
    merged = burden.merge(pheno, on="sample", how="left")
    assoc = (
        merged.groupby("chrom")[["cnv_burden", "trait"]]
        .apply(lambda d: np.corrcoef(d["cnv_burden"], d["trait"])[0, 1])
        .reset_index(name="corr")
    )
    assoc["p_value_proxy"] = np.exp(-np.abs(assoc["corr"]) * np.sqrt(80))
    assoc = assoc.sort_values("p_value_proxy")
    assoc.to_csv(out / "cnv_trait_association.csv", index=False)

    burden_chr = calls.groupby("chrom").size().reset_index(name="events")
    burden_chr.to_csv(out / "cnv_burden_by_chr.csv", index=False)

    plt.figure(figsize=(8, 4.4))
    plt.bar(burden_chr["chrom"], burden_chr["events"], color="#2ca02c")
    plt.title("CNV Event Burden by Chromosome")
    plt.xlabel("Chromosome")
    plt.ylabel("CNV events")
    plt.tight_layout()
    plt.savefig(out / "cnv_burden_plot.png", dpi=160)
    plt.close()

    plt.figure(figsize=(7, 4.4))
    plt.scatter(
        assoc["chrom"], -np.log10(assoc["p_value_proxy"] + 1e-12), s=80, c="#1f77b4"
    )
    plt.title("CNV-Phenotype Association Scan")
    plt.xlabel("Chromosome")
    plt.ylabel("-log10(p-value proxy)")
    plt.tight_layout()
    plt.savefig(out / "cnv_association_plot.png", dpi=160)
    plt.close()

    print("Saved CNV integration outputs and plots")


if __name__ == "__main__":
    main()
