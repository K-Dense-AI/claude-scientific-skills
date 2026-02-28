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

    rng = np.random.default_rng(42)
    n_samples, n_markers = 200, 1200
    maf = rng.uniform(0.05, 0.5, size=n_markers)
    geno = rng.binomial(2, maf, size=(n_samples, n_markers)).astype(float)

    missing_mask = rng.random((n_samples, n_markers)) < 0.02
    geno[missing_mask] = np.nan

    call_rate = 1 - np.mean(np.isnan(geno), axis=1)
    het = np.nanmean((geno == 1).astype(float), axis=1)
    sex = rng.choice(["F", "M"], size=n_samples)
    x_het = het + np.where(sex == "F", 0.03, -0.03)

    qc = pd.DataFrame(
        {
            "sample": [f"S{i + 1:03d}" for i in range(n_samples)],
            "sex_reported": sex,
            "call_rate": call_rate,
            "heterozygosity": het,
            "x_het_proxy": x_het,
        }
    )
    qc["call_rate_fail"] = qc["call_rate"] < 0.95
    h_mu, h_sd = qc["heterozygosity"].mean(), qc["heterozygosity"].std()
    qc["het_outlier"] = (qc["heterozygosity"] < h_mu - 3 * h_sd) | (
        qc["heterozygosity"] > h_mu + 3 * h_sd
    )

    qc.to_csv(out / "sample_qc_metrics.csv", index=False)

    plt.figure(figsize=(8, 5))
    plt.scatter(qc["call_rate"], qc["heterozygosity"], alpha=0.7)
    plt.axvline(0.95, color="red", linestyle="--", linewidth=1.5)
    plt.xlabel("Call rate")
    plt.ylabel("Heterozygosity")
    plt.title("Sample QC: Call rate vs heterozygosity")
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(out / "sample_qc_scatter.png", dpi=150)
    plt.close()

    sex_summary = qc.groupby("sex_reported")["x_het_proxy"].mean().reset_index()
    sex_summary.to_csv(out / "sex_check_summary.csv", index=False)
    print("Saved sample_qc_metrics.csv, sample_qc_scatter.png, sex_check_summary.csv")


if __name__ == "__main__":
    main()
