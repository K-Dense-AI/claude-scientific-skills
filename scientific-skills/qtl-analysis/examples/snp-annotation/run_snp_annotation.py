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

    rng = np.random.default_rng(11)
    n = 300
    effects = rng.choice(
        [
            "missense_variant",
            "synonymous_variant",
            "intron_variant",
            "intergenic_variant",
            "stop_gained",
        ],
        size=n,
        p=[0.25, 0.2, 0.3, 0.2, 0.05],
    )
    impacts = np.array(
        [
            "HIGH"
            if e in {"stop_gained"}
            else "MODERATE"
            if e in {"missense_variant"}
            else "LOW"
            for e in effects
        ]
    )
    genes = [f"Gene{rng.integers(1, 80):03d}" for _ in range(n)]

    df = pd.DataFrame(
        {
            "variant_id": [f"var_{i + 1:05d}" for i in range(n)],
            "chrom": rng.integers(1, 11, size=n),
            "pos": rng.integers(1, 5_000_000, size=n),
            "gene": genes,
            "effect": effects,
            "impact": impacts,
        }
    )
    df.to_csv(out / "annotated_variants.csv", index=False)

    effect_counts = df["effect"].value_counts().sort_values(ascending=False)
    effect_counts.to_csv(out / "effect_summary.csv", header=["count"])

    plt.figure(figsize=(8, 4.5))
    effect_counts.plot(kind="bar", color="#3498db")
    plt.ylabel("Count")
    plt.title("SNP Effect Annotation Summary")
    plt.tight_layout()
    plt.savefig(out / "effect_summary_plot.png", dpi=150)
    plt.close()
    print("Saved annotated_variants.csv, effect_summary.csv, effect_summary_plot.png")


if __name__ == "__main__":
    main()
