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

    n = 50
    geno_ids = [f"G{i + 1:03d}" for i in range(n)]

    trait_data = pd.DataFrame(
        {
            "genotype": geno_ids,
            "yield_t_ha": rng.normal(6.5, 0.8, size=n),
            "protein_pct": rng.normal(12.0, 1.2, size=n),
            "disease_score": rng.normal(4.0, 1.0, size=n),
        }
    )

    economic_weights = pd.Series(
        {"yield_t_ha": 0.6, "protein_pct": 0.3, "disease_score": -0.4}
    )

    x = trait_data[economic_weights.index]
    x_std = (x - x.mean()) / x.std(ddof=0)
    trait_data["selection_index"] = x_std.mul(economic_weights, axis=1).sum(axis=1)
    trait_data = trait_data.sort_values("selection_index", ascending=False).reset_index(
        drop=True
    )

    top = trait_data.head(12).copy()
    top["selection_decision"] = "selected"

    trait_data.to_csv(out / "selection_index_scores.csv", index=False)
    top.to_csv(out / "selection_index_top12.csv", index=False)

    weights = economic_weights.rename("weight").reset_index()
    weights.columns = ["trait", "weight"]
    weights.to_csv(out / "selection_index_weights.csv", index=False)

    plot_df = trait_data.head(15).iloc[::-1]
    plt.figure(figsize=(8, 6))
    plt.barh(plot_df["genotype"], plot_df["selection_index"], color="#2b8cbe")
    plt.xlabel("Selection Index")
    plt.title("Top 15 Selection Index Scores")
    plt.tight_layout()
    plt.savefig(out / "selection_index_top15.png", dpi=150)
    plt.close()

    print("Saved selection index scores, top selections, weights, and plot")


if __name__ == "__main__":
    main()
