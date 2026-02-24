#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from pathlib import Path
import numpy as np
import pandas as pd


def main():
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)
    rng = np.random.default_rng(2026)

    n = 60
    df = pd.DataFrame(
        {
            "genotype": [f"G{i + 1:03d}" for i in range(n)],
            "gebv": rng.normal(0.0, 1.0, size=n),
            "se": rng.uniform(0.15, 0.35, size=n),
        }
    )

    z = 1.96
    df["ci_lower"] = df["gebv"] - z * df["se"]
    df["ci_upper"] = df["gebv"] + z * df["se"]
    df = df.sort_values("gebv", ascending=False).reset_index(drop=True)
    df["rank"] = np.arange(1, n + 1)

    selected = df.head(15).copy()
    selected["selection_decision"] = "advance"

    df.to_csv(out / "gebv_ranking_with_ci.csv", index=False)
    selected.to_csv(out / "selection_list_top15.csv", index=False)

    print("Saved GEBV ranking with confidence intervals and selection list")


if __name__ == "__main__":
    main()
