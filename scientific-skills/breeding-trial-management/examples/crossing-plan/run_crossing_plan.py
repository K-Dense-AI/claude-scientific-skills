#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)
    rng = np.random.default_rng(99)

    n = 20
    parents = [f"P{i + 1:02d}" for i in range(n)]
    gebv = pd.Series(rng.normal(0.0, 1.0, n), index=parents, name="gebv")

    raw = rng.uniform(0.0, 0.35, size=(n, n))
    coancestry = (raw + raw.T) / 2
    np.fill_diagonal(coancestry, 0.5)
    coa = pd.DataFrame(coancestry, index=parents, columns=parents)
    coa.to_csv(out / "coancestry_matrix.csv", index=True)

    rows = []
    penalty = 1.5
    for p1, p2 in combinations(parents, 2):
        mean_parent_gebv = (gebv[p1] + gebv[p2]) / 2
        f_est = float(coa.loc[p1, p2])
        cross_score = mean_parent_gebv - penalty * f_est
        rows.append(
            {
                "parent_1": p1,
                "parent_2": p2,
                "mean_parent_gebv": round(mean_parent_gebv, 4),
                "coancestry": round(f_est, 4),
                "cross_score": round(cross_score, 4),
            }
        )

    plan = pd.DataFrame(rows).sort_values("cross_score", ascending=False)

    selected = []
    used = set()
    for _, row in plan.iterrows():
        if row["parent_1"] in used or row["parent_2"] in used:
            continue
        selected.append(row)
        used.add(row["parent_1"])
        used.add(row["parent_2"])
        if len(selected) == 8:
            break

    selected_df = pd.DataFrame(selected)
    plan.to_csv(out / "all_candidate_crosses.csv", index=False)
    selected_df.to_csv(out / "optimal_crossing_plan.csv", index=False)

    print("Saved candidate crosses, selected crossing plan, and coancestry matrix")


if __name__ == "__main__":
    main()
