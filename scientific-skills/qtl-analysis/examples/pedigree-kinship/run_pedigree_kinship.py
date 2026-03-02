#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def build_simple_nrm(pedigree: pd.DataFrame) -> pd.DataFrame:
    ids = pedigree["id"].tolist()
    idx = {a: i for i, a in enumerate(ids)}
    n = len(ids)
    a = np.zeros((n, n), dtype=float)

    for i, animal in enumerate(ids):
        sire = pedigree.loc[pedigree["id"] == animal, "sire"].iloc[0]
        dam = pedigree.loc[pedigree["id"] == animal, "dam"].iloc[0]
        if sire == "0" and dam == "0":
            a[i, i] = 1.0
            continue

        s = idx.get(sire)
        d = idx.get(dam)
        for j in range(i + 1):
            v = 0.0
            if s is not None:
                v += 0.5 * a[s, j]
            if d is not None:
                v += 0.5 * a[d, j]
            a[i, j] = v
            a[j, i] = v
        f_sd = 0.0 if (s is None or d is None) else a[s, d]
        a[i, i] = 1.0 + 0.5 * f_sd

    return pd.DataFrame(a, index=ids, columns=ids)


def main():
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)

    ped = pd.DataFrame(
        {
            "id": ["A", "B", "C", "D", "E", "F", "G", "H"],
            "sire": ["0", "0", "A", "A", "C", "C", "E", "F"],
            "dam": ["0", "0", "B", "B", "D", "D", "0", "G"],
        }
    )
    nrm = build_simple_nrm(ped)
    nrm.to_csv(out / "pedigree_nrm_matrix.csv")

    plt.figure(figsize=(6, 5))
    plt.imshow(nrm.to_numpy(), cmap="viridis")
    plt.colorbar(label="Relationship")
    plt.xticks(range(len(nrm.columns)), nrm.columns)
    plt.yticks(range(len(nrm.index)), nrm.index)
    plt.title("Pedigree NRM")
    plt.tight_layout()
    plt.savefig(out / "pedigree_nrm_heatmap.png", dpi=150)
    plt.close()
    print("Saved pedigree_nrm_matrix.csv and pedigree_nrm_heatmap.png")


if __name__ == "__main__":
    main()
