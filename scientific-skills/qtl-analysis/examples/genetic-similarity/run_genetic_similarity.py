#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform


def ibs_similarity(x: np.ndarray):
    n = x.shape[0]
    sim = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i, n):
            d = np.abs(x[i] - x[j])
            s = 1.0 - np.mean(d / 2.0)
            sim[i, j] = s
            sim[j, i] = s
    return sim


def main():
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)

    rng = np.random.default_rng(7)
    n_ind, n_markers = 50, 800
    p = rng.uniform(0.05, 0.45, size=n_markers)
    x = rng.binomial(2, p, size=(n_ind, n_markers))
    sim = ibs_similarity(x)
    dist = 1 - sim

    ids = [f"ID_{i + 1:03d}" for i in range(n_ind)]
    pd.DataFrame(sim, index=ids, columns=ids).to_csv(out / "ibs_similarity_matrix.csv")
    pd.DataFrame(dist, index=ids, columns=ids).to_csv(out / "ibs_distance_matrix.csv")

    condensed = squareform(dist, checks=False)
    z = linkage(condensed, method="average")
    plt.figure(figsize=(10, 4))
    dendrogram(z, labels=ids, leaf_rotation=90, leaf_font_size=6)
    plt.title("Genetic Similarity Dendrogram (IBS distance)")
    plt.tight_layout()
    plt.savefig(out / "genetic_similarity_dendrogram.png", dpi=150)
    plt.close()
    print("Saved IBS similarity/distance matrices and dendrogram")


if __name__ == "__main__":
    main()
