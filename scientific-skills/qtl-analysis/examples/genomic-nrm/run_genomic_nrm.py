#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def genomic_nrm_02_scale(x: np.ndarray, p: np.ndarray):
    z = x - 2 * p
    denom = 2 * np.sum(p * (1 - p))
    g = (z @ z.T) / max(denom, 1e-12)
    return g


def main():
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)

    rng = np.random.default_rng(42)
    n_ind, n_markers = 120, 600
    p = rng.uniform(0.05, 0.5, size=n_markers)
    x = rng.binomial(2, p, size=(n_ind, n_markers))
    g = genomic_nrm_02_scale(x, p)

    ids = [f"ID_{i + 1:03d}" for i in range(n_ind)]
    pd.DataFrame(g, index=ids, columns=ids).to_csv(out / "genomic_nrm_0_2.csv")

    plt.figure(figsize=(7, 6))
    plt.imshow(g, cmap="magma")
    plt.colorbar(label="Genomic relationship")
    plt.title("Genomic NRM (0-2 scale)")
    plt.tight_layout()
    plt.savefig(out / "genomic_nrm_heatmap.png", dpi=150)
    plt.close()
    print("Saved genomic_nrm_0_2.csv and genomic_nrm_heatmap.png")


if __name__ == "__main__":
    main()
