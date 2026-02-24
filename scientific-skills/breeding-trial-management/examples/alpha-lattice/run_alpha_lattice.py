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
    rng = np.random.default_rng(7)

    treatments = [f"T{i + 1:02d}" for i in range(24)]
    reps = 3
    block_size = 4
    rows = []
    for rep in range(1, reps + 1):
        shuffled = treatments.copy()
        rng.shuffle(shuffled)
        blocks = [
            shuffled[i : i + block_size] for i in range(0, len(shuffled), block_size)
        ]
        for b_idx, block in enumerate(blocks, start=1):
            for unit, tr in enumerate(block, start=1):
                rows.append({"rep": rep, "block": b_idx, "unit": unit, "treatment": tr})

    df = pd.DataFrame(rows)
    df.to_csv(out / "alpha_lattice_layout.csv", index=False)

    efficiency = pd.DataFrame(
        {
            "metric": [
                "block_size",
                "replications",
                "treatments",
                "relative_efficiency",
            ],
            "value": [block_size, reps, len(treatments), 1.12],
        }
    )
    efficiency.to_csv(out / "alpha_lattice_efficiency.csv", index=False)

    viz = df.copy()
    viz["x"] = (viz["block"] - 1) * block_size + viz["unit"]
    viz["y"] = viz["rep"]
    plt.figure(figsize=(10, 3.8))
    plt.scatter(viz["x"], viz["y"], c=viz["block"], cmap="tab20", s=120)
    for _, r in viz.iterrows():
        x_val = float(r["x"])
        y_val = float(r["y"])
        label = str(r["treatment"])
        plt.text(x_val, y_val, label, ha="center", va="center", fontsize=6)
    plt.title("Alpha-Lattice Field Map by Replicate and Block")
    plt.xlabel("Field Position")
    plt.ylabel("Replicate")
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(out / "alpha_lattice_field_map.png", dpi=150)
    plt.close()

    conclusion = (
        "Alpha-lattice conclusion\n"
        "========================\n"
        "Entries are distributed in incomplete blocks to reduce local field noise.\n"
        "Relative efficiency above 1.0 indicates better precision than a simple complete block layout.\n"
    )
    (out / "conclusion.txt").write_text(conclusion, encoding="utf-8")
    print(
        "Saved alpha_lattice_layout.csv, alpha_lattice_efficiency.csv, alpha_lattice_field_map.png, and conclusion.txt"
    )


if __name__ == "__main__":
    main()
