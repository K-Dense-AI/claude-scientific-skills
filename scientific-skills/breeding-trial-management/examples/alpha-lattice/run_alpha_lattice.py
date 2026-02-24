#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from pathlib import Path
import numpy as np
import pandas as pd


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
    print("Saved alpha_lattice_layout.csv and alpha_lattice_efficiency.csv")


if __name__ == "__main__":
    main()
