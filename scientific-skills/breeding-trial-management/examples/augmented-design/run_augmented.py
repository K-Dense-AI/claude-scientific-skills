#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from pathlib import Path
import numpy as np
import pandas as pd


def main():
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)
    rng = np.random.default_rng(9)

    checks = [f"C{i + 1:02d}" for i in range(4)]
    entries = [f"E{i + 1:03d}" for i in range(30)]
    blocks = 6
    rows = []
    e_idx = 0
    for b in range(1, blocks + 1):
        for c in checks:
            rows.append({"block": b, "plot_type": "check", "entry": c})
        for _ in range(5):
            rows.append({"block": b, "plot_type": "unrep", "entry": entries[e_idx]})
            e_idx += 1
    df = pd.DataFrame(rows)
    df["plot"] = np.arange(1, len(df) + 1)
    df.to_csv(out / "augmented_layout.csv", index=False)
    df.groupby(["block", "plot_type"]).size().reset_index(name="count").to_csv(
        out / "augmented_summary.csv", index=False
    )
    print("Saved augmented_layout.csv and augmented_summary.csv")


if __name__ == "__main__":
    main()
