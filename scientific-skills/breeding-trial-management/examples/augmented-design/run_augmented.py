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
    summary = df.groupby(["block", "plot_type"]).size().to_frame("count").reset_index()
    summary.to_csv(out / "augmented_summary.csv", index=False)

    viz = df.copy()
    viz["x"] = viz.groupby("block").cumcount() + 1
    viz["y"] = viz["block"]
    colors = np.where(viz["plot_type"] == "check", "#1f77b4", "#ff7f0e")
    plt.figure(figsize=(10, 4))
    plt.scatter(viz["x"], viz["y"], c=colors, s=110)
    for _, r in viz.iterrows():
        plt.text(
            float(r["x"]),
            float(r["y"]),
            str(r["entry"]),
            ha="center",
            va="center",
            fontsize=6,
        )
    plt.title("Augmented Trial Field Map (Checks vs Unreplicated Entries)")
    plt.xlabel("Plot Position Within Block")
    plt.ylabel("Block")
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(out / "augmented_field_map.png", dpi=150)
    plt.close()

    conclusion = (
        "Augmented design conclusion\n"
        "===========================\n"
        "Checks are repeated in each block to anchor comparisons for many unreplicated entries.\n"
        "This layout supports grower-style early screening with limited seed while preserving comparability.\n"
    )
    (out / "conclusion.txt").write_text(conclusion, encoding="utf-8")
    print(
        "Saved augmented_layout.csv, augmented_summary.csv, augmented_field_map.png, and conclusion.txt"
    )


if __name__ == "__main__":
    main()
