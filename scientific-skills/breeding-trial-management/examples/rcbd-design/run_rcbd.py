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

    genotypes = [f"G{i + 1:02d}" for i in range(12)]
    blocks = [f"B{i + 1}" for i in range(4)]
    rows = []
    for b in blocks:
        order = genotypes.copy()
        rng.shuffle(order)
        for p, g in enumerate(order, start=1):
            rows.append({"block": b, "plot": p, "genotype": g})

    df = pd.DataFrame(rows)
    df.to_csv(out / "rcbd_layout.csv", index=False)

    pivot = df.pivot(index="block", columns="plot", values="genotype")
    plt.figure(figsize=(10, 3))
    plt.imshow(np.arange(pivot.size).reshape(pivot.shape), cmap="tab20")
    plt.title("RCBD Plot Layout")
    plt.yticks(range(len(pivot.index)), list(pivot.index))
    plt.xticks(range(len(pivot.columns)), list(pivot.columns))
    plt.tight_layout()
    plt.savefig(out / "rcbd_layout.png", dpi=150)
    plt.close()

    conclusion = (
        "RCBD conclusion\n"
        "===============\n"
        "Blocking distributes each genotype across four blocks, which helps reduce field-position bias.\n"
        "This layout is suitable for early-stage comparisons where growers need fair, repeatable plot placement.\n"
    )
    (out / "conclusion.txt").write_text(conclusion, encoding="utf-8")
    print("Saved rcbd_layout.csv, rcbd_layout.png, and conclusion.txt")


if __name__ == "__main__":
    main()
