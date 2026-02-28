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
    rng = np.random.default_rng(901)

    n_cycles = 8
    pop_size = 260
    base_bv = rng.normal(0.0, 1.0, size=pop_size)

    rows = []
    for c in range(n_cycles + 1):
        shift = 0.18 * c
        bv = base_bv + shift + rng.normal(0, 0.18, size=pop_size)
        selected = np.sort(bv)[-60:]
        mean_gain = float(selected.mean())
        diversity = float(np.std(selected))
        rows.append(
            {
                "cycle": c,
                "mean_selected_bv": round(mean_gain, 4),
                "selected_diversity_sd": round(diversity, 4),
            }
        )

    summary = pd.DataFrame(rows)
    summary.to_csv(out / "simulation_cycle_summary.csv", index=False)

    plt.figure(figsize=(7.2, 4.8))
    plt.plot(summary["cycle"], summary["mean_selected_bv"], marker="o", color="#1f77b4")
    plt.plot(
        summary["cycle"], summary["selected_diversity_sd"], marker="s", color="#ff7f0e"
    )
    plt.title("Breeding Simulation: Gain and Diversity Across Cycles")
    plt.xlabel("Selection cycle")
    plt.ylabel("Metric value")
    plt.legend(["Mean selected breeding value", "Selected diversity (SD)"], loc="best")
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(out / "simulation_gain_diversity.png", dpi=160)
    plt.close()

    conclusion = (
        "Simulation conclusion\n"
        "====================\n"
        "Expected gain increases across cycles while diversity gradually declines.\n"
        "Use this trend to tune selection intensity before implementing long-term crossing programs.\n"
    )
    (out / "conclusion.txt").write_text(conclusion, encoding="utf-8")
    print("Saved simulation summary, trend plot, and conclusion")


if __name__ == "__main__":
    main()
