#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def main():
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)

    rows = []
    for r in range(1, 11):
        for c in range(1, 7):
            plot_id = f"P{r:02d}{c:02d}"
            rows.append(
                {
                    "plot_id": plot_id,
                    "row": r,
                    "col": c,
                    "entry": f"G{((r - 1) * 6 + c) % 24 + 1:02d}",
                    "qr_label": f"QR::{plot_id}",
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(out / "fieldbook.csv", index=False)
    df[["plot_id", "qr_label"]].to_csv(out / "qr_labels.csv", index=False)

    plt.figure(figsize=(8, 6))
    plt.scatter(
        df["col"], df["row"], c=df["row"], cmap="YlGn", s=130, edgecolor="black"
    )
    for _, r in df.iloc[::4].iterrows():
        plt.text(float(r["col"]) + 0.08, float(r["row"]), str(r["entry"]), fontsize=6)
    plt.gca().invert_yaxis()
    plt.title("Field Book Plot Map (Grower View)")
    plt.xlabel("Column")
    plt.ylabel("Row")
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(out / "fieldbook_plot_map.png", dpi=150)
    plt.close()

    conclusion = (
        "Field book conclusion\n"
        "====================\n"
        "The map aligns plot IDs, genotype entries, and label payloads for field crews.\n"
        "This supports consistent scoring and sampling in grower-scale trial operations.\n"
    )
    (out / "conclusion.txt").write_text(conclusion, encoding="utf-8")
    print(
        "Saved fieldbook.csv, qr_labels.csv, fieldbook_plot_map.png, and conclusion.txt"
    )


if __name__ == "__main__":
    main()
