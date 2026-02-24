#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from pathlib import Path
import pandas as pd


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
    print("Saved fieldbook.csv and qr_labels.csv")


if __name__ == "__main__":
    main()
