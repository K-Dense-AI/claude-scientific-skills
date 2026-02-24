#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def main():
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)

    catalog = pd.DataFrame(
        {
            "dataset": ["Ames panel", "Maize NAM", "Rice 3K"],
            "source": [
                "https://www.maizegenetics.net/",
                "https://www.panzea.org/",
                "https://snp-seek.irri.org/",
            ],
            "status": ["optional", "optional", "optional"],
            "notes": [
                "Use for larger-scale GWAS benchmarking",
                "Use for NAM-style QTL workflows",
                "Use for rice population analyses",
            ],
        }
    )
    catalog.to_csv(out / "real_dataset_catalog.csv", index=False)
    (out / "download_instructions.txt").write_text(
        "Optional real datasets are listed in real_dataset_catalog.csv.\n"
        "Download manually to comply with each source's terms and license.\n",
        encoding="utf-8",
    )

    plt.figure(figsize=(7, 4))
    y = range(len(catalog))
    plt.barh(y, [1] * len(catalog), color="#4c78a8")
    plt.yticks(y, list(catalog["dataset"]))
    plt.xticks([])
    for i, note in enumerate(catalog["notes"]):
        plt.text(0.02, i, note, va="center", fontsize=8)
    plt.title("Optional Real Dataset Catalog")
    plt.tight_layout()
    plt.savefig(out / "real_dataset_catalog_overview.png", dpi=150)
    plt.close()

    print(
        "Saved real_dataset_catalog.csv, download_instructions.txt, and real_dataset_catalog_overview.png"
    )


if __name__ == "__main__":
    main()
