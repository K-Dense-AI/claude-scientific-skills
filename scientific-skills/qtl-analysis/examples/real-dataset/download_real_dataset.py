#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from pathlib import Path
import pandas as pd


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
    print("Saved real_dataset_catalog.csv and download_instructions.txt")


if __name__ == "__main__":
    main()
