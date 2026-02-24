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

    n = 36
    raw = pd.DataFrame(
        {
            "PlotID": [f"P{i + 1:03d}" for i in range(n)],
            "Geno Name": [f"G{i % 12 + 1:02d}" for i in range(n)],
            "Env": rng.choice(["E1", "E2", "E3"], size=n),
            "Yield(kg/ha)": rng.normal(5400, 520, size=n).round(1),
            "Moisture%": rng.normal(13.5, 1.1, size=n).round(2),
            "Rep": rng.choice([1, 2, 3], size=n),
        }
    )

    raw.loc[[5, 18], "Yield(kg/ha)"] = np.nan
    raw.to_csv(out / "raw_phenotypes.csv", index=False)

    standardized = raw.rename(
        columns={
            "PlotID": "plot_id",
            "Geno Name": "genotype",
            "Env": "environment",
            "Yield(kg/ha)": "yield_kg_ha",
            "Moisture%": "moisture_pct",
            "Rep": "replicate",
        }
    )

    standardized["yield_kg_ha"] = standardized["yield_kg_ha"].fillna(
        standardized["yield_kg_ha"].median()
    )
    standardized["replicate"] = standardized["replicate"].astype(int)
    standardized = standardized[
        [
            "plot_id",
            "genotype",
            "environment",
            "replicate",
            "yield_kg_ha",
            "moisture_pct",
        ]
    ]

    standardized.to_csv(out / "standardized_phenotypes.csv", index=False)

    report_lines = [
        "Validation Report",
        "=================",
        f"Rows imported: {len(raw)}",
        f"Missing yield values imputed: {raw['Yield(kg/ha)'].isna().sum()}",
        f"Unique genotypes: {standardized['genotype'].nunique()}",
        f"Unique environments: {standardized['environment'].nunique()}",
        f"Replicate values valid: {standardized['replicate'].isin([1, 2, 3]).all()}",
    ]
    (out / "validation_report.txt").write_text("\n".join(report_lines) + "\n")

    print("Saved raw import, standardized phenotype file, and validation report")


if __name__ == "__main__":
    main()
