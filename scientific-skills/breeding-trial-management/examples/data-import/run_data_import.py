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
    standardized = pd.DataFrame(
        standardized[
            [
                "plot_id",
                "genotype",
                "environment",
                "replicate",
                "yield_kg_ha",
                "moisture_pct",
            ]
        ]
    )

    standardized.to_csv(out / "standardized_phenotypes.csv", index=False)

    env_centers = {"E1": (-97.0, 40.8), "E2": (-96.3, 41.0), "E3": (-95.8, 40.6)}
    geo = standardized.copy()
    env_list = [str(e) for e in geo["environment"].tolist()]
    geo["lon"] = [env_centers[e][0] for e in env_list] + rng.normal(
        0.0, 0.03, size=len(geo)
    )
    geo["lat"] = [env_centers[e][1] for e in env_list] + rng.normal(
        0.0, 0.03, size=len(geo)
    )
    geo.to_csv(out / "standardized_sites.csv", index=False)

    plt.figure(figsize=(7, 4.8))
    s = plt.scatter(
        geo["lon"],
        geo["lat"],
        c=geo["yield_kg_ha"],
        cmap="YlOrRd",
        s=90,
        edgecolor="black",
        linewidth=0.3,
    )
    plt.colorbar(s, label="Yield (kg/ha)")
    plt.title("Imported Phenotype Sites (Mock Geospatial View)")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(out / "standardized_sites_map.png", dpi=150)
    plt.close()

    report_lines = [
        "Validation Report",
        "=================",
        f"Rows imported: {len(raw)}",
        f"Missing yield values imputed: {raw['Yield(kg/ha)'].isna().sum()}",
        f"Unique genotypes: {standardized['genotype'].nunique()}",
        f"Unique environments: {standardized['environment'].nunique()}",
        f"Replicate values valid: {standardized['replicate'].isin([1, 2, 3]).all()}",
        "Conclusion: Data are standardized and ready for model fitting or breeding decisions.",
    ]
    (out / "validation_report.txt").write_text("\n".join(report_lines) + "\n")

    print("Saved raw import, standardized files, site map, and validation report")


if __name__ == "__main__":
    main()
