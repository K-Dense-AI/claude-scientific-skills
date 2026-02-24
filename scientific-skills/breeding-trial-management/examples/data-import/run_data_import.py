#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def draw_light_geopolitical_context(ax, lon, lat):
    min_lon, max_lon = float(min(lon)), float(max(lon))
    min_lat, max_lat = float(min(lat)), float(max(lat))
    dx = max(0.35, (max_lon - min_lon) * 0.6)
    dy = max(0.25, (max_lat - min_lat) * 0.6)
    ax.set_facecolor("#f7fbff")
    ax.set_xlim(min_lon - dx, max_lon + dx)
    ax.set_ylim(min_lat - dy, max_lat + dy)
    for x in [min_lon - 0.15, (min_lon + max_lon) / 2, max_lon + 0.15]:
        ax.axvline(x, color="#d9e2ec", linewidth=0.8, linestyle="--", zorder=0)
    for y in [min_lat - 0.1, (min_lat + max_lat) / 2, max_lat + 0.1]:
        ax.axhline(y, color="#d9e2ec", linewidth=0.8, linestyle="--", zorder=0)
    ax.text(
        min_lon - dx + 0.05,
        max_lat + dy - 0.08,
        "Regional context",
        fontsize=8,
        color="#6b7280",
    )


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

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    draw_light_geopolitical_context(ax, geo["lon"], geo["lat"])
    s = ax.scatter(
        geo["lon"],
        geo["lat"],
        c=geo["yield_kg_ha"],
        cmap="YlOrRd",
        s=90,
        edgecolor="black",
        linewidth=0.3,
        zorder=2,
    )
    fig.colorbar(s, ax=ax, label="Yield (kg/ha)")
    ax.set_title("Imported Phenotype Sites")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(alpha=0.18)
    fig.tight_layout()
    fig.savefig(out / "standardized_sites_map.png", dpi=150)
    plt.close(fig)

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
