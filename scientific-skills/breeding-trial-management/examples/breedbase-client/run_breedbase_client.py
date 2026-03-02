#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from pathlib import Path
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

    accessions = pd.DataFrame(
        {
            "accession_id": ["BB001", "BB002", "BB003"],
            "name": ["Elite-A", "Elite-B", "Donor-X"],
            "status": ["active", "active", "active"],
        }
    )
    accessions.to_csv(out / "breedbase_accessions_read.csv", index=False)

    payload = pd.DataFrame(
        {
            "accession_id": ["BB004", "BB005"],
            "name": ["Line-4", "Line-5"],
            "write_result": ["mock_created", "mock_created"],
        }
    )
    payload.to_csv(out / "breedbase_accessions_write_result.csv", index=False)

    site_map = pd.DataFrame(
        {
            "accession_id": accessions["accession_id"],
            "site_lon": [-96.7, -96.2, -95.9],
            "site_lat": [41.3, 40.9, 41.0],
        }
    )
    site_map.to_csv(out / "breedbase_accession_sites.csv", index=False)
    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    draw_light_geopolitical_context(ax, site_map["site_lon"], site_map["site_lat"])
    ax.scatter(site_map["site_lon"], site_map["site_lat"], c="#1f77b4", s=120, zorder=2)
    for _, r in site_map.iterrows():
        ax.annotate(
            str(r["accession_id"]),
            (float(r["site_lon"]), float(r["site_lat"])),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8,
        )
    ax.set_title("Breedbase Accessions by Trial Site")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(alpha=0.18)
    fig.tight_layout()
    fig.savefig(out / "breedbase_accession_site_map.png", dpi=150)
    plt.close(fig)

    conclusion = (
        "Breedbase client conclusion\n"
        "===========================\n"
        "Read/write accession workflows are represented alongside site placement metadata.\n"
        "This helps teams confirm germplasm records and deployment locations before planting.\n"
    )
    (out / "conclusion.txt").write_text(conclusion, encoding="utf-8")
    print("Saved Breedbase mock outputs, site map, and conclusion")


if __name__ == "__main__":
    main()
