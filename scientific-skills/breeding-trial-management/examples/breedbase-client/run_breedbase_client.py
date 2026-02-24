#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


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
    plt.figure(figsize=(6.5, 4.5))
    plt.scatter(site_map["site_lon"], site_map["site_lat"], c="#1f77b4", s=120)
    for _, r in site_map.iterrows():
        plt.annotate(
            str(r["accession_id"]),
            (float(r["site_lon"]), float(r["site_lat"])),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8,
        )
    plt.title("Breedbase Accessions by Trial Site (Mock Geospatial)")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(out / "breedbase_accession_site_map.png", dpi=150)
    plt.close()

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
