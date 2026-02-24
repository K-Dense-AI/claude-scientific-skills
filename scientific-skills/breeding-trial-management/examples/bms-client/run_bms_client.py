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

    trials = pd.DataFrame(
        {
            "trial_id": ["BMS-T01", "BMS-T02"],
            "location": ["Site-A", "Site-B"],
            "season": ["2026A", "2026A"],
            "status": ["retrieved_mock", "retrieved_mock"],
        }
    )
    trials.to_csv(out / "bms_trials_read.csv", index=False)

    update = pd.DataFrame(
        {
            "trial_id": ["BMS-T03"],
            "operation": ["create_mock"],
            "result": ["success"],
        }
    )
    update.to_csv(out / "bms_trials_write_result.csv", index=False)

    sites = pd.DataFrame(
        {
            "trial_id": ["BMS-T01", "BMS-T02", "BMS-T03"],
            "lon": [-97.1, -96.4, -95.7],
            "lat": [40.6, 40.8, 41.2],
            "status": ["retrieved", "retrieved", "created"],
        }
    )
    sites.to_csv(out / "bms_trial_sites.csv", index=False)

    colors = np.where(sites["status"] == "created", "#ff7f0e", "#2ca02c")
    plt.figure(figsize=(6.5, 4.5))
    plt.scatter(sites["lon"], sites["lat"], c=colors, s=130)
    for _, r in sites.iterrows():
        plt.annotate(
            str(r["trial_id"]),
            (float(r["lon"]), float(r["lat"])),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8,
        )
    plt.title("BMS Trial Footprint (Mock Geospatial)")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(out / "bms_trial_site_map.png", dpi=150)
    plt.close()

    conclusion = (
        "BMS client conclusion\n"
        "=====================\n"
        "Mock BMS trial transactions are paired with site footprints for operational review.\n"
        "This format supports quick checks by breeders, agronomists, and data managers.\n"
    )
    (out / "conclusion.txt").write_text(conclusion, encoding="utf-8")
    print("Saved BMS mock read/write outputs, site map, and conclusion")


if __name__ == "__main__":
    main()
