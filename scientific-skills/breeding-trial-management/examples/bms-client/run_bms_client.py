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
    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    draw_light_geopolitical_context(ax, sites["lon"], sites["lat"])
    ax.scatter(sites["lon"], sites["lat"], c=colors, s=130, zorder=2)
    for _, r in sites.iterrows():
        ax.annotate(
            str(r["trial_id"]),
            (float(r["lon"]), float(r["lat"])),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8,
        )
    ax.set_title("BMS Trial Footprint")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(alpha=0.18)
    fig.tight_layout()
    fig.savefig(out / "bms_trial_site_map.png", dpi=150)
    plt.close(fig)

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
