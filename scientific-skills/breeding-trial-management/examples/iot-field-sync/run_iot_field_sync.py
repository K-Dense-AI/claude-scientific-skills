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
    dx = max(0.12, (max_lon - min_lon) * 0.55)
    dy = max(0.08, (max_lat - min_lat) * 0.55)
    ax.set_facecolor("#f7fbff")
    ax.set_xlim(min_lon - dx, max_lon + dx)
    ax.set_ylim(min_lat - dy, max_lat + dy)
    for x in [min_lon - 0.05, (min_lon + max_lon) / 2, max_lon + 0.05]:
        ax.axvline(x, color="#d9e2ec", linewidth=0.8, linestyle="--", zorder=0)
    for y in [min_lat - 0.04, (min_lat + max_lat) / 2, max_lat + 0.04]:
        ax.axhline(y, color="#d9e2ec", linewidth=0.8, linestyle="--", zorder=0)
    ax.text(
        min_lon - dx + 0.01,
        max_lat + dy - 0.015,
        "Regional context",
        fontsize=8,
        color="#6b7280",
    )


def main():
    out = Path(__file__).parent / "output"
    out.mkdir(exist_ok=True)
    rng = np.random.default_rng(88)

    plots = [f"P{i + 1:03d}" for i in range(64)]
    records = []
    for t in range(12):
        ts = f"2026-06-{t + 1:02d}T08:00:00Z"
        for p in plots:
            records.append(
                {
                    "timestamp": ts,
                    "plot_id": p,
                    "soil_moisture": round(float(rng.normal(28, 4)), 2),
                    "canopy_temp_c": round(float(rng.normal(27, 2.2)), 2),
                    "ndvi": round(float(rng.normal(0.64, 0.06)), 3),
                }
            )
    stream = pd.DataFrame(records)
    stream.to_csv(out / "iot_sensor_stream.csv", index=False)

    daily = (
        stream.groupby(["timestamp", "plot_id"])
        .agg({"soil_moisture": "mean", "canopy_temp_c": "mean", "ndvi": "mean"})
        .reset_index()
    )
    daily["sync_target"] = np.where(np.arange(len(daily)) % 2 == 0, "breedbase", "bms")
    daily["sync_status"] = "synced"
    daily.to_csv(out / "iot_sync_log.csv", index=False)

    latest = daily[daily["timestamp"] == daily["timestamp"].max()].copy()
    latest["field_col"] = [int(p[1:]) % 8 for p in latest["plot_id"]]
    latest["field_row"] = [int(p[1:]) // 8 for p in latest["plot_id"]]
    latest["lon"] = -96.55 + latest["field_col"] * 0.012
    latest["lat"] = 40.72 + latest["field_row"] * 0.009

    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    draw_light_geopolitical_context(ax, latest["lon"], latest["lat"])
    sc = ax.scatter(
        latest["lon"],
        latest["lat"],
        c=latest["ndvi"],
        cmap="YlGn",
        s=150,
        edgecolor="black",
        linewidth=0.3,
        zorder=2,
    )
    for _, r in latest.iloc[::8].iterrows():
        ax.annotate(
            str(r["plot_id"]),
            (float(r["lon"]), float(r["lat"])),
            fontsize=7,
            xytext=(2, 2),
            textcoords="offset points",
        )
    fig.colorbar(sc, ax=ax, label="NDVI")
    ax.set_title("IoT Field Snapshot for Sync Dispatch")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.grid(alpha=0.18)
    fig.tight_layout()
    fig.savefig(out / "iot_field_status_map.png", dpi=160)
    plt.close(fig)

    conclusion = (
        "IoT sync conclusion\n"
        "===================\n"
        "Sensor observations are transformed into sync-ready records for Breedbase/BMS pipelines.\n"
        "This pattern can be wired to real MQTT/Kafka ingestion while preserving current offline demo behavior.\n"
    )
    (out / "conclusion.txt").write_text(conclusion, encoding="utf-8")
    print("Saved IoT stream, sync log, map, and conclusion")


if __name__ == "__main__":
    main()
