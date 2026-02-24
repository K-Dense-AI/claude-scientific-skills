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
    latest["x"] = [int(p[1:]) % 8 for p in latest["plot_id"]]
    latest["y"] = [int(p[1:]) // 8 for p in latest["plot_id"]]
    plt.figure(figsize=(7, 5.2))
    sc = plt.scatter(
        latest["x"],
        latest["y"],
        c=latest["ndvi"],
        cmap="YlGn",
        s=150,
        edgecolor="black",
        linewidth=0.3,
    )
    plt.colorbar(sc, label="NDVI")
    plt.title("IoT Field Snapshot for Sync Dispatch")
    plt.xlabel("Field column")
    plt.ylabel("Field row")
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(out / "iot_field_status_map.png", dpi=160)
    plt.close()

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
