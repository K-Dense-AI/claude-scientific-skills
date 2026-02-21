#!/usr/bin/env python3
"""
TimesFM Anomaly Detection Example

Demonstrates using TimesFM quantile forecasts as prediction intervals
for anomaly detection. Approach:
  1. Use 36 months of real data as context
  2. Create synthetic 12-month future (natural continuation of trend)
  3. Inject 3 clear anomalies into that future
  4. Forecast with quantile intervals → flag anomalies by severity

TimesFM has NO built-in anomaly detection. Quantile forecasts provide
natural prediction intervals — values outside them are statistically unusual.

Quantile index reference (index 0 = mean, 1-9 = q10-q90):
  80% PI = q10 (idx 1) to q90 (idx 9)
  60% PI = q20 (idx 2) to q80 (idx 8)
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import timesfm

# Configuration
HORIZON = 12  # Forecast horizon (months)
DATA_FILE = (
    Path(__file__).parent.parent / "global-temperature" / "temperature_anomaly.csv"
)
OUTPUT_DIR = Path(__file__).parent / "output"

# Anomaly thresholds using available quantile outputs
# 80% PI = q10-q90  →  "critical" if outside
# 60% PI = q20-q80  →  "warning" if outside
IDX_Q10, IDX_Q20, IDX_Q80, IDX_Q90 = 1, 2, 8, 9


def build_synthetic_future(
    context: np.ndarray, n: int, seed: int = 42
) -> tuple[np.ndarray, list[int]]:
    """Build synthetic future that looks like a natural continuation.

    Takes the mean/std of the last 6 context months as the baseline,
    then injects 3 clear anomalies (2 high, 1 low) at fixed positions.
    """
    rng = np.random.default_rng(seed)
    recent_mean = float(context[-6:].mean())
    recent_std = float(context[-6:].std())

    # Natural-looking continuation: small gaussian noise around recent mean
    future = recent_mean + rng.normal(0, recent_std * 0.4, n).astype(np.float32)

    # Inject 3 unmistakable anomalies
    anomaly_cfg = [
        (2, +0.55),  # month 3  — large spike up
        (7, -0.50),  # month 8  — large dip down
        (10, +0.48),  # month 11 — spike up
    ]
    anomaly_indices = []
    for idx, delta in anomaly_cfg:
        future[idx] = recent_mean + delta
        anomaly_indices.append(idx)

    return future, sorted(anomaly_indices)


def main() -> None:
    print("=" * 60)
    print("  TIMESFM ANOMALY DETECTION DEMO")
    print("=" * 60)

    OUTPUT_DIR.mkdir(exist_ok=True)

    # ── Load all 36 months as context ─────────────────────────────
    print("\n📊 Loading temperature data (all 36 months as context)...")
    df = pd.read_csv(DATA_FILE, parse_dates=["date"])
    df = df.sort_values("date").reset_index(drop=True)
    context_values = df["anomaly_c"].values.astype(np.float32)  # all 36 months
    context_dates = df["date"].tolist()

    print(
        f"   Context: {len(context_values)} months ({context_dates[0].strftime('%Y-%m')} → {context_dates[-1].strftime('%Y-%m')})"
    )

    # ── Build synthetic future with known anomalies ────────────────
    print("\n🔬 Building synthetic 12-month future with injected anomalies...")
    future_values, injected_at = build_synthetic_future(context_values, HORIZON)
    future_dates = pd.date_range(
        start=context_dates[-1] + pd.DateOffset(months=1),
        periods=HORIZON,
        freq="MS",
    )
    print(
        f"   Anomalies injected at months: {[future_dates[i].strftime('%Y-%m') for i in injected_at]}"
    )

    # ── Load TimesFM and forecast ──────────────────────────────────
    print("\n🤖 Loading TimesFM 1.0 (200M) PyTorch...")
    hparams = timesfm.TimesFmHparams(horizon_len=HORIZON)
    checkpoint = timesfm.TimesFmCheckpoint(
        huggingface_repo_id="google/timesfm-1.0-200m-pytorch"
    )
    model = timesfm.TimesFm(hparams=hparams, checkpoint=checkpoint)

    print("\n📈 Forecasting...")
    point_fc, quant_fc = model.forecast([context_values], freq=[0])

    # quantile_forecast shape: (1, horizon, 10)
    #   index 0 = mean, index 1 = q10, ..., index 9 = q90
    point = point_fc[0]  # shape (12,)
    q10 = quant_fc[0, :, IDX_Q10]  # 10th pct
    q20 = quant_fc[0, :, IDX_Q20]  # 20th pct
    q80 = quant_fc[0, :, IDX_Q80]  # 80th pct
    q90 = quant_fc[0, :, IDX_Q90]  # 90th pct

    print(f"   Forecast mean: {point.mean():.3f}°C")
    print(f"   80% PI width:  {(q90 - q10).mean():.3f}°C (avg)")

    # ── Detect anomalies ───────────────────────────────────────────
    print("\n🔍 Detecting anomalies...")
    records = []
    for i, (actual, fcast, lo60, hi60, lo80, hi80) in enumerate(
        zip(future_values, point, q20, q80, q10, q90)
    ):
        month = future_dates[i].strftime("%Y-%m")

        if actual < lo80 or actual > hi80:
            severity = "CRITICAL"  # outside 80% PI
        elif actual < lo60 or actual > hi60:
            severity = "WARNING"  # outside 60% PI
        else:
            severity = "NORMAL"

        records.append(
            {
                "month": month,
                "actual": round(float(actual), 4),
                "forecast": round(float(fcast), 4),
                "lower_60pi": round(float(lo60), 4),
                "upper_60pi": round(float(hi60), 4),
                "lower_80pi": round(float(lo80), 4),
                "upper_80pi": round(float(hi80), 4),
                "severity": severity,
                "injected": (i in injected_at),
            }
        )

        if severity != "NORMAL":
            dev = actual - fcast
            print(
                f"   [{severity}] {month}: actual={actual:.2f}  forecast={fcast:.2f}  Δ={dev:+.2f}°C"
            )

    # ── Visualise ─────────────────────────────────────────────────
    print("\n📊 Creating visualization...")

    fig, axes = plt.subplots(2, 1, figsize=(13, 9))

    clr = {"CRITICAL": "red", "WARNING": "orange", "NORMAL": "steelblue"}

    # — Panel 1: full series ———————————————————————————————————————
    ax = axes[0]
    ax.plot(
        context_dates,
        context_values,
        "b-",
        lw=2,
        marker="o",
        ms=4,
        label="Context (36 months)",
    )
    ax.fill_between(
        future_dates, q10, q90, alpha=0.18, color="tomato", label="80% PI (q10–q90)"
    )
    ax.fill_between(
        future_dates, q20, q80, alpha=0.28, color="tomato", label="60% PI (q20–q80)"
    )
    ax.plot(future_dates, point, "r-", lw=2, marker="s", ms=5, label="Forecast")
    ax.plot(
        future_dates,
        future_values,
        "k--",
        lw=1.3,
        alpha=0.5,
        label="Synthetic future (clean)",
    )

    # mark anomalies
    for rec in records:
        if rec["severity"] != "NORMAL":
            dt = pd.to_datetime(rec["month"])
            c = "red" if rec["severity"] == "CRITICAL" else "orange"
            mk = "X" if rec["severity"] == "CRITICAL" else "^"
            ax.scatter(
                [dt], [rec["actual"]], c=c, s=220, marker=mk, zorder=6, linewidths=2
            )

    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right")
    ax.set_ylabel("Temperature Anomaly (°C)", fontsize=11)
    ax.set_title(
        "TimesFM Anomaly Detection — Prediction Interval Method",
        fontsize=13,
        fontweight="bold",
    )
    ax.legend(loc="upper left", fontsize=9, ncol=2)
    ax.grid(True, alpha=0.25)
    ax.annotate(
        "X = Critical (outside 80% PI)\n▲ = Warning (outside 60% PI)",
        xy=(0.98, 0.04),
        xycoords="axes fraction",
        ha="right",
        fontsize=9,
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    # — Panel 2: deviation bars ———————————————————————————————————
    ax2 = axes[1]
    deviations = future_values - point
    lo80_dev = q10 - point
    hi80_dev = q90 - point
    lo60_dev = q20 - point
    hi60_dev = q80 - point
    x = np.arange(HORIZON)

    ax2.fill_between(x, lo80_dev, hi80_dev, alpha=0.15, color="tomato", label="80% PI")
    ax2.fill_between(x, lo60_dev, hi60_dev, alpha=0.25, color="tomato", label="60% PI")
    bar_colors = [clr[r["severity"]] for r in records]
    ax2.bar(x, deviations, color=bar_colors, alpha=0.75, edgecolor="black", lw=0.5)
    ax2.axhline(0, color="black", lw=1)

    ax2.set_xticks(x)
    ax2.set_xticklabels(
        [r["month"] for r in records], rotation=45, ha="right", fontsize=9
    )
    ax2.set_ylabel("Δ from Forecast (°C)", fontsize=11)
    ax2.set_title(
        "Deviation from Forecast with Anomaly Thresholds",
        fontsize=13,
        fontweight="bold",
    )
    ax2.legend(loc="upper right", fontsize=9)
    ax2.grid(True, alpha=0.25, axis="y")

    plt.tight_layout()
    png_path = OUTPUT_DIR / "anomaly_detection.png"
    plt.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"   Saved: {png_path}")

    # ── Save JSON results ──────────────────────────────────────────
    summary = {
        "total": len(records),
        "critical": sum(1 for r in records if r["severity"] == "CRITICAL"),
        "warning": sum(1 for r in records if r["severity"] == "WARNING"),
        "normal": sum(1 for r in records if r["severity"] == "NORMAL"),
    }
    out = {
        "method": "quantile_prediction_intervals",
        "description": (
            "Anomaly detection via TimesFM quantile forecasts. "
            "80% PI = q10–q90 (CRITICAL if violated). "
            "60% PI = q20–q80 (WARNING if violated)."
        ),
        "context": "36 months of real NOAA temperature anomaly data (2022-2024)",
        "future": "12 synthetic months with 3 injected anomalies",
        "quantile_indices": {"q10": 1, "q20": 2, "q80": 8, "q90": 9},
        "summary": summary,
        "detections": records,
    }
    json_path = OUTPUT_DIR / "anomaly_detection.json"
    with open(json_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"   Saved: {json_path}")

    # ── Summary ────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("  ✅ ANOMALY DETECTION COMPLETE")
    print("=" * 60)
    print(f"\n   Total future points : {summary['total']}")
    print(f"   Critical (80% PI)   : {summary['critical']}")
    print(f"   Warning  (60% PI)   : {summary['warning']}")
    print(f"   Normal              : {summary['normal']}")


if __name__ == "__main__":
    main()
