#!/usr/bin/env python3
"""
TimesFM Anomaly Detection Example

This example demonstrates how to use TimesFM's quantile forecasts for
anomaly detection. The approach:
1. Forecast with quantile intervals (10th-90th percentiles)
2. Compare actual values against prediction intervals
3. Flag values outside intervals as anomalies

TimesFM does NOT have built-in anomaly detection, but the quantile
forecasts provide natural anomaly detection via prediction intervals.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import timesfm

# Configuration
HORIZON = 12  # Forecast horizon
ANOMALY_THRESHOLD_WARNING = 0.80  # Outside 80% CI = warning
ANOMALY_THRESHOLD_CRITICAL = 0.90  # Outside 90% CI = critical

EXAMPLE_DIR = Path(__file__).parent
DATA_FILE = (
    Path(__file__).parent.parent / "global-temperature" / "temperature_anomaly.csv"
)
OUTPUT_DIR = EXAMPLE_DIR / "output"


def inject_anomalies(
    values: np.ndarray, n_anomalies: int = 3, seed: int = 42
) -> tuple[np.ndarray, list[int]]:
    """Inject synthetic anomalies into the data for demonstration."""
    rng = np.random.default_rng(seed)
    anomaly_indices = rng.choice(len(values), size=n_anomalies, replace=False).tolist()

    anomalous_values = values.copy()
    for idx in anomaly_indices:
        # Inject spike or dip (±40-60% of value)
        multiplier = rng.choice([0.4, 0.6]) * rng.choice([1, -1])
        anomalous_values[idx] = values[idx] * (1 + multiplier)

    return anomalous_values, sorted(anomaly_indices)


def main() -> None:
    print("=" * 60)
    print("  TIMESFM ANOMALY DETECTION DEMO")
    print("=" * 60)

    OUTPUT_DIR.mkdir(exist_ok=True)

    # Load temperature data
    print("\n📊 Loading temperature anomaly data...")
    df = pd.read_csv(DATA_FILE, parse_dates=["date"])
    df = df.sort_values("date").reset_index(drop=True)

    # Split into context (first 24 months) and test (last 12 months)
    context_values = df["anomaly_c"].values[:24].astype(np.float32)
    actual_future = df["anomaly_c"].values[24:36].astype(np.float32)
    dates_future = df["date"].values[24:36]

    print(f"   Context: 24 months (2022-01 to 2023-12)")
    print(f"   Test: 12 months (2024-01 to 2024-12)")

    # Inject anomalies into test data for demonstration
    print("\n🔬 Injecting synthetic anomalies for demonstration...")
    test_values_with_anomalies, anomaly_indices = inject_anomalies(
        actual_future, n_anomalies=3
    )
    print(f"   Injected anomalies at months: {anomaly_indices}")

    # Load TimesFM
    print("\n🤖 Loading TimesFM 1.0 (200M) PyTorch...")
    hparams = timesfm.TimesFmHparams(horizon_len=HORIZON)
    checkpoint = timesfm.TimesFmCheckpoint(
        huggingface_repo_id="google/timesfm-1.0-200m-pytorch"
    )
    model = timesfm.TimesFm(hparams=hparams, checkpoint=checkpoint)

    # Forecast with quantiles
    print("\n📈 Forecasting with quantile intervals...")
    point_forecast, quantile_forecast = model.forecast(
        [context_values],
        freq=[0],
    )

    # Extract quantiles
    # quantile_forecast shape: (1, 12, 10) - [mean, q10, q20, ..., q90]
    point = point_forecast[0]
    q10 = quantile_forecast[0, :, 0]  # 10th percentile
    q20 = quantile_forecast[0, :, 1]  # 20th percentile
    q50 = quantile_forecast[0, :, 4]  # 50th percentile (median)
    q80 = quantile_forecast[0, :, 7]  # 80th percentile
    q90 = quantile_forecast[0, :, 8]  # 90th percentile

    print(f"   Forecast mean: {point.mean():.3f}°C")
    print(f"   90% CI width: {(q90 - q10).mean():.3f}°C (avg)")

    # Detect anomalies
    print("\n🔍 Detecting anomalies...")
    anomalies = []
    for i, (actual, lower_80, upper_80, lower_90, upper_90) in enumerate(
        zip(test_values_with_anomalies, q20, q80, q10, q90)
    ):
        month = dates_future[i]
        month_str = pd.to_datetime(month).strftime("%Y-%m")

        if actual < lower_90 or actual > upper_90:
            severity = "CRITICAL"
            threshold = "90% CI"
            color = "red"
        elif actual < lower_80 or actual > upper_80:
            severity = "WARNING"
            threshold = "80% CI"
            color = "orange"
        else:
            severity = "NORMAL"
            threshold = "within bounds"
            color = "green"

        anomalies.append(
            {
                "month": month_str,
                "actual": float(actual),
                "forecast": float(point[i]),
                "lower_80": float(lower_80),
                "upper_80": float(upper_80),
                "lower_90": float(lower_90),
                "upper_90": float(upper_90),
                "severity": severity,
                "threshold": threshold,
                "color": color,
            }
        )

        if severity != "NORMAL":
            deviation = abs(actual - point[i])
            print(
                f"   [{severity}] {month_str}: {actual:.2f}°C (forecast: {point[i]:.2f}°C, deviation: {deviation:.2f}°C)"
            )

    # Create visualization
    print("\n📊 Creating anomaly visualization...")

    fig, axes = plt.subplots(2, 1, figsize=(14, 10))

    # Plot 1: Full time series with forecast and anomalies
    ax1 = axes[0]

    # Historical data
    historical_dates = df["date"].values[:24]
    ax1.plot(
        historical_dates,
        context_values,
        "b-",
        linewidth=2,
        label="Historical Data",
        marker="o",
        markersize=4,
    )

    # Actual future (with anomalies)
    ax1.plot(
        dates_future,
        actual_future,
        "g--",
        linewidth=1.5,
        label="Actual (clean)",
        alpha=0.5,
    )
    ax1.plot(
        dates_future,
        test_values_with_anomalies,
        "ko",
        markersize=8,
        label="Actual (with anomalies)",
        alpha=0.7,
    )

    # Forecast
    ax1.plot(
        dates_future,
        point,
        "r-",
        linewidth=2,
        label="Forecast (median)",
        marker="s",
        markersize=6,
    )

    # 90% CI
    ax1.fill_between(dates_future, q10, q90, alpha=0.2, color="red", label="90% CI")

    # 80% CI
    ax1.fill_between(dates_future, q20, q80, alpha=0.3, color="red", label="80% CI")

    # Highlight anomalies
    for anomaly in anomalies:
        if anomaly["severity"] != "NORMAL":
            idx = [pd.to_datetime(d).strftime("%Y-%m") for d in dates_future].index(
                anomaly["month"]
            )
            ax1.scatter(
                [dates_future[idx]],
                [test_values_with_anomalies[idx]],
                c=anomaly["color"],
                s=200,
                marker="x" if anomaly["severity"] == "CRITICAL" else "^",
                linewidths=3,
                zorder=5,
            )

    ax1.set_xlabel("Date", fontsize=12)
    ax1.set_ylabel("Temperature Anomaly (°C)", fontsize=12)
    ax1.set_title(
        "TimesFM Anomaly Detection: Forecast Intervals Method",
        fontsize=14,
        fontweight="bold",
    )
    ax1.legend(loc="upper left", fontsize=10)
    ax1.grid(True, alpha=0.3)

    # Add annotation for anomalies
    ax1.annotate(
        "× = Critical (outside 90% CI)\n▲ = Warning (outside 80% CI)",
        xy=(0.98, 0.02),
        xycoords="axes fraction",
        ha="right",
        va="bottom",
        fontsize=10,
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.8),
    )

    # Plot 2: Deviation from forecast with thresholds
    ax2 = axes[1]

    deviation = test_values_with_anomalies - point
    lower_90_dev = q10 - point
    upper_90_dev = q90 - point
    lower_80_dev = q20 - point
    upper_80_dev = q80 - point

    months = [pd.to_datetime(d).strftime("%Y-%m") for d in dates_future]
    x = np.arange(len(months))

    # Threshold bands
    ax2.fill_between(
        x, lower_90_dev, upper_90_dev, alpha=0.2, color="red", label="90% CI bounds"
    )
    ax2.fill_between(
        x, lower_80_dev, upper_80_dev, alpha=0.3, color="red", label="80% CI bounds"
    )

    # Deviation bars
    colors = [
        "red"
        if d < lower_90_dev[i] or d > upper_90_dev[i]
        else "orange"
        if d < lower_80_dev[i] or d > upper_80_dev[i]
        else "green"
        for i, d in enumerate(deviation)
    ]
    ax2.bar(x, deviation, color=colors, alpha=0.7, edgecolor="black", linewidth=0.5)

    # Zero line
    ax2.axhline(y=0, color="black", linestyle="-", linewidth=1)

    ax2.set_xlabel("Month", fontsize=12)
    ax2.set_ylabel("Deviation from Forecast (°C)", fontsize=12)
    ax2.set_title(
        "Deviation from Forecast with Anomaly Thresholds",
        fontsize=14,
        fontweight="bold",
    )
    ax2.set_xticks(x)
    ax2.set_xticklabels(months, rotation=45, ha="right")
    ax2.legend(loc="upper right", fontsize=10)
    ax2.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()

    output_path = OUTPUT_DIR / "anomaly_detection.png"
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"   Saved: {output_path}")
    plt.close()

    # Save results
    results = {
        "method": "quantile_intervals",
        "description": "Anomaly detection using TimesFM quantile forecasts as prediction intervals",
        "thresholds": {
            "warning": f"Outside {ANOMALY_THRESHOLD_WARNING * 100:.0f}% CI (q20-q80)",
            "critical": f"Outside {ANOMALY_THRESHOLD_CRITICAL * 100:.0f}% CI (q10-q90)",
        },
        "anomalies": anomalies,
        "summary": {
            "total_points": len(anomalies),
            "critical": sum(1 for a in anomalies if a["severity"] == "CRITICAL"),
            "warning": sum(1 for a in anomalies if a["severity"] == "WARNING"),
            "normal": sum(1 for a in anomalies if a["severity"] == "NORMAL"),
        },
    }

    results_path = OUTPUT_DIR / "anomaly_detection.json"
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"   Saved: {results_path}")

    # Print summary
    print("\n" + "=" * 60)
    print("  ✅ ANOMALY DETECTION COMPLETE")
    print("=" * 60)
    print(f"\n📊 Summary:")
    print(f"   Total test points: {results['summary']['total_points']}")
    print(f"   Critical anomalies: {results['summary']['critical']} (outside 90% CI)")
    print(f"   Warnings: {results['summary']['warning']} (outside 80% CI)")
    print(f"   Normal: {results['summary']['normal']}")

    print("\n💡 How It Works:")
    print("   1. TimesFM forecasts with quantile intervals (q10, q20, ..., q90)")
    print("   2. If actual value falls outside 90% CI → CRITICAL anomaly")
    print("   3. If actual value falls outside 80% CI → WARNING")
    print("   4. Otherwise → NORMAL")

    print("\n📁 Output Files:")
    print(f"   {output_path}")
    print(f"   {results_path}")


if __name__ == "__main__":
    main()
