#!/usr/bin/env python3
"""
Generate animated GIF showing forecast evolution.

Creates a GIF animation showing how the TimesFM forecast changes
as more historical data points are added.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pandas as pd
from PIL import Image

# Configuration
EXAMPLE_DIR = Path(__file__).parent
DATA_FILE = EXAMPLE_DIR / "animation_data.json"
OUTPUT_FILE = EXAMPLE_DIR / "forecast_animation.gif"
DURATION_MS = 500  # Time per frame in milliseconds


def create_frame(ax, step_data: dict, actual_data: dict, total_steps: int) -> None:
    """Create a single frame of the animation."""
    ax.clear()

    # Parse dates
    historical_dates = pd.to_datetime(step_data["historical_dates"])
    forecast_dates = pd.to_datetime(step_data["forecast_dates"])

    # Plot historical data
    ax.plot(
        historical_dates,
        step_data["historical_values"],
        color="#3b82f6",
        linewidth=2,
        marker="o",
        markersize=4,
        label="Historical",
    )

    # Plot 90% CI (outer)
    ax.fill_between(
        forecast_dates,
        step_data["q10"],
        step_data["q90"],
        alpha=0.1,
        color="#ef4444",
        label="90% CI",
    )

    # Plot 80% CI (inner)
    ax.fill_between(
        forecast_dates,
        step_data["q20"],
        step_data["q80"],
        alpha=0.2,
        color="#ef4444",
        label="80% CI",
    )

    # Plot forecast
    ax.plot(
        forecast_dates,
        step_data["point_forecast"],
        color="#ef4444",
        linewidth=2,
        marker="s",
        markersize=4,
        label="Forecast",
    )

    # Plot actual future data if available
    actual_dates = pd.to_datetime(actual_data["dates"])
    actual_values = actual_data["values"]

    # Find which actual points fall in forecast period
    forecast_start = forecast_dates[0]
    forecast_end = forecast_dates[-1]

    future_mask = (actual_dates >= forecast_start) & (actual_dates <= forecast_end)
    future_dates = actual_dates[future_mask]
    future_values = np.array(actual_values)[future_mask]

    if len(future_dates) > 0:
        ax.plot(
            future_dates,
            future_values,
            color="#10b981",
            linewidth=1,
            linestyle="--",
            marker="o",
            markersize=3,
            alpha=0.7,
            label="Actual (future)",
        )

    # Add vertical line at forecast boundary
    ax.axvline(
        x=historical_dates[-1],
        color="#6b7280",
        linestyle="--",
        linewidth=1,
        alpha=0.5,
    )

    # Formatting
    ax.set_xlabel("Date", fontsize=11)
    ax.set_ylabel("Temperature Anomaly (°C)", fontsize=11)
    ax.set_title(
        f"TimesFM Forecast Evolution\n"
        f"Step {step_data['step']}/{total_steps}: {step_data['n_points']} points → {step_data['last_historical_date']}",
        fontsize=13,
        fontweight="bold",
    )

    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=9)
    ax.set_ylim(0.5, 1.6)

    # Format x-axis
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=6))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right")


def main() -> None:
    print("=" * 60)
    print("  GENERATING ANIMATED GIF")
    print("=" * 60)

    # Load data
    with open(DATA_FILE) as f:
        data = json.load(f)

    total_steps = len(data["animation_steps"])
    print(f"\n📊 Total frames: {total_steps}")

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 6))

    # Generate frames
    frames = []

    for i, step in enumerate(data["animation_steps"]):
        print(f"   Frame {i + 1}/{total_steps}...")

        create_frame(ax, step, data["actual_data"], total_steps)

        # Save frame to buffer
        fig.canvas.draw()

        # Convert to PIL Image
        buf = fig.canvas.buffer_rgba()
        width, height = fig.canvas.get_width_height()
        img = Image.frombytes("RGBA", (width, height), buf)
        frames.append(img.convert("RGB"))

    plt.close()

    # Save as GIF
    print(f"\n💾 Saving GIF: {OUTPUT_FILE}")
    frames[0].save(
        OUTPUT_FILE,
        save_all=True,
        append_images=frames[1:],
        duration=DURATION_MS,
        loop=0,  # Loop forever
    )

    # Get file size
    size_kb = OUTPUT_FILE.stat().st_size / 1024
    print(f"   File size: {size_kb:.1f} KB")
    print(f"\n✅ Done!")


if __name__ == "__main__":
    main()
