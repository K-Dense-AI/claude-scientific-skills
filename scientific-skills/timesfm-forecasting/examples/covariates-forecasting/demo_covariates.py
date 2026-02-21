#!/usr/bin/env python3
"""
TimesFM Covariates (XReg) Example

This example demonstrates TimesFM's exogenous variable support through the
forecast_with_covariates() API. This requires `timesfm[xreg]` installation.

Covariate Types Supported:
- Dynamic Numerical: Time-varying numeric features (e.g., price, temperature)
- Dynamic Categorical: Time-varying categorical features (e.g., holiday, day_of_week)
- Static Numerical: Per-series numeric features (e.g., store_size)
- Static Categorical: Per-series categorical features (e.g., store_type, region)

Note: TimesFM 1.0 (used here) does NOT support forecast_with_covariates().
This example uses TimesFM 2.5 which requires a different API. We'll demonstrate
the concept with synthetic data and show the API signature.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Note: TimesFM 1.0 does not support forecast_with_covariates
# This example demonstrates the API with TimesFM 2.5
# Installation: pip install timesfm[xreg]

EXAMPLE_DIR = Path(__file__).parent
OUTPUT_DIR = EXAMPLE_DIR / "output"

# Synthetic data configuration
N_STORES = 3
CONTEXT_LEN = 48  # 48 weeks of history
HORIZON_LEN = 12  # 12 weeks forecast
TOTAL_LEN = CONTEXT_LEN + HORIZON_LEN


def generate_sales_data() -> dict:
    """Generate synthetic retail sales data with covariates."""
    rng = np.random.default_rng(42)

    # Store configurations
    stores = {
        "store_A": {"type": "premium", "region": "urban", "base_sales": 1000},
        "store_B": {"type": "standard", "region": "suburban", "base_sales": 750},
        "store_C": {"type": "discount", "region": "rural", "base_sales": 500},
    }

    data = {"stores": {}, "covariates": {}}

    for store_id, config in stores.items():
        # Base sales with trend
        weeks = np.arange(TOTAL_LEN)
        trend = config["base_sales"] * (1 + 0.005 * weeks)

        # Seasonality (yearly pattern)
        seasonality = 100 * np.sin(2 * np.pi * weeks / 52)

        # Noise
        noise = rng.normal(0, 50, TOTAL_LEN)

        # Price (affects sales negatively)
        price = 10 + rng.uniform(-1, 1, TOTAL_LEN)
        price_effect = -20 * (price - 10)

        # Holidays (boost sales)
        holidays = np.zeros(TOTAL_LEN)
        holiday_weeks = [0, 11, 23, 35, 47, 51]  # Major holidays
        for hw in holiday_weeks:
            if hw < TOTAL_LEN:
                holidays[hw] = 1

        holiday_effect = 200 * holidays

        # Promotion (boost sales)
        promotion = rng.choice([0, 1], TOTAL_LEN, p=[0.8, 0.2])
        promo_effect = 150 * promotion

        # Final sales
        sales = (
            trend + seasonality + noise + price_effect + holiday_effect + promo_effect
        )
        sales = np.maximum(sales, 50)  # Ensure positive

        # Day of week effect (0=Mon, 6=Sun) - simplified to weekly
        day_of_week = np.tile(np.arange(7), TOTAL_LEN // 7 + 1)[:TOTAL_LEN]

        data["stores"][store_id] = {
            "sales": sales.astype(np.float32),
            "config": config,
        }

        # Covariates (same structure for all stores, different values)
        if store_id == "store_A":
            data["covariates"] = {
                "price": {store_id: price.astype(np.float32) for store_id in stores},
                "promotion": {
                    store_id: promotion.astype(np.float32) for store_id in stores
                },
                "holiday": {
                    store_id: holidays.astype(np.float32) for store_id in stores
                },
                "day_of_week": {
                    store_id: day_of_week.astype(np.int32) for store_id in stores
                },
                "store_type": {store_id: config["type"] for store_id in stores},
                "region": {store_id: config["region"] for store_id in stores},
            }

    return data


def demonstrate_api() -> None:
    """Show the forecast_with_covariates API structure."""

    print("\n" + "=" * 70)
    print("  TIMESFM COVARIATES API (TimesFM 2.5)")
    print("=" * 70)

    api_code = """
# Installation
pip install timesfm[xreg]

# Import
import timesfm

# Load TimesFM 2.5 (supports covariates)
hparams = timesfm.TimesFmHparams(
    backend="cpu",  # or "gpu"
    per_core_batch_size=32,
    horizon_len=12,
)
checkpoint = timesfm.TimesFmCheckpoint(
    huggingface_repo_id="google/timesfm-2.5-200m-pytorch"
)
model = timesfm.TimesFm(hparams=hparams, checkpoint=checkpoint)

# Prepare inputs
inputs = [sales_store_a, sales_store_b, sales_store_c]  # List of historical sales

# Dynamic numerical covariates (context + horizon values per series)
dynamic_numerical_covariates = {
    "price": [
        price_history_store_a,  # Shape: (context_len + horizon_len,)
        price_history_store_b,
        price_history_store_c,
    ],
    "promotion": [promo_a, promo_b, promo_c],
}

# Dynamic categorical covariates
dynamic_categorical_covariates = {
    "holiday": [holiday_a, holiday_b, holiday_c],  # 0 or 1 flags
    "day_of_week": [dow_a, dow_b, dow_c],  # 0-6 integer values
}

# Static categorical covariates (one value per series)
static_categorical_covariates = {
    "store_type": ["premium", "standard", "discount"],
    "region": ["urban", "suburban", "rural"],
}

# Forecast with covariates
point_forecast, quantile_forecast = model.forecast_with_covariates(
    inputs=inputs,
    dynamic_numerical_covariates=dynamic_numerical_covariates,
    dynamic_categorical_covariates=dynamic_categorical_covariates,
    static_categorical_covariates=static_categorical_covariates,
    xreg_mode="xreg + timesfm",  # or "timesfm + xreg"
    ridge=0.0,                   # Ridge regularization
    normalize_xreg_target_per_input=True,
)

# Output shapes
# point_forecast: (num_series, horizon_len)
# quantile_forecast: (num_series, horizon_len, 10)
"""
    print(api_code)


def explain_xreg_modes() -> None:
    """Explain the two XReg modes."""

    print("\n" + "=" * 70)
    print("  XREG MODES EXPLAINED")
    print("=" * 70)

    print("""
┌─────────────────────────────────────────────────────────────────────┐
│  Mode 1: "xreg + timesfm" (DEFAULT)                                │
├─────────────────────────────────────────────────────────────────────┤
│  1. TimesFM makes baseline forecast (ignoring covariates)          │
│  2. Calculate residuals: actual - baseline                         │
│  3. Fit linear regression: residuals ~ covariates                  │
│  4. Final forecast = TimesFM baseline + XReg adjustment            │
│                                                                     │
│  Best for: Covariates capture residual patterns                    │
│           (e.g., promotions affecting baseline sales)              │
└─────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────┐
│  Mode 2: "timesfm + xreg"                                          │
├─────────────────────────────────────────────────────────────────────┤
│  1. Fit linear regression: target ~ covariates                     │
│  2. Calculate residuals: actual - regression_prediction            │
│  3. TimesFM forecasts residuals                                     │
│  4. Final forecast = XReg prediction + TimesFM residual forecast   │
│                                                                     │
│  Best for: Covariates explain main signal                          │
│           (e.g., temperature driving ice cream sales)              │
└─────────────────────────────────────────────────────────────────────┘
""")


def create_visualization(data: dict) -> None:
    """Create visualization of sales data with covariates."""

    OUTPUT_DIR.mkdir(exist_ok=True)

    fig, axes = plt.subplots(3, 2, figsize=(16, 12))

    weeks = np.arange(TOTAL_LEN)
    context_weeks = weeks[:CONTEXT_LEN]
    horizon_weeks = weeks[CONTEXT_LEN:]

    # Plot 1: Sales by store
    ax = axes[0, 0]
    for store_id, store_data in data["stores"].items():
        ax.plot(
            context_weeks,
            store_data["sales"][:CONTEXT_LEN],
            label=f"{store_id} ({store_data['config']['type']})",
            linewidth=2,
        )
    ax.axvline(x=CONTEXT_LEN, color="red", linestyle="--", label="Forecast Start")
    ax.set_xlabel("Week")
    ax.set_ylabel("Sales")
    ax.set_title("Historical Sales by Store")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: Price covariate
    ax = axes[0, 1]
    for store_id in data["stores"]:
        ax.plot(weeks, data["covariates"]["price"][store_id], label=store_id, alpha=0.7)
    ax.axvline(x=CONTEXT_LEN, color="red", linestyle="--")
    ax.set_xlabel("Week")
    ax.set_ylabel("Price ($)")
    ax.set_title("Dynamic Numerical Covariate: Price")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 3: Holiday covariate
    ax = axes[1, 0]
    holidays = data["covariates"]["holiday"]["store_A"]
    ax.bar(weeks, holidays, alpha=0.7, color="orange")
    ax.axvline(x=CONTEXT_LEN, color="red", linestyle="--")
    ax.set_xlabel("Week")
    ax.set_ylabel("Holiday Flag")
    ax.set_title("Dynamic Categorical Covariate: Holiday")
    ax.grid(True, alpha=0.3)

    # Plot 4: Promotion covariate
    ax = axes[1, 1]
    promotions = data["covariates"]["promotion"]["store_A"]
    ax.bar(weeks, promotions, alpha=0.7, color="green")
    ax.axvline(x=CONTEXT_LEN, color="red", linestyle="--")
    ax.set_xlabel("Week")
    ax.set_ylabel("Promotion Flag")
    ax.set_title("Dynamic Categorical Covariate: Promotion")
    ax.grid(True, alpha=0.3)

    # Plot 5: Store type (static)
    ax = axes[2, 0]
    store_types = [data["covariates"]["store_type"][s] for s in data["stores"]]
    store_ids = list(data["stores"].keys())
    colors = {"premium": "gold", "standard": "silver", "discount": "brown"}
    ax.bar(store_ids, [1, 1, 1], color=[colors[t] for t in store_types])
    ax.set_ylabel("Store Type")
    ax.set_title("Static Categorical Covariate: Store Type")
    ax.set_yticks([])
    for i, (sid, t) in enumerate(zip(store_ids, store_types)):
        ax.text(i, 0.5, t, ha="center", va="center", fontweight="bold")

    # Plot 6: Data structure summary
    ax = axes[2, 1]
    ax.axis("off")

    summary_text = """
    COVARIATE DATA STRUCTURE
    ─────────────────────────
    
    Dynamic Numerical Covariates:
    • price: np.ndarray[context_len + horizon_len] per series
    • promotion: np.ndarray[context_len + horizon_len] per series
    
    Dynamic Categorical Covariates:
    • holiday: np.ndarray[context_len + horizon_len] per series
    • day_of_week: np.ndarray[context_len + horizon_len] per series
    
    Static Categorical Covariates:
    • store_type: ["premium", "standard", "discount"]
    • region: ["urban", "suburban", "rural"]
    
    Note: Future covariate values must be known!
    (Price, promotion schedule, holidays are planned in advance)
    """
    ax.text(
        0.1,
        0.5,
        summary_text,
        transform=ax.transAxes,
        fontfamily="monospace",
        fontsize=10,
        verticalalignment="center",
    )

    plt.tight_layout()

    output_path = OUTPUT_DIR / "covariates_data.png"
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"\n📊 Saved visualization: {output_path}")
    plt.close()


def main() -> None:
    print("=" * 70)
    print("  TIMESFM COVARIATES (XREG) EXAMPLE")
    print("=" * 70)

    # Generate synthetic data
    print("\n📊 Generating synthetic retail sales data...")
    data = generate_sales_data()

    print(f"   Stores: {list(data['stores'].keys())}")
    print(f"   Context length: {CONTEXT_LEN} weeks")
    print(f"   Horizon length: {HORIZON_LEN} weeks")
    print(f"   Covariates: {list(data['covariates'].keys())}")

    # Show API
    demonstrate_api()

    # Explain modes
    explain_xreg_modes()

    # Create visualization
    print("\n📊 Creating data visualization...")
    create_visualization(data)

    # Save data
    print("\n💾 Saving synthetic data...")

    # Convert to DataFrame for CSV export
    records = []
    for store_id, store_data in data["stores"].items():
        for i, week in enumerate(range(TOTAL_LEN)):
            records.append(
                {
                    "store_id": store_id,
                    "week": week,
                    "sales": store_data["sales"][i],
                    "price": data["covariates"]["price"][store_id][i],
                    "promotion": data["covariates"]["promotion"][store_id][i],
                    "holiday": int(data["covariates"]["holiday"][store_id][i]),
                    "day_of_week": int(data["covariates"]["day_of_week"][store_id][i]),
                    "store_type": data["covariates"]["store_type"][store_id],
                    "region": data["covariates"]["region"][store_id],
                }
            )

    df = pd.DataFrame(records)
    csv_path = OUTPUT_DIR / "sales_with_covariates.csv"
    df.to_csv(csv_path, index=False)
    print(f"   Saved: {csv_path}")

    # Save metadata
    metadata = {
        "description": "Synthetic retail sales data with covariates for TimesFM XReg demo",
        "stores": {sid: sdata["config"] for sid, sdata in data["stores"].items()},
        "dimensions": {
            "context_length": CONTEXT_LEN,
            "horizon_length": HORIZON_LEN,
            "total_length": TOTAL_LEN,
        },
        "covariates": {
            "dynamic_numerical": ["price", "promotion"],
            "dynamic_categorical": ["holiday", "day_of_week"],
            "static_categorical": ["store_type", "region"],
        },
        "xreg_modes": {
            "xreg + timesfm": "Fit regression on residuals after TimesFM forecast",
            "timesfm + xreg": "TimesFM forecasts residuals after regression fit",
        },
    }

    meta_path = OUTPUT_DIR / "covariates_metadata.json"
    with open(meta_path, "w") as f:
        json.dump(metadata, f, indent=2)
    print(f"   Saved: {meta_path}")

    # Summary
    print("\n" + "=" * 70)
    print("  ✅ COVARIATES EXAMPLE COMPLETE")
    print("=" * 70)

    print("""
💡 Key Points:

1. INSTALLATION: Requires timesfm[xreg] extra
   pip install timesfm[xreg]

2. COVARIATE TYPES:
   • Dynamic: Changes over time (price, promotion, holiday)
   • Static: Fixed per series (store type, region)

3. DATA REQUIREMENTS:
   • Dynamic covariates need values for context + horizon
   • Future values must be known (e.g., planned prices, scheduled holidays)

4. XREG MODES:
   • "xreg + timesfm" (default): Regression on residuals
   • "timesfm + xreg": TimesFM on residuals after regression

5. LIMITATIONS:
   • String categorical values work but slower (use int encoding)
   • Requires TimesFM 2.5+ (v1.0 does not support XReg)

📁 Output Files:
   • output/covariates_data.png - Data visualization
   • output/sales_with_covariates.csv - Sample data
   • output/covariates_metadata.json - Metadata
""")


if __name__ == "__main__":
    main()
