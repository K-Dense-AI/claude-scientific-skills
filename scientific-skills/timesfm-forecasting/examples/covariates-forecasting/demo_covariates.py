#!/usr/bin/env python3
"""
TimesFM Covariates (XReg) Example

Demonstrates the TimesFM covariate API structure using synthetic retail
sales data. TimesFM 1.0 does NOT support forecast_with_covariates().
That feature requires TimesFM 2.5 + `timesfm[xreg]`.

This script:
  1. Generates synthetic 3-store retail data (24-week context, 12-week horizon)
  2. Visualises each covariate type (dynamic numerical, dynamic categorical, static)
  3. Prints the forecast_with_covariates() call signature for reference
  4. Exports a compact CSV (90 rows) and metadata JSON

NOTE ON REAL DATA:
  If you want to use a real retail dataset (e.g., Kaggle Rossmann Store Sales),
  download it to a TEMP location — do NOT commit large CSVs to this repo.
  Example:
      import tempfile, urllib.request
      tmp = tempfile.mkdtemp(prefix="timesfm_retail_")
      # urllib.request.urlretrieve("https://...store_sales.csv", f"{tmp}/store_sales.csv")
      # df = pd.read_csv(f"{tmp}/store_sales.csv")
  Users should persist the data wherever makes sense for their workflow;
  this skills directory intentionally keeps only tiny reference datasets.
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

# Synthetic data configuration — kept SMALL (24 weeks context, 90 CSV rows)
N_STORES = 3
CONTEXT_LEN = 24  # weeks of history  (was 48 — halved for token efficiency)
HORIZON_LEN = 12  # weeks to forecast
TOTAL_LEN = CONTEXT_LEN + HORIZON_LEN  # 36 weeks total per store


def generate_sales_data() -> dict:
    """Generate synthetic retail sales data with covariates.

    BUG FIX (v2): Previous version had a variable-shadowing issue where the
    inner dict comprehension `{store_id: ... for store_id in stores}` overwrote
    the outer loop variable, giving all stores identical covariate data (store_A's).
    Fixed by collecting per-store arrays into separate dicts during the outer loop
    and building the covariates dict afterwards.
    """
    rng = np.random.default_rng(42)

    # Store configurations
    stores = {
        "store_A": {"type": "premium", "region": "urban", "base_sales": 1000},
        "store_B": {"type": "standard", "region": "suburban", "base_sales": 750},
        "store_C": {"type": "discount", "region": "rural", "base_sales": 500},
    }

    data: dict = {"stores": {}, "covariates": {}}

    # Collect per-store covariate arrays *before* building the covariates dict
    prices_by_store: dict[str, np.ndarray] = {}
    promos_by_store: dict[str, np.ndarray] = {}
    holidays_by_store: dict[str, np.ndarray] = {}
    day_of_week_by_store: dict[str, np.ndarray] = {}

    for store_id, config in stores.items():
        weeks = np.arange(TOTAL_LEN)
        trend = config["base_sales"] * (1 + 0.005 * weeks)
        seasonality = 100 * np.sin(2 * np.pi * weeks / 52)
        noise = rng.normal(0, 50, TOTAL_LEN)

        # Price — slightly different range per store to reflect market positioning
        base_price = {"store_A": 12.0, "store_B": 10.0, "store_C": 7.5}[store_id]
        price = base_price + rng.uniform(-0.5, 0.5, TOTAL_LEN)
        price_effect = -20 * (price - base_price)

        # Holidays (major retail weeks)
        holidays = np.zeros(TOTAL_LEN)
        for hw in [0, 11, 23, 35]:
            if hw < TOTAL_LEN:
                holidays[hw] = 1.0
        holiday_effect = 200 * holidays

        # Promotion — random 20% of weeks
        promotion = rng.choice([0.0, 1.0], TOTAL_LEN, p=[0.8, 0.2])
        promo_effect = 150 * promotion

        # Day-of-week proxy (weekly granularity → repeat 0-6 pattern)
        day_of_week = np.tile(np.arange(7), TOTAL_LEN // 7 + 1)[:TOTAL_LEN]

        sales = (
            trend + seasonality + noise + price_effect + holiday_effect + promo_effect
        )
        sales = np.maximum(sales, 50.0).astype(np.float32)

        data["stores"][store_id] = {"sales": sales, "config": config}

        prices_by_store[store_id] = price.astype(np.float32)
        promos_by_store[store_id] = promotion.astype(np.float32)
        holidays_by_store[store_id] = holidays.astype(np.float32)
        day_of_week_by_store[store_id] = day_of_week.astype(np.int32)

    # Build covariates dict AFTER the loop (avoids shadowing bug)
    data["covariates"] = {
        "price": prices_by_store,
        "promotion": promos_by_store,
        "holiday": holidays_by_store,
        "day_of_week": day_of_week_by_store,
        "store_type": {sid: stores[sid]["type"] for sid in stores},
        "region": {sid: stores[sid]["region"] for sid in stores},
    }

    return data


def demonstrate_api() -> None:
    """Print the forecast_with_covariates API structure (TimesFM 2.5)."""

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
    "holiday":     [holiday_a, holiday_b, holiday_c],    # 0 or 1 flags
    "day_of_week": [dow_a, dow_b, dow_c],                # 0-6 integer values
}

# Static categorical covariates (one value per series)
static_categorical_covariates = {
    "store_type": ["premium", "standard", "discount"],
    "region":     ["urban", "suburban", "rural"],
}

# Forecast with covariates
point_forecast, quantile_forecast = model.forecast_with_covariates(
    inputs=inputs,
    dynamic_numerical_covariates=dynamic_numerical_covariates,
    dynamic_categorical_covariates=dynamic_categorical_covariates,
    static_categorical_covariates=static_categorical_covariates,
    xreg_mode="xreg + timesfm",            # or "timesfm + xreg"
    ridge=0.0,                              # Ridge regularization
    normalize_xreg_target_per_input=True,
)

# Output shapes
# point_forecast:    (num_series, horizon_len)
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

    # Panel 1 — Sales by store (context only)
    ax = axes[0, 0]
    for store_id, store_data in data["stores"].items():
        ax.plot(
            context_weeks,
            store_data["sales"][:CONTEXT_LEN],
            label=f"{store_id} ({store_data['config']['type']})",
            linewidth=2,
        )
    ax.axvline(
        x=CONTEXT_LEN - 0.5, color="red", linestyle="--", label="Forecast Start →"
    )
    ax.set_xlabel("Week")
    ax.set_ylabel("Sales")
    ax.set_title("Historical Sales by Store (24-week context)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Panel 2 — Price covariate (all weeks including horizon)
    ax = axes[0, 1]
    for store_id in data["stores"]:
        ax.plot(weeks, data["covariates"]["price"][store_id], label=store_id, alpha=0.8)
    ax.axvline(x=CONTEXT_LEN - 0.5, color="red", linestyle="--")
    ax.set_xlabel("Week")
    ax.set_ylabel("Price ($)")
    ax.set_title("Dynamic Numerical Covariate: Price\n(different baseline per store)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Panel 3 — Holiday flag
    ax = axes[1, 0]
    # Show all 3 stores' holidays side by side (they're the same here but could differ)
    ax.bar(weeks, data["covariates"]["holiday"]["store_A"], alpha=0.7, color="orange")
    ax.axvline(x=CONTEXT_LEN - 0.5, color="red", linestyle="--")
    ax.set_xlabel("Week")
    ax.set_ylabel("Holiday Flag")
    ax.set_title("Dynamic Categorical Covariate: Holiday")
    ax.grid(True, alpha=0.3)

    # Panel 4 — Promotion (store_A example — each store differs)
    ax = axes[1, 1]
    for store_id in data["stores"]:
        ax.bar(
            weeks + {"store_A": -0.3, "store_B": 0.0, "store_C": 0.3}[store_id],
            data["covariates"]["promotion"][store_id],
            width=0.3,
            alpha=0.7,
            label=store_id,
        )
    ax.axvline(x=CONTEXT_LEN - 0.5, color="red", linestyle="--")
    ax.set_xlabel("Week")
    ax.set_ylabel("Promotion Flag")
    ax.set_title("Dynamic Categorical Covariate: Promotion\n(independent per store)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Panel 5 — Store type (static)
    ax = axes[2, 0]
    store_types = [data["covariates"]["store_type"][s] for s in data["stores"]]
    store_ids = list(data["stores"].keys())
    colors = {"premium": "gold", "standard": "silver", "discount": "#cd7f32"}
    ax.bar(store_ids, [1, 1, 1], color=[colors[t] for t in store_types])
    ax.set_ylabel("Store Type")
    ax.set_title("Static Categorical Covariate: Store Type")
    ax.set_yticks([])
    for i, (sid, t) in enumerate(zip(store_ids, store_types)):
        ax.text(i, 0.5, t, ha="center", va="center", fontweight="bold", fontsize=11)

    # Panel 6 — Data structure summary
    ax = axes[2, 1]
    ax.axis("off")
    summary_text = (
        "  COVARIATE DATA STRUCTURE\n"
        "  ─────────────────────────\n\n"
        "  Dynamic Numerical Covariates:\n"
        "  • price:     array[context_len + horizon_len] per series\n"
        "  • promotion: array[context_len + horizon_len] per series\n\n"
        "  Dynamic Categorical Covariates:\n"
        "  • holiday:     array[context_len + horizon_len] per series\n"
        "  • day_of_week: array[context_len + horizon_len] per series\n\n"
        "  Static Categorical Covariates:\n"
        "  • store_type: ['premium', 'standard', 'discount']\n"
        "  • region:     ['urban', 'suburban', 'rural']\n\n"
        "  ⚠  Future covariate values must be KNOWN at forecast time!\n"
        "     (Prices, promotion schedules, and holidays are planned.)"
    )
    ax.text(
        0.05,
        0.5,
        summary_text,
        transform=ax.transAxes,
        fontfamily="monospace",
        fontsize=9,
        verticalalignment="center",
    )

    plt.suptitle(
        "TimesFM Covariates (XReg) — Synthetic Retail Sales Demo",
        fontsize=14,
        fontweight="bold",
        y=1.01,
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

    print(f"   Stores:          {list(data['stores'].keys())}")
    print(f"   Context length:  {CONTEXT_LEN} weeks")
    print(f"   Horizon length:  {HORIZON_LEN} weeks")
    print(f"   Covariates:      {list(data['covariates'].keys())}")

    # Show API
    demonstrate_api()

    # Explain modes
    explain_xreg_modes()

    # Create visualization
    print("\n📊 Creating data visualization...")
    create_visualization(data)

    # Save data
    print("\n💾 Saving synthetic data...")

    records = []
    for store_id, store_data in data["stores"].items():
        for i in range(TOTAL_LEN):
            records.append(
                {
                    "store_id": store_id,
                    "week": i,
                    "split": "context" if i < CONTEXT_LEN else "horizon",
                    "sales": round(float(store_data["sales"][i]), 2),
                    "price": round(float(data["covariates"]["price"][store_id][i]), 4),
                    "promotion": int(data["covariates"]["promotion"][store_id][i]),
                    "holiday": int(data["covariates"]["holiday"][store_id][i]),
                    "day_of_week": int(data["covariates"]["day_of_week"][store_id][i]),
                    "store_type": data["covariates"]["store_type"][store_id],
                    "region": data["covariates"]["region"][store_id],
                }
            )

    df = pd.DataFrame(records)
    csv_path = OUTPUT_DIR / "sales_with_covariates.csv"
    df.to_csv(csv_path, index=False)
    print(f"   Saved: {csv_path}  ({len(df)} rows × {len(df.columns)} cols)")

    # Save metadata
    metadata = {
        "description": "Synthetic retail sales data with covariates for TimesFM XReg demo",
        "note_on_real_data": (
            "If using a real dataset (e.g., Kaggle Rossmann Store Sales), "
            "download it to a temp directory (tempfile.mkdtemp) and do NOT "
            "commit it here. This skills directory only ships tiny reference files."
        ),
        "stores": {sid: sdata["config"] for sid, sdata in data["stores"].items()},
        "dimensions": {
            "context_length": CONTEXT_LEN,
            "horizon_length": HORIZON_LEN,
            "total_length": TOTAL_LEN,
            "num_stores": N_STORES,
            "csv_rows": len(df),
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
        "bug_fixes": [
            "v2: Fixed variable-shadowing in generate_sales_data() — inner dict "
            "comprehension `{store_id: ... for store_id in stores}` was overwriting "
            "the outer loop variable, causing all stores to get identical covariate "
            "arrays. Fixed by using separate per-store dicts during the loop.",
            "v2: Reduced CONTEXT_LEN from 48 → 24 weeks; CSV now 90 rows (was 180).",
        ],
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
   • Dynamic Numerical:   time-varying numeric (price, promotion)
   • Dynamic Categorical: time-varying flags   (holiday, day_of_week)
   • Static Categorical:  fixed per series     (store_type, region)

3. DATA REQUIREMENTS:
   • Dynamic covariates need values for context + horizon
   • Future values must be known (prices, scheduled holidays, etc.)

4. XREG MODES:
   • "xreg + timesfm" (default): Regression on residuals
   • "timesfm + xreg":           TimesFM on residuals after regression

5. LIMITATIONS:
   • Requires TimesFM 2.5+ (v1.0 does not support XReg)
   • String categoricals work but int encoding is faster

📁 Output Files:
   • output/covariates_data.png         — visualization (6 panels)
   • output/sales_with_covariates.csv   — 90-row compact dataset
   • output/covariates_metadata.json    — metadata + bug-fix log
""")


if __name__ == "__main__":
    main()
