# Bayesian Genomic Prediction

## What This Example Does

This example demonstrates **Bayesian Genomic Prediction** — a method for predicting breeding values using genome-wide markers with Bayesian statistical models.

## Why It Matters

Traditional genomic prediction (GBLUP) treats all markers equally. Bayesian methods allow markers to have different effect sizes, which can improve prediction accuracy when:
- Only a subset of markers contribute to trait variation
- Marker effects follow specific distributions (e.g., many small effects, few large effects)

## Methods Compared

| Method | Key Feature | Best When |
|--------|-------------|-----------|
| **BayesA** | Marker-specific variance (t prior) | Many QTL with varying effects |
| **BayesB** | Spike-slab (some markers = 0) | Sparse genetic architecture |
| **BayesCpi** | Unknown mixture proportion | Unclear if sparse or dense |
| **GBLUP** | Equal variance for all markers | Many small effects (BLUP-like) |

## Running the Example

```bash
cd scientific-skills/qtl-analysis/examples/bayesian-gp
python run_bayesian.py
```

## Input → Process → Output

### Input
- Synthetic genotype data: 300 individuals × 500 markers
- 30 true QTL with random effects
- Heritability: 0.5

### Process
- 5-fold cross-validation for each method
- Compare prediction accuracy (correlation, RMSE, R²)

### Output
- `output/bayesian_comparison.png` — Bar chart comparing methods
- `output/method_summary.csv` — Accuracy metrics by method
- `output/fold_results.csv` — Per-fold results

## QTLmax Equivalent

See QTLmax guide: https://open.qtlmax.com/guide/index.php/2025/07/12/bayesian-gp/

## Notes

This implementation uses sklearn proxies:
- Ridge (α=0.01) → BayesA
- LassoCV → BayesB  
- ElasticNetCV → BayesCpi

For production use, consider BGLR (R package) for true Bayesian implementations.
