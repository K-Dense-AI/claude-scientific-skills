<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Elastic Net Cross-Validation for SNP Selection

## What This Example Does

This example demonstrates **elastic net regularization** for selecting markers associated with a trait in genomic prediction.

## Why It Matters

Elastic net combines L1 (Lasso) and L2 (Ridge) regularization:
- **Lasso** → Sparsity (sets many coefficients to zero)
- **Ridge** → Shrinks large coefficients, handles correlated predictors

In GWAS, elastic net helps identify the most important markers while handling linkage disequilibrium (correlated markers).

## Running the Example

```bash
cd scientific-skills/qtl-analysis/examples/elastic-net-cv
python run_elastic_net.py
```

## Input → Process → Output

### Input
- Synthetic genotype: 300 individuals × 500 markers
- 20 true QTL, heritability 0.5

### Process
- ElasticNetCV with 5-fold CV
- Automatic lambda (alpha) selection from 20 values
- L1 ratio = 0.5 (equal L1/L2)

### Output
- `output/lambda_optimization.png` — Alpha selection and accuracy per fold
- `output/selected_markers.png` — Distribution of selected markers
- `output/cv_results.csv` — Per-fold results
- `output/selected_snps.csv` — List of selected markers with effects

## Key Metrics

| Metric | Description |
|--------|-------------|
| Alpha (λ) | Regularization strength — larger = more shrinkage |
| N Selected | Number of markers with non-zero coefficient |
| Correlation | Prediction accuracy vs true TBV |

## When to Use Elastic Net

- **Dense architecture**: Many small-effect markers
- **Correlated markers**: LD blocks handled better than pure Lasso
- **Feature selection**: When you need interpretable marker lists
