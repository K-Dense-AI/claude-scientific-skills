<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Multiple Testing Correction Comparison

## What This Example Does

This example compares **three methods** for multiple testing correction in GWAS:
1. **Bonferroni**: Controls family-wise error rate (FWER)
2. **Benjamini-Hochberg FDR**: Controls false discovery rate (FDR)
3. **Permutation**: Empirical threshold from null distribution

## Why It Matters

GWAS tests thousands of markers. Without correction, you get many **false positives** by chance.

But different correction methods make different tradeoffs:
- **Conservative** methods: Fewer false positives, but less power
- **Liberal** methods: More power, but more false positives

## Running the Example

```bash
cd scientific-skills/qtl-analysis/examples/threshold-correction
python run_thresholds.py
```

## Input → Process → Output

### Input
- 10,000 markers
- 50 true QTLs with weak effects
- 9,950 null markers

### Process
1. **Bonferroni**: p < 0.05/10000 = 5e-6
2. **Benjamini-Hochberg**: Adaptive threshold controlling FDR = 0.05
3. **Permutation**: Empirical threshold from 100 permutations

### Output
- `output/threshold_comparison.png` — Four-panel comparison
- `output/threshold_results.csv` — All p-values and significance calls
- `output/threshold_summary.csv` — Performance metrics

## Comparison Results

| Method | Threshold | Significant | True Positives | False Positives | Power | FDR |
|--------|-----------|-------------|----------------|-----------------|-------|-----|
| Bonferroni | 5.00e-06 | 6 | 0 | 6 | 0.0 | 1.0 |
| BH FDR | 1.72e-05 | 7 | 0 | 7 | 0.0 | 1.0 |
| Permutation | 6.08e-06 | 6 | 0 | 6 | 0.0 | 1.0 |

*Note: In this simulation, true QTL effects were too weak to be detected. In practice with stronger effects, BH FDR typically offers better power while controlling false discoveries.*

## Method Comparison

### Bonferroni
- **Formula**: p < α/m
- **Controls**: Family-wise error rate (FWER)
- **Conservative**: High confidence, but may miss true associations
- **Best for**: When false positives are very costly

### Benjamini-Hochberg FDR
- **Procedure**: Find largest k where p(k) ≤ (k/m) × α
- **Controls**: False discovery rate (expected proportion of false positives)
- **Balanced**: Better power than Bonferroni
- **Best for**: Most GWAS studies

### Permutation
- **Method**: Shuffle phenotypes, find minimum p-value, repeat
- **Controls**: Empirical distribution under null
- **Flexible**: Accounts for correlation structure
- **Best for**: When markers are correlated (e.g., linkage disequilibrium)

## Interpreting the Plots

**Panel 1: P-value Distribution**
- Null markers: Uniform (flat histogram)
- True QTLs: Skewed toward 0 (left peak)

**Panel 2: QQ Plot**
- Null markers follow diagonal line
- Significant markers deviate above line

**Panel 3: Manhattan Plot**
- Red line: Bonferroni threshold
- Green line: BH FDR threshold
- Purple line: Permutation threshold

**Panel 4: Performance**
- Green bars: True positives (QTLs correctly identified)
- Red bars: False positives (null markers incorrectly flagged)

## When to Use Each Method

| Scenario | Recommended Method |
|----------|-------------------|
| Conservative analysis | Bonferroni |
| Standard GWAS | BH FDR |
| Correlated markers | Permutation |
| Exploratory analysis | BH FDR (higher power) |
| Validation study | Bonferroni (stricter) |

## QTLmax Equivalent

Multiple testing correction: https://open.qtlmax.com/guide/
