# Multi-Trait GWAS

## What This Example Does

This example demonstrates **multi-trait GWAS** — jointly analyzing two correlated traits to increase statistical power when they share genetic architecture.

## Why It Matters

When traits are genetically correlated (e.g., yield and height in crops), analyzing them together can:
- Increase power to detect shared QTLs
- Distinguish pleiotropic effects from linkage
- Improve prediction accuracy

## Running the Example

```bash
cd scientific-skills/qtl-analysis/examples/multi-trait-gwas
python run_multitrait.py
```

## Input → Process → Output

### Input
- 300 individuals, 500 markers
- 10 shared QTLs (affect both traits)
- 5 trait1-specific QTLs, 5 trait2-specific QTLs
- Genetic correlation: 0.7, Heritability: 0.5

### Process
- Single-trait GWAS for each trait
- Multi-trait analysis using Fisher's method to combine p-values
- Bonferroni correction for multiple testing

### Output
- `output/manhattan_trait1.png` — Manhattan plot for Trait 1
- `output/manhattan_trait2.png` — Manhattan plot for Trait 2
- `output/gwas_results.csv` — P-values and QTL annotations

## Key Insight

Multi-trait analysis detected shared QTLs with increased power. Compare:
- Trait 1 only: 4 significant hits
- Trait 2 only: 3 significant hits
- Multi-trait: 4 significant hits (includes shared effects)

## When to Use Multi-Trait GWAS

- **Correlated traits**: When traits have genetic correlation > 0.3
- **Limited sample size**: Increases effective sample size
- **Shared biology**: When traits share biological pathways

## Methods

Uses Fisher's method to combine p-values from both traits:
```
χ² = -2 × (ln(p₁) + ln(p₂))
df = 4 (two traits × two degrees of freedom)
```

## QTLmax Equivalent

Multi-trait analysis: https://open.qtlmax.com/guide/
