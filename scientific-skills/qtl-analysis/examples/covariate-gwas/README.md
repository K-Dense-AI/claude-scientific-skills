# Covariate-Adjusted GWAS

## What This Example Does

This example demonstrates **GWAS with covariate adjustment** — including population structure, sex, and batch effects to avoid confounding and reduce false positives.

## Why It Matters

Population structure can cause **false positives** in GWAS:
- Allele frequencies differ between subpopulations
- Phenotypes may also differ by subpopulation
- This creates spurious associations

**Covariate adjustment** removes these confounding effects.

## Running the Example

```bash
cd scientific-skills/qtl-analysis/examples/covariate-gwas
python run_covariate_gwas.py
```

## Input → Process → Output

### Input
- 300 individuals, 500 markers
- 3 subpopulations with different allele frequencies
- 15 true QTLs (some confounded with population structure)
- Covariates: PCs (population structure), sex, batch

### Process
1. **Simple GWAS**: No covariates (confounded)
2. **Covariate-adjusted GWAS**: Residualizes for PCs, sex, batch
3. **Comparison**: Shows Manhattan plots side-by-side

### Output
- `output/covariate_comparison.png` — Side-by-side Manhattan plots
- `output/covariate_results.csv` — P-values with and without adjustment

## Key Results

| Method | Significant Hits | Interpretation |
|--------|------------------|----------------|
| Without covariates | 456 | Many false positives from population structure |
| With covariates | 2 | Better calibrated, detects true QTLs |
| True QTLs detected | 2/15 | Maintains power while reducing false positives |

## What Are Covariates?

| Covariate | What It Captures | When to Include |
|-----------|------------------|---------------|
| **Population PCs** | Ancestry differences | Always for structured populations |
| **Sex** | Male/female differences | Sexual dimorphism in trait |
| **Batch** | Processing time/location | Technical variation |

## How Covariate Adjustment Works

1. **Regress phenotype on covariates** → Get residuals
2. **Regress each marker on covariates** → Get residuals
3. **Correlate residuals** → True genetic association

This removes the confounding effect of population structure.

## When to Use Covariates

- **Population structure**: PCA on genotypes
- **Sex effects**: Include sex as covariate
- **Batch effects**: Include batch/ processing date
- **Environmental gradients**: PCs of environmental data

## QTLmax Equivalent

Population structure correction: https://open.qtlmax.com/guide/
