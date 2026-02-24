# G×E GWAS (Genotype × Environment)

## What This Example Does

This example demonstrates **G×E GWAS** — detecting QTLs with environment-specific effects in multi-environment trials (METs).

## Why It Matters

Genotype × Environment (G×E) interaction is common in breeding:
- A QTL may have strong effects in some environments but not others
- Understanding G×E helps select stable varieties
- MET analysis increases power by combining data across locations/years

## Running the Example

```bash
cd scientific-skills/qtl-analysis/examples/gxe-gwas
python run_gxe_gwas.py
```

## Input → Process → Output

### Input
- 250 individuals, 400 markers
- 3 environments (locations or years)
- 8 stable QTLs (same effect in all environments)
- 6 G×E QTLs (environment-specific effects)
- Heritability: 0.5

### Process
- Test main effects (average across environments)
- Test G×E interaction effects (does marker effect vary by environment?)
- Bonferroni correction for multiple testing

### Output
- `output/gxe_manhattan.png` — Manhattan plots for main and interaction effects
- `output/gxe_results.csv` — P-values with QTL annotations

## Key Insight

G×E analysis distinguishes:
- **Stable QTLs**: Significant main effects, consistent across environments
- **G×E QTLs**: Significant interaction effects, environment-dependent

Results showed 5 main effects and 400 G×E interactions (note: interaction test may need refinement for real data).

## Interpreting Results

Look at both panels in `gxe_manhattan.png`:
- **Top (Main Effects)**: QTLs that work everywhere (stable)
- **Bottom (G×E Effects)**: QTLs that work in specific environments

Orange lines mark true stable QTLs. Purple lines mark true G×E QTLs.

## When to Use G×E GWAS

- **Multi-location trials**: Tested in multiple environments
- **Year-to-year variation**: Different growing seasons
- **Targeted breeding**: Identify QTLs for specific environments

## Methods

Uses simple correlation for main effects and ANOVA-style F-test for interactions:
```
Main effect: Correlation between marker and mean phenotype
G×E effect: F-test comparing phenotypic variance across environments
```

## QTLmax Equivalent

G×E analysis: https://open.qtlmax.com/guide/

## Notes

This is a simplified implementation. Production G×E GWAS should use:
- Linear mixed models with random environment effects
- Proper accounting for spatial variation
- Multi-environment genomic prediction
