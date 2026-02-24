# Genomic Control

## What This Example Does
Runs a GWAS inflation check, computes genomic inflation factor lambda, applies genomic control, and compares QQ plots before and after adjustment.

## Input -> Process -> Output

### Input
- Simulated GWAS chi-square statistics for 5,000 markers
- Intentional inflation (population structure effect)

### Process
1. Convert chi-square statistics to p-values
2. Compute lambda: median(chi-square) / 0.456
3. Adjust statistics by dividing by lambda
4. Recompute p-values and generate QQ comparison

### Output
- `output/qq_before_after_gc.png`
- `output/genomic_control_results.csv`
- `output/lambda_summary.csv`

## Run
```bash
python run_genomic_control.py
```

## Interpretation
- Lambda > 1 indicates inflation from confounding.
- After genomic control, lambda should move closer to 1.
- Use this as a diagnostic step before claiming GWAS hits.
