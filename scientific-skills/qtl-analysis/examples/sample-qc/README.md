<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Sample QC

## Input -> Process -> Output

### Input
- Simulated sample-by-marker genotype matrix with missingness

### Process
1. Compute per-sample call rate
2. Compute heterozygosity
3. Flag call-rate failures and heterozygosity outliers
4. Generate QC scatter and sex proxy summary

### Output
- `output/sample_qc_metrics.csv`
- `output/sample_qc_scatter.png`
- `output/sex_check_summary.csv`

## Run
```bash
python run_sample_qc.py
```
