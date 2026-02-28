<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# SNP Annotation

## Input -> Process -> Output

### Input
- Simulated variant records (chr, pos, gene)

### Process
1. Assign variant effects (missense/synonymous/intron/intergenic/stop-gained)
2. Map effects to impact levels
3. Export annotated table and effect distribution

### Output
- `output/annotated_variants.csv`
- `output/effect_summary.csv`
- `output/effect_summary_plot.png`

## Run
```bash
python run_snp_annotation.py
```
