<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# CNV Integration for QTL Workflow

Input:
- Synthetic CNV call set (sample, chromosome, interval, copy-number state)
- Synthetic quantitative phenotype values

Process:
- Aggregate CNV burden per sample and chromosome
- Run chromosome-level CNV-trait association proxy
- Export burden and association visuals for downstream reporting

Output:
- output/cnv_calls.csv
- output/cnv_burden_by_chr.csv
- output/cnv_trait_association.csv
- output/cnv_burden_plot.png
- output/cnv_association_plot.png

Run:
```bash
python run_cnv_integration.py
```
