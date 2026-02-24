<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Breeding Value Ranking

Input:
- Synthetic genomic estimated breeding values (GEBV) and standard errors for 60 genotypes

Process:
- Calculate 95% confidence intervals for each GEBV
- Rank genotypes by GEBV from highest to lowest
- Generate a top-15 advancement list

Output:
- output/gebv_ranking_with_ci.csv
- output/selection_list_top15.csv

Run:
```bash
python run_breeding_value_ranking.py
```
