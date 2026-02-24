<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Crossing Plan

Input:
- Synthetic parent GEBV values and synthetic pairwise coancestry matrix

Process:
- Enumerate all parent pairs
- Compute cross score = mean parent GEBV minus inbreeding penalty
- Select top non-overlapping crosses to reduce relatedness risk

Output:
- output/coancestry_matrix.csv
- output/all_candidate_crosses.csv
- output/optimal_crossing_plan.csv

Run:
```bash
python run_crossing_plan.py
```
