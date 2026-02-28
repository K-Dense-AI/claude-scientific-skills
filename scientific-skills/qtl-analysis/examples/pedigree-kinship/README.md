<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Pedigree Kinship (NRM)

## Input -> Process -> Output

### Input
- Toy pedigree table (ID, sire, dam)

### Process
1. Build numerator relationship matrix (NRM)
2. Export full kinship matrix
3. Visualize as heatmap

### Output
- `output/pedigree_nrm_matrix.csv`
- `output/pedigree_nrm_heatmap.png`

## Run
```bash
python run_pedigree_kinship.py
```
