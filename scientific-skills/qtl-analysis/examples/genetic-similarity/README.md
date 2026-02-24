<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Genetic Similarity (IBS)

## Input -> Process -> Output

### Input
- Simulated genotype matrix for 50 individuals

### Process
1. Compute pairwise IBS similarity
2. Convert to IBS distance
3. Cluster individuals and draw dendrogram

### Output
- `output/ibs_similarity_matrix.csv`
- `output/ibs_distance_matrix.csv`
- `output/genetic_similarity_dendrogram.png`

## Run
```bash
python run_genetic_similarity.py
```
