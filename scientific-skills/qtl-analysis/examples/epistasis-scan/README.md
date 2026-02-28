<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Genome-wide Epistasis Scan

Input:
- Synthetic genotype matrix for 300 samples and 18 SNPs
- Synthetic trait with a planted pairwise interaction signal

Process:
- Compute pairwise SNP interaction terms across all SNP pairs
- Score each interaction against phenotype signal
- Rank candidate epistatic pairs and visualize interaction landscape

Output:
- output/epistasis_pair_scan.csv
- output/epistasis_top_pairs.csv
- output/epistasis_top_pairs_plot.png
- output/epistasis_signal_matrix.png

Run:
```bash
python run_epistasis_scan.py
```
