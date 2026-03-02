<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Pedigree Management

Input:
- Synthetic pedigree records with parent IDs

Process:
- Build lineage structure from parent-child links
- Compute inbreeding coefficients from relatedness paths
- Export both tables and lineage visualization

Output:
- output/pedigree_records.csv
- output/inbreeding_coefficients.csv
- output/pedigree_lineage.png
- output/conclusion.txt

Run:
```bash
python run_pedigree_management.py
```
