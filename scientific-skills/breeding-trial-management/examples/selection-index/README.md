<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Selection Index

Input:
- Synthetic trait values for yield, protein, and disease score for 50 genotypes

Process:
- Standardize traits to comparable scale
- Apply economic weights (yield and protein positive, disease score negative)
- Compute weighted selection index and rank genotypes

Output:
- output/selection_index_scores.csv
- output/selection_index_top12.csv
- output/selection_index_weights.csv

Run:
```bash
python run_selection_index.py
```
