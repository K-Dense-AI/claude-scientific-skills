<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Augmented Design

Input:
- Synthetic unreplicated test entries plus replicated check entries

Process:
- Assign checks across blocks
- Place unreplicated entries around checks
- Export layout and summary statistics

Output:
- output/augmented_layout.csv
- output/augmented_summary.csv
- output/augmented_field_map.png
- output/conclusion.txt

Run:
```bash
python run_augmented.py
```
