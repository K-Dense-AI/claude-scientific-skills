<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# RCBD Design

Input:
- Synthetic genotype list and block count for a randomized complete block design (RCBD)

Process:
- Randomize genotype order within each block
- Build plot layout table
- Export both table and quick layout image

Output:
- output/rcbd_layout.csv
- output/rcbd_layout.png

Run:
```bash
python run_rcbd.py
```
