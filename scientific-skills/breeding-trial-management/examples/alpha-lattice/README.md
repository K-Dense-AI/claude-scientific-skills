<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Alpha-Lattice Design

Input:
- Synthetic genotype list configured for incomplete blocks inside replicates

Process:
- Assign entries to alpha-lattice style blocks
- Export design table and simple efficiency summary

Output:
- output/alpha_lattice_layout.csv
- output/alpha_lattice_efficiency.csv
- output/alpha_lattice_field_map.png
- output/conclusion.txt

Run:
```bash
python run_alpha_lattice.py
```
