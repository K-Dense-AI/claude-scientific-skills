<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Simple Breeding Cycle Simulation

Input:
- Synthetic founder breeding values for one breeding population
- Configured number of recurrent selection cycles

Process:
- Simulate selection cycles with progressive gain shift
- Track top-selected mean breeding value and diversity per cycle
- Export trend chart and conclusion for planning decisions

Output:
- output/simulation_cycle_summary.csv
- output/simulation_gain_diversity.png
- output/conclusion.txt

Run:
```bash
python run_simple_cycle.py
```
