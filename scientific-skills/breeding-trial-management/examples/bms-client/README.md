<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# BMS Client (Mock)

Input:
- Synthetic trial records and a mock BMS write payload

Process:
- Simulate BMS trial read endpoint behavior
- Simulate trial write request/response output

Output:
- output/bms_trials_read.csv
- output/bms_trials_write_result.csv
- output/bms_trial_sites.csv
- output/bms_trial_site_map.png
- output/conclusion.txt

Run:
```bash
python run_bms_client.py
```
