<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# IoT Field Sync (Breedbase/BMS Pattern)

Input:
- Synthetic field sensor stream (soil moisture, canopy temperature, NDVI)
- Synthetic plot identifiers over daily snapshots

Process:
- Aggregate incoming IoT records into sync-ready daily payloads
- Alternate dispatch targets (Breedbase/BMS) to demonstrate routing
- Generate professional field status map for reviewer and grower reporting

Output:
- output/iot_sensor_stream.csv
- output/iot_sync_log.csv
- output/iot_field_status_map.png
- output/conclusion.txt

Run:
```bash
python run_iot_field_sync.py
```
