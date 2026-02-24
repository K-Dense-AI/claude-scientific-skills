<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Breeding Management System (BMS) API Notes

## Typical Resources

- Trials and studies
- Germplasm lists
- Observation units and measured traits

## Example Request Shapes

- Read trials: `GET /bmsapi/trials?program=<name>`
- Write trial: `POST /bmsapi/trials` with trial name, season, and location

## Mapping to Examples

- `examples/bms-client/run_bms_client.py` simulates read and write operations
- Output files represent common trial tables consumed by downstream analytics
