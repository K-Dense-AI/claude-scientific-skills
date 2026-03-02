<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Breedbase API Notes

## Typical Resources

- Accessions: create, read, update accession records
- Trials: list and query trial metadata
- Observations: trait values linked to plots and environments

## Example Request Shapes

- Read accessions: `GET /api/accessions?program=<name>&page=<n>`
- Write accession: `POST /api/accessions` with accession ID, name, and metadata fields

## Mapping to Examples

- `examples/breedbase-client/run_breedbase_client.py` simulates read and write flows
- Output tables mirror common fields used in breeding pipelines
