<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Breedbase Client (Mock)

Input:
- Synthetic accession records and a mock write payload

Process:
- Simulate Breedbase read endpoint behavior
- Simulate accession write request/response output

Output:
- output/breedbase_accessions_read.csv
- output/breedbase_accessions_write_result.csv
- output/breedbase_accession_sites.csv
- output/breedbase_accession_site_map.png
- output/conclusion.txt

Run:
```bash
python run_breedbase_client.py
```
