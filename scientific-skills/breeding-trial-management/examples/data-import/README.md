<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Data Import Pipeline

Input:
- Synthetic raw phenotype table that mimics field-device CSV naming patterns

Process:
- Normalize headers to a standard schema
- Impute missing yield values with median
- Validate row counts, replicate values, and key field coverage

Output:
- output/raw_phenotypes.csv
- output/standardized_phenotypes.csv
- output/validation_report.txt

Run:
```bash
python run_data_import.py
```
