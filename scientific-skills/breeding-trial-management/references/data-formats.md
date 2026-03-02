<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Data Format Conventions

## Core Tables

- Fieldbook table: `plot_id`, `block`, `row`, `col`, `genotype`
- Phenotype table: `plot_id`, `genotype`, `environment`, `replicate`, trait columns
- Germplasm table: `accession_id`, `accession_name`, source metadata
- Crossing table: `parent_1`, `parent_2`, expected score columns

## Naming Rules

- Use lowercase snake_case for standardized outputs
- Use unit suffixes where needed (`yield_kg_ha`, `moisture_pct`)
- Keep stable IDs (`plot_id`, `genotype`) for table joins

## Validation Rules

- No empty IDs
- Replicates should be positive integers
- Trait columns should be numeric after import normalization
