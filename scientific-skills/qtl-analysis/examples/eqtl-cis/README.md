<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# cis-eQTL Analysis Example

## Overview
This example demonstrates cis-eQTL mapping using tensorQTL for gene expression analysis.

## Input → Process → Output

### Input
| File | Description | Preview |
|------|-------------|---------|
| `expression_example.csv` | Gene expression values for 20 samples | Text preview below |
| `genotypes.bed` | PLINK genotype data (synthetic) | Binary format |

**Expression Preview:**
```csv
gene_id,sample,expression
ENSG000001234,Sample1,5.23
ENSG000001234,Sample2,4.89
...
```

### Process
1. **Load Expression**: Read phenotype BED format
2. **Load Genotypes**: Read PLINK binary format
3. **Map cis-eQTLs**: Test variants within 1Mb of gene TSS
4. **LocusZoom Plot**: Visualize significant associations

### Output
| File | Description |
|------|-------------|
| `eqtl_locuszoom_example.png` | Local association plot around gene |
| `tensorqtl_results.csv` | eQTL statistics (if tensorQTL available) |

**LocusZoom Plot:**
![LocusZoom](output/eqtl_locuszoom_example.png)

**Key Findings:**
- Lead SNP at position 1.0 Mb with p = 2e-15
- Cis-window: ±500kb from gene TSS
- Clear peak showing regulatory variant

## Running the Example

```bash
cd examples/eqtl-cis
python run_eqtl.py
```

## Expected Runtime
- Synthetic data generation: < 1 second
- tensorQTL mapping: ~10 seconds (GPU) or ~30 seconds (CPU)
- Plot generation: ~1 second

## Acceptance Criteria
- [x] Significant eQTL identified within cis-window
- [x] Lead SNP annotated with p-value
- [x] LocusZoom plot shows clear association peak

## Tools Used
- **tensorQTL**: Fast GPU-accelerated eQTL mapping
- **pandas**: Data manipulation
- **matplotlib**: Visualization
