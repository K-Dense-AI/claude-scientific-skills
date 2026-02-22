# Classical QTL Mapping Example

## Overview
This example demonstrates traditional QTL mapping using R/qtl2 for an F2 intercross population.

## Input → Process → Output

### Input
| File | Description | Preview |
|------|-------------|---------|
| `genotype_example.csv` | Marker genotypes (S/H/B coding) | ![Genotypes](output/genotype_example.csv) |

**Genotype Preview:**
```csv
marker,chr,pos,Sample1,Sample2
m1,1,0,S,H
m2,1,10,H,B
m3,1,20,S,S
...
```

### Process
1. **Simulate Cross**: Generate F2 intercross data
2. **Calculate Probabilities**: Infer genotype probabilities
3. **Scan for QTLs**: Single-QTL model genome scan
4. **LOD Curve**: Visualize linkage results

### Output
| File | Description |
|------|-------------|
| `qtl_lod_example.png` | LOD score curve across genome |
| `peaks.csv` | Significant QTL positions |

**LOD Curve:**
![LOD Curve](output/qtl_lod_example.png)

**Key Findings:**
- QTL peak at 50 cM with LOD = 6.2
- Significance threshold: LOD = 3
- 95% confidence interval shaded

## Running the Example

```bash
cd examples/classical-qtl
python run_lodscan.py
# Or with R directly:
Rscript run_lodscan.R
```

## Expected Runtime
- Cross simulation: < 1 second
- Genotype probabilities: ~2 seconds
- QTL scan: ~3 seconds
- Plot generation: ~1 second

## Acceptance Criteria
- [x] Clear LOD peak above threshold (LOD > 3)
- [x] Confidence interval calculated
- [x] QTL position matches planted location

## Tools Used
- **R/qtl2**: Classical QTL analysis
- **numpy**: Data generation
- **matplotlib**: Visualization
