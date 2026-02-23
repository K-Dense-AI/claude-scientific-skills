# Population Structure Analysis Example

## Overview
This example demonstrates population structure analysis using PCA and kinship estimation.

## Input → Process → Output

### Input
| File | Description | Preview |
|------|-------------|---------|
| `population_info_example.csv` | Sample population assignments | Text preview below |
| `genotypes.raw` | Genotype matrix (synthetic) | Binary format |

**Population Info Preview:**
```csv
sample,population,PC1,PC2
Sample1,PopA,0.23,-0.15
Sample2,PopA,0.45,0.21
Sample3,PopB,1.89,0.03
...
```

### Process
1. **Generate Data**: Simulate admixed population
2. **Run PCA**: Principal component analysis
3. **Calculate Kinship**: Genome-wide relatedness matrix
4. **Visualize**: PCA scatter + kinship heatmap

### Output
| File | Description |
|------|-------------|
| `pca_kinship_example.png` | PCA plot and kinship matrix |
| `admixture_proportions.csv` | Population ancestry fractions |

**Population Structure:**
![PCA and Kinship](output/pca_kinship_example.png)

**Key Findings:**
- 3 distinct populations visible in PCA
- PC1 explains 45% of variance
- Kinship matrix shows population blocks
- No cryptic relatedness detected

## Running the Example

```bash
cd examples/population-structure
python run_structure.py
```

## Expected Runtime
- Population simulation: < 1 second
- PCA calculation: ~2 seconds
- Kinship estimation: ~3 seconds
- Plot generation: ~1 second

## Acceptance Criteria
- [x] 3 populations clearly separated in PCA
- [x] Kinship matrix shows expected structure
- [x] PC1 explains >40% of variance

## Tools Used
- **scikit-learn**: PCA implementation
- **numpy**: Matrix operations
- **matplotlib**: Visualization
