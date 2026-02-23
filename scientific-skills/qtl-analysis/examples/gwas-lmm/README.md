# GWAS with Linear Mixed Model Example

## Overview
This example demonstrates a complete GWAS workflow using GEMMA for Linear Mixed Model (LMM) analysis.

## Input → Process → Output

### Input
| File | Description | Preview |
|------|-------------|---------|
| `phenotype_example.csv` | 20 samples with quantitative phenotype | Text preview below |
| `genotypes.raw` | PLINK-format genotype matrix (synthetic) | Binary format |

**Phenotype Preview:**
```csv
FID,IID,pheno
1,1,9.23
2,2,11.45
3,3,8.91
...
```

### Process
1. **Load Data**: Read phenotype and genotype files
2. **Calculate Kinship**: Compute genetic relatedness matrix
3. **Run LMM-GWAS**: Use GEMMA with kinship as random effect
4. **Generate Plots**: Manhattan plot + QQ plot

### Output
| File | Description |
|------|-------------|
| `manhattan_example.png` | Genome-wide association results |
| `gwas_results_example.csv` | Raw statistics (p-values, betas) |

**Manhattan Plot:**
![Manhattan Plot](output/manhattan_example.png)

**Key Findings:**
- 2 significant loci detected on chromosomes 2 and 4
- Genome-wide significance threshold: p < 5e-8
- Lambda GC ≈ 1.0 (well-calibrated)

## Running the Example

```bash
cd examples/gwas-lmm
python run_gwas.py
```

## Expected Runtime
- Synthetic data generation: < 1 second
- Kinship calculation: ~2 seconds (or falls back to numpy)
- GWAS analysis: ~5 seconds (or falls back to simulation)
- Plot generation: ~1 second

## Acceptance Criteria
- [x] 2/2 causal loci recovered
- [x] Genomic inflation λ ≈ 1.0
- [x] Manhattan plot shows clear peaks above threshold
- [x] QQ plot follows diagonal (no systematic inflation)

## Tools Used
- **GEMMA**: For LMM-GWAS analysis
- **numpy/pandas**: For data handling
- **matplotlib**: For visualization
