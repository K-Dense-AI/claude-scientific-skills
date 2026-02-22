# QTLmax vs Open-Source Tools Comparison

## Feature Matrix

| QTLmax Feature | Open-Source Equivalent | Language | Notes |
|---------------|------------------------|----------|-------|
| **LMM-GWAS** | GEMMA | C++ | Same algorithm, identical results |
| **GLM-GWAS** | PLINK 2 | C++ | Industry standard |
| **eQTL mapping** | tensorQTL | Python/GPU | 100x faster than QTLtools |
| **Manhattan/QQ plots** | pyQTL + matplotlib | Python | Fully customizable |
| **Kinship matrices** | GEMMA -gk or numpy | C++/Python | VanRaden method |
| **PCA** | PLINK --pca or scikit-learn | C++/Python | Fast, standard |
| **Admixture** | ADMIXTURE or sklearn | C++/Python | K-means alternative |
| **Classical QTL** | R/qtl2 | R | Only option for LOD scans |
| **Interval mapping** | R/qtl2 | R | Haley-Knott, EM, imp |
| **CIM/MQM** | R/qtl2 | R | Multiple QTL models |
| **LD analysis** | PLINK --ld | C++ | Pairwise r² |
| **Imputation** | Beagle | Java | Standard pipeline |
| **Genomic prediction** | rrBLUP or custom | R/Python | GBLUP implementation |
| **SNP annotation** | snpEff | Java | More databases |

## Performance Comparison

| Tool | Speed vs QTLmax | Memory | GPU | Notes |
|------|----------------|--------|-----|-------|
| tensorQTL | **100x faster** | Moderate | ✅ Required | For eQTL only |
| GEMMA | Same | High | ❌ | Standard LMM |
| PLINK 2 | 2-5x faster | Low | ❌ | For GWAS/PCA |
| R/qtl2 | Slower | Moderate | ❌ | Only for crosses |

## Cost Comparison

| Solution | Cost | License | Support |
|----------|------|---------|---------|
| QTLmax | $500-2000/year | Commercial | Vendor |
| Open-source stack | **Free** | GPL/BSD | Community |
| This skill | **Free** | Apache-2.0 | Community |

## When to Use What

### Use tensorQTL when:
- eQTL/molQTL mapping
- GPU available
- Large cohorts (>1000 samples)
- Need fast iteration

### Use GEMMA when:
- LMM-GWAS required
- Accounting for relatedness
- Small to medium cohorts
- Exact QTLmax compatibility needed

### Use R/qtl2 when:
- Experimental crosses (F2, BC, RIL, DO)
- Classical QTL mapping
- LOD scores and CIs needed
- No Python alternative exists

### Use PLINK when:
- Data manipulation
- QC filtering
- PCA/admixture
- Industry standard needed
