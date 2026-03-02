<!-- Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC) -->
<!-- Licensed under the Apache License, Version 2.0. -->

# Minimal Working Examples for QTL Analysis Tools

## Executive Summary

This document provides minimum viable data sizes, file format examples, and actual commands for running QTL analyses using tensorQTL, R/qtl2, and GEMMA.

---

## 1. TENSORQTL (GPU-Accelerated eQTL Mapping)

### Repository
- **GitHub**: https://github.com/broadinstitute/tensorqtl
- **Example Notebook**: https://github.com/broadinstitute/tensorqtl/blob/master/example/tensorqtl_examples.ipynb

### Minimum Viable Data Size

| Parameter | Minimum | Recommended | Notes |
|-----------|---------|-------------|-------|
| **Samples** | 50-100 | 100+ | GEUVADIS example uses 462 samples |
| **Genes** | 10-50 | 1000+ | For testing/validation |
| **SNPs** | 1,000-10,000 | 100,000+ | Per chromosome for cis-eQTL |
| **Memory** | ~50GB | 100GB+ | For full genome analysis |
| **GPU** | Required | CUDA-capable | 200x speedup vs CPU |

### File Formats

**Genotype Files (PLINK format):**
- `.bed` - Binary genotype data
- `.bim` - SNP map information
- `.fam` - Sample information

**Phenotype File (BED format):**
```
# chr start end gene_id sample1 sample2 sample3 ...
chr1 1000 2000 gene_A 5.2 3.1 4.5 ...
chr1 3000 4000 gene_B 2.1 6.7 3.2 ...
```

**Covariates File:**
- Tab-delimited with header
- First column: sample IDs
- Subsequent columns: covariate values

### Example Commands

```bash
# Install tensorQTL
pip install tensorqtl

# Cis-eQTL mapping (basic)
python -m tensorqtl \
    genotypes.bed \
    expression.bed \
    output_prefix \
    --mode cis

# Cis-eQTL with covariates
python -m tensorqtl \
    genotypes.bed \
    expression.bed \
    output_prefix \
    --mode cis \
    --covariates covariates.txt

# Trans-eQTL mapping
python -m tensorqtl \
    genotypes.bed \
    expression.bed \
    output_prefix \
    --mode trans

# Cis-eQTL with permutation (for empirical p-values)
python -m tensorqtl \
    genotypes.bed \
    expression.bed \
    output_prefix \
    --mode cis \
    --permutations 10000
```

### Expected Runtime
- **Cis-eQTL**: ~5-30 minutes for 20,000 genes × 10,000 SNPs (GPU)
- **Trans-eQTL**: Hours to days depending on genome-wide search space
- **GEUVADIS example** (chr18 only): ~10 minutes

### Output Format
```
gene_id variant_id tss_distance ma_samples ma_count maf pval_nominal slope slope_se pval_nominal_threshold min_pval_nominal pval_beta
```

---

## 2. R/QTL2 (Classical QTL Mapping)

### Repository
- **GitHub**: https://github.com/rqtl/qtl2
- **Example Data**: https://github.com/rqtl/qtl2data
- **Documentation**: https://kbroman.org/qtl2/

### Minimum Viable Data Size

| Parameter | Minimum | Recommended | Notes |
|-----------|---------|-------------|-------|
| **Samples** | 50-100 | 200+ | F2/BC populations |
| **RIL Lines** | 25-50 | 100+ | For RIL-based mapping |
| **Markers** | 100-500 | 1000+ | Genetic map density |
| **Phenotypes** | 1 | 10+ | Can map multiple traits |

### Smallest Cross Types (Easiest to Start)

1. **Backcross (BC)** - Simplest
   - 2 genotypes (AA, AB)
   - No cross_info needed
   - No sex covariate (unless X chr)

2. **F2 Intercross** - Most common
   - 3 genotypes (AA, AB, BB)
   - May need cross_info for X chr
   - Sex covariate for X chromosome

3. **RIL (Recombinant Inbred Lines)**
   - 2 genotypes (homozygous only)
   - No sex covariate needed
   - Good for replication

### File Format Structure

**Control File (YAML format):**
```yaml
# Data from Grant et al. (2006)
crosstype: f2
geno: iron_geno.csv
pheno: iron_pheno.csv
gmap: iron_gmap.csv
alleles:
  - S
  - B
genotypes:
  SS: 1
  SB: 2
  BB: 3
sex:
  covar: sex
  f: female
  m: male
x_chr: X
na.strings:
  - '-'
  - NA
```

**Genotype File (CSV):**
```csv
marker,sample1,sample2,sample3,...
rs1234567,SS,SB,BB,...
rs2345678,BB,SS,SB,...
```

**Phenotype File (CSV):**
```csv
id,sex,pheno1,pheno2
sample1,f,12.5,3.2
sample2,m,15.2,4.1
```

**Genetic Map File (CSV):**
```csv
marker,chr,pos
rs1234567,1,12.5
rs2345678,1,23.7
```

### Example R Commands

```r
# Install qtl2
install.packages("qtl2")
library(qtl2)

# Read data
grav2 <- read_cross2("~/my_data/grav2.yaml")

# Or read from zip
grav2 <- read_cross2("https://kbroman.org/qtl2/assets/sampledata/grav2/grav2.zip")

# Calculate genotype probabilities
pr <- calc_genoprob(grav2, map=grav2$gmap, error_prob=0.002)

# Perform genome scan
out <- scan1(pr, grav2$pheno)

# Find peaks
peaks <- find_peaks(out, grav2$gmap, threshold=3)

# For F2 with X chromosome (need cross_info)
# Example: iron dataset
iron <- read_cross2("iron.yaml")
pr <- calc_genoprob(iron, map=iron$gmap, error_prob=0.002)
out <- scan1(pr, iron$pheno[,1])  # First phenotype
```

### Sample Datasets (from qtl2data)

| Dataset | Type | Samples | Markers | Description |
|---------|------|---------|---------|-------------|
| **grav2** | RIL | 162 | 117 | Arabidopsis gravity response |
| **iron** | F2 | 284 | 66 | Mouse iron levels |
| **DOex** | DO | 500 | 1,000 | Diversity Outbred (subset) |
| **B6BTBR** | F2 | 450 | 1,000 | Mouse intercross |

### Expected Runtime
- **Genotype probabilities**: Seconds to minutes
- **Genome scan**: Seconds for simple crosses
- **Permutation testing**: Minutes to hours

---

## 3. GEMMA (GWAS with Linear Mixed Models)

### Repository
- **GitHub**: https://github.com/genetics-statistics/GEMMA
- **Manual**: https://www.xzlab.org/software/GEMMAmanual.pdf

### Minimum Viable Data Size

| Parameter | Minimum | Recommended | Notes |
|-----------|---------|-------------|-------|
| **Samples** | 100 | 500+ | For adequate power |
| **SNPs** | 10,000 | 100,000+ | For relatedness matrix |
| **Phenotypes** | 1 | 1+ | Univariate or multivariate |
| **Covariates** | 0 | 3-10 | Sex, age, PCs, etc. |

### File Formats

**PLINK Binary Format (Input):**
- `genotypes.bed` - Binary genotype data
- `genotypes.bim` - SNP information
- `genotypes.fam` - Sample information

**Phenotype File (optional, space/tab-delimited):**
```
# Sample format (no header):
FID IID pheno1 pheno2
1 sample1 1.5
2 sample2 2.3
```

**Covariates File:**
```
# First two columns: FID IID
# Subsequent columns: covariates
FID IID sex age PC1 PC2
1 sample1 1 25 0.1 0.2
2 sample2 2 30 -0.1 0.3
```

**Relatedness Matrix (K, optional):**
- Can be computed from GEMMA or provided

### Example Commands

```bash
# Install GEMMA
# Download from GitHub releases

# 1. Compute relatedness matrix (centered)
gemma -bfile genotypes \
      -gk 1 \
      -o output_prefix

# 2. Run LMM association (univariate)
gemma -bfile genotypes \
      -k output_prefix.cXX.txt \
      -p phenotypes.txt \
      -lmm 4 \
      -o gwas_results

# 3. Run LMM with covariates
gemma -bfile genotypes \
      -k output_prefix.cXX.txt \
      -p phenotypes.txt \
      -c covariates.txt \
      -lmm 4 \
      -o gwas_results_covar

# 4. Multivariate LMM (multiple phenotypes)
gemma -bfile genotypes \
      -k output_prefix.cXX.txt \
      -p phenotypes.txt \
      -n 3 \  # Number of phenotypes
      -lmm 4 \
      -o mv_gwas_results

# 5. Bayesian Sparse LMM (BSLMM)
gemma -bfile genotypes \
      -p phenotypes.txt \
      -bslmm 1 \
      -o bslmm_results
```

### LMM Test Options
- `-lmm 1`: Wald test
- `-lmm 2`: Likelihood ratio test
- `-lmm 3`: Score test
- `-lmm 4`: All tests (recommended)

### Expected Runtime
- **Relatedness matrix**: Minutes to hours (depends on N)
- **GWAS (LMM)**: Minutes to hours
- **BSLMM**: Hours to days (MCMC)

### Output Files

**Association Results:**
```
chr rs ps n_maf beta se p_wald p_lrt p_score
1 rs1234567 1000 0.25 0.5 0.1 0.001 0.002 0.001
```

**Relatedness Matrix:**
- `prefix.cXX.txt` - Centered relatedness matrix

---

## 4. MINIMUM SAMPLE SIZE GUIDELINES

### GWAS with LMM

| Effect Size | MAF | Power | Sample Size |
|-------------|-----|-------|-------------|
| 0.5 SD | 0.3 | 80% | ~100-200 |
| 0.3 SD | 0.2 | 80% | ~300-500 |
| 0.2 SD | 0.1 | 80% | ~1000+ |

**Rule of thumb**: Minimum 100 samples for proof-of-concept; 500+ for publication-quality results.

### eQTL Mapping

| Study Type | Minimum Samples | Notes |
|------------|-----------------|-------|
| **Cis-eQTL** | 50-100 | Per tissue type |
| **Trans-eQTL** | 200+ | Requires more power |
| **Single-cell eQTL** | 100+ | Cell type specific |

**GTEx standards**: 70+ samples per tissue for cis-eQTL; 200+ for trans-eQTL.

### Classical QTL (Experimental Crosses)

| Cross Type | Minimum | Recommended | Notes |
|------------|---------|-------------|-------|
| **F2** | 100 | 200-400 | Intercross |
| **Backcross** | 50 | 100-200 | Simpler genetics |
| **RIL** | 25 lines | 50-100 lines | Each line can have replicates |
| **DO/CC** | 100 | 200+ | Multi-parent populations |

### Population Structure / PCA

| Analysis | Minimum Samples | Minimum SNPs | Notes |
|----------|-----------------|--------------|-------|
| **PCA** | 10-20 per group | 1,000-10,000 | For population stratification |
| **Structure inference** | 30+ per population | 10,000+ | For ancestry analysis |
| **Relatedness** | 50+ | 10,000+ | For cryptic relatedness |

---

## 5. EXAMPLE DATA SOURCES

### Public Test Datasets

1. **GEUVADIS** (tensorQTL)
   - 462 samples
   - RNA-seq from lymphoblastoid cell lines
   - EUR and YRI populations
   - URL: https://www.ebi.ac.uk/Tools/geuvadis-das/

2. **R/qtl2 Sample Data**
   - grav2: Arabidopsis RIL (162 samples, 117 markers)
   - iron: Mouse F2 (284 samples, 66 markers)
   - DOex: Mouse DO subset (500 samples)
   - URL: https://kbroman.org/qtl2/pages/sampledata.html

3. **GEMMA Example Data**
   - PLINK-formatted test datasets
   - Available in repository

4. **PLINK Tutorial Data**
   - HapMap data
   - URL: https://www.cog-genomics.org/plink/1.9/resources

---

## 6. QUICK START WORKFLOW

### For Testing/Development

```bash
# 1. Create minimal PLINK dataset (100 samples, 1000 SNPs)
# Use PLINK to subset existing data
plink --bfile large_dataset \
      --thin 0.01 \
      --make-bed \
      --out minimal_test

# 2. Run GEMMA on minimal data
gemma -bfile minimal_test -gk 1 -o test_kinship
gemma -bfile minimal_test -k test_kinship.cXX.txt -p pheno.txt -lmm 4 -o test_gwas

# 3. Run tensorQTL on single chromosome
python -m tensorqtl genotypes_chr22.bed expression.bed test_eqtl --mode cis

# 4. Run R/qtl2 on sample data
R -e 'library(qtl2); grav2 <- read_cross2("grav2.zip"); pr <- calc_genoprob(grav2); out <- scan1(pr, grav2$pheno)'
```

---

## 7. FILE FORMAT CONVERSION

### VCF to PLINK
```bash
plink --vcf input.vcf --make-bed --output output_prefix
```

### PLINK to R/qtl2
```r
library(qtl2convert)
cross <- read_cross2("data.yaml")  # Use control file
```

### Text to BED (tensorQTL)
```python
import pandas as pd
# Format: chr start end gene_id sample1 sample2 ...
df = pd.read_csv("expression.txt", sep="\t")
df.to_parquet("expression.bed.parquet")  # Or use PLINK BED
```

---

## References

1. Taylor-Weiner et al. (2019). TensorQTL: GPU-based QTL mapper. Genome Biology.
2. Broman et al. (2019). R/qtl2: QTL mapping in experimental crosses. Genetics.
3. Zhou et al. (2012). GEMMA: Genome-wide Efficient Mixed Model Association.
4. GTEx Consortium (2020). GTEx Analysis V8.
5. Hong & Park (2012). Sample Size and Statistical Power in Genetic Association Studies.
