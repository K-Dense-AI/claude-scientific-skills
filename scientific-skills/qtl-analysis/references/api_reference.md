# API Reference

Complete API documentation for all tools in the QTL Analysis skill.

## tensorQTL (cis-eQTL/mapping)

Python library for fast eQTL mapping using PyTorch.

### Installation
```bash
pip install tensorqtl
```

### Python API

#### `tensorqtl.read_phenotype_bed(phenotype_bed)`
Load phenotype data from BED format.

**Parameters:**
- `phenotype_bed` (str): Path to phenotype BED file

**Returns:**
- `phenotype_df` (pandas.DataFrame): Phenotype expression matrix
- `phenotype_pos_df` (pandas.DataFrame): Phenotype positions

#### `tensorqtl.read_plink(plink_prefix, verbose=True)`
Load genotype data from PLINK format.

**Parameters:**
- `plink_prefix` (str): Prefix for PLINK files (.bed/.bim/.fam)
- `verbose` (bool): Print progress information

**Returns:**
- `genotype_df` (pandas.DataFrame): Genotype dosage matrix
- `variant_df` (pandas.DataFrame): Variant information

#### `tensorqtl.cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df=None, window=1000000, nperm=10000)`
Perform cis-eQTL mapping.

**Parameters:**
- `genotype_df` (DataFrame): Genotype dosage matrix (samples x variants)
- `variant_df` (DataFrame): Variant information with columns: chrom, pos
- `phenotype_df` (DataFrame): Phenotype expression matrix (samples x phenotypes)
- `phenotype_pos_df` (DataFrame): Phenotype positions with columns: chrom, start, end
- `covariates_df` (DataFrame, optional): Covariates matrix
- `window` (int): Cis window size in base pairs (default: 1Mb)
- `nperm` (int): Number of permutations for empirical p-values

**Returns:**
- `cis_df` (DataFrame): Results with columns: phenotype_id, variant_id, pval_nominal, slope, slope_se, etc.

**Example:**
```python
import tensorqtl
from tensorqtl import read_phenotype_bed, read_plink, cis

# Load data
phenotype_df, phenotype_pos_df = read_phenotype_bed('expression.bed.gz')
genotype_df, variant_df = read_plink('genotypes', verbose=True)

# Map cis-eQTLs
cis_df = cis.map_cis(
    genotype_df, variant_df, phenotype_df, phenotype_pos_df,
    window=500000, nperm=10000
)
```

#### `tensorqtl.cis.map_independent(genotype_df, variant_df, phenotype_df, phenotype_pos_df, cis_df, n_cis=5)`
Find independent secondary eQTLs.

**Parameters:**
- `cis_df` (DataFrame): Primary cis-eQTL results
- `n_cis` (int): Maximum number of independent signals per gene

**Returns:**
- `indep_df` (DataFrame): Independent eQTL results

---

## GEMMA (GWAS-LMM)

Linear Mixed Model for genome-wide association studies.

### CLI Usage

#### Association Analysis
```bash
gemma -bfile [prefix] -k [kinship] -lmm [num] -o [prefix]
```

**Required Flags:**
- `-bfile [prefix]` : Prefix for PLINK binary files (.bed/.bim/.fam)
- `-k [file]` : Kinship matrix file
- `-lmm [num]` : LMM test type (1=Wald, 2=likelihood ratio, 4=score)

**Optional Flags:**
- `-c [file]` : Covariates file
- `-p [file]` : Phenotype file (if not in .fam)
- `-o [prefix]` : Output prefix
- `-maf [float]` : Minor allele frequency threshold
- `-miss [float]` : Missingness threshold
- `-r2 [float]` : R² threshold for pruning
- `-n [int]` : Number of eigenvectors to use

**Example:**
```bash
# Generate kinship matrix first
gemma -bfile genotypes -gk 1 -o kinship

# Run GWAS with LMM
gemma -bfile genotypes -k output/kinship.cXX.txt -lmm 4 -c covariates.txt -o gwas_results
```

#### Kinship Matrix Generation
```bash
gemma -bfile [prefix] -gk [num] -o [prefix]
```

**Flag `-gk` values:**
- `1` : Centered kinship (GEMMA default)
- `2` : Standardized kinship
- `3` : Centered and scaled kinship
- `4` : Identity matrix

#### Output Files

| File | Description |
|------|-------------|
| `*.assoc.txt` | Association results |
| `*.log.txt` | Log file |
| `*.cXX.txt` | Kinship matrix |
| `*.eigenD.txt` | Eigenvalues |
| `*.eigenU.txt` | Eigenvectors |

#### Association Results Format

| Column | Description |
|--------|-------------|
| `chr` | Chromosome |
| `rs` | SNP ID |
| `ps` | Physical position (bp) |
| `n_mis` | Number of missing genotypes |
| `n_obs` | Number of observed genotypes |
| `allele1` | Reference allele |
| `allele0` | Alternative allele |
| `af` | Allele frequency |
| `beta` | Effect size estimate |
| `se` | Standard error |
| `logl_H1` | Log-likelihood under H1 |
| `l_remle` | REML estimate of lambda |
| `p_wald` | Wald test p-value |
| `p_lrt` | Likelihood ratio test p-value |
| `p_score` | Score test p-value |

---

## PLINK (genotype processing)

Swiss-army knife for genotype data manipulation.

### Data Conversion

#### VCF to PLINK
```bash
plink2 --vcf [file.vcf] --make-bed --out [prefix]
```

#### Text to Binary
```bash
plink --file [prefix] --make-bed --out [prefix]
```

### Quality Control

#### MAF and Missingness Filtering
```bash
plink2 --bfile [prefix] \
       --maf 0.05 \
       --geno 0.02 \
       --hwe 1e-6 \
       --make-bed \
       --out [prefix]_filtered
```

**Flags:**
- `--maf [float]` : Minor allele frequency threshold
- `--geno [float]` : Maximum per-variant missingness
- `--mind [float]` : Maximum per-sample missingness
- `--hwe [float]` : Hardy-Weinberg p-value threshold
- `--max-allele-ct [int]` : Maximum allele count (biallelic: 2)

#### LD Pruning
```bash
plink2 --bfile [prefix] \
       --indep-pairwise 50 10 0.1 \
       --out [prefix]_pruned
```

### Association Testing

#### Basic Association
```bash
plink2 --bfile [prefix] \
       --pheno [pheno.txt] \
       --allow-no-sex \
       --glm allow-no-covars \
       --out [prefix]_assoc
```

#### Linear Regression with Covariates
```bash
plink2 --bfile [prefix] \
       --pheno [pheno.txt] \
       --covar [covar.txt] \
       --glm hide-covar \
       --out [prefix]_assoc
```

### PCA for Population Structure
```bash
plink2 --bfile [prefix] \
       --pca [count] \
       --out [prefix]_pca
```

---

## R/qtl2 (Classical QTL)

Next-generation QTL analysis for multi-parent populations.

### R API

#### Data Input

##### `read_cross2(file, overwrite=FALSE)`
Read cross data in qtl2 format.

**Parameters:**
- `file` (character): Path to YAML control file or zip file
- `overwrite` (logical): Overwrite existing data if TRUE

**Returns:**
- Object of class `"cross2"`

##### `read_pheno(file)`
Read phenotype data from CSV.

##### `read_geno(file)`
Read genotype data from CSV.

##### `read_map(file)`
Read genetic map from CSV.

#### Simulation

##### `sim_cross(map, n_ind, type=c("f2", "bc", "riself", "risib", "dh", "haploid"), ...)`
Simulate a cross.

**Parameters:**
- `map` (list): Genetic map (from `insert_pseudomarkers`)
- `n_ind` (int): Number of individuals
- `type` (character): Cross type

**Returns:**
- Object of class `"cross2"`

##### `sim_geno(cross, n_draws=1, step=0, off_end=0, error_prob=0.0001, map_function=c("haldane", "kosambi", "c-f", "morgan"))`
Simulate genotypes from inferred genotypes.

#### QTL Mapping

##### `scan1(genoprobs, pheno, kinship=NULL, addcovar=NULL, Xcovar=NULL, intcovar=NULL, weights=NULL, reml=TRUE, model=c("normal", "binary"), ...)`
Genome scan with a single-QTL model.

**Parameters:**
- `genoprobs` : Genotype probabilities (from `calc_genoprob`)
- `pheno` (matrix): Phenotypes (individuals x traits)
- `kinship` (list): Kinship matrices (optional)
- `addcovar` (matrix): Additive covariates
- `Xcovar` (matrix): Sex covariates for X chromosome
- `intcovar` (matrix): Interactive covariates
- `reml` (logical): Use REML (TRUE) or ML (FALSE)
- `model` (character): "normal" or "binary"

**Returns:**
- Object of class `"scan1"`

##### `scan1perm(genoprobs, pheno, n_perm=1000, ...)`
Permutation test for `scan1`.

**Parameters:**
- `n_perm` (int): Number of permutations

**Returns:**
- Matrix of permutation LOD scores

##### `find_peaks(scan1_output, map, threshold=3, peakdrop=Inf, prob=0.95, expandtomarkers=FALSE)`
Find peaks in scan output.

**Parameters:**
- `scan1_output` : Output from `scan1`
- `map` (list): Genetic map
- `threshold` (float): LOD threshold for peaks
- `peakdrop` (float): Drop in LOD required for distinct peaks
- `prob` (float): Probability for Bayes credible interval

**Returns:**
- Data frame with peak information

#### Multi-QTL Models

##### `fitqtl(cross, pheno.col=1, qtl, covar=NULL, formula, method=c("imp", "hk", "em"), model=c("normal", "binary"))`
Fit a QTL model with multiple QTLs.

##### `makeqtl(cross, chr, pos, qtl.name=NULL)`
Create a QTL object for model fitting.

##### `addqtl(cross, qtl, covar=NULL, formula, method=c("imp", "hk", "em"))`
Add one QTL to an existing model.

##### `addpair(cross, qtl, covar=NULL, formula, method=c("imp", "hk", "em"))`
Add a pair of QTLs (epistasis scan).

#### Genotype Probabilities

##### `calc_genoprob(cross, step=0, off_end=0, error_prob=0.0001, map_function=c("haldane", "kosambi", "c-f", "morgan"))`
Calculate conditional genotype probabilities.

##### `genoprob_to_alleleprob(genoprobs)`
Convert genotype probabilities to allele probabilities.

#### Kinship

##### `calc_kinship(cross, type=c("overall", "loco"), omit_x=FALSE, use_allele_probs=TRUE)`
Calculate kinship matrix from genotype probabilities.

**Parameters:**
- `type` (character): "overall" or "loco" (leave-one-chromosome-out)
- `omit_x` (logical): Omit X chromosome
- `use_allele_probs` (logical): Use allele probabilities

### Complete Example

```r
library(qtl2)

# Read data
iron <- read_cross2("https://raw.githubusercontent.com/rqtl/qtl2data/main/iron/iron.zip")

# Insert pseudomarkers
map <- insert_pseudomarkers(iron$gmap, step=1)

# Calculate genotype probabilities
pr <- calc_genoprob(iron, map, error_prob=0.002)

# Calculate kinship
kinship <- calc_kinship(pr)

# Scan for QTLs
out <- scan1(pr, iron$pheno[, "liver"], kinship)

# Find peaks
peaks <- find_peaks(out, map, threshold=3)

# Plot results
plot(out, map)
```

---

## Quick Command Reference

| Task | Tool | Command |
|------|------|---------|
| cis-eQTL mapping | tensorQTL | `tensorqtl.cis.map_cis()` |
| GWAS with LMM | GEMMA | `gemma -bfile geno -k kin -lmm 4` |
| QC and filtering | PLINK | `plink2 --bfile in --maf 0.05 --make-bed --out out` |
| Classical QTL scan | R/qtl2 | `scan1(genoprobs, pheno)` |
| Kinship matrix | GEMMA | `gemma -bfile geno -gk 1` |
| Kinship matrix | R/qtl2 | `calc_kinship(genoprobs)` |
| PCA | PLINK | `plink2 --bfile geno --pca 10` |
| Convert VCF→PLINK | PLINK | `plink2 --vcf file.vcf --make-bed --out prefix` |
| Permutation test | R/qtl2 | `scan1perm(genoprobs, pheno, n_perm=1000)` |
