# Data Formats Reference

Complete guide to file formats used in QTL analysis workflows.

## VCF (Variant Call Format)

Standard format for storing gene sequence variations.

### Format Structure

```
##fileformat=VCFv4.2
##fileDate=20240101
##source=Example
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample1 Sample2 Sample3
1       1000    rs1     A       G       99      PASS    NS=3    GT      0/0     0/1     1/1
1       2000    rs2     C       T       99      PASS    NS=3    GT      0/1     0/0     0/1
```

### Required Columns

| Column | Description | Example |
|--------|-------------|---------|
| `#CHROM` | Chromosome identifier | `1`, `chr1` |
| `POS` | Position (1-based) | `1000`, `1234567` |
| `ID` | Variant identifier | `rs12345`, `.` |
| `REF` | Reference allele | `A`, `CG` |
| `ALT` | Alternative allele(s) | `G`, `T,A` |
| `QUAL` | Phred-scaled quality | `99`, `.` |
| `FILTER` | Filter status | `PASS`, `LowQual` |
| `INFO` | Additional information | `NS=3;AF=0.5` |
| `FORMAT` | Format of genotype fields | `GT:DP:GQ` |
| `[Samples]` | Genotype data per sample | `0/1:30:99` |

### Genotype Encoding

| Code | Meaning |
|------|---------|
| `0/0` | Homozygous reference |
| `0/1` | Heterozygous |
| `1/1` | Homozygous alternate |
| `./.` | Missing |
| `0` | Haploid reference |
| `1` | Haploid alternate |

### Compression and Indexing
```bash
# Compress with bgzip
bgzip variants.vcf

# Index for random access
tabix -p vcf variants.vcf.gz

# Extract region
bcftools view variants.vcf.gz 1:1000000-2000000
```

---

## PLINK Binary Format (.bed/.bim/.fam)

Compact binary format for genotype data.

### .bed (Binary PED)

- **Binary**: SNP-major or individual-major order
- **Magic bytes**: 0x6c 0x1b (SNP-major) or 0x6c 0x1b 0x01 (individual-major)
- **Encoding**: 2 bits per genotype (00=homA1, 01=missing, 10=het, 11=homA2)

### .bim (Extended MAP file)

Tab-delimited, 6 columns:

| Column | Name | Description |
|--------|------|-------------|
| 1 | Chromosome | Chromosome code (1-22, X, Y, XY, MT, 0) |
| 2 | SNP ID | Variant identifier |
| 3 | CM position | Genetic distance (cM) or 0 |
| 4 | BP position | Base-pair position |
| 5 | Allele 1 | Minor allele (usually) |
| 6 | Allele 2 | Major allele (usually) |

### .fam (Family Information)

Tab/space-delimited, 6 columns:

| Column | Name | Description |
|--------|------|-------------|
| 1 | Family ID | Family or population ID |
| 2 | Individual ID | Sample identifier |
| 3 | Paternal ID | Father's ID (0 if unknown) |
| 4 | Maternal ID | Mother's ID (0 if unknown) |
| 5 | Sex | 1=male, 2=female, 0=unknown |
| 6 | Phenotype | -9/0=missing, 1=unaffected, 2=affected, or continuous value |

### PLINK 1.9 Text Format (.ped/.map)

#### .ped Format
Space/tab-delimited, 6+ columns:
- Cols 1-6: Same as .fam
- Cols 7+: Genotypes (two columns per SNP: allele1 allele2)

```
FID IID PID MID Sex Pheno rs1_A rs1_G rs2_C rs2_T
001 S01 0 0 1 1 A A C T
001 S02 0 0 2 2 A G C C
```

#### .map Format
Tab-delimited, 4 columns:
- Chromosome, SNP ID, Genetic position (cM), Physical position (bp)

```
1 rs1 0 1000
1 rs2 0 2000
```

---

## Phenotype Files

### Basic Format

Tab/space-delimited with header:

```
FID IID Pheno1 Pheno2 Covar1 Covar2
001 S01 10.5 20.3 1 0.5
001 S02 12.1 18.7 1 1.2
```

### Long Format (for eQTL)

```
sample_id    expression_A1BG    expression_A2M    ...
Sample1      5.23               7.11               ...
Sample2      4.89               6.92               ...
```

### BED Format for tensorQTL

Tab-delimited, first 4 columns fixed:

| Column | Description |
|--------|-------------|
| `#Chr` | Chromosome |
| `start` | Start position |
| `end` | End position |
| `phenotype_id` | Gene/probe ID |
| `Sample1` | Expression value |
| `Sample2` | Expression value |
| ... | ... |

```
#Chr start end phenotype_id Sample1 Sample2 Sample3
1 1000 2000 ENSG000001234 5.2 4.8 6.1
1 3000 4000 ENSG000002345 7.1 6.9 7.3
```

**Requirements:**
- First 4 columns must match above
- Coordinates are 0-based, half-open [start, end)
- Must be sorted by chromosome and position
- Can be gzip compressed

---

## Covariates Files

### Standard Format

Tab-delimited, first two columns are FID/IID:

```
FID IID Sex Age PCs1 PCs2 PCs3
001 S01 1 25 0.01 0.02 -0.01
001 S02 2 30 0.02 -0.01 0.03
```

### Format for tensorQTL

Samples as columns, covariates as rows:

```
id Sample1 Sample2 Sample3
Sex 1 2 1
Age 25 30 28
PC1 0.01 0.02 -0.01
```

### Principal Components Format

From PLINK `--pca`:

```
#FID IID PC1 PC2 PC3 PC4 PC5
001 S01 0.023 -0.015 0.008 0.001 -0.003
001 S02 -0.018 0.021 -0.007 0.004 0.001
```

---

## Kinship/Relationship Matrices

### GEMMA Format

Space-delimited square matrix:

```
Sample1 1.0 0.5 0.3 0.2
Sample2 0.5 1.0 0.4 0.1
Sample3 0.3 0.4 1.0 0.6
Sample4 0.2 0.1 0.6 1.0
```

- Header row: Sample IDs
- First column: Sample IDs
- Values: pairwise kinship coefficients

### GCTA GRM Format

Binary format:
- `.grm.bin`: Kinship values (N(N+1)/2 elements, 4-byte float)
- `.grm.id`: Sample IDs
- `.grm.N`: Number of SNPs used per pair

### R/qtl2 Kinship

R list with per-chromosome or overall matrices:

```r
$`1`
        Sample1 Sample2 Sample3
Sample1 1.00    0.45    0.32
Sample2 0.45    1.00    0.28
Sample3 0.32    0.28    1.00

$`2`
...
```

---

## Association Results

### GEMMA Output

Space-delimited with header:

```
chr rs ps n_mis n_obs allele1 allele0 af beta se logl_H1 l_remle p_wald p_lrt p_score
1 rs1 1000 0 100 A G 0.25 0.5 0.2 -50.2 0.85 0.01 0.008 0.012
```

### PLINK Output

Tab-delimited:

```
#CHROM POS ID REF ALT A1 TEST OBS_CT BETA SE T_STAT P
1 1000 rs1 A G G ADD 100 0.5 0.2 2.5 0.0123
```

### tensorQTL Output

Tab-delimited from `map_cis`:

| Column | Description |
|--------|-------------|
| `phenotype_id` | Gene/expression ID |
| `variant_id` | Lead variant ID |
| `tss_distance` | Distance from TSS (bp) |
| `ma_samples` | Number of samples with minor allele |
| `ma_count` | Minor allele count |
| `maf` | Minor allele frequency |
| `pval_nominal` | Nominal p-value |
| `slope` | Effect size (beta) |
| `slope_se` | Standard error |
| `pval_perm` | Empirical p-value (if permutations run) |
| `pval_beta` | Beta-approximated p-value |

---

## R/qtl2 Cross Data Format

YAML-based format for complex crosses.

### Directory Structure

```
my_cross/
├── cross.yaml          # Control file
├── geno.csv            # Genotype data
├── gmap.csv            # Genetic map
├── pmap.csv            # Physical map (optional)
├── pheno.csv           # Phenotypes
├── covar.csv           # Covariates (optional)
├── phenocovar.csv      # Phenotype covariates (optional)
└── founder_geno.csv    # Founder genotypes (optional)
```

### cross.yaml (Control File)

```yaml
crosstype: f2       # f2, bc, riself, risib, dh, haploid, do, etc.
geno: geno.csv
geno_transposed: true
pheno: pheno.csv
gmap: gmap.csv
alleles:
  - S
  - B
geno_codes:
  S: 1
  H: 2
  B: 3
  -: -
sex_covar: sex     # Column in covar.csv for sex
```

### geno.csv (Genotypes)

Markers as rows, individuals as columns:

```
marker,Sample1,Sample2,Sample3,...
rs1,S,H,B,...
rs2,H,S,S,...
...
```

Genotype codes:
- `S`, `B`: Parental homozygotes
- `H`: Heterozygote
- `-`: Missing

### gmap.csv (Genetic Map)

```
marker,chr,pos
rs1,1,0.0
rs2,1,2.5
rs3,1,5.1
...
```

### pheno.csv (Phenotypes)

Individuals as rows, phenotypes as columns:

```
id,pheno1,pheno2,...
Sample1,10.5,20.3,...
Sample2,12.1,18.7,...
...
```

---

## Format Conversion Quick Reference

| From | To | Tool | Command |
|------|-----|------|---------|
| VCF | PLINK | PLINK2 | `plink2 --vcf in.vcf --make-bed --out out` |
| PLINK | VCF | PLINK2 | `plink2 --bfile in --recode vcf --out out` |
| VCF | BGEN | QCTOOL | `qctool -g in.vcf -og out.bgen` |
| PLINK | BED | tensorQTL | `tensorqtl.read_plink(prefix)` |
| CSV | PLINK | Custom | R/Python scripts |
| Text | Binary | PLINK | `plink --file text --make-bed --out bin` |

---

## Validation Commands

```bash
# Check VCF
bcftools stats file.vcf.gz > vcf_stats.txt

# Check PLINK
plink2 --bfile prefix --validate

# Check for duplicates
plink2 --bfile prefix --list-duplicate-vars

# Check sex concordance
plink2 --bfile prefix --check-sex

# Check missingness
plink2 --bfile prefix --missing
```
