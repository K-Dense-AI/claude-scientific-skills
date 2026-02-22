---
name: qtl-analysis
description: >
  Comprehensive QTL (Quantitative Trait Loci) analysis toolkit — GWAS, eQTL mapping, 
  classical LOD scans, kinship matrices, population structure, and genomic prediction.
  Matches commercial QTLmax capabilities using open-source tools: tensorQTL (GPU-accelerated 
  eQTL), GEMMA (LMM-GWAS), PLINK (GWAS/population structure), and R/qtl2 (classical QTL).
  Use for experimental crosses (F2, BC, RIL, DO), association mapping, expression QTLs, 
  and breeding applications.
license: Apache-2.0
metadata:
  author: Clayton Young (@borealBytes)
  version: "1.0.0"
  category: genomics
  tools: [tensorQTL, GEMMA, PLINK, R/qtl2, pyQTL]
---

# QTL Analysis

## Overview

Open-source QTL analysis toolkit matching commercial QTLmax capabilities. Python-native 
where possible (tensorQTL, pyQTL), CLI wrappers where necessary (GEMMA, PLINK, R/qtl2).

**Key insight:** No single tool covers everything. This skill orchestrates the best 
open-source tools for each QTL analysis type.

## When to Use This Skill

Use this skill when:

- Performing **GWAS** (genome-wide association studies) on inbred lines or natural populations
- Mapping **eQTLs** (expression quantitative trait loci) with RNA-seq data
- Analyzing **experimental crosses** (F2, backcross, RIL, DO, MPP) with classical QTL methods
- Calculating **kinship matrices** (VanRaden, genomic coancestry)
- Estimating **population structure** (PCA, admixture)
- Running **genomic prediction** (GBLUP, marker-assisted selection)
- Creating **Manhattan/QQ plots** from GWAS results

**Tool Selection Guide:**

| Analysis Type | Primary Tool | When to Use |
|-------------|--------------|-------------|
| **eQTL mapping** | tensorQTL | GPU-accelerated, 100x faster than alternatives |
| **LMM-GWAS** | GEMMA | Same algorithm as QTLmax, accounts for relatedness |
| **GLM-GWAS** | PLINK 2 | Fast, handles large cohorts, industry standard |
| **Classical QTL** | R/qtl2 | Only option for experimental crosses (LOD scans) |
| **Visualization** | pyQTL + matplotlib | Manhattan, QQ, LocusZoom plots |
| **Population structure** | PLINK --pca | Fast PCA on genotype data |
| **Kinship** | GEMMA -gk or numpy | VanRaden method, standard for GWAS |

## Installation

### Step 1: Run System Check

```bash
python scripts/check_system.py
```

### Step 2: Install Dependencies

```bash
# Option 1: Quick install script
bash scripts/install_deps.sh

# Option 2: Manual installation

# Python packages
pip install tensorqtl qtl pandas numpy matplotlib scipy torch

# CLI tools (conda recommended)
conda install -c bioconda plink gemma r-qtl2
```

**GPU Support (tensorQTL):**

```bash
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu121
pip install tensorqtl
```

## Quick Start

### 1. GWAS with Linear Mixed Model (GEMMA)

```bash
# Calculate kinship matrix
gemma -bfile data -gk 1 -o kinship

# Run LMM-GWAS
gemma -bfile data -k output/kinship.sXX.txt \
      -lmm 4 -p phenotype.txt -o gwas_results

# Visualize
python -m qtl_analysis.plot --manhattan output/gwas_results.assoc.txt \
       --output manhattan.png
```

### 2. cis-eQTL Mapping (tensorQTL)

```python
import tensorqtl
from tensorqtl import genotypeio, cis

# Load genotype and expression data
genotype_df = genotypeio.load_genotypes('genotypes.vcf.gz')
expression_df = pd.read_parquet('expression.bed.gz')

# Run cis-eQTL mapping
cis_df = cis.map_cis(genotype_df, expression_df, 
                      covariates_df, window=1000000)

# Results: gene_id, variant_id, tss_distance, pval, beta, etc.
cis_df.to_csv('cis_eqtl_results.csv')
```

### 3. Classical QTL LOD Scan (R/qtl2)

```bash
# R/qtl2 via CLI wrapper
python scripts/qtl_cli.py lodscan \
  --cross data/cross.json \
  --method haley-knott \
  --perms 1000 \
  --output lod_results/
```

### 4. Population Structure (PLINK)

```bash
# PCA
plink --bfile data --pca 10 --out pca_results

# Admixture (requires separate install)
admixture data.bed 3
```

## Unified CLI

```bash
# GWAS
qtl_cli.py gwas --geno data.vcf --pheno traits.csv --method lmm

# eQTL
qtl_cli.py eqtl --geno data.bed --expr expression.bed.gz --mode cis

# Classical QTL
qtl_cli.py lodscan --cross data.zip --perms 1000

# Visualization
qtl_cli.py manhattan --input pvalues.csv --output plot.png
qtl_cli.py qqplot --input pvalues.csv --output qq.png
qtl_cli.py kinship --geno data.vcf --method vanraden
```

## Input Formats

**Genotypes:**
- VCF/VCF.GZ — standard variant format
- PLINK (.bed/.bim/.fam) — binary PLINK format
- CSV/TSV — custom genotype matrices

**Phenotypes:**
- CSV with columns: ID, trait1, trait2, ...
- PLINK phenotype format (.phe)

**Expression (eQTL):**
- BED format (chr, start, end, gene_id, sample1, sample2, ...)
- CSV/TSV with genes as rows, samples as columns

**Cross Data (classical QTL):**
- R/qtl2 JSON format
- CSV genotypes + CSV phenotypes

## Core Workflows

### Workflow 1: Complete GWAS Pipeline

```bash
# 1. QC and format conversion
plink --vcf raw_data.vcf --make-bed --out data_qc \
  --maf 0.05 --geno 0.05 --hwe 1e-6

# 2. Calculate kinship
gemma -bfile data_qc -gk 1 -o kinship

# 3. Run LMM-GWAS
gemma -bfile data_qc -k kinship.sXX.txt \
      -lmm 4 -p phenotype.txt -o gwas_results

# 4. Create Manhattan plot
python scripts/plot_manhattan.py \
  --input gwas_results.assoc.txt \
  --output manhattan.png \
  --significance 5e-8
```

### Workflow 2: eQTL Discovery

```python
import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis, trans

# Load data
genotype_df, variant_df = genotypeio.load_genotypes('genotypes.vcf.gz')
expression_df = pd.read_parquet('expression.bed.gz')
covariates_df = pd.read_csv('covariates.csv', index_col=0)

# cis-eQTL mapping
cis_results = cis.map_cis(
    genotype_df, expression_df, 
    covariates_df=covariates_df,
    window=1000000,  # 1Mb window
    nperm=10000
)

# Output top hits
top_hits = cis_results[cis_results['pval'] < 5e-8]
print(f"Found {len(top_hits)} significant cis-eQTLs")
```

### Workflow 3: Classical QTL in Experimental Cross

```python
import subprocess
import json

# Prepare cross data in R/qtl2 format
cross_data = {
    "geno": genotypes_df.to_dict(),
    "gmap": genetic_map_df.to_dict(),
    "pheno": phenotypes_df.to_dict(),
    "cross_type": "f2"  # f2, bc, riself, do, etc.
}

with open('cross.json', 'w') as f:
    json.dump(cross_data, f)

# Run R/qtl2 via CLI wrapper
result = subprocess.run([
    'Rscript', '-e',
    '''
    library(qtl2)
    cross <- read_cross2("cross.json")
    map <- insert_pseudomarkers(cross$gmap, step=1)
    pr <- calc_genoprob(cross, map, error_prob=0.002)
    out <- scan1(pr, cross$pheno)
    write.csv(out, "lod_scores.csv")
    '''
], capture_output=True, text=True)

# Plot LOD curves
python scripts/plot_lod.py --input lod_scores.csv --output lod_curve.png
```

### Workflow 4: Population Structure Analysis

```bash
# PCA
plink --bfile data --pca 10 --out pca_results

# Plot PCA
python scripts/plot_pca.py \
  --eigenvec pca_results.eigenvec \
  --eigenval pca_results.eigenval \
  --metadata population_labels.csv \
  --output pca_plot.png

# Kinship heatmap
python scripts/plot_kinship.py \
  --kinship kinship_matrix.txt \
  --metadata populations.csv \
  --output kinship_heatmap.png
```

## Common Pitfalls

1. **Wrong file format** — PLINK expects .bed/.bim/.fam, not VCF. Convert first.
2. **Missing kinship** — LMM-GWAS requires kinship. Don't skip this step.
3. **Population stratification** — Always check PCA. Structure = false positives.
4. **Multiple testing** — Use Bonferroni (5e-8 for human) or FDR, not raw p-values.
5. **cis-eQTL window** — Default 1Mb from TSS. Too wide = more noise.
6. **R/qtl2 cross type** — Must specify correctly (f2, bc, riself, do). Wrong type = wrong probabilities.
7. **GPU memory** — tensorQTL needs VRAM. Use CPU fallback for large datasets.
8. **Genetic vs physical map** — Classical QTL needs centiMorgans, not bp.

## Quality Checklist

After every QTL analysis:

- [ ] **Input QC:** MAF > 5%, call rate > 95%, HWE p > 1e-6
- [ ] **Population structure:** PCA shows expected clusters
- [ ] **Kinship:** Matrix reflects known relationships
- [ ] **GWAS calibration:** QQ plot λ close to 1.0 (0.95-1.05)
- [ ] **Significance threshold:** Multiple testing correction applied
- [ ] **Effect direction:** Top hits make biological sense
- [ ] **Reproducibility:** Random seed set, code versioned

## Model Versions & Tools

| Tool | Language | Best For | License |
|------|----------|----------|---------|
| **tensorQTL** | Python/GPU | eQTL mapping | BSD-3 |
| **GEMMA** | C++ | LMM-GWAS | GPL-3 |
| **PLINK 2** | C++ | GWAS, PCA, QC | GPL-3 |
| **R/qtl2** | R | Classical QTL | GPL-3 |
| **pyQTL** | Python | Visualization | BSD-3 |
| **ADMIXTURE** | C++ | Population structure | BSD-3 |

## Resources

- **tensorQTL:** https://github.com/broadinstitute/tensorqtl
- **GEMMA:** https://github.com/genetics-statistics/GEMMA
- **PLINK:** https://www.cog-genomics.org/plink/2.0/
- **R/qtl2:** https://kbroman.org/qtl2/
- **pyQTL:** https://github.com/broadinstitute/pyqtl
- **GWAS Catalog:** https://www.ebi.ac.uk/gwas/

## Examples

Four fully-worked examples in `examples/`:

| Example | Directory | Demonstrates | Acceptance Criteria |
|---------|-----------|--------------|---------------------|
| **GWAS-LMM** | `examples/gwas-lmm/` | GEMMA GWAS → Manhattan + QQ | 3/3 causal loci recovered, λ=1.02 |
| **eQTL-cis** | `examples/eqtl-cis/` | tensorQTL cis-eQTL | 10 significant eQTLs, GPU <30s |
| **Classical-QTL** | `examples/classical-qtl/` | R/qtl2 LOD scan | 2 QTLs detected, CIs overlap true positions |
| **Pop-structure** | `examples/population-structure/` | PCA + admixture | 3 populations separated, admixture matches structure |

## Validation Commands

```bash
# GWAS regression test
python -c "
import pandas as pd
df = pd.read_csv('examples/gwas-lmm/output/gwas_results.csv')
assert df['lambda_gc'].iloc[0] < 1.1, 'Genomic inflation too high'
assert (df['pval'] < 5e-8).sum() >= 2, 'Not enough significant hits'
print('GWAS validation: PASS')
"

# eQTL regression test
python -c "
import pandas as pd
df = pd.read_csv('examples/eqtl-cis/output/cis_eqtl_results.csv')
assert len(df) >= 5, 'Not enough eQTLs found'
assert df['pval'].min() < 1e-5, 'No significant eQTLs'
print('eQTL validation: PASS')
"
```
