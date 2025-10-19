# Biomni Task Examples

Comprehensive collection of biomedical task examples with code patterns and best practices.

## Table of Contents

1. [Single-Cell RNA-seq Analysis](#single-cell-rna-seq-analysis)
2. [CRISPR Screening](#crispr-screening)
3. [Genomic Analysis (GWAS, Variant Calling)](#genomic-analysis)
4. [Protein Structure and Function](#protein-structure-and-function)
5. [Drug Discovery and ADMET](#drug-discovery-and-admet)
6. [Pathway and Network Analysis](#pathway-and-network-analysis)
7. [Disease Classification](#disease-classification)
8. [Multi-Omics Integration](#multi-omics-integration)
9. [Proteomics Analysis](#proteomics-analysis)
10. [Biomarker Discovery](#biomarker-discovery)

---

## Single-Cell RNA-seq Analysis

### Basic scRNA-seq Pipeline

```python
from biomni.agent import A1

agent = A1(path='./data', llm='claude-sonnet-4-20250514')

agent.go("""
Analyze the 10X Genomics scRNA-seq dataset located at 'data/pbmc_10k.h5ad'.

Workflow:
1. Load the data and perform QC:
   - Filter cells with <200 genes or >5000 genes
   - Filter cells with >10% mitochondrial reads
   - Filter genes expressed in <3 cells

2. Normalize and identify highly variable genes:
   - Use SCTransform or standard log-normalization
   - Identify top 2000 HVGs

3. Dimensionality reduction:
   - PCA (50 components)
   - UMAP for visualization

4. Clustering:
   - Find neighbors (k=10)
   - Leiden clustering with resolution 0.5

5. Visualization:
   - UMAP colored by cluster
   - QC metrics on UMAP

Save processed data as 'results/pbmc_processed.h5ad'
""")
```

### Cell Type Annotation

```python
agent.go("""
Using the processed PBMC data at 'results/pbmc_processed.h5ad':

1. Find marker genes for each cluster:
   - Wilcoxon rank-sum test
   - Log fold change > 0.5
   - Adjusted p-value < 0.01
   - Present in >25% of cluster cells

2. Annotate cell types using markers:
   - T cells: CD3D, CD3E, CD3G
   - B cells: CD19, MS4A1 (CD20)
   - NK cells: GNLY, NKG7, NCAM1
   - Monocytes: CD14, LYZ, CD68
   - Dendritic cells: FCER1A, CD1C

3. Create visualization:
   - UMAP with cell type labels
   - Dotplot of marker genes by cell type
   - Proportion of cell types (bar plot)

4. Save annotated data with cell types
""")
```

### Differential Expression Between Conditions

```python
agent.go("""
Compare gene expression between stimulated and control conditions:

Data: 'data/immune_stim_experiment.h5ad' (contains 'condition' metadata)

Analysis:
1. Subset to T cells only (cell_type == 'T cell')

2. Differential expression between stim vs control:
   - Use pseudobulk approach (aggregate by donor + condition)
   - DESeq2 or edgeR for statistical testing
   - Filter: |log2FC| > 1, padj < 0.05

3. Pathway enrichment on DEGs:
   - Use GO biological processes
   - Use KEGG pathways
   - Run enrichment analysis with gprofiler or enrichr

4. Visualization:
   - Volcano plot of DEGs
   - Heatmap of top 50 DEGs
   - Bar plot of top enriched pathways

5. Export results table with gene symbols, log2FC, p-values, and pathway annotations
""")
```

### Trajectory Analysis

```python
agent.go("""
Perform pseudotime trajectory analysis on hematopoietic differentiation data:

Data: 'data/hematopoiesis.h5ad'

Steps:
1. Subset to progenitor and mature cell types:
   - HSC, MPP, GMP, Monocytes, Neutrophils

2. Run trajectory inference:
   - Use PAGA or Monocle3
   - Set HSC as root cell type

3. Calculate pseudotime for all cells

4. Identify trajectory-associated genes:
   - Genes that change along pseudotime
   - Statistical test with FDR < 0.05
   - Cluster genes by expression pattern (early, middle, late)

5. Visualization:
   - UMAP colored by pseudotime
   - Heatmap of trajectory genes
   - Gene expression along pseudotime for key TFs

6. Functional analysis:
   - GO enrichment for early/middle/late gene clusters
""")
```

### Integration of Multiple Datasets

```python
agent.go("""
Integrate three scRNA-seq datasets from different batches:

Data files:
- 'data/batch1_pbmc.h5ad'
- 'data/batch2_pbmc.h5ad'
- 'data/batch3_pbmc.h5ad'

Integration workflow:
1. Load all three datasets

2. Perform individual QC on each batch:
   - Same filters as standard QC
   - Note batch-specific statistics

3. Integration using Harmony or Scanorama:
   - Concatenate datasets
   - Identify HVGs on combined data
   - Run batch correction
   - Verify batch mixing with LISI score

4. Re-cluster integrated data:
   - Use corrected embeddings
   - Leiden clustering

5. Cell type annotation on integrated data

6. Visualization:
   - UMAP split by batch (before/after correction)
   - UMAP colored by cell type
   - Batch mixing statistics

7. Save integrated dataset
""")
```

---

## CRISPR Screening

### Guide RNA Design

```python
agent.go("""
Design guide RNAs for CRISPR knockout screening of cell cycle genes:

Target genes:
- CDK1, CDK2, CDK4, CDK6
- CCNA2, CCNB1, CCND1, CCNE1
- TP53, RB1, MYC

Requirements:
1. Design 4-6 guides per gene targeting early exons

2. For each guide, evaluate:
   - On-target efficiency score (Doench 2016)
   - Off-target potential (CFD score < 0.3)
   - Avoid common SNPs (1000 Genomes)

3. Add control guides:
   - 100 non-targeting controls
   - 20 positive controls (essential genes)

4. Output:
   - Table with: gene, guide_sequence, PAM, position, on_target_score, off_target_count
   - Sequences in format for oligonucleotide ordering
   - Visual summary of guide distribution per gene

Reference genome: hg38
""")
```

### CRISPR Screen Analysis

```python
agent.go("""
Analyze data from a genome-wide CRISPR knockout screen:

Data: 'data/crispr_screen_counts.csv'
- Columns: guide_id, gene, sample_T0, sample_T15, replicate
- ~80,000 guides targeting ~18,000 genes

Analysis:
1. Quality control:
   - Guide representation (reads per guide)
   - Sample correlation
   - Remove guides with <30 reads in T0

2. Normalize counts:
   - Reads per million (RPM)
   - Log2 fold change (T15 vs T0)

3. Statistical analysis using MAGeCK:
   - Identify significantly depleted/enriched genes
   - FDR < 0.05
   - Rank genes by robust rank aggregation (RRA)

4. Functional analysis:
   - Pathway enrichment of hit genes
   - Known vs novel essential genes
   - Correlation with Cancer Dependency Map

5. Visualization:
   - Scatterplot: log2FC vs -log10(FDR)
   - Heatmap: top 50 depleted genes across replicates
   - Network: PPI network of hit genes

6. Export:
   - Ranked gene list with statistics
   - Enriched pathways table
""")
```

### Pooled Optical Screening Analysis

```python
agent.go("""
Analyze pooled CRISPR screen with imaging readout (e.g., Cell Painting):

Data structure:
- 'data/guide_assignments.csv': cell_id, guide_id, gene
- 'data/morphology_features.csv': cell_id, feature_1...feature_500

Analysis:
1. Feature preprocessing:
   - Remove low-variance features
   - Normalize features (z-score per plate)
   - PCA for dimensionality reduction

2. Associate phenotypes with perturbations:
   - Aggregate cells by guide (mean/median)
   - Calculate morphological distance from controls
   - Statistical test for phenotype change

3. Identify phenotype-altering genes:
   - Mahalanobis distance from control distribution
   - Bonferroni correction for multiple testing
   - Effect size threshold

4. Cluster genes by phenotype similarity:
   - Hierarchical clustering of gene profiles
   - Identify phenotype classes

5. Validation and interpretation:
   - Compare to known gene functions
   - Pathway enrichment per phenotype cluster

6. Visualization:
   - UMAP of all perturbations
   - Heatmap of gene clusters × morphology features
   - Representative images for each cluster
""")
```

---

## Genomic Analysis

### GWAS Analysis

```python
agent.go("""
Perform genome-wide association study for Type 2 Diabetes:

Data:
- 'data/genotypes.bed' (PLINK format, 500K SNPs, 5000 cases, 5000 controls)
- 'data/phenotypes.txt' (sample_id, T2D_status, age, sex, BMI, ancestry_PCs)

Workflow:
1. Quality control:
   - SNP QC: MAF > 0.01, HWE p > 1e-6, genotyping rate > 0.95
   - Sample QC: genotyping rate > 0.95, heterozygosity check
   - Remove related individuals (kinship > 0.125)

2. Association testing:
   - Logistic regression: T2D ~ SNP + age + sex + BMI + PC1-10
   - Genome-wide significance threshold: p < 5e-8
   - Suggestive threshold: p < 1e-5

3. Post-GWAS analysis:
   - LD clumping (r² > 0.1, 500kb window)
   - Annotate lead SNPs with nearby genes (±100kb)
   - Query GWAS Catalog for known associations

4. Functional annotation:
   - Overlap with regulatory elements (ENCODE)
   - eQTL colocalization (GTEx)
   - GWAS prioritization scores (PoPS, ABC)

5. Visualization:
   - Manhattan plot
   - QQ plot
   - Regional association plots for top loci
   - Locus zoom plots

6. Heritability and genetic correlation:
   - SNP heritability (LDSC)
   - Genetic correlation with related traits

Export summary statistics for meta-analysis
""")
```

### Whole Exome Sequencing Analysis

```python
agent.go("""
Analyze whole exome sequencing data for rare disease diagnosis:

Data: Family trio (proband, mother, father)
- 'data/proband.bam'
- 'data/mother.bam'
- 'data/father.bam'

Phenotype: Developmental delay, seizures, intellectual disability

Pipeline:
1. Variant calling:
   - GATK HaplotypeCaller on each sample
   - Joint genotyping across trio
   - VQSR filtering (SNPs and indels separately)

2. Variant annotation:
   - Functional consequence (VEP or ANNOVAR)
   - Population frequencies (gnomAD)
   - Pathogenicity predictions (CADD, REVEL, SpliceAI)
   - Disease databases (ClinVar, OMIM)

3. Inheritance analysis:
   - De novo variants (absent in both parents)
   - Compound heterozygous variants
   - Rare homozygous variants (autozygosity)
   - X-linked variants (if proband is male)

4. Filtering strategy:
   - Population AF < 0.001 (gnomAD)
   - High-quality variants (GQ > 20, DP > 10)
   - Loss-of-function or missense with CADD > 20
   - Match phenotype to gene function (HPO terms)

5. Prioritization:
   - Known disease genes for phenotype
   - De novo in intolerant genes (pLI > 0.9)
   - Protein-truncating variants

6. Report:
   - Top candidate variants with evidence
   - Gene function and disease association
   - Segregation analysis
   - Recommended validation (Sanger sequencing)
   - ACMG variant classification

Save VCF with annotations and prioritized candidate list
""")
```

### Variant Calling from RNA-seq

```python
agent.go("""
Identify expressed variants from RNA-seq data:

Data: Tumor RNA-seq BAM file
- 'data/tumor_RNAseq.bam'
- Reference: hg38

Purpose: Identify expressed somatic mutations for neoantigen prediction

Steps:
1. Pre-processing:
   - Mark duplicates (Picard)
   - Split reads at junctions (GATK SplitNCigarReads)
   - Base quality recalibration

2. Variant calling:
   - GATK HaplotypeCaller (RNA-seq mode)
   - Filter: DP > 10, AF > 0.05

3. Filtering artifacts:
   - Remove common SNPs (gnomAD AF > 0.001)
   - Filter intronic/intergenic variants
   - Remove known RNA editing sites (RADAR database)
   - Panel of normals (if available)

4. Annotation:
   - Functional impact (VEP)
   - Identify non-synonymous variants
   - Predict MHC binding (NetMHCpan for patient HLA type)

5. Prioritize neoantigens:
   - Strong MHC binding (IC50 < 500nM)
   - High expression (TPM > 5)
   - High variant allele frequency

6. Output:
   - Annotated VCF
   - Neoantigen candidates table
   - Peptide sequences for validation

This requires patient HLA typing data
""")
```

---

## Protein Structure and Function

### Protein Structure Prediction and Analysis

```python
agent.go("""
Predict and analyze structure for novel protein sequence:

Sequence (FASTA format):
>Novel_Kinase_Domain
MKLLVVDDDGVADYSKRDGAFMVAYCIEPGDG...

Tasks:
1. Structure prediction:
   - Use AlphaFold2 or ESMFold
   - Generate 5 models, rank by confidence

2. Quality assessment:
   - pLDDT scores (per-residue confidence)
   - pTM score (global confidence)
   - Identify low-confidence regions

3. Domain identification:
   - InterProScan for domain architecture
   - Pfam domain search
   - Identify catalytic residues

4. Functional site prediction:
   - Active site prediction
   - Substrate binding pocket identification
   - Post-translational modification sites

5. Structural alignment:
   - Search for similar structures (PDB)
   - Align to close homologs
   - Identify conserved structural motifs

6. Mutation analysis:
   - Known disease mutations in homologs
   - Predict impact on structure (Rosetta ddG)

7. Visualization and output:
   - PyMOL/Chimera visualization scripts
   - Structural alignment figures
   - Annotated PDB file with functional sites
   - Summary report with predictions
""")
```

### Protein-Protein Interaction Prediction

```python
agent.go("""
Predict and validate protein-protein interactions:

Target protein: BRCA1
Species: Human

Analysis:
1. Literature-based interactions:
   - Query BioGRID, STRING, IntAct databases
   - Extract high-confidence interactors (score > 0.7)

2. Structure-based prediction:
   - Predict BRCA1 structure (if not available)
   - Dock with known interactors (BRCA2, BARD1, etc.)
   - Score interfaces (PISA, PDBePISA)

3. Sequence-based prediction:
   - Coevolution analysis (EVcouplings)
   - Domain-domain interaction prediction
   - Linear motif search (ELM database)

4. Functional analysis of interactors:
   - GO enrichment analysis
   - KEGG pathway membership
   - Tissue/cell type expression patterns

5. Network analysis:
   - Build PPI network
   - Identify network modules
   - Central hub proteins

6. Experimental validation suggestions:
   - Prioritize interactions for validation
   - Suggest Co-IP or Y2H experiments
   - Identify commercially available antibodies

7. Output:
   - Ranked interaction list with evidence
   - PPI network visualization
   - Structural models of key interactions
""")
```

### Protein Engineering Design

```python
agent.go("""
Design improved enzyme variant with enhanced thermostability:

Target enzyme: TEM-1 β-lactamase
Goal: Increase melting temperature by >10°C while maintaining activity

Strategy:
1. Analyze current structure:
   - Load PDB structure (1BTL)
   - Identify flexible regions (B-factors)
   - Find potential disulfide bond sites

2. Computational design:
   - Rosetta design simulations
   - Identify stabilizing mutations (ΔΔG < -1.0 kcal/mol)
   - Avoid active site and substrate binding regions

3. Prioritize mutations:
   - Surface entropy reduction (SER)
   - Disulfide bond introduction
   - Salt bridge formation
   - Hydrophobic core packing

4. Check conservation:
   - Multiple sequence alignment of β-lactamases
   - Avoid highly conserved positions
   - Prefer positions with natural variation

5. Design library:
   - Rank top 20 single mutants
   - Design 5 combinatorial variants (2-3 mutations)
   - Ensure codon optimization for E. coli

6. Validation plan:
   - Expression and purification protocol
   - Thermal shift assay (DSF)
   - Activity assay (nitrocefin)
   - Recommended high-throughput screening

7. Output:
   - Ranked mutation list with predicted ΔΔG
   - Structural visualizations
   - Codon-optimized sequences
   - Cloning primers
   - Experimental validation protocol
""")
```

---

## Drug Discovery and ADMET

### Virtual Screening

```python
agent.go("""
Perform virtual screening for SARS-CoV-2 Mpro inhibitors:

Target: SARS-CoV-2 Main protease (Mpro)
Crystal structure: PDB 6LU7

Compound library: ZINC15 drug-like subset (~100K compounds)
File: 'data/zinc_druglike_100k.smi' (SMILES format)

Workflow:
1. Protein preparation:
   - Remove crystallographic waters (keep catalytic waters)
   - Add hydrogens, optimize H-bond network
   - Define binding site (residues within 5Å of native ligand)

2. Ligand preparation:
   - Generate 3D coordinates from SMILES
   - Enumerate tautomers and protonation states
   - Energy minimization

3. Molecular docking:
   - Dock all compounds (AutoDock Vina or Glide)
   - Generate top 3 poses per compound
   - Score binding affinity

4. Consensus scoring:
   - Combine multiple scoring functions
   - Rank compounds by consensus score

5. ADMET filtering:
   - Lipinski's rule of 5
   - BBB permeability (not needed for this target)
   - hERG liability (pIC50 > 5)
   - CYP450 inhibition prediction
   - Toxicity prediction (Tox21)

6. Visual inspection:
   - Top 100 compounds
   - Check key interactions (His41, Cys145 catalytic dyad)
   - Remove PAINS and frequent hitters

7. Final selection:
   - Top 20 compounds for experimental testing
   - Cluster by scaffold diversity

8. Output:
   - Ranked compound list with scores and ADMET properties
   - Docking poses (mol2 or PDB format)
   - 2D interaction diagrams
   - Purchase availability from vendors
""")
```

### ADMET Property Prediction

```python
agent.go("""
Predict ADMET properties for drug candidate series:

Input: 'data/compound_series.smi' (25 analogs, SMILES format)
Lead scaffold: Novel kinase inhibitor series

Properties to predict:
1. Absorption:
   - Caco-2 permeability
   - Human intestinal absorption (HIA)
   - P-glycoprotein substrate

2. Distribution:
   - Plasma protein binding (% bound)
   - Volume of distribution (VDss)
   - Blood-brain barrier permeability (LogBB)

3. Metabolism:
   - CYP450 substrate (1A2, 2C9, 2C19, 2D6, 3A4)
   - CYP450 inhibition (same isoforms)
   - Sites of metabolism (SOM prediction)

4. Excretion:
   - Clearance estimation
   - Half-life prediction
   - Renal excretion likelihood

5. Toxicity:
   - hERG inhibition (cardiotoxicity)
   - AMES mutagenicity
   - Hepatotoxicity
   - Skin sensitization
   - Rat acute toxicity (LD50)

6. Drug-likeness:
   - Lipinski's Ro5
   - QED score
   - Synthetic accessibility

Analysis:
- Compare all analogs in the series
- Structure-property relationships
- Identify best balanced compound
- Suggest modifications for improvement

Output:
- Comprehensive ADMET table
- Radar plots for each compound
- SAR analysis for each property
- Recommendations for next design iteration
""")
```

### Lead Optimization

```python
agent.go("""
Optimize lead compound balancing potency and selectivity:

Current lead:
- IC50 (target kinase): 50 nM
- IC50 (off-target kinases): 100-500 nM (poor selectivity)
- Microsomal stability: t1/2 = 20 min (too short)
- Solubility: 5 μM (low)

Goal: Maintain potency, improve selectivity (>100x), improve PK properties

Strategy:
1. Analyze current binding mode:
   - Docking to target and off-targets
   - Identify selectivity-determining residues
   - Map interaction hotspots

2. Design focused library:
   - Modifications to improve selectivity:
     * Target residues unique to on-target
     * Avoid conserved kinase regions
   - Modifications to improve solubility:
     * Add polar groups to solvent-exposed regions
     * Replace lipophilic groups
   - Modifications to improve metabolic stability:
     * Block metabolically labile positions
     * Replace metabolically unstable groups

3. Virtual enumeration:
   - Generate ~200 analogs
   - Predict binding affinity (docking)
   - Predict ADMET properties

4. Multi-parameter optimization:
   - Calculate MPO score (potency + selectivity + ADMET)
   - Pareto optimization
   - Select top 20 compounds

5. Clustering and diversity:
   - Ensure structural diversity
   - Test different modification strategies

6. Synthetic feasibility:
   - Retrosynthetic analysis
   - Flag difficult syntheses
   - Prioritize 10 compounds for synthesis

7. Deliverables:
   - Ranked compound designs
   - Predicted properties table
   - Binding mode visualizations
   - Synthetic routes
   - Recommended testing cascade
""")
```

---

## Pathway and Network Analysis

### Pathway Enrichment Analysis

```python
agent.go("""
Perform comprehensive pathway enrichment on differentially expressed genes:

Input: 'data/DEGs.csv'
Columns: gene_symbol, log2FC, padj
Significant DEGs: padj < 0.05, |log2FC| > 1
Total: 450 upregulated, 380 downregulated genes

Background: all detected genes in the experiment (~15,000)

Analysis:
1. GO enrichment (biological processes):
   - Test upregulated and downregulated genes separately
   - Use hypergeometric test
   - FDR correction (Benjamini-Hochberg)
   - Filter: padj < 0.05, fold enrichment > 2

2. KEGG pathway enrichment:
   - Same approach as GO
   - Focus on signaling and metabolic pathways

3. Reactome pathway enrichment:
   - More detailed pathway hierarchy

4. Disease association:
   - DisGeNET disease enrichment
   - Compare to disease gene signatures (MSigDB)

5. Transcription factor enrichment:
   - Predict upstream regulators (ChEA3)
   - ENCODE ChIP-seq enrichment

6. Drug/compound perturbations:
   - L1000 connectivity map
   - Identify drugs that reverse/mimic signature

7. Cross-pathway analysis:
   - Pathway crosstalk
   - Hierarchical clustering of pathways by gene overlap
   - Network visualization of enriched pathways

8. Visualization:
   - Dot plots (GO, KEGG, Reactome)
   - Enrichment map network
   - Chord diagram (genes-pathways)
   - Treemap of hierarchical GO terms

9. Export:
   - All enrichment tables
   - Pathway gene lists
   - Interactive HTML report
""")
```

### Protein-Protein Interaction Network

```python
agent.go("""
Build and analyze PPI network for Alzheimer's disease genes:

Seed genes: Known AD risk genes (APP, PSEN1, PSEN2, APOE, MAPT, etc.)
File: 'data/AD_seed_genes.txt'

Network construction:
1. Build network from seed genes:
   - Query STRING database (confidence > 0.7)
   - Include direct and second-degree interactors
   - Maximum network size: 500 proteins

2. Network enrichment:
   - Add disease associations (DisGeNET)
   - Add tissue expression (GTEx - prioritize brain)
   - Add functional annotations (GO, Reactome)

3. Network analysis:
   - Calculate centrality measures:
     * Degree centrality
     * Betweenness centrality
     * Eigenvector centrality
   - Identify hub proteins
   - Community detection (Louvain algorithm)

4. Module analysis:
   - Functional enrichment per community
   - Identify disease-relevant modules
   - Key bridge proteins between modules

5. Druggability analysis:
   - Identify druggable targets (DGIdb)
   - Known drugs targeting network proteins
   - Clinical trial status

6. Network perturbation:
   - Simulate gene knockout
   - Network robustness analysis
   - Identify critical nodes

7. Visualization:
   - Interactive network (Cytoscape format)
   - Layout by module membership
   - Color by centrality/expression
   - Size by degree

8. Prioritization:
   - Rank proteins by:
     * Network centrality
     * Brain expression
     * Druggability
     * Genetic evidence (GWAS)
   - Top therapeutic targets

Output:
- Network file (graphML, SIF)
- Module membership table
- Prioritized target list
- Druggable targets with existing compounds
""")
```

### Gene Regulatory Network Inference

```python
agent.go("""
Infer gene regulatory network from scRNA-seq data:

Data: 'data/development_timecourse.h5ad'
- Cells from 5 developmental timepoints
- 3000 HVGs quantified

Goal: Identify TF→target relationships during development

Methods:
1. Preprocessing:
   - Select TFs (from TF census list)
   - Select potential target genes (HVGs)
   - Normalize expression

2. GRN inference using multiple methods:
   - GENIE3 (random forest)
   - SCENIC (motif-based)
   - CellOracle (perturbation-based)
   - Pearson/Spearman correlation (baseline)

3. Integrate predictions:
   - Combine scores from multiple methods
   - Weight by motif evidence (JASPAR)
   - Filter low-confidence edges

4. Network refinement:
   - Remove indirect edges (transitive reduction)
   - Validate with ChIP-seq data (if available)
   - Literature validation (TRRUST database)

5. Dynamic network analysis:
   - TF activity per timepoint/cell state
   - Identify stage-specific regulators
   - Find regulatory switches

6. Downstream analysis:
   - Master regulators (high out-degree)
   - Regulatory cascades
   - Feed-forward loops
   - Coherent vs incoherent motifs

7. Experimental validation priorities:
   - Rank TF→target edges for validation
   - Suggest ChIP-seq or CUT&RUN experiments
   - Suggest perturbation experiments (knockout/CRISPRi)

8. Visualization:
   - Full GRN network (Cytoscape)
   - Key TF subnetworks
   - TF activity heatmap across development
   - Sankey diagram of regulatory flow

Output:
- Edge list with confidence scores
- TF activity matrix
- Validated vs novel interactions
- Prioritized validation experiments
""")
```

---

## Disease Classification

### Cancer Type Classification from Gene Expression

```python
agent.go("""
Build multi-class classifier for cancer type prediction:

Data: TCGA pan-cancer RNA-seq data
- Training: 8000 samples across 33 cancer types
- Expression: 'data/tcga_expression.csv' (samples × genes)
- Labels: 'data/tcga_labels.csv' (sample_id, cancer_type)

Task: Classify tumor samples by cancer type

Pipeline:
1. Data preprocessing:
   - Log2(TPM + 1) transformation
   - Remove low-variance genes (variance < 0.1)
   - Z-score normalization

2. Feature selection:
   - Variance filtering (top 5000 genes)
   - Univariate feature selection (ANOVA F-test)
   - Select top 500 features

3. Train-test split:
   - 80% train, 20% test
   - Stratified by cancer type

4. Model training (compare multiple algorithms):
   - Random Forest
   - Gradient Boosting (XGBoost)
   - Neural Network (MLP)
   - Elastic Net logistic regression

5. Model evaluation:
   - Accuracy, precision, recall per class
   - Confusion matrix
   - ROC curves (one-vs-rest)
   - Feature importance ranking

6. Model interpretation:
   - SHAP values for predictions
   - Top predictive genes per cancer type
   - Pathway enrichment of predictive features

7. Clinical validation:
   - Test on independent dataset (if available)
   - Analyze misclassifications
   - Identify hard-to-classify subtypes

8. Deliverables:
   - Trained model (pickle)
   - Performance metrics report
   - Feature importance table
   - Confusion matrix heatmap
   - Prediction script for new samples
""")
```

### Disease Risk Prediction from Multi-Omics

```python
agent.go("""
Develop integrative model predicting cardiovascular disease risk:

Data sources:
1. Genotypes: 'data/genotypes.csv' (500K SNPs, polygenic risk scores)
2. Clinical: 'data/clinical.csv' (age, sex, BMI, blood pressure, cholesterol)
3. Proteomics: 'data/proteomics.csv' (200 plasma proteins, Olink panel)
4. Metabolomics: 'data/metabolomics.csv' (150 metabolites)

Outcome: 10-year CVD incidence (binary)
- Cases: 800
- Controls: 3200

Approach:
1. Data preprocessing:
   - Impute missing values (missForest)
   - Transform skewed features (log/Box-Cox)
   - Normalize each omics layer separately

2. Feature engineering:
   - Calculate PRS from SNP data
   - Interaction terms (age × metabolites, etc.)
   - Metabolite ratios (known CVD markers)

3. Feature selection per omics:
   - Lasso for each data type
   - Select informative features

4. Integration strategies (compare):
   - Early integration: concatenate all features
   - Late integration: separate models, combine predictions
   - Intermediate integration: Multi-omics factor analysis (MOFA)

5. Model development:
   - Logistic regression (interpretable baseline)
   - Random Forest
   - Elastic Net
   - Neural network with omics-specific layers

6. Cross-validation:
   - 5-fold CV, stratified
   - Hyperparameter tuning
   - Calculate confidence intervals

7. Model evaluation:
   - AUC-ROC, AUC-PR
   - Calibration plots
   - Net reclassification improvement (NRI)
   - Compare to clinical models (Framingham, SCORE)

8. Interpretation:
   - Feature importance (permutation importance)
   - SHAP values for individuals
   - Identify most informative omics layer

9. Clinical utility:
   - Decision curve analysis
   - Risk stratification groups
   - Biomarker panel selection

Outputs:
- Model comparison table
- ROC curves all models
- Feature importance per omics
- Reclassification table
- Clinical implementation recommendations
""")
```

---

## Multi-Omics Integration

### Multi-Omics Data Integration

```python
agent.go("""
Integrate transcriptomics, proteomics, and metabolomics data:

Study: Drug response in cancer cell lines
Data:
- RNA-seq: 'data/transcriptomics.csv' (15000 genes × 50 cell lines)
- Proteomics: 'data/proteomics.csv' (3000 proteins × 50 cell lines)
- Metabolomics: 'data/metabolomics.csv' (200 metabolites × 50 cell lines)
- Drug response: 'data/drug_response.csv' (cell line, drug, IC50)

Goal: Identify multi-omics signatures of drug sensitivity

Analysis:
1. Data preprocessing:
   - Match samples across omics layers
   - Filter low-variance features per omics
   - Normalize each omics separately (z-score)

2. Integration methods (compare):

   **Method 1: MOFA (Multi-Omics Factor Analysis)**
   - Identify latent factors capturing variance across omics
   - Determine factor contributions per omics
   - Relate factors to drug response

   **Method 2: DIABLO (sparse PLS-DA)**
   - Supervised integration
   - Maximize covariance between omics and drug response
   - Select features from each omics layer

   **Method 3: Similarity Network Fusion (SNF)**
   - Build patient similarity networks per omics
   - Fuse networks
   - Cluster cell lines by integrated similarity

3. Association with drug response:
   - Correlation of factors/components with IC50
   - Identify drug-sensitive vs resistant groups
   - Multi-omics biomarkers

4. Network analysis:
   - Build multi-layer network:
     * Gene regulatory network (RNA)
     * Protein-protein interactions (proteins)
     * Gene-metabolite associations
   - Integrate layers
   - Find dysregulated pathways

5. Predictive modeling:
   - Train model predicting drug response from multi-omics
   - Compare: using all omics vs individual omics
   - Feature selection across omics

6. Biological interpretation:
   - Map features to pathways
   - Identify mechanism of drug action
   - Suggest combination therapies

7. Validation:
   - Leave-one-out cross-validation
   - Test in independent cell line panel

Outputs:
- Factor loadings per omics (MOFA)
- Multi-omics biomarker signature
- Integrated network visualization
- Predictive model of drug response
- Mechanistic hypotheses
""")
```

---

## Proteomics Analysis

### Label-Free Quantitative Proteomics

```python
agent.go("""
Analyze label-free proteomics data from mass spectrometry:

Study: Comparison of normal vs diseased tissue (n=6 per group)
Data: MaxQuant output
- 'data/proteinGroups.txt' (MaxQuant protein quantification)
- 'data/peptides.txt' (peptide-level data)

Experimental design:
- 6 normal samples
- 6 disease samples
- TMT-labeled, 3 fractions each

Analysis:
1. Data loading and QC:
   - Load proteinGroups.txt
   - Remove contaminants, reverse hits
   - Filter: valid values in ≥50% of samples per group
   - Check sample correlations and outliers
   - PCA for quality assessment

2. Imputation:
   - Impute missing values (MAR vs MNAR approach)
   - Use MinProb for low-abundance missing values
   - Use kNN for random missing values

3. Normalization:
   - Median normalization
   - Or: VSN (variance stabilizing normalization)

4. Differential expression:
   - Two-sample t-test (for each protein)
   - Moderated t-test (limma)
   - Filter: |log2FC| > 0.58 (~1.5-fold), adj.p < 0.05

5. Visualization:
   - Volcano plot
   - Heatmap of significant proteins
   - PCA colored by condition
   - Intensity distributions (before/after normalization)

6. Functional enrichment:
   - GO enrichment (up and down separately)
   - KEGG pathways
   - Reactome pathways
   - STRING PPI network of DEPs

7. PTM analysis (if available):
   - Phosphorylation site analysis
   - Kinase enrichment analysis (KEA3)

8. Orthogonal validation:
   - Compare to RNA-seq data (if available)
   - Protein-RNA correlation
   - Identify discordant genes

9. Biomarker candidates:
   - Rank proteins by fold-change and significance
   - Filter for secreted proteins (potential biomarkers)
   - Check if targetable (druggable)

Outputs:
- Differential abundance table
- QC report with plots
- Enrichment analysis results
- PPI network of DEPs
- Candidate biomarkers list
""")
```

---

## Biomarker Discovery

### Diagnostic Biomarker Discovery

```python
agent.go("""
Discover diagnostic biomarkers for early cancer detection:

Study: Plasma proteomics comparing early-stage cancer vs healthy controls
Data:
- 'data/proteomics.csv' (1000 proteins × 200 samples)
- 'data/metadata.csv' (sample_id, group [cancer/healthy], age, sex)

Groups:
- Early-stage cancer: 100 samples
- Healthy controls: 100 samples

Goal: Identify protein panel for early detection (target AUC > 0.90)

Workflow:
1. Exploratory analysis:
   - PCA, tSNE to visualize separation
   - Univariate differential abundance
   - Volcano plot

2. Feature selection:
   - Rank proteins by:
     * Fold change
     * Statistical significance (t-test, Mann-Whitney)
     * AUC (each protein individually)
   - Select proteins with AUC > 0.70

3. Biomarker panel construction:
   - Correlation analysis (remove redundant markers)
   - Forward selection:
     * Start with best single marker
     * Add markers improving panel performance
     * Stop when no improvement
   - Aim for 5-10 marker panel (practical for assay)

4. Model building:
   - Logistic regression on selected panel
   - Calculate combined risk score
   - Cross-validation (10-fold)

5. Performance evaluation:
   - AUC-ROC, AUC-PR
   - Sensitivity/specificity at different thresholds
   - Clinical decision threshold (e.g., 90% sensitivity)
   - Calibration plot

6. Biological validation:
   - Literature support for cancer association
   - Expression in tumor vs blood
   - Mechanism of release/shedding

7. Clinical utility:
   - Compare to existing biomarkers (CEA, CA19-9, etc.)
   - Cost-effectiveness consideration
   - Assay feasibility (ELISA, MRM, etc.)

8. Independent validation plan:
   - Power calculation for validation cohort
   - Suggested sample size
   - Pre-analytical variables to control

Outputs:
- Ranked protein list with individual performance
- Final biomarker panel
- Logistic regression model
- ROC curves (individual + panel)
- Clinical characteristics table
- Validation study protocol
""")
```

---

## Additional Advanced Examples

### Spatial Transcriptomics Analysis

```python
agent.go("""
Analyze Visium spatial transcriptomics data:

Data: 'data/visium_brain_tumor.h5ad'
- Contains spatial coordinates and gene expression
- Tissue: Brain tumor biopsy

Analysis:
1. Data QC and normalization:
   - Filter low-quality spots (total counts, detected genes)
   - Normalize, log-transform
   - Calculate spatial statistics

2. Spatial clustering:
   - Graph-based clustering considering spatial proximity
   - Identify tumor regions, stroma, necrosis, etc.

3. Spatially variable genes:
   - Test for spatial patterns (Moran's I, SpatialDE)
   - Identify genes with spatial gradients

4. Deconvolution:
   - Estimate cell type composition per spot
   - Use scRNA-seq reference (if available)
   - Methods: Cell2location, RCTD, SPOTlight

5. Niche analysis:
   - Define tissue niches by cell type composition
   - Identify tumor-stroma interface
   - Analyze cell-cell interactions

6. Spatial pathway analysis:
   - Map pathway activity onto tissue
   - Identify spatially localized processes

7. Visualization:
   - Spatial plots colored by cluster, gene expression
   - Cell type composition maps
   - Pathway activity maps

Output:
- Annotated spatial data object
- Spatially variable gene list
- Cell type composition per spot
- Niche definitions and cell-cell interactions
""")
```

---

## Tips for Effective Task Specification

### 1. Be Specific About Data Formats and Locations

✅ Good:
```python
agent.go("Analyze scRNA-seq data in AnnData format at 'data/experiment1.h5ad'")
```

❌ Vague:
```python
agent.go("Analyze my data")
```

### 2. Specify Analysis Parameters

✅ Good:
```python
agent.go("""
Cluster cells using Leiden algorithm with resolution 0.5,
k-neighbors=10, using PCA components 1-30
""")
```

❌ Vague:
```python
agent.go("Cluster the cells")
```

### 3. Request Specific Outputs

✅ Good:
```python
agent.go("""
... and save results as:
- CSV table with statistics
- PNG figures at 300 DPI
- Processed data as AnnData at 'results/processed.h5ad'
""")
```

❌ Vague:
```python
agent.go("... and save the results")
```

### 4. Provide Biological Context

✅ Good:
```python
agent.go("""
This is a drug treatment experiment. Compare vehicle vs treated groups
to identify drug-induced transcriptional changes. Focus on apoptosis and
cell cycle pathways.
""")
```

❌ Vague:
```python
agent.go("Compare the two groups")
```

### 5. Break Complex Analyses into Steps

✅ Good:
```python
# Step 1
agent.go("Load and QC the data, save QC metrics")

# Step 2
agent.go("Based on QC, normalize and find HVGs")

# Step 3
agent.go("Cluster and annotate cell types")
```

❌ Overwhelming:
```python
agent.go("Do a complete scRNA-seq analysis pipeline")
```
