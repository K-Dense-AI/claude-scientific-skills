# Biomni Use Cases and Examples

Comprehensive examples demonstrating biomni across biomedical research domains.

## Table of Contents

1. [CRISPR Screening and Gene Editing](#crispr-screening-and-gene-editing)
2. [Single-Cell RNA-seq Analysis](#single-cell-rna-seq-analysis)
3. [Drug Discovery and ADMET](#drug-discovery-and-admet)
4. [GWAS and Genetic Analysis](#gwas-and-genetic-analysis)
5. [Clinical Genomics and Diagnostics](#clinical-genomics-and-diagnostics)
6. [Protein Structure and Function](#protein-structure-and-function)
7. [Literature and Knowledge Synthesis](#literature-and-knowledge-synthesis)
8. [Multi-Omics Integration](#multi-omics-integration)

---

## CRISPR Screening and Gene Editing

### Example 1: Genome-Wide CRISPR Screen Design

**Task:** Design a CRISPR knockout screen to identify genes regulating autophagy.

```python
from biomni.agent import A1

agent = A1(path='./data', llm='claude-sonnet-4-20250514')

result = agent.go("""
Design a genome-wide CRISPR knockout screen to identify genes regulating
autophagy in HEK293 cells.

Requirements:
1. Generate comprehensive sgRNA library targeting all protein-coding genes
2. Design 4 sgRNAs per gene with optimal on-target and minimal off-target scores
3. Include positive controls (known autophagy regulators: ATG5, BECN1, ULK1)
4. Include negative controls (non-targeting sgRNAs)
5. Prioritize genes based on:
   - Existing autophagy pathway annotations
   - Protein-protein interactions with known autophagy factors
   - Expression levels in HEK293 cells
6. Output sgRNA sequences, scores, and gene prioritization rankings

Provide analysis as Python code and interpret results.
""")

agent.save_conversation_history("autophagy_screen_design.pdf")
```

**Expected Output:**
- sgRNA library with ~80,000 guides (4 per gene × ~20,000 genes)
- On-target and off-target scores for each sgRNA
- Prioritized gene list based on pathway enrichment
- Quality control metrics for library design

### Example 2: CRISPR Off-Target Prediction

```python
result = agent.go("""
Analyze potential off-target effects for this sgRNA sequence:
GCTGAAGATCCAGTTCGATG

Tasks:
1. Identify all genomic locations with ≤3 mismatches
2. Score each potential off-target site
3. Assess likelihood of cleavage at off-target sites
4. Recommend whether sgRNA is suitable for use
5. If unsuitable, suggest alternative sgRNAs for the same gene
""")
```

### Example 3: Screen Hit Analysis

```python
result = agent.go("""
Analyze CRISPR screen results from autophagy phenotype screen.

Input file: screen_results.csv
Columns: sgRNA_ID, gene, log2_fold_change, p_value, FDR

Tasks:
1. Identify significant hits (FDR < 0.05, |LFC| > 1.5)
2. Perform gene ontology enrichment on hit genes
3. Map hits to known autophagy pathways
4. Identify novel candidates not previously linked to autophagy
5. Predict functional relationships between hit genes
6. Generate visualization of hit genes in pathway context
""")
```

---

## Single-Cell RNA-seq Analysis

### Example 1: Cell Type Annotation

**Task:** Analyze single-cell RNA-seq data and annotate cell populations.

```python
agent = A1(path='./data', llm='claude-sonnet-4-20250514')

result = agent.go("""
Analyze single-cell RNA-seq dataset from human PBMC sample.

File: pbmc_data.h5ad (10X Genomics format)

Workflow:
1. Quality control:
   - Filter cells with <200 or >5000 detected genes
   - Remove cells with >20% mitochondrial content
   - Filter genes detected in <3 cells

2. Normalization and preprocessing:
   - Normalize to 10,000 reads per cell
   - Log-transform
   - Identify highly variable genes
   - Scale data

3. Dimensionality reduction:
   - PCA (50 components)
   - UMAP visualization

4. Clustering:
   - Leiden algorithm with resolution=0.8
   - Identify cluster markers (Wilcoxon rank-sum test)

5. Cell type annotation:
   - Annotate clusters using marker genes:
     * T cells (CD3D, CD3E)
     * B cells (CD79A, MS4A1)
     * NK cells (GNLY, NKG7)
     * Monocytes (CD14, LYZ)
     * Dendritic cells (FCER1A, CST3)

6. Generate UMAP plots with annotations and export results
""")

agent.save_conversation_history("pbmc_scrna_analysis.pdf")
```

### Example 2: Differential Expression Analysis

```python
result = agent.go("""
Perform differential expression analysis between conditions in scRNA-seq data.

Data: pbmc_treated_vs_control.h5ad
Conditions: treated (drug X) vs control

Tasks:
1. Identify differentially expressed genes for each cell type
2. Use statistical tests appropriate for scRNA-seq (MAST or Wilcoxon)
3. Apply multiple testing correction (Benjamini-Hochberg)
4. Threshold: |log2FC| > 0.5, adjusted p < 0.05
5. Perform pathway enrichment on DE genes per cell type
6. Identify cell-type-specific drug responses
7. Generate volcano plots and heatmaps
""")
```

### Example 3: Trajectory Analysis

```python
result = agent.go("""
Perform pseudotime trajectory analysis on differentiation dataset.

Data: hematopoiesis_scrna.h5ad
Starting population: Hematopoietic stem cells (HSCs)

Analysis:
1. Subset to hematopoietic lineages
2. Compute diffusion map or PAGA for trajectory inference
3. Order cells along pseudotime
4. Identify genes with dynamic expression along trajectory
5. Cluster genes by expression patterns
6. Map trajectories to known differentiation pathways
7. Visualize key transcription factors driving differentiation
""")
```

---

## Drug Discovery and ADMET

### Example 1: ADMET Property Prediction

**Task:** Predict ADMET properties for drug candidates.

```python
agent = A1(path='./data', llm='claude-sonnet-4-20250514')

result = agent.go("""
Predict ADMET properties for these drug candidates:

Compounds (SMILES format):
1. CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5
2. CN1CCN(CC1)C2=C(C=C3C(=C2)N=CN=C3NC4=CC=C(C=C4)F)OC
3. CC(C)(C)NC(=O)N(CC1=CC=CC=C1)C2CCN(CC2)C(=O)C3=CC4=C(C=C3)OCO4

For each compound, predict:

**Absorption:**
- Caco-2 permeability (cm/s)
- Human intestinal absorption (HIA %)
- Oral bioavailability

**Distribution:**
- Plasma protein binding (%)
- Blood-brain barrier penetration (BBB+/-)
- Volume of distribution (L/kg)

**Metabolism:**
- CYP450 substrate/inhibitor predictions (2D6, 3A4, 2C9, 2C19)
- Metabolic stability (T1/2)

**Excretion:**
- Clearance (mL/min/kg)
- Half-life (hours)

**Toxicity:**
- hERG IC50 (cardiotoxicity risk)
- Hepatotoxicity prediction
- Ames mutagenicity
- LD50 estimates

Provide predictions with confidence scores and flag any red flags.
""")

agent.save_conversation_history("admet_predictions.pdf")
```

### Example 2: Target Identification

```python
result = agent.go("""
Identify potential protein targets for Alzheimer's disease drug development.

Tasks:
1. Query GWAS data for Alzheimer's-associated genes
2. Identify genes with druggable domains (kinases, GPCRs, ion channels, etc.)
3. Check for brain expression patterns
4. Assess disease relevance via literature mining
5. Evaluate existing chemical probe availability
6. Rank targets by:
   - Genetic evidence strength
   - Druggability
   - Lack of existing therapies
7. Suggest target validation experiments
""")
```

### Example 3: Virtual Screening

```python
result = agent.go("""
Perform virtual screening for EGFR kinase inhibitors.

Database: ZINC15 lead-like subset (~6M compounds)
Target: EGFR kinase domain (PDB: 1M17)

Workflow:
1. Prepare protein structure (remove waters, add hydrogens)
2. Define binding pocket (based on erlotinib binding site)
3. Generate pharmacophore model from known EGFR inhibitors
4. Filter ZINC database by:
   - Molecular weight: 200-500 Da
   - LogP: 0-5
   - Lipinski's rule of five
   - Pharmacophore match
5. Dock top 10,000 compounds
6. Score by docking energy and predicted binding affinity
7. Select top 100 for further analysis
8. Predict ADMET properties for top hits
9. Recommend top 10 compounds for experimental validation
""")
```

---

## GWAS and Genetic Analysis

### Example 1: GWAS Summary Statistics Analysis

**Task:** Interpret GWAS results and identify causal genes.

```python
agent = A1(path='./data', llm='claude-sonnet-4-20250514')

result = agent.go("""
Analyze GWAS summary statistics for Type 2 Diabetes.

Input file: t2d_gwas_summary.txt
Columns: CHR, BP, SNP, P, OR, BETA, SE, A1, A2

Analysis steps:
1. Identify genome-wide significant variants (P < 5e-8)
2. Perform LD clumping to identify independent signals
3. Map variants to genes using:
   - Nearest gene
   - eQTL databases (GTEx)
   - Hi-C chromatin interactions
4. Prioritize causal genes using multiple evidence:
   - Fine-mapping scores
   - Coding variant consequences
   - Gene expression in relevant tissues (pancreas, liver, adipose)
   - Pathway enrichment
5. Identify druggable targets among causal genes
6. Compare with known T2D genes and highlight novel associations
7. Generate Manhattan plot, QQ plot, and gene prioritization table
""")

agent.save_conversation_history("t2d_gwas_analysis.pdf")
```

### Example 2: Polygenic Risk Score

```python
result = agent.go("""
Develop and validate polygenic risk score (PRS) for coronary artery disease (CAD).

Training GWAS: CAD_discovery_summary_stats.txt (N=180,000)
Validation cohort: CAD_validation_genotypes.vcf (N=50,000)

Tasks:
1. Select variants for PRS using p-value thresholding (P < 1e-5)
2. Perform LD clumping (r² < 0.1, 500kb window)
3. Calculate PRS weights from GWAS betas
4. Compute PRS for validation cohort individuals
5. Evaluate PRS performance:
   - AUC for CAD case/control discrimination
   - Odds ratios across PRS deciles
   - Compare to traditional risk factors (age, sex, BMI, smoking)
6. Assess PRS calibration and create risk stratification plot
7. Identify high-risk individuals (top 5% PRS)
""")
```

### Example 3: Variant Pathogenicity Prediction

```python
result = agent.go("""
Predict pathogenicity of rare coding variants in candidate disease genes.

Variants (VCF format):
- chr17:41234451:A>G (BRCA1 p.Arg1347Gly)
- chr2:179428448:C>T (TTN p.Trp13579*)
- chr7:117188679:G>A (CFTR p.Gly542Ser)

For each variant, assess:
1. In silico predictions (SIFT, PolyPhen2, CADD, REVEL)
2. Population frequency (gnomAD)
3. Evolutionary conservation (PhyloP, PhastCons)
4. Protein structure impact (using AlphaFold structures)
5. Functional domain location
6. ClinVar annotations (if available)
7. Literature evidence
8. ACMG/AMP classification criteria

Provide pathogenicity classification (benign, likely benign, VUS, likely pathogenic, pathogenic) with supporting evidence.
""")
```

---

## Clinical Genomics and Diagnostics

### Example 1: Rare Disease Diagnosis

**Task:** Diagnose rare genetic disease from whole exome sequencing.

```python
agent = A1(path='./data', llm='claude-sonnet-4-20250514')

result = agent.go("""
Analyze whole exome sequencing (WES) data for rare disease diagnosis.

Patient phenotypes (HPO terms):
- HP:0001250 (Seizures)
- HP:0001249 (Intellectual disability)
- HP:0001263 (Global developmental delay)
- HP:0001252 (Hypotonia)

VCF file: patient_trio.vcf (proband + parents)

Analysis workflow:
1. Variant filtering:
   - Quality filters (QUAL > 30, DP > 10, GQ > 20)
   - Frequency filters (gnomAD AF < 0.01)
   - Functional impact (missense, nonsense, frameshift, splice site)

2. Inheritance pattern analysis:
   - De novo variants
   - Autosomal recessive (compound het, homozygous)
   - X-linked

3. Phenotype-driven prioritization:
   - Match candidate genes to HPO terms
   - Use HPO-gene associations
   - Check gene expression in relevant tissues (brain)

4. Variant pathogenicity assessment:
   - In silico predictions
   - ACMG classification
   - Literature evidence

5. Generate diagnostic report with:
   - Top candidate variants
   - Supporting evidence
   - Functional validation suggestions
   - Genetic counseling recommendations
""")

agent.save_conversation_history("rare_disease_diagnosis.pdf")
```

### Example 2: Cancer Genomics Analysis

```python
result = agent.go("""
Analyze tumor-normal paired sequencing for cancer genomics.

Files:
- tumor_sample.vcf (somatic variants)
- tumor_rnaseq.bam (gene expression)
- tumor_cnv.seg (copy number variants)

Analysis:
1. Identify driver mutations:
   - Known cancer genes (COSMIC, OncoKB)
   - Recurrent hotspot mutations
   - Truncating mutations in tumor suppressors

2. Analyze mutational signatures:
   - Decompose signatures (COSMIC signatures)
   - Identify mutagenic processes

3. Copy number analysis:
   - Identify amplifications and deletions
   - Focal vs. arm-level events
   - Assess oncogene amplifications and TSG deletions

4. Gene expression analysis:
   - Identify outlier gene expression
   - Fusion transcript detection
   - Pathway dysregulation

5. Therapeutic implications:
   - Match alterations to FDA-approved therapies
   - Identify clinical trial opportunities
   - Predict response to targeted therapies

6. Generate precision oncology report
""")
```

### Example 3: Pharmacogenomics

```python
result = agent.go("""
Generate pharmacogenomics report for patient genotype data.

VCF file: patient_pgx.vcf

Analyze variants affecting drug metabolism:

**CYP450 genes:**
- CYP2D6 (affects ~25% of drugs)
- CYP2C19 (clopidogrel, PPIs, antidepressants)
- CYP2C9 (warfarin, NSAIDs)
- CYP3A5 (tacrolimus, immunosuppressants)

**Drug transporter genes:**
- SLCO1B1 (statin myopathy risk)
- ABCB1 (P-glycoprotein)

**Drug targets:**
- VKORC1 (warfarin dosing)
- DPYD (fluoropyrimidine toxicity)
- TPMT (thiopurine toxicity)

For each gene:
1. Determine diplotype (*1/*1, *1/*2, etc.)
2. Assign metabolizer phenotype (PM, IM, NM, RM, UM)
3. Provide dosing recommendations using CPIC/PharmGKB guidelines
4. Flag high-risk drug-gene interactions
5. Suggest alternative medications if needed

Generate patient-friendly report with actionable recommendations.
""")
```

---

## Protein Structure and Function

### Example 1: AlphaFold Structure Analysis

```python
agent = A1(path='./data', llm='claude-sonnet-4-20250514')

result = agent.go("""
Analyze AlphaFold structure prediction for novel protein.

Protein: Hypothetical protein ABC123 (UniProt: Q9XYZ1)

Tasks:
1. Retrieve AlphaFold structure from database
2. Assess prediction quality:
   - pLDDT scores per residue
   - Identify high-confidence regions (pLDDT > 90)
   - Flag low-confidence regions (pLDDT < 50)

3. Structural analysis:
   - Identify domains using structural alignment
   - Predict fold family
   - Identify secondary structure elements

4. Functional prediction:
   - Search for structural homologs in PDB
   - Identify conserved functional sites
   - Predict binding pockets
   - Suggest possible ligands/substrates

5. Variant impact analysis:
   - Map disease-associated variants to structure
   - Predict structural consequences
   - Identify variants affecting binding sites

6. Generate PyMOL visualization scripts highlighting key features
""")

agent.save_conversation_history("alphafold_analysis.pdf")
```

### Example 2: Protein-Protein Interaction Prediction

```python
result = agent.go("""
Predict and analyze protein-protein interactions for autophagy pathway.

Query proteins: ATG5, ATG12, ATG16L1

Analysis:
1. Retrieve known interactions from:
   - STRING database
   - BioGRID
   - IntAct
   - Literature mining

2. Predict novel interactions using:
   - Structural modeling (AlphaFold-Multimer)
   - Coexpression analysis
   - Phylogenetic profiling

3. Analyze interaction interfaces:
   - Identify binding residues
   - Assess interface properties (area, hydrophobicity)
   - Predict binding affinity

4. Functional analysis:
   - Map interactions to autophagy pathway steps
   - Identify regulatory interactions
   - Predict complex stoichiometry

5. Therapeutic implications:
   - Identify druggable interfaces
   - Suggest peptide inhibitors
   - Design disruption strategies

Generate network visualization and interaction details.
""")
```

---

## Literature and Knowledge Synthesis

### Example 1: Systematic Literature Review

```python
agent = A1(path='./data', llm='claude-sonnet-4-20250514')

result = agent.go("""
Perform systematic literature review on CRISPR base editing applications.

Search query: "CRISPR base editing" OR "base editor" OR "CBE" OR "ABE"
Date range: 2016-2025

Tasks:
1. Search PubMed and retrieve relevant abstracts
2. Filter for original research articles
3. Extract key information:
   - Base editor type (CBE, ABE, dual)
   - Target organism/cell type
   - Application (disease model, therapy, crop improvement)
   - Editing efficiency
   - Off-target assessment

4. Categorize applications:
   - Therapeutic applications (by disease)
   - Agricultural applications
   - Basic research

5. Analyze trends:
   - Publications over time
   - Most studied diseases
   - Evolution of base editor technology

6. Synthesize findings:
   - Clinical trial status
   - Remaining challenges
   - Future directions

Generate comprehensive review document with citation statistics.
""")

agent.save_conversation_history("crispr_base_editing_review.pdf")
```

### Example 2: Gene Function Synthesis

```python
result = agent.go("""
Synthesize knowledge about gene function from multiple sources.

Target gene: PARK7 (DJ-1)

Integrate information from:
1. **Genetic databases:**
   - NCBI Gene
   - UniProt
   - OMIM

2. **Expression data:**
   - GTEx tissue expression
   - Human Protein Atlas
   - Single-cell expression atlases

3. **Functional data:**
   - GO annotations
   - KEGG pathways
   - Reactome
   - Protein interactions (STRING)

4. **Disease associations:**
   - ClinVar variants
   - GWAS catalog
   - Disease databases (DisGeNET)

5. **Literature:**
   - PubMed abstracts
   - Key mechanistic studies
   - Review articles

Synthesize into comprehensive gene report:
- Molecular function
- Biological processes
- Cellular localization
- Tissue distribution
- Disease associations
- Known drug targets/inhibitors
- Unresolved questions

Generate structured summary suitable for research planning.
""")
```

---

## Multi-Omics Integration

### Example 1: Multi-Omics Disease Analysis

```python
agent = A1(path='./data', llm='claude-sonnet-4-20250514')

result = agent.go("""
Integrate multi-omics data to understand disease mechanism.

Disease: Alzheimer's disease
Data types:
- Genomics: GWAS summary statistics (gwas_ad.txt)
- Transcriptomics: Brain RNA-seq (controls vs AD, rnaseq_data.csv)
- Proteomics: CSF proteomics (proteomics_csf.csv)
- Metabolomics: Plasma metabolomics (metabolomics_plasma.csv)
- Epigenomics: Brain methylation array (methylation_data.csv)

Integration workflow:
1. Analyze each omics layer independently:
   - Identify significantly altered features
   - Perform pathway enrichment

2. Cross-omics correlation:
   - Correlate gene expression with protein levels
   - Link genetic variants to expression (eQTL)
   - Associate methylation with gene expression
   - Connect proteins to metabolites

3. Network analysis:
   - Build multi-omics network
   - Identify key hub genes/proteins
   - Detect disease modules

4. Causal inference:
   - Prioritize drivers vs. consequences
   - Identify therapeutic targets
   - Predict drug mechanisms

5. Generate integrative model of AD pathogenesis

Provide visualization and therapeutic target recommendations.
""")

agent.save_conversation_history("ad_multiomics_analysis.pdf")
```

### Example 2: Systems Biology Modeling

```python
result = agent.go("""
Build systems biology model of metabolic pathway.

Pathway: Glycolysis
Data sources:
- Enzyme kinetics (BRENDA database)
- Metabolite concentrations (literature)
- Gene expression (tissue-specific, GTEx)
- Flux measurements (C13 labeling studies)

Modeling tasks:
1. Construct pathway model:
   - Define reactions and stoichiometry
   - Parameterize enzyme kinetics (Km, Vmax, Ki)
   - Set initial metabolite concentrations

2. Simulate pathway dynamics:
   - Steady-state analysis
   - Time-course simulations
   - Sensitivity analysis

3. Constraint-based modeling:
   - Flux balance analysis (FBA)
   - Identify bottleneck reactions
   - Predict metabolic engineering strategies

4. Integrate with gene expression:
   - Tissue-specific model predictions
   - Disease vs. normal comparisons

5. Therapeutic predictions:
   - Enzyme inhibition effects
   - Metabolic rescue strategies
   - Drug target identification

Generate model in SBML format and simulation results.
""")
```

---

## Best Practices for Task Formulation

### 1. Be Specific and Detailed

**Poor:**
```python
agent.go("Analyze this RNA-seq data")
```

**Good:**
```python
agent.go("""
Analyze bulk RNA-seq data from cancer vs. normal samples.

Files: cancer_rnaseq.csv (TPM values, 50 cancer, 50 normal)

Tasks:
1. Differential expression (DESeq2, padj < 0.05, |log2FC| > 1)
2. Pathway enrichment (KEGG, Reactome)
3. Generate volcano plot and top DE gene heatmap
""")
```

### 2. Include File Paths and Formats

Always specify:
- Exact file paths
- File formats (VCF, BAM, CSV, H5AD, etc.)
- Data structure (columns, sample IDs)

### 3. Set Clear Success Criteria

Define thresholds and cutoffs:
- Statistical significance (P < 0.05, FDR < 0.1)
- Fold change thresholds
- Quality filters
- Expected outputs

### 4. Request Visualizations

Explicitly ask for plots:
- Volcano plots, MA plots
- Heatmaps, PCA plots
- Network diagrams
- Manhattan plots

### 5. Specify Biological Context

Include:
- Organism (human, mouse, etc.)
- Tissue/cell type
- Disease/condition
- Treatment details

### 6. Request Interpretations

Ask agent to:
- Interpret biological significance
- Suggest follow-up experiments
- Identify limitations
- Provide literature context

---

## Common Patterns

### Data Quality Control

```python
"""
Before analysis, perform quality control:
1. Check for missing values
2. Assess data distributions
3. Identify outliers
4. Generate QC report
Only proceed with analysis if data passes QC.
"""
```

### Iterative Refinement

```python
"""
Perform analysis in stages:
1. Initial exploratory analysis
2. Based on results, refine parameters
3. Focus on interesting findings
4. Generate final report

Show intermediate results for each stage.
"""
```

### Reproducibility

```python
"""
Ensure reproducibility:
1. Set random seeds where applicable
2. Log all parameters used
3. Save intermediate files
4. Export environment info (package versions)
5. Generate methods section for paper
"""
```

These examples demonstrate the breadth of biomedical tasks biomni can handle. Adapt the patterns to your specific research questions, and always include sufficient detail for the agent to execute autonomously.
