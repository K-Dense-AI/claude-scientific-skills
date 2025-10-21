# Claude Scientific Skills

A comprehensive collection of ready-to-use scientific skills for Claude, curated by the K-Dense team. These skills enable Claude to work with specialized scientific libraries and databases across bioinformatics, cheminformatics, machine learning, materials science, and data analysis. Using these set of skills with Claude Code allows you to create an 'AI Scientist' on your desktop! If you want substantially more advanced capabilties, compute infrastructure and enterprise ready offering check out https://k-dense.ai/.

## Getting Started

### Claude Code
You can register this repository as a Claude Code Plugin marketplace by running the following command in Claude Code:

```
/plugin marketplace add K-Dense-AI/claude-scientific-skills
```

Then, to install a specific set of skills:

1. Select Browse and install plugins
2. Select claude-scientific-skills
3. Select scientific-databases, scientific-packages, or scientific-thinking (includes document processing)
4. Select Install now

After installing the plugin, you can use the skill by just mentioning it. Additionally, in most case, Claude Code will figure out what to use based on the task.

## Available Skills

### Scientific Databases

- **AlphaFold DB** - AI-predicted protein structure database with 200M+ predictions, confidence metrics (pLDDT, PAE), and Google Cloud bulk access
- **ChEMBL** - Bioactive molecule database with drug-like properties (2M+ compounds, 19M+ activities, 13K+ targets)
- **ClinPGx** - Clinical pharmacogenomics database (successor to PharmGKB) providing gene-drug interactions, CPIC clinical guidelines, allele functions, drug labels, and pharmacogenomic annotations for precision medicine and personalized pharmacotherapy (consolidates PharmGKB, CPIC, and PharmCAT resources)
- **ClinVar** - NCBI's public archive of genomic variants and their clinical significance with standardized classifications (pathogenic, benign, VUS), E-utilities API access, and bulk FTP downloads for variant interpretation and precision medicine research
- **COSMIC** - Catalogue of Somatic Mutations in Cancer, the world's largest database of somatic cancer mutations (millions of mutations across thousands of cancer types, Cancer Gene Census, mutational signatures, structural variants, and drug resistance data)
- **ENA (European Nucleotide Archive)** - Comprehensive public repository for nucleotide sequence data and metadata with REST APIs for accessing sequences, assemblies, samples, studies, and reads; supports advanced search, taxonomy lookups, and bulk downloads via FTP/Aspera (rate limit: 50 req/sec)
- **Ensembl** - Genome browser and bioinformatics database providing genomic annotations, sequences, variants, and comparative genomics data for 250+ vertebrate species (Release 115, 2025) with comprehensive REST API for gene lookups, sequence retrieval, variant effect prediction (VEP), ortholog finding, assembly mapping (GRCh37/GRCh38), and region analysis
- **GEO (Gene Expression Omnibus)** - High-throughput gene expression and functional genomics data repository (264K+ studies, 8M+ samples) with microarray, RNA-seq, and expression profile access
- **GWAS Catalog** - NHGRI-EBI catalog of published genome-wide association studies with curated SNP-trait associations (thousands of studies, genome-wide significant associations p≤5×10⁻⁸), full summary statistics, REST API access for variant/trait/gene queries, and FTP downloads for genetic epidemiology and precision medicine research
- **HMDB (Human Metabolome Database)** - Comprehensive metabolomics resource with 220K+ metabolite entries, detailed chemical/biological data, concentration ranges, disease associations, pathways, and spectral data for metabolite identification and biomarker discovery
- **KEGG** - Kyoto Encyclopedia of Genes and Genomes for biological pathway analysis, gene-to-pathway mapping, compound searches, and molecular interaction networks (pathway enrichment, metabolic pathways, gene annotations, drug-drug interactions, ID conversion)
- **Metabolomics Workbench** - NIH Common Fund metabolomics data repository with 4,200+ processed studies, standardized nomenclature (RefMet), mass spectrometry searches, and comprehensive REST API for accessing metabolite structures, study metadata, experimental results, and gene/protein-metabolite associations
- **NCBI Gene** - Work with NCBI Gene database to search, retrieve, and analyze gene information including nomenclature, sequences, variations, phenotypes, and pathways using E-utilities and Datasets API
- **Protein Data Bank (PDB)** - Access 3D structural data of proteins, nucleic acids, and biological macromolecules (200K+ structures) with search, retrieval, and analysis capabilities
- **PubChem** - Access chemical compound data from the world's largest free chemical database (110M+ compounds, 270M+ bioactivities)
- **PubMed** - Access to PubMed literature database with advanced search capabilities
- **Reactome** - Curated pathway database for biological processes and molecular interactions (2,825+ human pathways, 16K+ reactions, 11K+ proteins) with pathway enrichment analysis, expression data analysis, and species comparison using Content Service and Analysis Service APIs
- **STRING** - Protein-protein interaction network database (5000+ genomes, 59.3M proteins, 20B+ interactions) with functional enrichment analysis, interaction partner discovery, and network visualization from experimental data, computational prediction, and text-mining
- **UniProt** - Universal Protein Resource for protein sequences, annotations, and functional information (UniProtKB/Swiss-Prot reviewed entries, TrEMBL unreviewed entries) with REST API access for search, retrieval, ID mapping, and batch operations across 200+ databases
- **ZINC** - Free database of commercially-available compounds for virtual screening and drug discovery (230M+ purchasable compounds in ready-to-dock 3D formats)

### Scientific Packages

**Bioinformatics & Genomics:**
- **AnnData** - Annotated data matrices for single-cell genomics and h5ad files
- **Arboreto** - Gene regulatory network inference using GRNBoost2 and GENIE3
- **BioPython** - Sequence manipulation, NCBI database access, BLAST searches, alignments, and phylogenetics
- **BioServices** - Programmatic access to 40+ biological web services (KEGG, UniProt, ChEBI, ChEMBL)
- **Cellxgene Census** - Query and analyze large-scale single-cell RNA-seq data
- **gget** - Efficient genomic database queries (Ensembl, UniProt, NCBI, PDB, COSMIC)
- **pysam** - Read, write, and manipulate genomic data files (SAM/BAM/CRAM alignments, VCF/BCF variants, FASTA/FASTQ sequences) with pileup analysis, coverage calculations, and bioinformatics workflows
- **PyDESeq2** - Differential gene expression analysis for bulk RNA-seq data
- **Scanpy** - Single-cell RNA-seq analysis with clustering, marker genes, and UMAP/t-SNE visualization

**Cheminformatics & Drug Discovery:**
- **Datamol** - Molecular manipulation and featurization with enhanced RDKit workflows
- **DeepChem** - Molecular machine learning, graph neural networks, and MoleculeNet benchmarks
- **DiffDock** - Diffusion-based molecular docking for protein-ligand binding prediction
- **MedChem** - Medicinal chemistry analysis, ADMET prediction, and drug-likeness assessment
- **Molfeat** - 100+ molecular featurizers including fingerprints, descriptors, and pretrained models
- **PyTDC** - Therapeutics Data Commons for drug discovery datasets and benchmarks
- **RDKit** - Cheminformatics toolkit for molecular I/O, descriptors, fingerprints, and SMARTS

**Proteomics & Mass Spectrometry:**
- **matchms** - Processing and similarity matching of mass spectrometry data with 40+ filters, spectral library matching (Cosine, Modified Cosine, Neutral Losses), metadata harmonization, molecular fingerprint comparison, and support for multiple file formats (MGF, MSP, mzML, JSON)
- **pyOpenMS** - Comprehensive mass spectrometry data analysis for proteomics and metabolomics (LC-MS/MS processing, peptide identification, feature detection, quantification, chemical calculations, and integration with search engines like Comet, Mascot, MSGF+)

**Machine Learning & Deep Learning:**
- **PyMC** - Bayesian statistical modeling and probabilistic programming
- **PyMOO** - Multi-objective optimization with evolutionary algorithms
- **PyTorch Lightning** - Structured PyTorch training with automatic optimization
- **scikit-learn** - Machine learning algorithms, preprocessing, and model selection
- **Torch Geometric** - Graph Neural Networks for molecular and geometric data
- **Transformers** - Hugging Face transformers for NLU, image classification, and generation
- **UMAP-learn** - Dimensionality reduction and manifold learning

**Materials Science & Chemistry:**
- **Astropy** - Astronomy and astrophysics (coordinates, cosmology, FITS files)
- **COBRApy** - Constraint-based metabolic modeling and flux balance analysis
- **Pymatgen** - Materials structure analysis, phase diagrams, and electronic structure

**Data Analysis & Visualization:**
- **Dask** - Parallel computing for larger-than-memory datasets with distributed DataFrames, Arrays, Bags, and Futures
- **Matplotlib** - Publication-quality plotting and visualization
- **Polars** - High-performance DataFrame operations with lazy evaluation
- **Seaborn** - Statistical data visualization
- **ReportLab** - Programmatic PDF generation for reports and documents

**Phylogenetics & Trees:**
- **ETE Toolkit** - Phylogenetic tree manipulation, visualization, and analysis

**Genomics Tools:**
- **deepTools** - NGS data analysis (ChIP-seq, RNA-seq, ATAC-seq) with BAM/bigWig files
- **FlowIO** - Flow Cytometry Standard (FCS) file reading and manipulation
- **scikit-bio** - Bioinformatics sequence analysis and diversity metrics
- **Zarr** - Chunked, compressed N-dimensional array storage

**Multi-omics & Integration:**
- **BIOMNI** - Multi-omics data integration with LLM-powered analysis

### Scientific Thinking & Analysis

**Analysis & Methodology:**
- **Exploratory Data Analysis** - Comprehensive EDA toolkit with automated statistics, visualizations, and insights for any tabular dataset
- **Hypothesis Generation** - Structured frameworks for generating and evaluating scientific hypotheses
- **Peer Review** - Comprehensive toolkit for conducting high-quality scientific peer review with structured evaluation of methodology, statistics, reproducibility, ethics, and presentation across all scientific disciplines
- **Scientific Brainstorming** - Conversational brainstorming partner for generating novel research ideas, exploring connections, challenging assumptions, and developing creative approaches through structured ideation workflows
- **Scientific Critical Thinking** - Tools and approaches for rigorous scientific reasoning and evaluation
- **Scientific Visualization** - Best practices and templates for creating publication-quality scientific figures
- **Statistical Analysis** - Comprehensive statistical testing, power analysis, and experimental design

**Document Processing:**
- **DOCX** - Comprehensive document creation, editing, and analysis with support for tracked changes, comments, formatting preservation, and text extraction
- **PDF** - PDF manipulation toolkit for extracting text and tables, creating new PDFs, merging/splitting documents, and handling forms
- **PPTX** - Presentation creation, editing, and analysis with support for layouts, comments, and speaker notes
- **XLSX** - Spreadsheet creation, editing, and analysis with support for formulas, formatting, data analysis, and visualization

## TODO: Future Scientific Capabilities

### Scientific Databases
- **ArrayExpress** - EMBL-EBI gene expression database with functional genomics experiments
- **BioGRID** - Biological General Repository for Interaction Datasets (protein, genetic, and chemical interactions)
- **DAVID** - Database for Annotation, Visualization and Integrated Discovery for functional enrichment analysis
- **dbSNP** - NCBI's database of single nucleotide polymorphisms and short genetic variations
- **GenBank** - NIH genetic sequence database (part of NCBI but with specific access patterns)
- **InterPro** - Protein sequence analysis and classification with functional annotations
- **MetaboLights** - EMBL-EBI metabolomics database with experimental data and metadata
- **OMIM** - Online Mendelian Inheritance in Man for genetic disorders and genes
- **Pfam** - Protein families database with multiple sequence alignments and HMMs
- **RefSeq** - NCBI's non-redundant reference sequence database
- **UCSC Genome Browser** - Genomic data visualization and custom track integration
- **WikiPathways** - Community-curated biological pathway database

### Bioinformatics & Genomics
- **pybedtools** - Wrapper for BEDTools genome arithmetic operations
- **mygene** - Python client for MyGene.Info gene query service
- **pyensembl** - Python interface to Ensembl reference genome metadata
- **nglview** - IPython/Jupyter widget for molecular visualization
- **pyvcf** - Variant Call Format (VCF) file parser
- **pyfaidx** - Efficient FASTA file indexing and retrieval
- **kipoiseq** - Genomic sequence data loading for ML models
- **genomepy** - Download and manage genome assemblies
- **MACS2/3** - Peak calling for ChIP-seq data

### Cheminformatics & Drug Discovery
- **Open Babel** - Chemical file format conversion and molecular mechanics
- **ChemPy** - Chemistry and thermodynamics calculations
- **Psi4** - Quantum chemistry software for ab initio calculations
- **pmapper** - Pharmacophore modeling and fingerprinting
- **ODDT** - Open Drug Discovery Toolkit for structure-based drug design
- **ProLIF** - Protein-ligand interaction fingerprints
- **Mordred** - Molecular descriptor calculator (1800+ descriptors)
- **ProteinMPNN** - Deep learning for protein sequence design
- **ESM** - Evolutionary Scale Modeling for protein language models
- **OpenMM** - Molecular dynamics simulation toolkit

### Proteomics & Mass Spectrometry
- **pyteomics** - Mass spectrometry data analysis
- **MSstats** - Statistical analysis of quantitative proteomics

### Systems Biology & Networks
- **NetworkX** - Complex network analysis and graph algorithms
- **igraph** - Fast network analysis library
- **PyBioNetFit** - Biological network modeling and fitting
- **PINT** - Pathway integration analysis
- **GEMEditor** - Graphical tool for genome-scale metabolic models

### Structural Biology
- **MDAnalysis** - Molecular dynamics trajectory analysis
- **ProDy** - Protein dynamics and structure analysis
- **PyMOL** - Molecular visualization scripting
- **Chimera/ChimeraX** - UCSF molecular visualization
- **FreeSASA** - Solvent accessible surface area calculations
- **DSSP** - Secondary structure assignment

### Machine Learning for Science
- **DGL-LifeSci** - Deep Graph Library for life sciences
- **ChemBERTa** - Transformer models for chemistry
- **TorchDrug** - PyTorch library for drug discovery
- **GraNNField** - Graph neural networks for force fields
- **SchNet/DimeNet** - Continuous-filter convolutional networks for molecules
- **MoleculeNet** - Benchmark datasets for molecular machine learning
- **TorchMD** - Molecular dynamics with PyTorch
- **jax-md** - Differentiable molecular dynamics in JAX

### Imaging & Microscopy
- **scikit-image** - Image processing algorithms
- **CellProfiler** - Cell image analysis
- **Napari** - Multi-dimensional image viewer
- **Fiji/ImageJ** - Image processing scripting
- **StarDist** - Cell/nucleus detection with deep learning
- **Cellpose** - Generalist cell segmentation

### Phylogenetics & Evolution
- **DendroPy** - Phylogenetic computing library
- **PyCogent** - Comparative genomics toolkit
- **TreeTime** - Phylodynamic analysis and molecular clock inference

### Metabolomics
- **PyCytoData** - Cytometry data processing
- **MS-DIAL** - Data-independent MS/MS deconvolution
- **XCMS** - LC/MS and GC/MS data processing

### Climate & Environmental Science
- **xarray** - N-dimensional labeled arrays and datasets
- **Iris** - Climate and weather data analysis
- **MetPy** - Meteorological data toolkit
- **climlab** - Climate modeling and analysis

### Statistics & Experimental Design
- **statsmodels** - Statistical models and hypothesis testing
- **pingouin** - Statistical tests with clear output
- **PyDOE2** - Design of experiments
- **scipy.stats** - Statistical functions and distributions

### Data Management & Processing
- **Parquet** - Columnar storage format for big data
- **DuckDB** - Analytical SQL database
- **SQLAlchemy** - SQL toolkit and ORM

### Visualization
- **Plotly** - Interactive graphing library
- **Bokeh** - Interactive visualization for web browsers
- **Altair** - Declarative statistical visualization
- **PyVista** - 3D plotting and mesh analysis
