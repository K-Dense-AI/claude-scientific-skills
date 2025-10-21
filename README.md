# Claude Scientific Skills

A comprehensive collection of ready-to-use scientific skills for Claude, curated by the K-Dense team. These skills enable Claude to work with specialized scientific libraries and databases across bioinformatics, cheminformatics, machine learning, materials science, and data analysis. Using these set of skills with Claude Code allows you to create an 'AI Scientist' on your desktop! If you want substantially more advanced capabilties, compute infrastructure and enterprise ready offering check out https://k-dense.ai/.

This repository provides access to **21 scientific databases**, **44 scientific packages**, **4 scientific integrations**, and **103 unique workflows** covering a wide range of scientific computing tasks.

## Getting Started

### Claude Code
You can register this repository as a Claude Code Plugin marketplace by running the following command in Claude Code:

```
/plugin marketplace add K-Dense-AI/claude-scientific-skills
```

Then, to install a specific set of skills:

1. Select Browse and install plugins
2. Select claude-scientific-skills
3. Select scientific-databases, scientific-packages, scientific-thinking (includes document processing), or scientific-integrations
4. Select Install now

After installing the plugin, you can use the skill by just mentioning it. Additionally, in most case, Claude Code will figure out what to use based on the task.

## Available Skills

### Scientific Databases

- **AlphaFold DB** - AI-predicted protein structure database with 200M+ predictions, confidence metrics (pLDDT, PAE), and Google Cloud bulk access
- **ChEMBL** - Bioactive molecule database with drug-like properties (2M+ compounds, 19M+ activities, 13K+ targets)
- **ClinPGx** - Clinical pharmacogenomics database (successor to PharmGKB) providing gene-drug interactions, CPIC clinical guidelines, allele functions, drug labels, and pharmacogenomic annotations for precision medicine and personalized pharmacotherapy (consolidates PharmGKB, CPIC, and PharmCAT resources)
- **ClinVar** - NCBI's public archive of genomic variants and their clinical significance with standardized classifications (pathogenic, benign, VUS), E-utilities API access, and bulk FTP downloads for variant interpretation and precision medicine research
- **ClinicalTrials.gov** - Comprehensive registry of clinical studies conducted worldwide (maintained by U.S. National Library of Medicine) with API v2 access for searching trials by condition, intervention, location, sponsor, study status, and phase; retrieve detailed trial information including eligibility criteria, outcomes, contacts, and locations; export to CSV/JSON formats for analysis (public API, no authentication required, ~50 req/min rate limit)
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
- **statsmodels** - Statistical modeling and econometrics (OLS, GLM, logit/probit, ARIMA, time series forecasting, hypothesis testing, diagnostics)
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

### Scientific Integrations

**Laboratory Information Management Systems (LIMS) & R&D Platforms:**
- **Benchling Integration** - Toolkit for integrating with Benchling's R&D platform, providing programmatic access to laboratory data management including registry entities (DNA sequences, proteins), inventory systems (samples, containers, locations), electronic lab notebooks (entries, protocols), workflows (tasks, automation), and data exports using Python SDK and REST API

**Cloud Platforms for Genomics & Biomedical Data:**
- **DNAnexus Integration** - Comprehensive toolkit for working with the DNAnexus cloud platform for genomics and biomedical data analysis. Covers building and deploying apps/applets (Python/Bash), managing data objects (files, records, databases), running analyses and workflows, using the dxpy Python SDK, and configuring app metadata and dependencies (dxapp.json setup, system packages, Docker, assets). Enables processing of FASTQ/BAM/VCF files, bioinformatics pipelines, job execution, workflow orchestration, and platform operations including project management and permissions

**Laboratory Automation:**
- **Opentrons Integration** - Toolkit for creating, editing, and debugging Opentrons Python Protocol API v2 protocols for laboratory automation using Flex and OT-2 robots. Enables automated liquid handling, pipetting workflows, hardware module control (thermocycler, temperature, magnetic, heater-shaker, absorbance plate reader), labware management, and complex protocol development for biological and chemical experiments

**Electronic Lab Notebooks (ELN):**
- **LabArchives Integration** - Toolkit for interacting with LabArchives Electronic Lab Notebook (ELN) REST API. Provides programmatic access to notebooks (backup, retrieval, management), entries (creation, comments, attachments), user authentication, site reports and analytics, and third-party integrations (Protocols.io, GraphPad Prism, SnapGene, Geneious, Jupyter, REDCap). Includes Python scripts for configuration setup, notebook operations, and entry management. Supports multi-regional API endpoints (US, UK, Australia) and OAuth authentication

**Microscopy & Bio-image Data:**
- **OMERO Integration** - Toolkit for interacting with OMERO microscopy data management systems using Python. Provides comprehensive access to microscopy images stored in OMERO servers, including dataset and screening data retrieval, pixel data analysis, annotation and metadata management, regions of interest (ROIs) creation and analysis, batch processing, OMERO.scripts development, and OMERO.tables for structured data storage. Essential for researchers working with high-content screening data, multi-dimensional microscopy datasets, or collaborative image repositories

## TODO: Future Scientific Capabilities

### Scientific Integrations
- **PerkinElmer Signals** - Scientific data management and ELN platform integration
- **CDD Vault** - Collaborative Drug Discovery platform integration for chemical registration and bioassay data
- **Geneious** - Molecular biology and NGS analysis software integration
- **SnapGene** - Molecular cloning and DNA visualization platform integration
- **Synthego ICE** - CRISPR editing analysis platform integration
- **TeselaGen** - Synthetic biology design and automation platform integration
- **Galaxy** - Web-based bioinformatics workflow platform integration
- **Nextflow/nf-core** - Workflow management system integration for reproducible pipelines
- **Seven Bridges** - Genomics analysis platform and workspace integration
- **DNAnexus** - Cloud-based genome sequencing analysis platform integration
- **BaseSpace** - Illumina genomics data analysis and management platform integration

### Scientific Databases
- **BioGRID** - Biological General Repository for Interaction Datasets (protein, genetic, and chemical interactions)
- **dbSNP** - NCBI's database of single nucleotide polymorphisms and short genetic variations
- **InterPro** - Protein sequence analysis and classification with functional annotations
- **OMIM** - Online Mendelian Inheritance in Man for genetic disorders and genes
- **Pfam** - Protein families database with multiple sequence alignments and HMMs
- **RefSeq** - NCBI's non-redundant reference sequence database
- **UCSC Genome Browser** - Genomic data visualization and custom track integration
- **WikiPathways** - Community-curated biological pathway database
- **MetaboLights** - EMBL-EBI metabolomics database with experimental data and metadata

### Bioinformatics & Genomics
- **pybedtools** - Wrapper for BEDTools genome arithmetic operations
- **mygene** - Python client for MyGene.Info gene query service
- **nglview** - IPython/Jupyter widget for molecular visualization
- **pyfaidx** - Efficient FASTA file indexing and retrieval
- **MACS2/3** - Peak calling for ChIP-seq data

### Cheminformatics & Drug Discovery
- **Open Babel** - Chemical file format conversion and molecular mechanics
- **Psi4** - Quantum chemistry software for ab initio calculations
- **ProteinMPNN** - Deep learning for protein sequence design
- **ESM (Evolutionary Scale Modeling)** - Protein language models for structure and function prediction
- **OpenMM** - Molecular dynamics simulation toolkit

### Proteomics & Mass Spectrometry
- **pyteomics** - Mass spectrometry data analysis and peptide/protein identification

### Systems Biology & Networks
- **NetworkX** - Complex network analysis and graph algorithms
- **igraph** - Fast network analysis library with efficient algorithms

### Structural Biology
- **MDAnalysis** - Molecular dynamics trajectory analysis
- **ProDy** - Protein dynamics and structure analysis
- **PyMOL** - Molecular visualization scripting

### Machine Learning for Science
- **DGL-LifeSci** - Deep Graph Library for life sciences
- **ChemBERTa** - Transformer models for chemistry
- **TorchDrug** - PyTorch library for drug discovery
- **SchNet/DimeNet** - Continuous-filter convolutional networks for molecules

### Imaging & Microscopy
- **scikit-image** - Image processing algorithms
- **Napari** - Multi-dimensional image viewer
- **CellProfiler** - Cell image analysis
- **Cellpose** - Generalist cell segmentation
- **StarDist** - Cell/nucleus detection with deep learning

### Phylogenetics & Evolution
- **DendroPy** - Phylogenetic computing library

### Climate & Environmental Science
- **xarray** - N-dimensional labeled arrays and datasets for scientific computing

### Statistics & Experimental Design
- **pingouin** - Statistical tests with clear output and effect sizes
- **scipy.stats** - Statistical functions and distributions

### Data Management & Processing
- **DuckDB** - Analytical SQL database for in-process analytics
- **Parquet** - Columnar storage format for big data

### Visualization
- **Plotly** - Interactive graphing library for web-based visualizations
- **Altair** - Declarative statistical visualization
- **PyVista** - 3D plotting and mesh analysis

## Contributing

We welcome contributions to expand and improve this scientific skills repository! There are several ways you can contribute:

### Improving Existing Skills
- Enhance documentation with more examples and use cases
- Add new workflows and reference materials
- Improve code examples and scripts
- Fix bugs or update outdated information

### How to Contribute
1. Fork the repository
2. Create a feature branch for your contribution
3. Follow the existing directory structure and documentation patterns
4. Ensure all new skills include comprehensive `SKILL.md` files
5. Test your examples and workflows
6. Submit a pull request with a clear description of your changes

### Guidelines
- Maintain consistency with existing skill documentation format
- Include practical, working examples in all contributions
- Ensure all code examples are tested and functional
- Follow scientific best practices in examples and workflows
- Update relevant sections of this README when adding new capabilities

Your contributions help make scientific computing more accessible and enable researchers to leverage AI tools more effectively!
