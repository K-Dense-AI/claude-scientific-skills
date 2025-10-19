# Claude Scientific Skills

A comprehensive collection of ready-to-use scientific skills for Claude, curated by the K-Dense team. These skills enable Claude to work with specialized scientific libraries and databases across bioinformatics, cheminformatics, machine learning, materials science, and data analysis.

## Available Skills

### Scientific Databases

- **PubChem** - Access chemical compound data from the world's largest free chemical database (110M+ compounds, 270M+ bioactivities)

### Scientific Packages

**Bioinformatics & Genomics:**
- **AnnData** - Annotated data matrices for single-cell genomics and h5ad files
- **Arboreto** - Gene regulatory network inference using GRNBoost2 and GENIE3
- **BioPython** - Sequence manipulation, NCBI database access, BLAST searches, alignments, and phylogenetics
- **BioServices** - Programmatic access to 40+ biological web services (KEGG, UniProt, ChEBI, ChEMBL)
- **Cellxgene Census** - Query and analyze large-scale single-cell RNA-seq data
- **gget** - Efficient genomic database queries (Ensembl, UniProt, NCBI, PDB, COSMIC)
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
- **BioMNI** - Multi-omics data integration with LLM-powered analysis

## Try in Claude Code, Claude.ai, and the API

### Claude Code
You can register this repository as a Claude Code Plugin marketplace by running the following command in Claude Code:

```
/plugin marketplace add K-Dense-AI/claude-scientific-skills
```

Then, to install a specific set of skills:

1. Select Browse and install plugins
2. Select claude-scientific-skills
3. Select scientific-databases or scientific-packages
4. Select Install now


After installing the plugin, you can use the skill by just mentioning it. Additionally, in most case, Claude Code will figure out what to use based on the task.

### Claude.ai
These example skills are all already available to paid plans in Claude.ai.

To use any skill from this repository or upload custom skills, follow the instructions in [Using skills in Claude](https://docs.anthropic.com/claude/skills).

### Claude API
You can use Anthropic's pre-built skills, and upload custom skills, via the Claude API. See the [Skills API Quickstart](https://docs.anthropic.com/claude/skills-api-quickstart) for more.
