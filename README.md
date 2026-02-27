# Claude Scientific Skills

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE.md)
[![Skills](https://img.shields.io/badge/Skills-148-brightgreen.svg)](#whats-included)
[![Databases](https://img.shields.io/badge/Databases-250%2B-orange.svg)](#whats-included)
[![Agent Skills](https://img.shields.io/badge/Standard-Agent_Skills-blueviolet.svg)](https://agentskills.io/)
[![Works with](https://img.shields.io/badge/Works_with-Cursor_|_Claude_Code_|_Codex-blue.svg)](#getting-started)

A comprehensive collection of **148+ ready-to-use scientific and research skills** (now including financial/SEC research, U.S. Treasury fiscal data, OFR Hedge Fund Monitor, and Alpha Vantage market data) for any AI agent that supports the open [Agent Skills](https://agentskills.io/) standard, created by [K-Dense](https://k-dense.ai). Works with **Cursor, Claude Code, Codex, and more**. Transform your AI agent into a research assistant capable of executing complex multi-step scientific workflows across biology, chemistry, medicine, and beyond.

**Looking for the full AI co-scientist experience?** Try [K-Dense Web](https://k-dense.ai) for 200+ skills, cloud compute, and publication-ready outputs.

<p align="center">
  <a href="https://k-dense.ai">
    <img src="docs/k-dense-web.gif" alt="K-Dense Web Demo" width="800"/>
  </a>
</p>

---

## K-Dense Web - The Full Experience

Want 10x the power with zero setup? **[K-Dense Web](https://k-dense.ai)** is the complete AI co-scientist platform—everything in this repo, plus:

| Feature | This Repo | K-Dense Web |
|---------|-----------|-------------|
| Scientific Skills | 148 skills | **200+ skills** (exclusive access) |
| Setup Required | Manual installation | **Zero setup** — works instantly |
| Compute | Your machine | **Cloud GPUs & HPC** included |
| Workflows | Basic prompts | **End-to-end research pipelines** |
| Outputs | Code & analysis | **Publication-ready** figures, reports & papers |
| Integrations | Local tools | **Lab systems, ELNs, cloud storage** |

**Researchers at Stanford, MIT, and leading pharma companies use K-Dense Web to accelerate discoveries.**

**Get $50 in free credits** — no credit card required.

<a href="https://k-dense.ai"><img src="https://img.shields.io/badge/Try_K--Dense_Web-Start_Free-blue?style=for-the-badge&logo=data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHdpZHRoPSIyNCIgaGVpZ2h0PSIyNCIgdmlld0JveD0iMCAwIDI0IDI0IiBmaWxsPSJub25lIiBzdHJva2U9IndoaXRlIiBzdHJva2Utd2lkdGg9IjIiIHN0cm9rZS1saW5lY2FwPSJyb3VuZCIgc3Ryb2tlLWxpbmVqb2luPSJyb3VuZCI+PHBhdGggZD0iTTUgMTJoMTQiLz48cGF0aCBkPSJtMTIgNSA3IDctNyA3Ii8+PC9zdmc+" alt="Try K-Dense Web"></a>

*Learn more at [k-dense.ai](https://k-dense.ai)* | *[Read our detailed comparison →](https://k-dense.ai/blog/k-dense-web-vs-claude-scientific-skills)*

---

These skills enable your AI agent to seamlessly work with specialized scientific libraries, databases, and tools across multiple scientific domains. While the agent can use any Python package or API on its own, these explicitly defined skills provide curated documentation and examples that make it significantly stronger and more reliable for the workflows below:
- 🧬 Bioinformatics & Genomics - Sequence analysis, single-cell RNA-seq, gene regulatory networks, variant annotation, phylogenetic analysis
- 🧪 Cheminformatics & Drug Discovery - Molecular property prediction, virtual screening, ADMET analysis, molecular docking, lead optimization
- 🔬 Proteomics & Mass Spectrometry - LC-MS/MS processing, peptide identification, spectral matching, protein quantification
- 🏥 Clinical Research & Precision Medicine - Clinical trials, pharmacogenomics, variant interpretation, drug safety, clinical decision support, treatment planning
- 🧠 Healthcare AI & Clinical ML - EHR analysis, physiological signal processing, medical imaging, clinical prediction models
- 🖼️ Medical Imaging & Digital Pathology - DICOM processing, whole slide image analysis, computational pathology, radiology workflows
- 🤖 Machine Learning & AI - Deep learning, reinforcement learning, time series analysis, model interpretability, Bayesian methods
- 🔮 Materials Science & Chemistry - Crystal structure analysis, phase diagrams, metabolic modeling, computational chemistry
- 🌌 Physics & Astronomy - Astronomical data analysis, coordinate transformations, cosmological calculations, symbolic mathematics, physics computations
- ⚙️ Engineering & Simulation - Discrete-event simulation, multi-objective optimization, metabolic engineering, systems modeling, process optimization
- 📊 Data Analysis & Visualization - Statistical analysis, network analysis, time series, publication-quality figures, large-scale data processing, EDA
- 🧪 Laboratory Automation - Liquid handling protocols, lab equipment control, workflow automation, LIMS integration
- 📚 Scientific Communication - Literature review, peer review, scientific writing, document processing, posters, slides, schematics, citation management
- 🔬 Multi-omics & Systems Biology - Multi-modal data integration, pathway analysis, network biology, systems-level insights
- 🧬 Protein Engineering & Design - Protein language models, structure prediction, sequence design, function annotation
- 🎓 Research Methodology - Hypothesis generation, scientific brainstorming, critical thinking, grant writing, scholar evaluation

**Transform your AI coding agent into an 'AI Scientist' on your desktop!**

> ⭐ **If you find this repository useful**, please consider giving it a star! It helps others discover these tools and encourages us to continue maintaining and expanding this collection.

> 🎬 **New to Claude Scientific Skills?** Watch our [Getting Started with Claude Scientific Skills](https://youtu.be/ZxbnDaD_FVg) video for a quick walkthrough.

---

## 📦 What's Included

This repository provides **148 scientific and research skills** organized into the following categories:

- **250+ Scientific & Financial Databases** - Collectively, these skills provide access to over 250 databases and data sources. Dedicated skills cover PubMed, ChEMBL, UniProt, COSMIC, ClinicalTrials.gov, SEC EDGAR, Alpha Vantage, and more; multi-database packages like BioServices (~40 bioinformatics services + 30+ PSICQUIC interaction databases), BioPython (38 NCBI sub-databases via Entrez), and gget (20+ genomics databases) account for the rest
- **55+ Optimized Python Package Skills** - Explicitly defined skills for RDKit, Scanpy, PyTorch Lightning, scikit-learn, BioPython, pyzotero, BioServices, PennyLane, Qiskit, and others — with curated documentation, examples, and best practices. Note: the agent can write code using *any* Python package, not just these; these skills simply provide stronger, more reliable performance for the packages listed
- **15+ Scientific Integration Skills** - Explicitly defined skills for Benchling, DNAnexus, LatchBio, OMERO, Protocols.io, and more. Again, the agent is not limited to these — any API or platform reachable from Python is fair game; these skills are the optimized, pre-documented paths
- **30+ Analysis & Communication Tools** - Literature review, scientific writing, peer review, document processing, posters, slides, schematics, and more
- **10+ Research & Clinical Tools** - Hypothesis generation, grant writing, clinical decision support, treatment plans, regulatory compliance

Each skill includes:
- ✅ Comprehensive documentation (`SKILL.md`)
- ✅ Practical code examples
- ✅ Use cases and best practices
- ✅ Integration guides
- ✅ Reference materials

---

## 📋 Table of Contents

- [What's Included](#whats-included)
- [Why Use This?](#why-use-this)
- [Getting Started](#getting-started)
- [Support Open Source](#-support-the-open-source-community)
- [Prerequisites](#prerequisites)
- [Quick Examples](#quick-examples)
- [Use Cases](#use-cases)
- [Available Skills](#available-skills)
- [Contributing](#contributing)
- [Troubleshooting](#troubleshooting)
- [FAQ](#faq)
- [Support](#support)
- [Join Our Community](#join-our-community)
- [Citation](#citation)
- [License](#license)

---

## 🚀 Why Use This?

### ⚡ **Accelerate Your Research**
- **Save Days of Work** - Skip API documentation research and integration setup
- **Production-Ready Code** - Tested, validated examples following scientific best practices
- **Multi-Step Workflows** - Execute complex pipelines with a single prompt

### 🎯 **Comprehensive Coverage**
- **148 Skills** - Extensive coverage across all major scientific domains
- **250+ Databases** - Collective access to 250+ databases and data sources spanning genomics, chemistry, clinical, financial, and more — through dedicated database skills and multi-database packages like BioServices, BioPython, and gget
- **55+ Optimized Python Package Skills** - RDKit, Scanpy, PyTorch Lightning, scikit-learn, BioServices, PennyLane, Qiskit, and others (the agent can use any Python package; these are the pre-documented, higher-performing paths)

### 🔧 **Easy Integration**
- **Simple Setup** - Copy skills to your skills directory and start working
- **Automatic Discovery** - Your agent automatically finds and uses relevant skills
- **Well Documented** - Each skill includes examples, use cases, and best practices

### 🌟 **Maintained & Supported**
- **Regular Updates** - Continuously maintained and expanded by K-Dense team
- **Community Driven** - Open source with active community contributions
- **Enterprise Ready** - Commercial support available for advanced needs

---

## 🎯 Getting Started

Claude Scientific Skills follows the open [Agent Skills](https://agentskills.io/) standard. Simply copy the skill folders into your skills directory and your AI agent will automatically discover and use them.

### Step 1: Clone the Repository

```bash
git clone https://github.com/K-Dense-AI/claude-scientific-skills.git
```

### Step 2: Copy Skills to Your Skills Directory

Copy the individual skill folders from `scientific-skills/` to one of the supported skill directories below. You can install skills **globally** (available across all projects) or **per-project** (available only in that project).

**Global installation** (recommended — skills available everywhere):

| Tool | Directory |
|------|-----------|
| Cursor | `~/.cursor/skills/` |
| Claude Code | `~/.claude/skills/` |
| Codex | `~/.codex/skills/` |

**Project-level installation** (skills scoped to a single project):

| Tool | Directory |
|------|-----------|
| Cursor | `.cursor/skills/` (in your project root) |
| Claude Code | `.claude/skills/` (in your project root) |
| Codex | `.codex/skills/` (in your project root) |

> **Note:** Cursor also reads from `.claude/skills/` and `.codex/skills/` directories, and vice versa, so skills are cross-compatible between tools.

**Example — global install for Cursor:**
```bash
cp -r claude-scientific-skills/scientific-skills/* ~/.cursor/skills/
```

**Example — global install for Claude Code:**
```bash
cp -r claude-scientific-skills/scientific-skills/* ~/.claude/skills/
```

**Example — project-level install:**
```bash
mkdir -p .cursor/skills
cp -r /path/to/claude-scientific-skills/scientific-skills/* .cursor/skills/
```

**That's it!** Your AI agent will automatically discover the skills and use them when relevant to your scientific tasks. You can also invoke any skill manually by mentioning the skill name in your prompt.

---

## ❤️ Support the Open Source Community

Claude Scientific Skills is powered by **50+ incredible open source projects** maintained by dedicated developers and research communities worldwide. Projects like Biopython, Scanpy, RDKit, scikit-learn, PyTorch Lightning, and many others form the foundation of these skills.

**If you find value in this repository, please consider supporting the projects that make it possible:**

- ⭐ **Star their repositories** on GitHub
- 💰 **Sponsor maintainers** via GitHub Sponsors or NumFOCUS
- 📝 **Cite projects** in your publications
- 💻 **Contribute** code, docs, or bug reports

👉 **[View the full list of projects to support](docs/open-source-sponsors.md)**

---

## ⚙️ Prerequisites

- **Python**: 3.9+ (3.12+ recommended for best compatibility)
- **uv**: Python package manager (required for installing skill dependencies)
- **Client**: Any agent that supports the [Agent Skills](https://agentskills.io/) standard (Cursor, Claude Code, Codex, etc.)
- **System**: macOS, Linux, or Windows with WSL2
- **Dependencies**: Automatically handled by individual skills (check `SKILL.md` files for specific requirements)

### Installing uv

The skills use `uv` as the package manager for installing Python dependencies. Install it using the instructions for your operating system:

**macOS and Linux:**
```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**Windows:**
```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

**Alternative (via pip):**
```bash
pip install uv
```

After installation, verify it works by running:
```bash
uv --version
```

For more installation options and details, visit the [official uv documentation](https://docs.astral.sh/uv/).

---

## 💡 Quick Examples

Once you've installed the skills, you can ask your AI agent to execute complex multi-step scientific workflows. Here are some example prompts:

### 🧪 Drug Discovery Pipeline
**Goal**: Find novel EGFR inhibitors for lung cancer treatment

**Prompt**:
```
Use available skills you have access to whenever possible. Query ChEMBL for EGFR inhibitors (IC50 < 50nM), analyze structure-activity relationships 
with RDKit, generate improved analogs with datamol, perform virtual screening with DiffDock 
against AlphaFold EGFR structure, search PubMed for resistance mechanisms, check COSMIC for 
mutations, and create visualizations and a comprehensive report.
```

**Skills Used**: ChEMBL, RDKit, datamol, DiffDock, AlphaFold DB, PubMed, COSMIC, scientific visualization

---

### 🔬 Single-Cell RNA-seq Analysis
**Goal**: Comprehensive analysis of 10X Genomics data with public data integration

**Prompt**:
```
Use available skills you have access to whenever possible. Load 10X dataset with Scanpy, perform QC and doublet removal, integrate with Cellxgene 
Census data, identify cell types using NCBI Gene markers, run differential expression with 
PyDESeq2, infer gene regulatory networks with Arboreto, enrich pathways via Reactome/KEGG, 
and identify therapeutic targets with Open Targets.
```

**Skills Used**: Scanpy, Cellxgene Census, NCBI Gene, PyDESeq2, Arboreto, Reactome, KEGG, Open Targets

---

### 🧬 Multi-Omics Biomarker Discovery
**Goal**: Integrate RNA-seq, proteomics, and metabolomics to predict patient outcomes

**Prompt**:
```
Use available skills you have access to whenever possible. Analyze RNA-seq with PyDESeq2, process mass spec with pyOpenMS, integrate metabolites from 
HMDB/Metabolomics Workbench, map proteins to pathways (UniProt/KEGG), find interactions via 
STRING, correlate omics layers with statsmodels, build predictive model with scikit-learn, 
and search ClinicalTrials.gov for relevant trials.
```

**Skills Used**: PyDESeq2, pyOpenMS, HMDB, Metabolomics Workbench, UniProt, KEGG, STRING, statsmodels, scikit-learn, ClinicalTrials.gov

---

### 🎯 Virtual Screening Campaign
**Goal**: Discover allosteric modulators for protein-protein interactions

**Prompt**:
```
Use available skills you have access to whenever possible. Retrieve AlphaFold structures, identify interaction interface with BioPython, search ZINC 
for allosteric candidates (MW 300-500, logP 2-4), filter with RDKit, dock with DiffDock, 
rank with DeepChem, check PubChem suppliers, search USPTO patents, and optimize leads with 
MedChem/molfeat.
```

**Skills Used**: AlphaFold DB, BioPython, ZINC, RDKit, DiffDock, DeepChem, PubChem, USPTO, MedChem, molfeat

---

### 🏥 Clinical Variant Interpretation
**Goal**: Analyze VCF file for hereditary cancer risk assessment

**Prompt**:
```
Use available skills you have access to whenever possible. Parse VCF with pysam, annotate variants with Ensembl VEP, query ClinVar for pathogenicity, 
check COSMIC for cancer mutations, retrieve gene info from NCBI Gene, analyze protein impact 
with UniProt, search PubMed for case reports, check ClinPGx for pharmacogenomics, generate 
clinical report with document processing tools, and find matching trials on ClinicalTrials.gov.
```

**Skills Used**: pysam, Ensembl, ClinVar, COSMIC, NCBI Gene, UniProt, PubMed, ClinPGx, Document Skills, ClinicalTrials.gov

---

### 🌐 Systems Biology Network Analysis
**Goal**: Analyze gene regulatory networks from RNA-seq data

**Prompt**:
```
Use available skills you have access to whenever possible. Query NCBI Gene for annotations, retrieve sequences from UniProt, identify interactions via 
STRING, map to Reactome/KEGG pathways, analyze topology with Torch Geometric, reconstruct 
GRNs with Arboreto, assess druggability with Open Targets, model with PyMC, visualize 
networks, and search GEO for similar patterns.
```

**Skills Used**: NCBI Gene, UniProt, STRING, Reactome, KEGG, Torch Geometric, Arboreto, Open Targets, PyMC, GEO

> 📖 **Want more examples?** Check out [docs/examples.md](docs/examples.md) for comprehensive workflow examples and detailed use cases across all scientific domains.

---

## 🔬 Use Cases

### 🧪 Drug Discovery & Medicinal Chemistry
- **Virtual Screening**: Screen millions of compounds from PubChem/ZINC against protein targets
- **Lead Optimization**: Analyze structure-activity relationships with RDKit, generate analogs with datamol
- **ADMET Prediction**: Predict absorption, distribution, metabolism, excretion, and toxicity with DeepChem
- **Molecular Docking**: Predict binding poses and affinities with DiffDock
- **Bioactivity Mining**: Query ChEMBL for known inhibitors and analyze SAR patterns

### 🧬 Bioinformatics & Genomics
- **Sequence Analysis**: Process DNA/RNA/protein sequences with BioPython and pysam
- **Single-Cell Analysis**: Analyze 10X Genomics data with Scanpy, identify cell types, infer GRNs with Arboreto
- **Variant Annotation**: Annotate VCF files with Ensembl VEP, query ClinVar for pathogenicity
- **Variant Database Management**: Build scalable VCF databases with TileDB-VCF for incremental sample addition, efficient population-scale queries, and compressed storage of genomic variant data
- **Gene Discovery**: Query NCBI Gene, UniProt, and Ensembl for comprehensive gene information
- **Network Analysis**: Identify protein-protein interactions via STRING, map to pathways (KEGG, Reactome)

### 🏥 Clinical Research & Precision Medicine
- **Clinical Trials**: Search ClinicalTrials.gov for relevant studies, analyze eligibility criteria
- **Variant Interpretation**: Annotate variants with ClinVar, COSMIC, and ClinPGx for pharmacogenomics
- **Drug Safety**: Query FDA databases for adverse events, drug interactions, and recalls
- **Precision Therapeutics**: Match patient variants to targeted therapies and clinical trials

### 🔬 Multi-Omics & Systems Biology
- **Multi-Omics Integration**: Combine RNA-seq, proteomics, and metabolomics data
- **Pathway Analysis**: Enrich differentially expressed genes in KEGG/Reactome pathways
- **Network Biology**: Reconstruct gene regulatory networks, identify hub genes
- **Biomarker Discovery**: Integrate multi-omics layers to predict patient outcomes

### 📊 Data Analysis & Visualization
- **Statistical Analysis**: Perform hypothesis testing, power analysis, and experimental design
- **Publication Figures**: Create publication-quality visualizations with matplotlib and seaborn
- **Network Visualization**: Visualize biological networks with NetworkX
- **Report Generation**: Generate comprehensive PDF reports with Document Skills

### 🧪 Laboratory Automation
- **Protocol Design**: Create Opentrons protocols for automated liquid handling
- **LIMS Integration**: Integrate with Benchling and LabArchives for data management
- **Workflow Automation**: Automate multi-step laboratory workflows

---

## 📚 Available Skills

This repository contains **148 scientific and research skills** organized across multiple domains. Each skill provides comprehensive documentation, code examples, and best practices for working with scientific libraries, databases, and tools.

### Skill Categories

#### 🧬 **Bioinformatics & Genomics** (16+ skills)
- Sequence analysis: BioPython, pysam, scikit-bio, BioServices
- Single-cell analysis: Scanpy, AnnData, scvi-tools, Arboreto, Cellxgene Census
- Genomic tools: gget, geniml, gtars, deepTools, FlowIO, Zarr, TileDB-VCF
- Phylogenetics: ETE Toolkit

#### 🧪 **Cheminformatics & Drug Discovery** (11+ skills)
- Molecular manipulation: RDKit, Datamol, Molfeat
- Deep learning: DeepChem, TorchDrug
- Docking & screening: DiffDock
- Cloud quantum chemistry: Rowan (pKa, docking, cofolding)
- Drug-likeness: MedChem
- Benchmarks: PyTDC

#### 🔬 **Proteomics & Mass Spectrometry** (2 skills)
- Spectral processing: matchms, pyOpenMS

#### 🏥 **Clinical Research & Precision Medicine** (12+ skills)
- Clinical databases: ClinicalTrials.gov, ClinVar, ClinPGx, COSMIC, FDA Databases
- Healthcare AI: PyHealth, NeuroKit2, Clinical Decision Support
- Clinical documentation: Clinical Reports, Treatment Plans
- Variant analysis: Ensembl, NCBI Gene

#### 🖼️ **Medical Imaging & Digital Pathology** (3 skills)
- DICOM processing: pydicom
- Whole slide imaging: histolab, PathML

#### 🧠 **Neuroscience & Electrophysiology** (1 skill)
- Neural recordings: Neuropixels-Analysis (extracellular spikes, silicon probes, spike sorting)

#### 🤖 **Machine Learning & AI** (15+ skills)
- Deep learning: PyTorch Lightning, Transformers, Stable Baselines3, PufferLib
- Classical ML: scikit-learn, scikit-survival, SHAP
- Time series: aeon
- Bayesian methods: PyMC
- Optimization: PyMOO
- Graph ML: Torch Geometric
- Dimensionality reduction: UMAP-learn
- Statistical modeling: statsmodels

#### 🔮 **Materials Science, Chemistry & Physics** (7 skills)
- Materials: Pymatgen
- Metabolic modeling: COBRApy
- Astronomy: Astropy
- Quantum computing: Cirq, PennyLane, Qiskit, QuTiP

#### ⚙️ **Engineering & Simulation** (4 skills)
- Numerical computing: MATLAB/Octave
- Computational fluid dynamics: FluidSim
- Discrete-event simulation: SimPy
- Data processing: Dask, Polars, Vaex

#### 📊 **Data Analysis & Visualization** (14+ skills)
- Visualization: Matplotlib, Seaborn, Plotly, Scientific Visualization
- Geospatial analysis: GeoPandas
- Network analysis: NetworkX
- Symbolic math: SymPy
- Document processing: Document Skills (PDF, DOCX, PPTX, XLSX)
- Data access: Data Commons
- Exploratory data analysis: EDA workflows
- Statistical analysis: Statistical Analysis workflows

#### 🧪 **Laboratory Automation** (3 skills)
- Liquid handling: PyLabRobot
- Protocol management: Protocols.io
- LIMS integration: Benchling, LabArchives

#### 🔬 **Multi-omics & Systems Biology** (5+ skills)
- Pathway analysis: KEGG, Reactome, STRING
- Multi-omics: Denario, HypoGeniC
- Data management: LaminDB

#### 🧬 **Protein Engineering & Design** (2 skills)
- Protein language models: ESM
- Cloud laboratory platform: Adaptyv (automated protein testing and validation)

#### 📚 **Scientific Communication** (20+ skills)
- Literature: OpenAlex, PubMed, bioRxiv, Literature Review
- Web search: Perplexity Search (AI-powered search with real-time information)
- Writing: Scientific Writing, Peer Review
- Document processing: XLSX, MarkItDown, Document Skills
- Publishing: Paper-2-Web, Venue Templates
- Presentations: Scientific Slides, LaTeX Posters, PPTX Posters
- Diagrams: Scientific Schematics
- Citations: Citation Management
- Illustration: Generate Image (AI image generation with FLUX.2 Pro and Gemini 3 Pro (Nano Banana Pro))

#### 🔬 **Scientific Databases** (28+ skills)
- Protein: UniProt, PDB, AlphaFold DB
- Chemical: PubChem, ChEMBL, DrugBank, ZINC, HMDB
- Genomic: Ensembl, NCBI Gene, GEO, ENA, GWAS Catalog
- Literature: bioRxiv (preprints)
- Clinical: ClinVar, COSMIC, ClinicalTrials.gov, ClinPGx, FDA Databases
- Pathways: KEGG, Reactome, STRING
- Targets: Open Targets
- Metabolomics: Metabolomics Workbench
- Enzymes: BRENDA
- Patents: USPTO

#### 🔧 **Infrastructure & Platforms** (6+ skills)
- Cloud compute: Modal
- Genomics platforms: DNAnexus, LatchBio
- Microscopy: OMERO
- Automation: Opentrons
- Resource detection: Get Available Resources

#### 🎓 **Research Methodology & Planning** (8+ skills)
- Ideation: Scientific Brainstorming, Hypothesis Generation
- Critical analysis: Scientific Critical Thinking, Scholar Evaluation
- Funding: Research Grants
- Discovery: Research Lookup
- Market analysis: Market Research Reports

#### ⚖️ **Regulatory & Standards** (1 skill)
- Medical device standards: ISO 13485 Certification

#### 💹 **Financial & SEC Research** (4 skills)
- SEC filings & financial data: edgartools (10-K, 10-Q, 8-K, 13F, Form 4, XBRL, insider trading, institutional holdings)
- U.S. federal fiscal data: usfiscaldata (national debt, Daily/Monthly Treasury Statements, Treasury auctions, interest rates, exchange rates, savings bonds)
- Hedge fund systemic risk: hedgefundmonitor (OFR Hedge Fund Monitor API — Form PF aggregated stats, CFTC futures positioning, FICC sponsored repo, SCOOS dealer financing)
- Global market data: alpha-vantage (real-time & historical stocks, options, forex, crypto, commodities, economic indicators, 50+ technical indicators via Alpha Vantage API)

> 📖 **For complete details on all skills**, see [docs/scientific-skills.md](docs/scientific-skills.md)

> 💡 **Looking for practical examples?** Check out [docs/examples.md](docs/examples.md) for comprehensive workflow examples across all scientific domains.

---

## 🤝 Contributing

We welcome contributions to expand and improve this scientific skills repository!

### Ways to Contribute

✨ **Add New Skills**
- Create skills for additional scientific packages or databases
- Add integrations for scientific platforms and tools

📚 **Improve Existing Skills**
- Enhance documentation with more examples and use cases
- Add new workflows and reference materials
- Improve code examples and scripts
- Fix bugs or update outdated information

🐛 **Report Issues**
- Submit bug reports with detailed reproduction steps
- Suggest improvements or new features

### How to Contribute

1. **Fork** the repository
2. **Create** a feature branch (`git checkout -b feature/amazing-skill`)
3. **Follow** the existing directory structure and documentation patterns
4. **Ensure** all new skills include comprehensive `SKILL.md` files
5. **Test** your examples and workflows thoroughly
6. **Commit** your changes (`git commit -m 'Add amazing skill'`)
7. **Push** to your branch (`git push origin feature/amazing-skill`)
8. **Submit** a pull request with a clear description of your changes

### Contribution Guidelines

✅ **Adhere to the [Agent Skills Specification](https://agentskills.io/specification)** — Every skill must follow the official spec (valid `SKILL.md` frontmatter, naming conventions, directory structure)  
✅ Maintain consistency with existing skill documentation format  
✅ Ensure all code examples are tested and functional  
✅ Follow scientific best practices in examples and workflows  
✅ Update relevant documentation when adding new capabilities  
✅ Provide clear comments and docstrings in code  
✅ Include references to official documentation

### Recognition

Contributors are recognized in our community and may be featured in:
- Repository contributors list
- Special mentions in release notes
- K-Dense community highlights

Your contributions help make scientific computing more accessible and enable researchers to leverage AI tools more effectively!

### Support Open Source

This project builds on 50+ amazing open source projects. If you find value in these skills, please consider [supporting the projects we depend on](docs/open-source-sponsors.md).

---

## 🔧 Troubleshooting

### Common Issues

**Problem: Skills not loading**
- Verify skill folders are in the correct directory (see [Getting Started](#getting-started))
- Each skill folder must contain a `SKILL.md` file
- Restart your agent/IDE after copying skills
- In Cursor, check Settings → Rules to confirm skills are discovered

**Problem: Missing Python dependencies**
- Solution: Check the specific `SKILL.md` file for required packages
- Install dependencies: `uv pip install package-name`

**Problem: API rate limits**
- Solution: Many databases have rate limits. Review the specific database documentation
- Consider implementing caching or batch requests

**Problem: Authentication errors**
- Solution: Some services require API keys. Check the `SKILL.md` for authentication setup
- Verify your credentials and permissions

**Problem: Outdated examples**
- Solution: Report the issue via GitHub Issues
- Check the official package documentation for updated syntax

---

## ❓ FAQ

### General Questions

**Q: Is this free to use?**  
A: Yes! This repository is MIT licensed. However, each individual skill has its own license specified in the `license` metadata field within its `SKILL.md` file—be sure to review and comply with those terms.

**Q: Why are all skills grouped together instead of separate packages?**  
A: We believe good science in the age of AI is inherently interdisciplinary. Bundling all skills together makes it trivial for you (and your agent) to bridge across fields—e.g., combining genomics, cheminformatics, clinical data, and machine learning in one workflow—without worrying about which individual skills to install or wire together.

**Q: Can I use this for commercial projects?**  
A: The repository itself is MIT licensed, which allows commercial use. However, individual skills may have different licenses—check the `license` field in each skill's `SKILL.md` file to ensure compliance with your intended use.

**Q: Do all skills have the same license?**  
A: No. Each skill has its own license specified in the `license` metadata field within its `SKILL.md` file. These licenses may differ from the repository's MIT License. Users are responsible for reviewing and adhering to the license terms of each individual skill they use.

**Q: How often is this updated?**  
A: We regularly update skills to reflect the latest versions of packages and APIs. Major updates are announced in release notes.

**Q: Can I use this with other AI models?**  
A: The skills follow the open [Agent Skills](https://agentskills.io/) standard and work with any compatible agent, including Cursor, Claude Code, and Codex.

### Installation & Setup

**Q: Do I need all the Python packages installed?**  
A: No! Only install the packages you need. Each skill specifies its requirements in its `SKILL.md` file.

**Q: What if a skill doesn't work?**  
A: First check the [Troubleshooting](#troubleshooting) section. If the issue persists, file an issue on GitHub with detailed reproduction steps.

**Q: Do the skills work offline?**  
A: Database skills require internet access to query APIs. Package skills work offline once Python dependencies are installed.

### Contributing

**Q: Can I contribute my own skills?**  
A: Absolutely! We welcome contributions. See the [Contributing](#contributing) section for guidelines and best practices.

**Q: How do I report bugs or suggest features?**  
A: Open an issue on GitHub with a clear description. For bugs, include reproduction steps and expected vs actual behavior.

---

## 💬 Support

Need help? Here's how to get support:

- 📖 **Documentation**: Check the relevant `SKILL.md` and `references/` folders
- 🐛 **Bug Reports**: [Open an issue](https://github.com/K-Dense-AI/claude-scientific-skills/issues)
- 💡 **Feature Requests**: [Submit a feature request](https://github.com/K-Dense-AI/claude-scientific-skills/issues/new)
- 💼 **Enterprise Support**: Contact [K-Dense](https://k-dense.ai/) for commercial support
- 🌐 **Community**: [Join our Slack](https://join.slack.com/t/k-densecommunity/shared_invite/zt-3iajtyls1-EwmkwIZk0g_o74311Tkf5g)

---

## 🎉 Join Our Community!

**We'd love to have you join us!** 🚀

Connect with other scientists, researchers, and AI enthusiasts using AI agents for scientific computing. Share your discoveries, ask questions, get help with your projects, and collaborate with the community!

🌟 **[Join our Slack Community](https://join.slack.com/t/k-densecommunity/shared_invite/zt-3iajtyls1-EwmkwIZk0g_o74311Tkf5g)** 🌟

Whether you're just getting started or you're a power user, our community is here to support you. We share tips, troubleshoot issues together, showcase cool projects, and discuss the latest developments in AI-powered scientific research.

**See you there!** 💬

---

## 📖 Citation

If you use Claude Scientific Skills in your research or project, please cite it as:

### BibTeX
```bibtex
@software{claude_scientific_skills_2026,
  author = {{K-Dense Inc.}},
  title = {Claude Scientific Skills: A Comprehensive Collection of Scientific Tools for Claude AI},
  year = {2026},
  url = {https://github.com/K-Dense-AI/claude-scientific-skills},
  note = {skills covering databases, packages, integrations, and analysis tools}
}
```

### APA
```
K-Dense Inc. (2026). Claude Scientific Skills: A comprehensive collection of scientific tools for Claude AI [Computer software]. https://github.com/K-Dense-AI/claude-scientific-skills
```

### MLA
```
K-Dense Inc. Claude Scientific Skills: A Comprehensive Collection of Scientific Tools for Claude AI. 2026, github.com/K-Dense-AI/claude-scientific-skills.
```

### Plain Text
```
Claude Scientific Skills by K-Dense Inc. (2026)
Available at: https://github.com/K-Dense-AI/claude-scientific-skills
```

We appreciate acknowledgment in publications, presentations, or projects that benefit from these skills!

---

## 📄 License

This project is licensed under the **MIT License**.

**Copyright © 2026 K-Dense Inc.** ([k-dense.ai](https://k-dense.ai/))

### Key Points:
- ✅ **Free for any use** (commercial and noncommercial)
- ✅ **Open source** - modify, distribute, and use freely
- ✅ **Permissive** - minimal restrictions on reuse
- ⚠️ **No warranty** - provided "as is" without warranty of any kind

See [LICENSE.md](LICENSE.md) for full terms.

### Individual Skill Licenses

> ⚠️ **Important**: Each skill has its own license specified in the `license` metadata field within its `SKILL.md` file. These licenses may differ from the repository's MIT License and may include additional terms or restrictions. **Users are responsible for reviewing and adhering to the license terms of each individual skill they use.**

## Star History

[![Star History Chart](https://api.star-history.com/svg?repos=K-Dense-AI/claude-scientific-skills&type=date&legend=top-left)](https://www.star-history.com/#K-Dense-AI/claude-scientific-skills&type=date&legend=top-left)
