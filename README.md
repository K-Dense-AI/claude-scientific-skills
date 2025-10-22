# Claude Scientific Skills

[![License: PolyForm Noncommercial 1.0.0](https://img.shields.io/badge/License-PolyForm%20Noncommercial-blue.svg)](LICENSE.md)
[![GitHub Stars](https://img.shields.io/github/stars/K-Dense-AI/claude-scientific-skills?style=social)](https://github.com/K-Dense-AI/claude-scientific-skills)
[![Skills](https://img.shields.io/badge/Skills-73%2B-brightgreen.svg)](#what-s-included)
[![Workflows](https://img.shields.io/badge/Workflows-122-orange.svg)](#what-s-included)

A comprehensive collection of ready-to-use scientific skills for Claude, curated by the K-Dense team.

These skills enable Claude to work with specialized scientific libraries and databases across multiple scientific domains:
- 🧬 Bioinformatics & Genomics
- 🧪 Cheminformatics & Drug Discovery
- 🔬 Proteomics & Mass Spectrometry
- 🤖 Machine Learning & AI
- 🔮 Materials Science & Chemistry
- 📊 Data Analysis & Visualization

**Transform Claude Code into an 'AI Scientist' on your desktop!**

> 💼 For substantially more advanced capabilities, compute infrastructure, and enterprise-ready offerings, check out [k-dense.ai](https://k-dense.ai/).

---

## 📋 Table of Contents

- [What's Included](#what-s-included)
- [Why Use This?](#why-use-this)
- [Getting Started](#getting-started)
  - [Claude Code](#claude-code)
  - [Any MCP Client](#any-mcp-client-including-chatgpt-cursor-google-adk-openai-agent-sdk-etc)
- [Prerequisites](#prerequisites)
- [Quick Examples](#quick-examples)
- [Use Cases](#use-cases)
- [Available Skills](#available-skills)
- [Contributing](#contributing)
- [Troubleshooting](#troubleshooting)
- [FAQ](#faq)
- [Support](#support)
- [License](#license)

---

## 📦 What's Included

| Category | Count | Description |
|----------|-------|-------------|
| 📊 **Scientific Databases** | 24 | PubMed, PubChem, UniProt, ChEMBL, COSMIC, AlphaFold DB, and more |
| 🔬 **Scientific Packages** | 43 | BioPython, RDKit, PyTorch, Scanpy, and specialized tools |
| 🔌 **Scientific Integrations** | 6 | Benchling, DNAnexus, Opentrons, LabArchives, LatchBio, OMERO |
| 🎯 **Context Initialization** | 1 | Auto-invoked skill to ensure Claude uses existing skills effectively |
| 📚 **Documented Workflows** | 122 | Ready-to-use examples and reference materials |

---

## 🚀 Why Use This?

✅ **Save Time** - Skip days of API documentation research and integration work  
✅ **Best Practices** - Curated workflows following scientific computing standards  
✅ **Production Ready** - Tested and validated code examples  
✅ **Regular Updates** - Maintained and expanded by K-Dense team  
✅ **Comprehensive** - Coverage across major scientific domains  
✅ **Enterprise Support** - Commercial offerings available for advanced needs

---

## 🎯 Getting Started

### Claude Code
Register this repository as a Claude Code Plugin marketplace by running:

```bash
/plugin marketplace add K-Dense-AI/claude-scientific-skills
```

Then, to install a specific set of skills:

1. Select **Browse and install plugins**
2. Select **claude-scientific-skills**
3. Choose from:
   - `scientific-databases` - Access to 24 scientific databases
   - `scientific-packages` - 43 specialized Python packages
   - `scientific-thinking` - Analysis tools and document processing
   - `scientific-integrations` - Lab automation and platform integrations
   - `scientific-context-initialization` - Ensures Claude searches for and uses existing skills
4. Select **Install now**

After installation, simply mention the skill or describe your task - Claude Code will automatically use the appropriate skills!

> 💡 **Tip**: If you find that Claude isn't utilizing the installed skills as much as you'd like, install the `scientific-context-initialization` skill. It automatically creates/updates an `AGENT.md` file in your workspace that instructs Claude to always search for and use existing skills before attempting any scientific task. This ensures Claude leverages documented patterns, authentication methods, working examples, and best practices from the repository.

### Any MCP Client (including ChatGPT, Cursor, Google ADK, OpenAI Agent SDK, etc.)
Use our newly released MCP server that allows you to use any Claude Skill in any client!

🔗 **[claude-skills-mcp](https://github.com/K-Dense-AI/claude-skills-mcp)**

---

## ⚙️ Prerequisites

- **Python**: 3.8+ (3.10+ recommended for best compatibility)
- **Claude Code**: Latest version or any MCP-compatible client
- **System**: macOS, Linux, or Windows with WSL2
- **Dependencies**: Automatically handled by individual skills (check `SKILL.md` files for specific requirements)

---

## 💡 Quick Examples

Once you've installed the skills, you can ask Claude to execute complex multi-step scientific workflows:

### End-to-End Drug Discovery Pipeline
```
"I need to find novel EGFR inhibitors for lung cancer treatment. Query ChEMBL for existing 
EGFR inhibitors with IC50 < 50nM, analyze their structure-activity relationships using RDKit, 
generate similar molecules with improved properties using datamol, perform virtual screening 
with DiffDock against the AlphaFold-predicted EGFR structure, and search PubMed for recent 
papers on resistance mechanisms to prioritize scaffolds. Finally, check COSMIC for common 
EGFR mutations and assess how our candidates might interact with mutant forms."
```

### Comprehensive Single-Cell Analysis Workflow
```
"Load this 10X Genomics dataset using Scanpy, perform quality control and doublet removal, 
integrate with public data from Cellxgene Census for the same tissue type, identify cell 
populations using known markers from NCBI Gene, perform differential expression analysis 
with PyDESeq2, run gene regulatory network inference with Arboreto, query Reactome and 
KEGG for pathway enrichment, and create publication-quality visualizations with matplotlib. 
Then cross-reference top dysregulated genes with Open Targets to identify potential 
therapeutic targets."
```

### Multi-Omics Integration for Biomarker Discovery
```
"I have RNA-seq, proteomics, and metabolomics data from cancer patients. Use PyDESeq2 for 
differential expression, pyOpenMS to analyze mass spec data, and integrate metabolite 
information from HMDB and Metabolomics Workbench. Map proteins to pathways using UniProt 
and KEGG, identify protein-protein interactions via STRING, correlate multi-omics layers 
using statsmodels, and build a machine learning model with scikit-learn to predict patient 
outcomes. Search ClinicalTrials.gov for ongoing trials targeting the top candidates."
```

### Structure-Based Virtual Screening Campaign
```
"I want to discover allosteric modulators for a protein-protein interaction. Retrieve the 
AlphaFold structure for both proteins, identify the interaction interface using BioPython, 
search ZINC15 for molecules with suitable properties for allosteric binding (MW 300-500, 
logP 2-4), filter for drug-likeness using RDKit, perform molecular docking with DiffDock 
to identify potential allosteric sites, rank candidates using DeepChem's property prediction 
models, check PubChem for suppliers, and search USPTO patents to assess freedom to operate. 
Finally, generate analogs with MedChem and molfeat for lead optimization."
```

### Clinical Genomics Variant Interpretation Pipeline
```
"Analyze this VCF file from a patient with suspected hereditary cancer. Use pysam to parse 
variants, annotate with Ensembl for functional consequences, query ClinVar for known 
pathogenic variants, check COSMIC for somatic mutations in cancer, retrieve gene information 
from NCBI Gene, analyze protein impact using UniProt, search PubMed for case reports of 
similar variants, query ClinPGx for pharmacogenomic implications, and generate a clinical 
report with ReportLab. Then search ClinicalTrials.gov for precision medicine trials matching 
the patient's profile."
```

### Systems Biology Network Analysis
```
"Starting with a list of differentially expressed genes from my RNA-seq experiment, query 
NCBI Gene for detailed annotations, retrieve protein sequences from UniProt, identify 
protein-protein interactions using STRING, map to biological pathways in Reactome and KEGG, 
analyze network topology with Torch Geometric, identify hub genes and bottleneck proteins, 
perform gene regulatory network reconstruction with Arboreto, integrate with Open Targets 
for druggability assessment, use PyMC for Bayesian network modeling, and create interactive 
network visualizations. Finally, search GEO for similar expression patterns across diseases."
```

---

## 🔬 Use Cases

### Drug Discovery Research
- Screen compound libraries from PubChem and ZINC
- Analyze bioactivity data from ChEMBL
- Predict molecular properties with RDKit and DeepChem
- Perform molecular docking with DiffDock

### Bioinformatics Analysis
- Process genomic sequences with BioPython
- Analyze single-cell RNA-seq data with Scanpy
- Query gene information from Ensembl and NCBI Gene
- Identify protein-protein interactions via STRING

### Materials Science
- Analyze crystal structures with Pymatgen
- Predict material properties
- Design novel compounds and materials

### Clinical Research
- Search clinical trials on ClinicalTrials.gov
- Analyze genetic variants in ClinVar
- Review pharmacogenomic data from ClinPGx
- Access cancer mutations from COSMIC

### Academic Research
- Literature searches via PubMed
- Patent landscape analysis using USPTO
- Data visualization for publications
- Statistical analysis and hypothesis testing

---

## 📚 Available Skills

### 🗄️ Scientific Databases
**24 comprehensive databases** including PubMed, PubChem, UniProt, ChEMBL, AlphaFold DB, COSMIC, Ensembl, KEGG, and more.

📖 **[Full Database Documentation →](docs/scientific-databases.md)**

<details>
<summary><strong>View all databases</strong></summary>

- **AlphaFold DB** - AI-predicted protein structures (200M+ predictions)
- **ChEMBL** - Bioactive molecules and drug-like properties
- **ClinPGx** - Clinical pharmacogenomics and gene-drug interactions
- **ClinVar** - Genomic variants and clinical significance
- **ClinicalTrials.gov** - Global clinical studies registry
- **COSMIC** - Somatic cancer mutations database
- **ENA** - European Nucleotide Archive
- **Ensembl** - Genome browser and annotations
- **FDA Databases** - Drug approvals, adverse events, recalls
- **GEO** - Gene expression and functional genomics
- **GWAS Catalog** - Genome-wide association studies
- **HMDB** - Human metabolome database
- **KEGG** - Biological pathways and molecular interactions
- **Metabolomics Workbench** - NIH metabolomics data
- **NCBI Gene** - Gene information and annotations
- **Open Targets** - Therapeutic target identification
- **PDB** - Protein structure database
- **PubChem** - Chemical compound data (110M+ compounds)
- **PubMed** - Biomedical literature database
- **Reactome** - Curated biological pathways
- **STRING** - Protein-protein interaction networks
- **UniProt** - Protein sequences and annotations
- **USPTO** - Patent and trademark data
- **ZINC** - Commercially-available compounds for screening

</details>

---

### 🔬 Scientific Packages
**43 specialized Python packages** organized by domain.

📖 **[Full Package Documentation →](docs/scientific-packages.md)**

<details>
<summary><strong>Bioinformatics & Genomics (11 packages)</strong></summary>

- AnnData, Arboreto, BioPython, BioServices, Cellxgene Census
- deepTools, FlowIO, gget, pysam, PyDESeq2, Scanpy

</details>

<details>
<summary><strong>Cheminformatics & Drug Discovery (8 packages)</strong></summary>

- Datamol, DeepChem, DiffDock, MedChem, Molfeat, PyTDC, RDKit, TorchDrug

</details>

<details>
<summary><strong>Proteomics & Mass Spectrometry (2 packages)</strong></summary>

- matchms, pyOpenMS

</details>

<details>
<summary><strong>Machine Learning & Deep Learning (8 packages)</strong></summary>

- PyMC, PyMOO, PyTorch Lightning, scikit-learn, statsmodels
- Torch Geometric, Transformers, UMAP-learn

</details>

<details>
<summary><strong>Materials Science & Chemistry (3 packages)</strong></summary>

- Astropy, COBRApy, Pymatgen

</details>

<details>
<summary><strong>Data Analysis & Visualization (5 packages)</strong></summary>

- Dask, Matplotlib, Polars, ReportLab, Seaborn

</details>

<details>
<summary><strong>Additional Packages (6 packages)</strong></summary>

- BIOMNI (Multi-omics), ETE Toolkit (Phylogenetics)
- Paper-2-Web (Academic paper dissemination and presentation)
- scikit-bio (Sequence analysis), ToolUniverse (600+ scientific tool ecosystem)
- Zarr (Array storage)

</details>

---

### 🧠 Scientific Thinking & Analysis
**Comprehensive analysis tools** and document processing capabilities.

📖 **[Full Thinking & Analysis Documentation →](docs/scientific-thinking.md)**

**Analysis & Methodology:**
- Exploratory Data Analysis (automated statistics and insights)
- Hypothesis Generation (structured frameworks)
- Peer Review (comprehensive evaluation toolkit)
- Scientific Brainstorming (ideation workflows)
- Scientific Critical Thinking (rigorous reasoning)
- Scientific Visualization (publication-quality figures)
- Scientific Writing (IMRAD format, citation styles)
- Statistical Analysis (testing and experimental design)

**Document Processing:**
- DOCX, PDF, PPTX, XLSX manipulation and analysis
- Tracked changes, comments, and formatting preservation
- Text extraction, table parsing, and data analysis

---

### 🔌 Scientific Integrations
**6 platform integrations** for lab automation and workflow management.

📖 **[Full Integration Documentation →](docs/scientific-integrations.md)**

- **Benchling** - R&D platform and LIMS integration
- **DNAnexus** - Cloud genomics and biomedical data analysis
- **LabArchives** - Electronic Lab Notebook (ELN) integration
- **LatchBio** - Workflow platform and cloud execution
- **OMERO** - Microscopy and bio-image data management
- **Opentrons** - Laboratory automation protocols

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

✅ Maintain consistency with existing skill documentation format  
✅ Include practical, working examples in all contributions  
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

📖 **[Contributing Guidelines →](CONTRIBUTING.md)** *(coming soon)*

---

## 🔧 Troubleshooting

### Common Issues

**Problem: Claude not using installed skills**
- Solution: Install the `scientific-context-initialization` skill
- This creates an `AGENT.md` file that instructs Claude to search for and use existing skills before attempting tasks
- After installation, Claude will automatically leverage documented patterns, examples, and best practices

**Problem: Skills not loading in Claude Code**
- Solution: Ensure you've installed the latest version of Claude Code
- Try reinstalling the plugin: `/plugin marketplace add K-Dense-AI/claude-scientific-skills`

**Problem: Missing Python dependencies**
- Solution: Check the specific `SKILL.md` file for required packages
- Install dependencies: `pip install package-name`

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

**Q: Is this free to use?**  
A: Yes, for noncommercial use. See the [License](#license) section for details.

**Q: Do I need all the Python packages installed?**  
A: No, only install the packages you need. Each skill specifies its requirements.

**Q: Can I use this with other AI models?**  
A: The skills are designed for Claude but can be adapted for other models with MCP support.

**Q: How often is this updated?**  
A: We regularly update skills to reflect the latest versions of packages and APIs.

**Q: Can I use this for commercial projects?**  
A: For commercial use, please visit [K-Dense](https://k-dense.ai/) for enterprise licensing.

**Q: What if a skill doesn't work?**  
A: First check the troubleshooting section, then file an issue on GitHub with details.

**Q: Can I contribute my own skills?**  
A: Absolutely! See the [Contributing](#contributing) section for guidelines.

**Q: Do the skills work offline?**  
A: Database skills require internet access. Package skills work offline once dependencies are installed.

---

## 💬 Support

Need help? Here's how to get support:

- 📖 **Documentation**: Check the relevant `SKILL.md` and `references/` folders
- 🐛 **Bug Reports**: [Open an issue](https://github.com/K-Dense-AI/claude-scientific-skills/issues)
- 💡 **Feature Requests**: [Submit a feature request](https://github.com/K-Dense-AI/claude-scientific-skills/issues/new)
- 💼 **Enterprise Support**: Contact [K-Dense](https://k-dense.ai/) for commercial support
- 🌐 **MCP Support**: Visit the [claude-skills-mcp](https://github.com/K-Dense-AI/claude-skills-mcp) repository

---

## 📄 License

This project is licensed under the **PolyForm Noncommercial License 1.0.0**.

**Copyright © K-Dense Inc.** ([k-dense.ai](https://k-dense.ai/))

### Key Points:
- ✅ **Free for noncommercial use** (research, education, personal projects)
- ✅ **Free for noncommercial organizations** (universities, research institutions)
- ❌ **Commercial use requires separate license** (contact K-Dense)

See [LICENSE.md](LICENSE.md) for full terms.
