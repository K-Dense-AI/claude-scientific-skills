---
name: biomni
description: "Use this skill for autonomous biomedical research execution across genomics, proteomics, drug discovery, and computational biology. Biomni is an AI agent that combines LLM reasoning with retrieval-augmented planning and code generation to autonomously execute complex biomedical tasks. Use when you need: CRISPR guide RNA design and screening experiments, single-cell RNA-seq analysis workflows, molecular ADMET property prediction, GWAS analysis, protein structure prediction, disease classification from multi-omics data, pathway analysis, drug repurposing, biomarker discovery, variant interpretation, cell type annotation, or any biomedical computational task requiring automated code generation, data analysis, and scientific reasoning. The agent autonomously decomposes tasks, retrieves relevant biomedical knowledge from its 11GB knowledge base, generates and executes analysis code, and provides comprehensive results. Ideal for researchers needing automated execution of complex biomedical workflows without manual coding."
---

# Biomni

## Overview

Biomni is a general-purpose biomedical AI agent that autonomously executes research tasks across diverse biomedical subfields. It combines large language model reasoning with retrieval-augmented planning and code-based execution to enhance scientific productivity and hypothesis generation. The system operates with an ~11GB biomedical knowledge base covering molecular, genomic, and clinical domains.

## Quick Start

Initialize and use the Biomni agent with these basic steps:

```python
from biomni.agent import A1

# Initialize agent with data path and LLM model
agent = A1(path='./data', llm='claude-sonnet-4-20250514')

# Execute a biomedical research task
agent.go("Your biomedical task description")
```

The agent will autonomously decompose the task, retrieve relevant biomedical knowledge, generate and execute code, and provide results.

## Installation and Setup

### Environment Preparation

1. **Set up the conda environment:**
   - Follow instructions in `biomni_env/README.md` from the repository
   - Activate the environment: `conda activate biomni_e1`

2. **Install the package:**
   ```bash
   pip install biomni --upgrade
   ```

   Or install from source:
   ```bash
   git clone https://github.com/snap-stanford/biomni.git
   cd biomni
   pip install -e .
   ```

3. **Configure API keys:**

   Set up credentials via environment variables or `.env` file:
   ```bash
   export ANTHROPIC_API_KEY="your-key-here"
   export OPENAI_API_KEY="your-key-here"  # Optional
   ```

4. **Data initialization:**

   On first use, the agent will automatically download the ~11GB biomedical knowledge base.

### LLM Provider Configuration

Biomni supports multiple LLM providers. Configure the default provider using:

```python
from biomni.config import default_config

# Set the default LLM model
default_config.llm = "claude-sonnet-4-20250514"  # Anthropic
# default_config.llm = "gpt-4"  # OpenAI
# default_config.llm = "azure/gpt-4"  # Azure OpenAI
# default_config.llm = "gemini/gemini-pro"  # Google Gemini

# Set timeout (optional)
default_config.timeout_seconds = 1200

# Set data path (optional)
default_config.data_path = "./custom/data/path"
```

Refer to `references/llm_providers.md` for detailed configuration options for each provider.

## Core Biomedical Research Tasks

### 1. CRISPR Screening and Design

Execute CRISPR screening tasks including guide RNA design, off-target analysis, and screening experiment planning:

```python
agent.go("Design a CRISPR screening experiment to identify genes involved in cancer cell resistance to drug X")
```

The agent will:
- Retrieve relevant gene databases
- Design guide RNAs with specificity analysis
- Plan experimental controls and readout strategies
- Generate analysis code for screening results

### 2. Single-Cell RNA-seq Analysis

Perform comprehensive scRNA-seq analysis workflows:

```python
agent.go("Analyze this 10X Genomics scRNA-seq dataset, identify cell types, and find differentially expressed genes between clusters")
```

Capabilities include:
- Quality control and preprocessing
- Dimensionality reduction and clustering
- Cell type annotation using marker databases
- Differential expression analysis
- Pathway enrichment analysis

### 3. Molecular Property Prediction (ADMET)

Predict absorption, distribution, metabolism, excretion, and toxicity properties:

```python
agent.go("Predict ADMET properties for these drug candidates: [SMILES strings]")
```

The agent handles:
- Molecular descriptor calculation
- Property prediction using integrated models
- Toxicity screening
- Drug-likeness assessment

### 4. Genomic Analysis

Execute genomic data analysis tasks:

```python
agent.go("Perform GWAS analysis to identify SNPs associated with disease phenotype in this cohort")
```

Supports:
- Genome-wide association studies (GWAS)
- Variant calling and annotation
- Population genetics analysis
- Functional genomics integration

### 5. Protein Structure and Function

Analyze protein sequences and structures:

```python
agent.go("Predict the structure of this protein sequence and identify potential binding sites")
```

Capabilities:
- Sequence analysis and domain identification
- Structure prediction integration
- Binding site prediction
- Protein-protein interaction analysis

### 6. Disease Diagnosis and Classification

Perform disease classification from multi-omics data:

```python
agent.go("Build a classifier to diagnose disease X from patient RNA-seq and clinical data")
```

### 7. Systems Biology and Pathway Analysis

Analyze biological pathways and networks:

```python
agent.go("Identify dysregulated pathways in this differential expression dataset")
```

### 8. Drug Discovery and Repurposing

Support drug discovery workflows:

```python
agent.go("Identify FDA-approved drugs that could be repurposed for treating disease Y based on mechanism of action")
```

## Advanced Features

### Custom Configuration per Agent

Override global configuration for specific agent instances:

```python
agent = A1(
    path='./project_data',
    llm='gpt-4o',
    timeout=1800
)
```

### Conversation History and Reporting

Save execution traces as formatted PDF reports:

```python
# After executing tasks
agent.save_conversation_history(
    output_path='./reports/experiment_log.pdf',
    format='pdf'
)
```

Requires one of: WeasyPrint, markdown2pdf, or Pandoc.

### Model Context Protocol (MCP) Integration

Extend agent capabilities with external tools:

```python
# Add MCP-compatible tools
agent.add_mcp(config_path='./mcp_config.json')
```

MCP enables integration with:
- Laboratory information management systems (LIMS)
- Specialized bioinformatics databases
- Custom analysis pipelines
- External computational resources

### Using Biomni-R0 (Specialized Reasoning Model)

Deploy the 32B parameter Biomni-R0 model for enhanced biological reasoning:

```bash
# Install SGLang
pip install "sglang[all]"

# Deploy Biomni-R0
python -m sglang.launch_server \
    --model-path snap-stanford/biomni-r0 \
    --port 30000 \
    --trust-remote-code
```

Then configure the agent:

```python
from biomni.config import default_config

default_config.llm = "openai/biomni-r0"
default_config.api_base = "http://localhost:30000/v1"
```

Biomni-R0 provides specialized reasoning for:
- Complex multi-step biological workflows
- Hypothesis generation and evaluation
- Experimental design optimization
- Literature-informed analysis

## Best Practices

### Task Specification

Provide clear, specific task descriptions:

✅ **Good:** "Analyze this scRNA-seq dataset (file: data.h5ad) to identify T cell subtypes, then perform differential expression analysis comparing activated vs. resting T cells"

❌ **Vague:** "Analyze my RNA-seq data"

### Data Organization

Structure data directories for efficient retrieval:

```
project/
├── data/              # Biomni knowledge base
├── raw_data/          # Your experimental data
├── results/           # Analysis outputs
└── reports/           # Generated reports
```

### Iterative Refinement

Use iterative task execution for complex analyses:

```python
# Step 1: Exploratory analysis
agent.go("Load and perform initial QC on the proteomics dataset")

# Step 2: Based on results, refine analysis
agent.go("Based on the QC results, remove low-quality samples and normalize using method X")

# Step 3: Downstream analysis
agent.go("Perform differential abundance analysis with adjusted parameters")
```

### Security Considerations

**CRITICAL:** Biomni executes LLM-generated code with full system privileges. For production use:

1. **Use sandboxed environments:** Deploy in Docker containers or VMs with restricted permissions
2. **Validate sensitive operations:** Review code before execution for file access, network calls, or credential usage
3. **Limit data access:** Restrict agent access to only necessary data directories
4. **Monitor execution:** Log all executed code for audit trails

Never run Biomni with:
- Unrestricted file system access
- Direct access to sensitive credentials
- Network access to production systems
- Elevated system privileges

### Model Selection Guidelines

Choose models based on task complexity:

- **Claude Sonnet 4:** Recommended for most biomedical tasks, excellent biological reasoning
- **GPT-4/GPT-4o:** Strong general capabilities, good for diverse tasks
- **Biomni-R0:** Specialized for complex biological reasoning, multi-step workflows
- **Smaller models:** Use for simple, well-defined tasks to reduce cost

## Evaluation and Benchmarking

Biomni-Eval1 benchmark contains 433 evaluation instances across 10 biological tasks:

- GWAS analysis
- Disease diagnosis
- Gene detection and classification
- Molecular property prediction
- Pathway analysis
- Protein function prediction
- Drug response prediction
- Variant interpretation
- Cell type annotation
- Biomarker discovery

Use the benchmark to:
- Evaluate custom agent configurations
- Compare LLM providers for specific tasks
- Validate analysis pipelines

## Troubleshooting

### Common Issues

**Issue:** Data download fails or times out
**Solution:** Manually download the knowledge base or increase timeout settings

**Issue:** Package dependency conflicts
**Solution:** Some optional dependencies cannot be installed by default due to conflicts. Install specific packages manually and uncomment relevant code sections as documented in the repository

**Issue:** LLM API errors
**Solution:** Verify API key configuration, check rate limits, ensure sufficient credits

**Issue:** Memory errors with large datasets
**Solution:** Process data in chunks, use data subsampling, or deploy on higher-memory instances

### Getting Help

For detailed troubleshooting:
- Review the Biomni GitHub repository issues
- Check `references/api_reference.md` for detailed API documentation
- Consult `references/task_examples.md` for comprehensive task patterns

## Resources

### references/
Detailed reference documentation for advanced usage:

- **api_reference.md:** Complete API documentation for A1 agent, configuration objects, and utility functions
- **llm_providers.md:** Comprehensive guide for configuring all supported LLM providers (Anthropic, OpenAI, Azure, Gemini, Groq, Ollama, AWS Bedrock)
- **task_examples.md:** Extensive collection of biomedical task examples with code patterns

### scripts/
Helper scripts for common operations:

- **setup_environment.py:** Automated environment setup and validation
- **generate_report.py:** Enhanced PDF report generation with custom formatting

Load reference documentation as needed:
```python
# Claude can read reference files when needed for detailed information
# Example: "Check references/llm_providers.md for Azure OpenAI configuration"
```
