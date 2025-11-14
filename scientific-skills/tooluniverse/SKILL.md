---
name: tooluniverse
description: Use this skill when working with scientific research tools and workflows across bioinformatics, cheminformatics, genomics, structural biology, proteomics, and drug discovery. This skill provides access to 600+ scientific tools including machine learning models, datasets, APIs, and analysis packages. Use when searching for scientific tools, executing computational biology workflows, composing multi-step research pipelines, accessing databases like OpenTargets/PubChem/UniProt/PDB/ChEMBL, performing tool discovery for research tasks, or integrating scientific computational resources into LLM workflows.
---

# ToolUniverse

## Overview

ToolUniverse is a unified ecosystem that enables AI agents to function as research scientists by providing standardized access to 600+ scientific resources. Use this skill to discover, execute, and compose scientific tools across multiple research domains including bioinformatics, cheminformatics, genomics, structural biology, proteomics, and drug discovery.

**Key Capabilities:**
- Access 600+ scientific tools, models, datasets, and APIs
- Discover tools using natural language, semantic search, or keywords
- Execute tools through standardized AI-Tool Interaction Protocol
- Compose multi-step workflows for complex research problems
- Integration with Claude Desktop/Code via Model Context Protocol (MCP)

## When to Use This Skill

Use this skill when:
- Searching for scientific tools by function or domain (e.g., "find protein structure prediction tools")
- Executing computational biology workflows (e.g., disease target identification, drug discovery, genomics analysis)
- Accessing scientific databases (OpenTargets, PubChem, UniProt, PDB, ChEMBL, KEGG, etc.)
- Composing multi-step research pipelines (e.g., target discovery → structure prediction → virtual screening)
- Working with bioinformatics, cheminformatics, or structural biology tasks
- Analyzing gene expression, protein sequences, molecular structures, or clinical data
- Performing literature searches, pathway enrichment, or variant annotation
- Building automated scientific research workflows

## Quick Start

### Basic Setup
```python
from tooluniverse import ToolUniverse

# Initialize and load tools
tu = ToolUniverse()
tu.load_tools()  # Loads 600+ scientific tools

# Discover tools
tools = tu.run({
    "name": "Tool_Finder_Keyword",
    "arguments": {
        "description": "disease target associations",
        "limit": 10
    }
})

# Execute a tool
result = tu.run({
    "name": "OpenTargets_get_associated_targets_by_disease_efoId",
    "arguments": {"efoId": "EFO_0000537"}  # Hypertension
})
```

### Model Context Protocol (MCP)
For Claude Desktop/Code integration:
```bash
tooluniverse-smcp
```

## Core Workflows

### 1. Tool Discovery

Find relevant tools for your research task:

**Three discovery methods:**
- `Tool_Finder` - Embedding-based semantic search (requires GPU)
- `Tool_Finder_LLM` - LLM-based semantic search (no GPU required)
- `Tool_Finder_Keyword` - Fast keyword search

**Example:**
```python
# Search by natural language description
tools = tu.run({
    "name": "Tool_Finder_LLM",
    "arguments": {
        "description": "Find tools for RNA sequencing differential expression analysis",
        "limit": 10
    }
})

# Review available tools
for tool in tools:
    print(f"{tool['name']}: {tool['description']}")
```

**See `references/tool-discovery.md` for:**
- Detailed discovery methods and search strategies
- Domain-specific keyword suggestions
- Best practices for finding tools

### 2. Tool Execution

Execute individual tools through the standardized interface:

**Example:**
```python
# Execute disease-target lookup
targets = tu.run({
    "name": "OpenTargets_get_associated_targets_by_disease_efoId",
    "arguments": {"efoId": "EFO_0000616"}  # Breast cancer
})

# Get protein structure
structure = tu.run({
    "name": "AlphaFold_get_structure",
    "arguments": {"uniprot_id": "P12345"}
})

# Calculate molecular properties
properties = tu.run({
    "name": "RDKit_calculate_descriptors",
    "arguments": {"smiles": "CCO"}  # Ethanol
})
```

**See `references/tool-execution.md` for:**
- Real-world execution examples across domains
- Tool parameter handling and validation
- Result processing and error handling
- Best practices for production use

### 3. Tool Composition and Workflows

Compose multiple tools for complex research workflows:

**Drug Discovery Example:**
```python
# 1. Find disease targets
targets = tu.run({
    "name": "OpenTargets_get_associated_targets_by_disease_efoId",
    "arguments": {"efoId": "EFO_0000616"}
})

# 2. Get protein structures
structures = []
for target in targets[:5]:
    structure = tu.run({
        "name": "AlphaFold_get_structure",
        "arguments": {"uniprot_id": target['uniprot_id']}
    })
    structures.append(structure)

# 3. Screen compounds
hits = []
for structure in structures:
    compounds = tu.run({
        "name": "ZINC_virtual_screening",
        "arguments": {
            "structure": structure,
            "library": "lead-like",
            "top_n": 100
        }
    })
    hits.extend(compounds)

# 4. Evaluate drug-likeness
drug_candidates = []
for compound in hits:
    props = tu.run({
        "name": "RDKit_calculate_drug_properties",
        "arguments": {"smiles": compound['smiles']}
    })
    if props['lipinski_pass']:
        drug_candidates.append(compound)
```

**See `references/tool-composition.md` for:**
- Complete workflow examples (drug discovery, genomics, clinical)
- Sequential and parallel tool composition patterns
- Output processing hooks
- Workflow best practices

## Scientific Domains

ToolUniverse supports 600+ tools across major scientific domains:

**Bioinformatics:**
- Sequence analysis, alignment, BLAST
- Gene expression (RNA-seq, DESeq2)
- Pathway enrichment (KEGG, Reactome, GO)
- Variant annotation (VEP, ClinVar)

**Cheminformatics:**
- Molecular descriptors and fingerprints
- Drug discovery and virtual screening
- ADMET prediction and drug-likeness
- Chemical databases (PubChem, ChEMBL, ZINC)

**Structural Biology:**
- Protein structure prediction (AlphaFold)
- Structure retrieval (PDB)
- Binding site detection
- Protein-protein interactions

**Proteomics:**
- Mass spectrometry analysis
- Protein databases (UniProt, STRING)
- Post-translational modifications

**Genomics:**
- Genome assembly and annotation
- Copy number variation
- Clinical genomics workflows

**Medical/Clinical:**
- Disease databases (OpenTargets, OMIM)
- Clinical trials and FDA data
- Variant classification

**See `references/domains.md` for:**
- Complete domain categorization
- Tool examples by discipline
- Cross-domain applications
- Search strategies by domain

## Reference Documentation

This skill includes comprehensive reference files that provide detailed information for specific aspects:

- **`references/installation.md`** - Installation, setup, MCP configuration, platform integration
- **`references/tool-discovery.md`** - Discovery methods, search strategies, listing tools
- **`references/tool-execution.md`** - Execution patterns, real-world examples, error handling
- **`references/tool-composition.md`** - Workflow composition, complex pipelines, parallel execution
- **`references/domains.md`** - Tool categorization by domain, use case examples
- **`references/api_reference.md`** - Python API documentation, hooks, protocols

**Workflow:** When helping with specific tasks, reference the appropriate file for detailed instructions. For example, if searching for tools, consult `references/tool-discovery.md` for search strategies.

## Example Scripts

Two executable example scripts demonstrate common use cases:

**`scripts/example_tool_search.py`** - Demonstrates all three discovery methods:
- Keyword-based search
- LLM-based search
- Domain-specific searches
- Getting detailed tool information

**`scripts/example_workflow.py`** - Complete workflow examples:
- Drug discovery pipeline (disease → targets → structures → screening → candidates)
- Genomics analysis (expression data → differential analysis → pathways)

Run examples to understand typical usage patterns and workflow composition.

## Best Practices

1. **Tool Discovery:**
   - Start with broad searches, then refine based on results
   - Use `Tool_Finder_Keyword` for fast searches with known terms
   - Use `Tool_Finder_LLM` for complex semantic queries
   - Set appropriate `limit` parameter (default: 10)

2. **Tool Execution:**
   - Always verify tool parameters before execution
   - Implement error handling for production workflows
   - Validate input data formats (SMILES, UniProt IDs, gene symbols)
   - Check result types and structures

3. **Workflow Composition:**
   - Test each step individually before composing full workflows
   - Implement checkpointing for long workflows
   - Consider rate limits for remote APIs
   - Use parallel execution when tools are independent

4. **Integration:**
   - Initialize ToolUniverse once and reuse the instance
   - Call `load_tools()` once at startup
   - Cache frequently used tool information
   - Enable logging for debugging

## Key Terminology

- **Tool**: A scientific resource (model, dataset, API, package) accessible through ToolUniverse
- **Tool Discovery**: Finding relevant tools using search methods (Finder, LLM, Keyword)
- **Tool Execution**: Running a tool with specific arguments via `tu.run()`
- **Tool Composition**: Chaining multiple tools for multi-step workflows
- **MCP**: Model Context Protocol for integration with Claude Desktop/Code
- **AI-Tool Interaction Protocol**: Standardized interface for LLM-tool communication

## Resources

- **Official Website**: https://aiscientist.tools
- **GitHub**: https://github.com/mims-harvard/ToolUniverse
- **Documentation**: https://zitniklab.hms.harvard.edu/ToolUniverse/
- **Installation**: `uv pip install tooluniverse`
- **MCP Server**: `tooluniverse-smcp`
