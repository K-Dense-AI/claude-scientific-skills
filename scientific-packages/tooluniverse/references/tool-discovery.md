# Tool Discovery in ToolUniverse

## Overview

ToolUniverse provides multiple methods to discover and search through 600+ scientific tools using natural language, keywords, or embeddings.

## Discovery Methods

### 1. Tool_Finder (Embedding-Based Search)

Uses semantic embeddings to find relevant tools. **Requires GPU** for optimal performance.

```python
from tooluniverse import ToolUniverse

tu = ToolUniverse()
tu.load_tools()

# Search by natural language description
tools = tu.run({
    "name": "Tool_Finder",
    "arguments": {
        "description": "protein structure prediction",
        "limit": 10
    }
})

print(tools)
```

**When to use:**
- Natural language queries
- Semantic similarity search
- When GPU is available

### 2. Tool_Finder_LLM (LLM-Based Search)

Alternative to embedding-based search that uses LLM reasoning. **No GPU required**.

```python
tools = tu.run({
    "name": "Tool_Finder_LLM",
    "arguments": {
        "description": "Find tools for analyzing gene expression data",
        "limit": 10
    }
})
```

**When to use:**
- When GPU is not available
- Complex queries requiring reasoning
- Semantic understanding needed

### 3. Tool_Finder_Keyword (Keyword Search)

Fast keyword-based search through tool names and descriptions.

```python
tools = tu.run({
    "name": "Tool_Finder_Keyword",
    "arguments": {
        "description": "disease target associations",
        "limit": 10
    }
})
```

**When to use:**
- Fast searches
- Known keywords
- Exact term matching

## Listing Available Tools

### List All Tools
```python
all_tools = tu.list_tools()
print(f"Total tools available: {len(all_tools)}")
```

### List Tools with Limit
```python
tools = tu.list_tools(limit=20)
for tool in tools:
    print(f"{tool['name']}: {tool['description']}")
```

## Tool Information

### Get Tool Details
```python
# After finding a tool, inspect its details
tool_info = tu.get_tool_info("OpenTargets_get_associated_targets_by_disease_efoId")
print(tool_info)
```

## Search Strategies

### By Domain
Use domain-specific keywords:
- Bioinformatics: "sequence alignment", "genomics", "RNA-seq"
- Cheminformatics: "molecular dynamics", "drug design", "SMILES"
- Machine Learning: "classification", "prediction", "neural network"
- Structural Biology: "protein structure", "PDB", "crystallography"

### By Functionality
Search by what you want to accomplish:
- "Find disease-gene associations"
- "Predict protein interactions"
- "Analyze clinical trial data"
- "Generate molecular descriptors"

### By Data Source
Search for specific databases or APIs:
- "OpenTargets", "PubChem", "UniProt"
- "AlphaFold", "ChEMBL", "PDB"
- "KEGG", "Reactome", "STRING"

## Best Practices

1. **Start Broad**: Begin with general terms, then refine
2. **Use Multiple Methods**: Try different discovery methods if results aren't satisfactory
3. **Set Appropriate Limits**: Use `limit` parameter to control result size (default: 10)
4. **Check Tool Descriptions**: Review returned tool descriptions to verify relevance
5. **Iterate**: Refine search terms based on initial results
