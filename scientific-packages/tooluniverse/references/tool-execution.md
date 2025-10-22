# Tool Execution in ToolUniverse

## Overview

Execute individual tools through ToolUniverse's standardized interface using the `run()` method.

## Basic Tool Execution

### Standard Pattern
```python
from tooluniverse import ToolUniverse

tu = ToolUniverse()
tu.load_tools()

# Execute a tool
result = tu.run({
    "name": "tool_name_here",
    "arguments": {
        "param1": "value1",
        "param2": "value2"
    }
})

print(result)
```

## Real-World Examples

### Example 1: Disease-Target Associations (OpenTargets)
```python
# Find targets associated with hypertension
result = tu.run({
    "name": "OpenTargets_get_associated_targets_by_disease_efoId",
    "arguments": {
        "efoId": "EFO_0000537"  # Hypertension
    }
})

print(f"Found {len(result)} targets associated with hypertension")
```

### Example 2: Protein Structure Prediction
```python
# Get AlphaFold structure prediction
result = tu.run({
    "name": "AlphaFold_get_structure",
    "arguments": {
        "uniprot_id": "P12345"
    }
})
```

### Example 3: Chemical Property Calculation
```python
# Calculate molecular descriptors
result = tu.run({
    "name": "RDKit_calculate_descriptors",
    "arguments": {
        "smiles": "CCO"  # Ethanol
    }
})
```

### Example 4: Gene Expression Analysis
```python
# Analyze differential gene expression
result = tu.run({
    "name": "GeneExpression_differential_analysis",
    "arguments": {
        "dataset_id": "GSE12345",
        "condition1": "control",
        "condition2": "treatment"
    }
})
```

## Tool Execution Workflow

### 1. Discover the Tool
```python
# Find relevant tools
tools = tu.run({
    "name": "Tool_Finder_Keyword",
    "arguments": {
        "description": "pathway enrichment",
        "limit": 5
    }
})

# Review available tools
for tool in tools:
    print(f"Name: {tool['name']}")
    print(f"Description: {tool['description']}")
    print(f"Parameters: {tool['parameters']}")
    print("---")
```

### 2. Check Tool Parameters
```python
# Get detailed tool information
tool_info = tu.get_tool_info("KEGG_pathway_enrichment")
print(tool_info['parameters'])
```

### 3. Execute with Proper Arguments
```python
# Execute the tool
result = tu.run({
    "name": "KEGG_pathway_enrichment",
    "arguments": {
        "gene_list": ["TP53", "BRCA1", "EGFR"],
        "organism": "hsa"  # Homo sapiens
    }
})
```

## Handling Tool Results

### Check Result Type
```python
result = tu.run({
    "name": "some_tool",
    "arguments": {"param": "value"}
})

# Results can be various types
if isinstance(result, dict):
    print("Dictionary result")
elif isinstance(result, list):
    print(f"List with {len(result)} items")
elif isinstance(result, str):
    print("String result")
```

### Process Results
```python
# Example: Processing multiple results
results = tu.run({
    "name": "PubMed_search",
    "arguments": {
        "query": "cancer immunotherapy",
        "max_results": 10
    }
})

for idx, paper in enumerate(results, 1):
    print(f"{idx}. {paper['title']}")
    print(f"   PMID: {paper['pmid']}")
    print(f"   Authors: {', '.join(paper['authors'][:3])}")
    print()
```

## Error Handling

```python
try:
    result = tu.run({
        "name": "some_tool",
        "arguments": {"param": "value"}
    })
except Exception as e:
    print(f"Tool execution failed: {e}")
    # Check if tool exists
    # Verify parameter names and types
    # Review tool documentation
```

## Best Practices

1. **Verify Tool Parameters**: Always check required parameters before execution
2. **Start Simple**: Test with simple cases before complex workflows
3. **Handle Results Appropriately**: Check result type and structure
4. **Error Recovery**: Implement try-except blocks for production code
5. **Documentation**: Review tool descriptions for parameter requirements and output formats
6. **Rate Limiting**: Be aware of API rate limits for remote tools
7. **Data Validation**: Validate input data format (e.g., SMILES, UniProt IDs, gene symbols)
