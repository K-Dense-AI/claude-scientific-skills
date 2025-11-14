# Tool Composition and Workflows in ToolUniverse

## Overview

ToolUniverse enables chaining multiple tools together to create complex scientific workflows. Tools can be composed sequentially or in parallel to solve multi-step research problems.

## Sequential Tool Composition

Execute tools in sequence where each tool's output feeds into the next tool.

### Basic Pattern
```python
from tooluniverse import ToolUniverse

tu = ToolUniverse()
tu.load_tools()

# Step 1: Get disease-associated targets
targets = tu.run({
    "name": "OpenTargets_get_associated_targets_by_disease_efoId",
    "arguments": {"efoId": "EFO_0000537"}  # Hypertension
})

# Step 2: For each target, get protein structure
structures = []
for target in targets[:5]:  # First 5 targets
    structure = tu.run({
        "name": "AlphaFold_get_structure",
        "arguments": {"uniprot_id": target['uniprot_id']}
    })
    structures.append(structure)

# Step 3: Analyze structures
for structure in structures:
    analysis = tu.run({
        "name": "ProteinAnalysis_calculate_properties",
        "arguments": {"structure": structure}
    })
```

## Complex Workflow Examples

### Drug Discovery Workflow

Complete workflow from disease to drug candidates:

```python
# 1. Find disease-associated targets
print("Finding disease targets...")
targets = tu.run({
    "name": "OpenTargets_get_associated_targets_by_disease_efoId",
    "arguments": {"efoId": "EFO_0000616"}  # Breast cancer
})

# 2. Get target protein sequences
print("Retrieving protein sequences...")
sequences = []
for target in targets[:10]:
    seq = tu.run({
        "name": "UniProt_get_sequence",
        "arguments": {"uniprot_id": target['uniprot_id']}
    })
    sequences.append(seq)

# 3. Predict protein structures
print("Predicting structures...")
structures = []
for seq in sequences:
    structure = tu.run({
        "name": "AlphaFold_get_structure",
        "arguments": {"sequence": seq}
    })
    structures.append(structure)

# 4. Find binding sites
print("Identifying binding sites...")
binding_sites = []
for structure in structures:
    sites = tu.run({
        "name": "Fpocket_find_binding_sites",
        "arguments": {"structure": structure}
    })
    binding_sites.append(sites)

# 5. Screen compound libraries
print("Screening compounds...")
hits = []
for site in binding_sites:
    compounds = tu.run({
        "name": "ZINC_virtual_screening",
        "arguments": {
            "binding_site": site,
            "library": "lead-like",
            "top_n": 100
        }
    })
    hits.extend(compounds)

# 6. Calculate drug-likeness
print("Evaluating drug-likeness...")
drug_candidates = []
for compound in hits:
    properties = tu.run({
        "name": "RDKit_calculate_drug_properties",
        "arguments": {"smiles": compound['smiles']}
    })
    if properties['lipinski_pass']:
        drug_candidates.append(compound)

print(f"Found {len(drug_candidates)} drug candidates")
```

### Genomics Analysis Workflow

```python
# 1. Download gene expression data
expression_data = tu.run({
    "name": "GEO_download_dataset",
    "arguments": {"geo_id": "GSE12345"}
})

# 2. Perform differential expression analysis
de_genes = tu.run({
    "name": "DESeq2_differential_expression",
    "arguments": {
        "data": expression_data,
        "condition1": "control",
        "condition2": "treated"
    }
})

# 3. Pathway enrichment analysis
pathways = tu.run({
    "name": "KEGG_pathway_enrichment",
    "arguments": {
        "gene_list": de_genes['significant_genes'],
        "organism": "hsa"
    }
})

# 4. Find relevant literature
papers = tu.run({
    "name": "PubMed_search",
    "arguments": {
        "query": f"{pathways[0]['pathway_name']} AND cancer",
        "max_results": 20
    }
})

# 5. Summarize findings
summary = tu.run({
    "name": "LLM_summarize",
    "arguments": {
        "text": papers,
        "focus": "therapeutic implications"
    }
})
```

### Clinical Genomics Workflow

```python
# 1. Load patient variants
variants = tu.run({
    "name": "VCF_parse",
    "arguments": {"vcf_file": "patient_001.vcf"}
})

# 2. Annotate variants
annotated = tu.run({
    "name": "VEP_annotate_variants",
    "arguments": {"variants": variants}
})

# 3. Filter pathogenic variants
pathogenic = tu.run({
    "name": "ClinVar_filter_pathogenic",
    "arguments": {"variants": annotated}
})

# 4. Find disease associations
diseases = tu.run({
    "name": "OMIM_disease_lookup",
    "arguments": {"genes": pathogenic['affected_genes']}
})

# 5. Generate clinical report
report = tu.run({
    "name": "Report_generator",
    "arguments": {
        "variants": pathogenic,
        "diseases": diseases,
        "format": "clinical"
    }
})
```

## Parallel Tool Execution

Execute multiple tools simultaneously when they don't depend on each other:

```python
import concurrent.futures

def run_tool(tu, tool_config):
    return tu.run(tool_config)

# Define parallel tasks
tasks = [
    {"name": "PubMed_search", "arguments": {"query": "cancer", "max_results": 10}},
    {"name": "OpenTargets_get_diseases", "arguments": {"therapeutic_area": "oncology"}},
    {"name": "ChEMBL_search_compounds", "arguments": {"target": "EGFR"}}
]

# Execute in parallel
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
    futures = [executor.submit(run_tool, tu, task) for task in tasks]
    results = [future.result() for future in concurrent.futures.as_completed(futures)]
```

## Output Processing Hooks

ToolUniverse supports post-processing hooks for:
- Summarization
- File saving
- Data transformation
- Visualization

```python
# Example: Save results to file
result = tu.run({
    "name": "some_tool",
    "arguments": {"param": "value"}
},
hooks={
    "save_to_file": {"filename": "results.json"},
    "summarize": {"format": "brief"}
})
```

## Best Practices

1. **Error Handling**: Implement try-except blocks for each tool in workflow
2. **Data Validation**: Verify output from each step before passing to next tool
3. **Checkpointing**: Save intermediate results for long workflows
4. **Logging**: Track progress through complex workflows
5. **Resource Management**: Consider rate limits and computational resources
6. **Modularity**: Break complex workflows into reusable functions
7. **Testing**: Test each step individually before composing full workflow
