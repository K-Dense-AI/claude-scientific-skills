# Biomni API Reference

This document provides comprehensive API documentation for the Biomni biomedical AI agent system.

## Core Classes

### A1 Agent

The primary agent class for executing biomedical research tasks.

#### Initialization

```python
from biomni.agent import A1

agent = A1(
    path='./data',              # Path to biomedical knowledge base
    llm='claude-sonnet-4-20250514',  # LLM model identifier
    timeout=None,               # Optional timeout in seconds
    verbose=True               # Enable detailed logging
)
```

**Parameters:**

- `path` (str, required): Directory path where the biomedical knowledge base is stored or will be downloaded. First-time initialization will download ~11GB of data.
- `llm` (str, optional): LLM model identifier. Defaults to the value in `default_config.llm`. Supports multiple providers (see LLM Providers section).
- `timeout` (int, optional): Maximum execution time in seconds for agent operations. Overrides `default_config.timeout_seconds`.
- `verbose` (bool, optional): Enable verbose logging for debugging. Default: True.

**Returns:** A1 agent instance ready for task execution.

#### Methods

##### `go(task_description: str) -> None`

Execute a biomedical research task autonomously.

```python
agent.go("Analyze this scRNA-seq dataset and identify cell types")
```

**Parameters:**
- `task_description` (str, required): Natural language description of the biomedical task to execute. Be specific about:
  - Data location and format
  - Desired analysis or output
  - Any specific methods or parameters
  - Expected results format

**Behavior:**
1. Decomposes the task into executable steps
2. Retrieves relevant biomedical knowledge from the data lake
3. Generates and executes Python/R code
4. Provides results and visualizations
5. Handles errors and retries with refinement

**Notes:**
- Executes code with system privileges - use in sandboxed environments
- Long-running tasks may require timeout adjustments
- Intermediate results are displayed during execution

##### `save_conversation_history(output_path: str, format: str = 'pdf') -> None`

Export conversation history and execution trace as a formatted report.

```python
agent.save_conversation_history(
    output_path='./reports/analysis_log.pdf',
    format='pdf'
)
```

**Parameters:**
- `output_path` (str, required): File path for the output report
- `format` (str, optional): Output format. Options: 'pdf', 'markdown'. Default: 'pdf'

**Requirements:**
- For PDF: Install one of: WeasyPrint, markdown2pdf, or Pandoc
  ```bash
  pip install weasyprint  # Recommended
  # or
  pip install markdown2pdf
  # or install Pandoc system-wide
  ```

**Report Contents:**
- Task description and parameters
- Retrieved biomedical knowledge
- Generated code with execution traces
- Results, visualizations, and outputs
- Timestamps and execution metadata

##### `add_mcp(config_path: str) -> None`

Add Model Context Protocol (MCP) tools to extend agent capabilities.

```python
agent.add_mcp(config_path='./mcp_tools_config.json')
```

**Parameters:**
- `config_path` (str, required): Path to MCP configuration JSON file

**MCP Configuration Format:**
```json
{
  "tools": [
    {
      "name": "tool_name",
      "endpoint": "http://localhost:8000/tool",
      "description": "Tool description for LLM",
      "parameters": {
        "param1": "string",
        "param2": "integer"
      }
    }
  ]
}
```

**Use Cases:**
- Connect to laboratory information systems
- Integrate proprietary databases
- Access specialized computational resources
- Link to institutional data repositories

## Configuration

### default_config

Global configuration object for Biomni settings.

```python
from biomni.config import default_config
```

#### Attributes

##### `llm: str`

Default LLM model identifier for all agent instances.

```python
default_config.llm = "claude-sonnet-4-20250514"
```

**Supported Models:**

**Anthropic:**
- `claude-sonnet-4-20250514` (Recommended)
- `claude-opus-4-20250514`
- `claude-3-5-sonnet-20241022`
- `claude-3-opus-20240229`

**OpenAI:**
- `gpt-4o`
- `gpt-4`
- `gpt-4-turbo`
- `gpt-3.5-turbo`

**Azure OpenAI:**
- `azure/gpt-4`
- `azure/<deployment-name>`

**Google Gemini:**
- `gemini/gemini-pro`
- `gemini/gemini-1.5-pro`

**Groq:**
- `groq/llama-3.1-70b-versatile`
- `groq/mixtral-8x7b-32768`

**Ollama (Local):**
- `ollama/llama3`
- `ollama/mistral`
- `ollama/<model-name>`

**AWS Bedrock:**
- `bedrock/anthropic.claude-v2`
- `bedrock/anthropic.claude-3-sonnet`

**Custom/Biomni-R0:**
- `openai/biomni-r0` (requires local SGLang deployment)

##### `timeout_seconds: int`

Default timeout for agent operations in seconds.

```python
default_config.timeout_seconds = 1200  # 20 minutes
```

**Recommended Values:**
- Simple tasks (QC, basic analysis): 300-600 seconds
- Medium tasks (differential expression, clustering): 600-1200 seconds
- Complex tasks (full pipelines, ML models): 1200-3600 seconds
- Very complex tasks: 3600+ seconds

##### `data_path: str`

Default path to biomedical knowledge base.

```python
default_config.data_path = "/path/to/biomni/data"
```

**Storage Requirements:**
- Initial download: ~11GB
- Extracted size: ~15GB
- Additional working space: ~5-10GB recommended

##### `api_base: str`

Custom API endpoint for LLM providers (advanced usage).

```python
# For local Biomni-R0 deployment
default_config.api_base = "http://localhost:30000/v1"

# For custom OpenAI-compatible endpoints
default_config.api_base = "https://your-endpoint.com/v1"
```

##### `max_retries: int`

Number of retry attempts for failed operations.

```python
default_config.max_retries = 3
```

#### Methods

##### `reset() -> None`

Reset all configuration values to system defaults.

```python
default_config.reset()
```

## Database Query System

Biomni includes a retrieval-augmented generation (RAG) system for querying the biomedical knowledge base.

### Query Functions

#### `query_genes(query: str, top_k: int = 10) -> List[Dict]`

Query gene information from integrated databases.

```python
from biomni.database import query_genes

results = query_genes(
    query="genes involved in p53 pathway",
    top_k=20
)
```

**Parameters:**
- `query` (str): Natural language or gene identifier query
- `top_k` (int): Number of results to return

**Returns:** List of dictionaries containing:
- `gene_symbol`: Official gene symbol
- `gene_name`: Full gene name
- `description`: Functional description
- `pathways`: Associated biological pathways
- `go_terms`: Gene Ontology annotations
- `diseases`: Associated diseases
- `similarity_score`: Relevance score (0-1)

#### `query_proteins(query: str, top_k: int = 10) -> List[Dict]`

Query protein information from UniProt and other sources.

```python
from biomni.database import query_proteins

results = query_proteins(
    query="kinase proteins in cell cycle",
    top_k=15
)
```

**Returns:** List of dictionaries with protein metadata:
- `uniprot_id`: UniProt accession
- `protein_name`: Protein name
- `function`: Functional annotation
- `domains`: Protein domains
- `subcellular_location`: Cellular localization
- `similarity_score`: Relevance score

#### `query_drugs(query: str, top_k: int = 10) -> List[Dict]`

Query drug and compound information.

```python
from biomni.database import query_drugs

results = query_drugs(
    query="FDA approved cancer drugs targeting EGFR",
    top_k=10
)
```

**Returns:** Drug information including:
- `drug_name`: Common name
- `drugbank_id`: DrugBank identifier
- `indication`: Therapeutic indication
- `mechanism`: Mechanism of action
- `targets`: Molecular targets
- `approval_status`: Regulatory status
- `smiles`: Chemical structure (SMILES notation)

#### `query_diseases(query: str, top_k: int = 10) -> List[Dict]`

Query disease information from clinical databases.

```python
from biomni.database import query_diseases

results = query_diseases(
    query="autoimmune diseases affecting joints",
    top_k=10
)
```

**Returns:** Disease data:
- `disease_name`: Standard disease name
- `disease_id`: Ontology identifier
- `symptoms`: Clinical manifestations
- `associated_genes`: Genetic associations
- `prevalence`: Epidemiological data

#### `query_pathways(query: str, top_k: int = 10) -> List[Dict]`

Query biological pathways from KEGG, Reactome, and other sources.

```python
from biomni.database import query_pathways

results = query_pathways(
    query="immune response signaling pathways",
    top_k=15
)
```

**Returns:** Pathway information:
- `pathway_name`: Pathway name
- `pathway_id`: Database identifier
- `genes`: Genes in pathway
- `description`: Functional description
- `source`: Database source (KEGG, Reactome, etc.)

## Data Structures

### TaskResult

Result object returned by complex agent operations.

```python
class TaskResult:
    success: bool           # Whether task completed successfully
    output: Any            # Task output (varies by task)
    code: str             # Generated code
    execution_time: float # Execution time in seconds
    error: Optional[str]  # Error message if failed
    metadata: Dict        # Additional metadata
```

### BiomedicalEntity

Base class for biomedical entities in the knowledge base.

```python
class BiomedicalEntity:
    entity_id: str        # Unique identifier
    entity_type: str      # Type (gene, protein, drug, etc.)
    name: str            # Entity name
    description: str     # Description
    attributes: Dict     # Additional attributes
    references: List[str] # Literature references
```

## Utility Functions

### `download_data(path: str, force: bool = False) -> None`

Manually download or update the biomedical knowledge base.

```python
from biomni.utils import download_data

download_data(
    path='./data',
    force=True  # Force re-download
)
```

### `validate_environment() -> Dict[str, bool]`

Check if the environment is properly configured.

```python
from biomni.utils import validate_environment

status = validate_environment()
# Returns: {
#   'conda_env': True,
#   'api_keys': True,
#   'data_available': True,
#   'dependencies': True
# }
```

### `list_available_models() -> List[str]`

Get a list of available LLM models based on configured API keys.

```python
from biomni.utils import list_available_models

models = list_available_models()
# Returns: ['claude-sonnet-4-20250514', 'gpt-4o', ...]
```

## Error Handling

### Common Exceptions

#### `BiomniConfigError`

Raised when configuration is invalid or incomplete.

```python
from biomni.exceptions import BiomniConfigError

try:
    agent = A1(path='./data')
except BiomniConfigError as e:
    print(f"Configuration error: {e}")
```

#### `BiomniExecutionError`

Raised when code generation or execution fails.

```python
from biomni.exceptions import BiomniExecutionError

try:
    agent.go("invalid task")
except BiomniExecutionError as e:
    print(f"Execution failed: {e}")
    # Access failed code: e.code
    # Access error details: e.details
```

#### `BiomniDataError`

Raised when knowledge base or data access fails.

```python
from biomni.exceptions import BiomniDataError

try:
    results = query_genes("unknown query format")
except BiomniDataError as e:
    print(f"Data access error: {e}")
```

#### `BiomniTimeoutError`

Raised when operations exceed timeout limit.

```python
from biomni.exceptions import BiomniTimeoutError

try:
    agent.go("very complex long-running task")
except BiomniTimeoutError as e:
    print(f"Task timed out after {e.duration} seconds")
    # Partial results may be available: e.partial_results
```

## Best Practices

### Efficient Knowledge Retrieval

Pre-query databases for relevant context before complex tasks:

```python
from biomni.database import query_genes, query_pathways

# Gather relevant biological context first
genes = query_genes("cell cycle genes", top_k=50)
pathways = query_pathways("cell cycle regulation", top_k=20)

# Then execute task with enriched context
agent.go(f"""
Analyze the cell cycle progression in this dataset.
Focus on these genes: {[g['gene_symbol'] for g in genes]}
Consider these pathways: {[p['pathway_name'] for p in pathways]}
""")
```

### Error Recovery

Implement robust error handling for production workflows:

```python
from biomni.exceptions import BiomniExecutionError, BiomniTimeoutError

max_attempts = 3
for attempt in range(max_attempts):
    try:
        agent.go("complex biomedical task")
        break
    except BiomniTimeoutError:
        # Increase timeout and retry
        default_config.timeout_seconds *= 2
        print(f"Timeout, retrying with {default_config.timeout_seconds}s timeout")
    except BiomniExecutionError as e:
        # Refine task based on error
        print(f"Execution failed: {e}, refining task...")
        # Optionally modify task description
    else:
        print("Task failed after max attempts")
```

### Memory Management

For large-scale analyses, manage memory explicitly:

```python
import gc

# Process datasets in chunks
for chunk_id in range(num_chunks):
    agent.go(f"Process data chunk {chunk_id} located at data/chunk_{chunk_id}.h5ad")

    # Force garbage collection between chunks
    gc.collect()

    # Save intermediate results
    agent.save_conversation_history(f"./reports/chunk_{chunk_id}.pdf")
```

### Reproducibility

Ensure reproducible analyses by:

1. **Fixing random seeds:**
```python
agent.go("Set random seed to 42 for all analyses, then perform clustering...")
```

2. **Logging configuration:**
```python
import json
config_log = {
    'llm': default_config.llm,
    'timeout': default_config.timeout_seconds,
    'data_path': default_config.data_path,
    'timestamp': datetime.now().isoformat()
}
with open('config_log.json', 'w') as f:
    json.dump(config_log, f, indent=2)
```

3. **Saving execution traces:**
```python
# Always save detailed reports
agent.save_conversation_history('./reports/full_analysis.pdf')
```

## Performance Optimization

### Model Selection Strategy

Choose models based on task characteristics:

```python
# For exploratory, simple tasks
default_config.llm = "gpt-3.5-turbo"  # Fast, cost-effective

# For standard biomedical analyses
default_config.llm = "claude-sonnet-4-20250514"  # Recommended

# For complex reasoning and hypothesis generation
default_config.llm = "claude-opus-4-20250514"  # Highest quality

# For specialized biological reasoning
default_config.llm = "openai/biomni-r0"  # Requires local deployment
```

### Timeout Tuning

Set appropriate timeouts based on task complexity:

```python
# Quick queries and simple analyses
agent = A1(path='./data', timeout=300)

# Standard workflows
agent = A1(path='./data', timeout=1200)

# Full pipelines with ML training
agent = A1(path='./data', timeout=3600)
```

### Caching and Reuse

Reuse agent instances for multiple related tasks:

```python
# Create agent once
agent = A1(path='./data', llm='claude-sonnet-4-20250514')

# Execute multiple related tasks
tasks = [
    "Load and QC the scRNA-seq dataset",
    "Perform clustering with resolution 0.5",
    "Identify marker genes for each cluster",
    "Annotate cell types based on markers"
]

for task in tasks:
    agent.go(task)

# Save complete workflow
agent.save_conversation_history('./reports/full_workflow.pdf')
```
