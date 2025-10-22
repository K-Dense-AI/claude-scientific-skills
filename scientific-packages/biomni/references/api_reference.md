# Biomni API Reference

Comprehensive API documentation for the biomni framework.

## A1 Agent Class

The A1 class is the primary interface for interacting with biomni.

### Initialization

```python
from biomni.agent import A1

agent = A1(
    path: str,              # Path to data lake directory
    llm: str,               # LLM model identifier
    verbose: bool = True,   # Enable verbose logging
    mcp_config: str = None  # Path to MCP server configuration
)
```

**Parameters:**

- **`path`** (str, required) - Directory path for biomni data lake (~11GB). Data is automatically downloaded on first use if not present.

- **`llm`** (str, required) - LLM model identifier. Options include:
  - `'claude-sonnet-4-20250514'` - Recommended for balanced performance
  - `'claude-opus-4-20250514'` - Maximum capability
  - `'gpt-4'`, `'gpt-4-turbo'` - OpenAI models
  - `'gemini-2.0-flash-exp'` - Google Gemini
  - `'llama-3.3-70b-versatile'` - Via Groq
  - Custom model endpoints via provider configuration

- **`verbose`** (bool, optional, default=True) - Enable detailed logging of agent reasoning, tool use, and code execution.

- **`mcp_config`** (str, optional) - Path to MCP (Model Context Protocol) server configuration file for external tool integration.

**Example:**
```python
# Basic initialization
agent = A1(path='./biomni_data', llm='claude-sonnet-4-20250514')

# With MCP integration
agent = A1(
    path='./biomni_data',
    llm='claude-sonnet-4-20250514',
    mcp_config='./.biomni/mcp_config.json'
)
```

### Core Methods

#### `go(query: str) -> str`

Execute a biomedical research task autonomously.

```python
result = agent.go(query: str)
```

**Parameters:**
- **`query`** (str) - Natural language description of the biomedical task to execute

**Returns:**
- **`str`** - Final answer or analysis result from the agent

**Behavior:**
1. Decomposes query into executable sub-tasks
2. Retrieves relevant knowledge from integrated databases
3. Generates and executes Python code for analysis
4. Iterates on results until task completion
5. Returns final synthesized answer

**Example:**
```python
result = agent.go("""
Identify genes associated with Alzheimer's disease from GWAS data.
Perform pathway enrichment analysis on top hits.
""")
print(result)
```

#### `save_conversation_history(output_path: str, format: str = 'pdf')`

Save complete conversation history including task, reasoning, code, and results.

```python
agent.save_conversation_history(
    output_path: str,
    format: str = 'pdf'
)
```

**Parameters:**
- **`output_path`** (str) - File path for saved report
- **`format`** (str, optional, default='pdf') - Output format: `'pdf'`, `'html'`, or `'markdown'`

**Example:**
```python
agent.save_conversation_history('reports/alzheimers_gwas_analysis.pdf')
```

#### `reset()`

Reset agent state and clear conversation history.

```python
agent.reset()
```

Use when starting a new independent task to clear previous context.

**Example:**
```python
# Task 1
agent.go("Analyze dataset A")
agent.save_conversation_history("task1.pdf")

# Reset for fresh context
agent.reset()

# Task 2 - independent of Task 1
agent.go("Analyze dataset B")
```

### Configuration via default_config

Global configuration parameters accessible via `biomni.config.default_config`.

```python
from biomni.config import default_config

# LLM Configuration
default_config.llm = "claude-sonnet-4-20250514"
default_config.llm_temperature = 0.7

# Execution Parameters
default_config.timeout_seconds = 1200  # 20 minutes
default_config.max_iterations = 50     # Max reasoning loops
default_config.max_tokens = 4096       # Max tokens per LLM call

# Code Execution
default_config.enable_code_execution = True
default_config.sandbox_mode = False    # Enable for restricted execution

# Data and Caching
default_config.data_cache_dir = "./biomni_cache"
default_config.enable_caching = True
```

**Key Parameters:**

- **`timeout_seconds`** (int, default=1200) - Maximum time for task execution. Increase for complex analyses.

- **`max_iterations`** (int, default=50) - Maximum agent reasoning loops. Prevents infinite loops.

- **`enable_code_execution`** (bool, default=True) - Allow agent to execute generated code. Disable for code generation only.

- **`sandbox_mode`** (bool, default=False) - Enable sandboxed code execution (requires additional setup).

## BiomniEval1 Evaluation Framework

Framework for benchmarking agent performance on biomedical tasks.

### Initialization

```python
from biomni.eval import BiomniEval1

evaluator = BiomniEval1(
    dataset_path: str = None,  # Path to evaluation dataset
    metrics: list = None        # Evaluation metrics to compute
)
```

**Example:**
```python
evaluator = BiomniEval1()
```

### Methods

#### `evaluate(task_type: str, instance_id: str, answer: str) -> float`

Evaluate agent answer against ground truth.

```python
score = evaluator.evaluate(
    task_type: str,     # Task category
    instance_id: str,   # Specific task instance
    answer: str         # Agent-generated answer
)
```

**Parameters:**
- **`task_type`** (str) - Task category: `'crispr_design'`, `'scrna_analysis'`, `'gwas_interpretation'`, `'drug_admet'`, `'clinical_diagnosis'`
- **`instance_id`** (str) - Unique identifier for task instance from dataset
- **`answer`** (str) - Agent's answer to evaluate

**Returns:**
- **`float`** - Evaluation score (0.0 to 1.0)

**Example:**
```python
# Generate answer
result = agent.go("Design CRISPR screen for autophagy genes")

# Evaluate
score = evaluator.evaluate(
    task_type='crispr_design',
    instance_id='autophagy_001',
    answer=result
)
print(f"Score: {score:.2f}")
```

#### `load_dataset() -> dict`

Load the Biomni-Eval1 benchmark dataset.

```python
dataset = evaluator.load_dataset()
```

**Returns:**
- **`dict`** - Dictionary with task instances organized by task type

**Example:**
```python
dataset = evaluator.load_dataset()

for task_type, instances in dataset.items():
    print(f"{task_type}: {len(instances)} instances")
```

#### `run_benchmark(agent: A1, task_types: list = None) -> dict`

Run full benchmark evaluation on agent.

```python
results = evaluator.run_benchmark(
    agent: A1,
    task_types: list = None  # Specific task types or None for all
)
```

**Returns:**
- **`dict`** - Results with scores, timing, and detailed metrics per task

**Example:**
```python
results = evaluator.run_benchmark(
    agent=agent,
    task_types=['crispr_design', 'scrna_analysis']
)

print(f"Overall accuracy: {results['mean_score']:.2f}")
print(f"Average time: {results['mean_time']:.1f}s")
```

## Data Lake API

Access integrated biomedical databases programmatically.

### Gene Database Queries

```python
from biomni.data import GeneDB

gene_db = GeneDB(path='./biomni_data')

# Query gene information
gene_info = gene_db.get_gene('BRCA1')
# Returns: {'symbol': 'BRCA1', 'name': '...', 'function': '...', ...}

# Search genes by pathway
pathway_genes = gene_db.search_by_pathway('DNA repair')
# Returns: List of gene symbols in pathway

# Get gene interactions
interactions = gene_db.get_interactions('TP53')
# Returns: List of interacting genes with interaction types
```

### Protein Structure Access

```python
from biomni.data import ProteinDB

protein_db = ProteinDB(path='./biomni_data')

# Get AlphaFold structure
structure = protein_db.get_structure('P38398')  # BRCA1 UniProt ID
# Returns: Path to PDB file or structure object

# Search PDB database
pdb_entries = protein_db.search_pdb('kinase', resolution_max=2.5)
# Returns: List of PDB IDs matching criteria
```

### Clinical Data Access

```python
from biomni.data import ClinicalDB

clinical_db = ClinicalDB(path='./biomni_data')

# Query ClinVar variants
variant_info = clinical_db.get_variant('rs429358')  # APOE4 variant
# Returns: {'significance': '...', 'disease': '...', 'frequency': ...}

# Search OMIM for disease
disease_info = clinical_db.search_omim('Alzheimer')
# Returns: List of OMIM entries with gene associations
```

### Literature Search

```python
from biomni.data import LiteratureDB

lit_db = LiteratureDB(path='./biomni_data')

# Search PubMed abstracts
papers = lit_db.search('CRISPR screening cancer', max_results=10)
# Returns: List of paper dictionaries with titles, abstracts, PMIDs

# Get citations for paper
citations = lit_db.get_citations('PMID:12345678')
# Returns: List of citing papers
```

## MCP Server Integration

Extend biomni with external tools via Model Context Protocol.

### Configuration Format

Create `.biomni/mcp_config.json`:

```json
{
  "servers": {
    "fda-drugs": {
      "command": "python",
      "args": ["-m", "mcp_server_fda"],
      "env": {
        "FDA_API_KEY": "${FDA_API_KEY}"
      }
    },
    "web-search": {
      "command": "npx",
      "args": ["-y", "@modelcontextprotocol/server-brave-search"],
      "env": {
        "BRAVE_API_KEY": "${BRAVE_API_KEY}"
      }
    }
  }
}
```

### Using MCP Tools in Tasks

```python
# Initialize with MCP config
agent = A1(
    path='./data',
    llm='claude-sonnet-4-20250514',
    mcp_config='./.biomni/mcp_config.json'
)

# Agent can now use MCP tools automatically
result = agent.go("""
Search for FDA-approved drugs targeting EGFR.
Get their approval dates and indications.
""")
# Agent uses fda-drugs MCP server automatically
```

## Error Handling

Common exceptions and handling strategies:

```python
from biomni.exceptions import (
    BiomniException,
    LLMError,
    CodeExecutionError,
    DataNotFoundError,
    TimeoutError
)

try:
    result = agent.go("Complex biomedical task")
except TimeoutError:
    # Task exceeded timeout_seconds
    print("Task timed out. Consider increasing timeout.")
    default_config.timeout_seconds = 3600
except CodeExecutionError as e:
    # Generated code failed to execute
    print(f"Code execution error: {e}")
    # Review generated code in conversation history
except DataNotFoundError:
    # Required data not in data lake
    print("Data not found. Ensure data lake is downloaded.")
except LLMError as e:
    # LLM API error
    print(f"LLM error: {e}")
    # Check API keys and rate limits
```

## Best Practices

### Efficient API Usage

1. **Reuse agent instances** for related tasks to maintain context
2. **Set appropriate timeouts** based on task complexity
3. **Use caching** to avoid redundant data downloads
4. **Monitor iterations** to detect reasoning loops early

### Production Deployment

```python
from biomni.agent import A1
from biomni.config import default_config
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)

# Production settings
default_config.timeout_seconds = 3600
default_config.max_iterations = 100
default_config.sandbox_mode = True  # Enable sandboxing

# Initialize with error handling
try:
    agent = A1(path='/data/biomni', llm='claude-sonnet-4-20250514')
    result = agent.go(task_query)
    agent.save_conversation_history(f'reports/{task_id}.pdf')
except Exception as e:
    logging.error(f"Task {task_id} failed: {e}")
    # Handle failure appropriately
```

### Memory Management

For large-scale analyses:

```python
# Process datasets in chunks
chunk_results = []
for chunk in dataset_chunks:
    agent.reset()  # Clear memory between chunks
    result = agent.go(f"Analyze chunk: {chunk}")
    chunk_results.append(result)

# Combine results
final_result = combine_results(chunk_results)
```
