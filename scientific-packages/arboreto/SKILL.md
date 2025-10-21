---
name: arboreto
description: "Gene regulatory network inference with GRNBoost2/GENIE3 algorithms. Infer TF-target relationships from expression data, scalable with Dask, for scRNA-seq and GRN analysis."
---

# Arboreto - Gene Regulatory Network Inference

## Overview

Arboreto is a Python library for inferring gene regulatory networks (GRNs) from gene expression data using machine learning algorithms. It enables scalable GRN inference from single machines to multi-node clusters using Dask for distributed computing. The skill provides comprehensive support for both GRNBoost2 (fast gradient boosting) and GENIE3 (Random Forest) algorithms.

## When to Use This Skill

This skill should be used when:
- Inferring regulatory relationships between genes from expression data
- Analyzing single-cell or bulk RNA-seq data to identify transcription factor targets
- Building the GRN inference component of a pySCENIC pipeline
- Comparing GRNBoost2 and GENIE3 algorithm performance
- Setting up distributed computing for large-scale genomic analyses
- Troubleshooting arboreto installation or runtime issues

## Core Capabilities

### 1. Basic GRN Inference

For standard gene regulatory network inference tasks:

**Key considerations:**
- Expression data format: Rows = observations (cells/samples), Columns = genes
- If data has genes as rows, transpose it first: `expression_df.T`
- Always include `seed` parameter for reproducible results
- Transcription factor list is optional but recommended for focused analysis

**Typical workflow:**
```python
import pandas as pd
from arboreto.algo import grnboost2
from arboreto.utils import load_tf_names

# Load expression data (ensure correct orientation)
expression_data = pd.read_csv('expression_data.tsv', sep='\t', index_col=0)

# Optional: Load TF names
tf_names = load_tf_names('transcription_factors.txt')

# Run inference
network = grnboost2(
    expression_data=expression_data,
    tf_names=tf_names,
    seed=42  # For reproducibility
)

# Save results
network.to_csv('network_output.tsv', sep='\t', index=False)
```

**Output format:**
- DataFrame with columns: `['TF', 'target', 'importance']`
- Higher importance scores indicate stronger predicted regulatory relationships
- Typically sorted by importance (descending)

**Multiprocessing requirement:**
All arboreto code must include `if __name__ == '__main__':` protection due to Dask's multiprocessing requirements:

```python
if __name__ == '__main__':
    # Arboreto code goes here
    network = grnboost2(expression_data=expr_data, seed=42)
```

### 2. Algorithm Selection

**GRNBoost2 (Recommended for most cases):**
- ~10-100x faster than GENIE3
- Uses stochastic gradient boosting with early-stopping
- Best for: Large datasets (>10k observations), time-sensitive analyses
- Function: `arboreto.algo.grnboost2()`

**GENIE3:**
- Uses Random Forest regression
- More established, classical approach
- Best for: Small datasets, methodological comparisons, reproducing published results
- Function: `arboreto.algo.genie3()`

**When to compare both algorithms:**
Use the provided `compare_algorithms.py` script when:
- Validating results for critical analyses
- Benchmarking performance on new datasets
- Publishing research requiring methodological comparisons

### 3. Distributed Computing

**Local execution (default):**
Arboreto automatically creates a local Dask client. No configuration needed:
```python
network = grnboost2(expression_data=expr_data)
```

**Custom local cluster (recommended for better control):**
```python
from dask.distributed import Client, LocalCluster

# Configure cluster
cluster = LocalCluster(
    n_workers=4,
    threads_per_worker=2,
    memory_limit='4GB',
    diagnostics_port=8787  # Dashboard at http://localhost:8787
)
client = Client(cluster)

# Run inference
network = grnboost2(
    expression_data=expr_data,
    client_or_address=client
)

# Clean up
client.close()
cluster.close()
```

**Distributed cluster (multi-node):**
On scheduler node:
```bash
dask-scheduler --no-bokeh
```

On worker nodes:
```bash
dask-worker scheduler-address:8786 --local-dir /tmp
```

In Python:
```python
from dask.distributed import Client

client = Client('scheduler-address:8786')
network = grnboost2(expression_data=expr_data, client_or_address=client)
```

### 4. Data Preparation

**Common data format issues:**

1. **Transposed data** (genes as rows instead of columns):
```python
# If genes are rows, transpose
expression_data = pd.read_csv('data.tsv', sep='\t', index_col=0).T
```

2. **Missing gene names:**
```python
# Provide gene names if using numpy array
network = grnboost2(
    expression_data=expr_array,
    gene_names=['Gene1', 'Gene2', 'Gene3', ...],
    seed=42
)
```

3. **Transcription factor specification:**
```python
# Option 1: Python list
tf_names = ['Sox2', 'Oct4', 'Nanog', 'Klf4']

# Option 2: Load from file (one TF per line)
from arboreto.utils import load_tf_names
tf_names = load_tf_names('tf_names.txt')
```

### 5. Reproducibility

Always specify a seed for consistent results:
```python
network = grnboost2(expression_data=expr_data, seed=42)
```

Without a seed, results will vary between runs due to algorithm randomness.

### 6. Result Interpretation

**Understanding the output:**
- `TF`: Transcription factor (regulator) gene
- `target`: Target gene being regulated
- `importance`: Strength of predicted regulatory relationship

**Typical post-processing:**
```python
# Filter by importance threshold
high_confidence = network[network['importance'] > 10]

# Get top N predictions
top_predictions = network.head(1000)

# Find all targets of a specific TF
sox2_targets = network[network['TF'] == 'Sox2']

# Count regulations per TF
tf_counts = network['TF'].value_counts()
```

## Installation

**Recommended (via conda):**
```bash
conda install -c bioconda arboreto
```

**Via pip:**
```bash
pip install arboreto
```

**From source:**
```bash
git clone https://github.com/tmoerman/arboreto.git
cd arboreto
pip install .
```

**Dependencies:**
- pandas
- numpy
- scikit-learn
- scipy
- dask
- distributed

## Troubleshooting

### Issue: Bokeh error when launching Dask scheduler

**Error:** `TypeError: got an unexpected keyword argument 'host'`

**Solutions:**
- Use `dask-scheduler --no-bokeh` to disable Bokeh
- Upgrade to Dask distributed >= 0.20.0

### Issue: Workers not connecting to scheduler

**Symptoms:** Worker processes start but fail to establish connections

**Solutions:**
- Remove `dask-worker-space` directory before restarting workers
- Specify adequate `local_dir` when creating cluster:
```python
cluster = LocalCluster(
    worker_kwargs={'local_dir': '/tmp'}
)
```

### Issue: Memory errors with large datasets

**Solutions:**
- Increase worker memory limits: `memory_limit='8GB'`
- Distribute across more nodes
- Reduce dataset size through preprocessing (e.g., feature selection)
- Ensure expression matrix fits in available RAM

### Issue: Inconsistent results across runs

**Solution:** Always specify a `seed` parameter:
```python
network = grnboost2(expression_data=expr_data, seed=42)
```

### Issue: Import errors or missing dependencies

**Solution:** Use conda installation to handle numerical library dependencies:
```bash
conda create --name arboreto-env
conda activate arboreto-env
conda install -c bioconda arboreto
```

## Provided Scripts

This skill includes ready-to-use scripts for common workflows:

### scripts/basic_grn_inference.py

Command-line tool for standard GRN inference workflow.

**Usage:**
```bash
python scripts/basic_grn_inference.py expression_data.tsv \
    -t tf_names.txt \
    -o network.tsv \
    -s 42 \
    --transpose  # if genes are rows
```

**Features:**
- Automatic data loading and validation
- Optional TF list specification
- Configurable output format
- Data transposition support
- Summary statistics

### scripts/distributed_inference.py

GRN inference with custom Dask cluster configuration.

**Usage:**
```bash
python scripts/distributed_inference.py expression_data.tsv \
    -t tf_names.txt \
    -w 8 \
    -m 4GB \
    --threads 2 \
    --dashboard-port 8787
```

**Features:**
- Configurable worker count and memory limits
- Dask dashboard integration
- Thread configuration
- Resource monitoring

### scripts/compare_algorithms.py

Compare GRNBoost2 and GENIE3 side-by-side.

**Usage:**
```bash
python scripts/compare_algorithms.py expression_data.tsv \
    -t tf_names.txt \
    --top-n 100
```

**Features:**
- Runtime comparison
- Network statistics
- Prediction overlap analysis
- Top prediction comparison

## Reference Documentation

Detailed API documentation is available in [references/api_reference.md](references/api_reference.md), including:
- Complete parameter descriptions for all functions
- Data format specifications
- Distributed computing configuration
- Performance optimization tips
- Integration with pySCENIC
- Comprehensive examples

Load this reference when:
- Working with advanced Dask configurations
- Troubleshooting complex deployment scenarios
- Understanding algorithm internals
- Optimizing performance for specific use cases

## Integration with pySCENIC

Arboreto is the first step in the pySCENIC single-cell analysis pipeline:

1. **GRN Inference (arboreto)** ← This skill
   - Input: Expression matrix
   - Output: Regulatory network

2. **Regulon Prediction (pySCENIC)**
   - Input: Network from arboreto
   - Output: Refined regulons

3. **Cell Type Identification (pySCENIC)**
   - Input: Regulons
   - Output: Cell type scores

When working with pySCENIC, use arboreto to generate the initial network, then pass results to the pySCENIC pipeline.

## Best Practices

1. **Always use seed parameter** for reproducible research
2. **Validate data orientation** (rows = observations, columns = genes)
3. **Specify TF list** when known to focus inference and improve speed
4. **Monitor with Dask dashboard** for distributed computing
5. **Save intermediate results** to avoid re-running long computations
6. **Filter results** by importance threshold for downstream analysis
7. **Use GRNBoost2 by default** unless specifically requiring GENIE3
8. **Include multiprocessing guard** (`if __name__ == '__main__':`) in all scripts

## Quick Reference

**Basic inference:**
```python
from arboreto.algo import grnboost2
network = grnboost2(expression_data=expr_df, seed=42)
```

**With TF specification:**
```python
network = grnboost2(expression_data=expr_df, tf_names=tf_list, seed=42)
```

**With custom Dask client:**
```python
from dask.distributed import Client, LocalCluster
cluster = LocalCluster(n_workers=4)
client = Client(cluster)
network = grnboost2(expression_data=expr_df, client_or_address=client, seed=42)
client.close()
cluster.close()
```

**Load TF names:**
```python
from arboreto.utils import load_tf_names
tf_names = load_tf_names('transcription_factors.txt')
```

**Transpose data:**
```python
expression_df = pd.read_csv('data.tsv', sep='\t', index_col=0).T
```
