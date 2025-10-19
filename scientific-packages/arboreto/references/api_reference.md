# Arboreto API Reference

This document provides comprehensive API documentation for the arboreto package, a Python library for gene regulatory network (GRN) inference.

## Overview

Arboreto enables inference of gene regulatory networks from expression data using machine learning algorithms. It supports distributed computing via Dask for scalability from single machines to multi-node clusters.

**Current Version:** 0.1.5
**GitHub:** https://github.com/tmoerman/arboreto
**License:** BSD 3-Clause

## Core Algorithms

### GRNBoost2

The flagship algorithm for fast gene regulatory network inference using stochastic gradient boosting.

**Function:** `arboreto.algo.grnboost2()`

**Parameters:**
- `expression_data` (pandas.DataFrame or numpy.ndarray): Expression matrix where rows are observations (cells/samples) and columns are genes. Required.
- `gene_names` (list, optional): List of gene names matching column order. If None, uses DataFrame column names.
- `tf_names` (list, optional): List of transcription factor names to consider as regulators. If None, all genes are considered potential regulators.
- `seed` (int, optional): Random seed for reproducibility. Recommended when consistent results are needed across runs.
- `client_or_address` (dask.distributed.Client or str, optional): Custom Dask client or scheduler address for distributed computing. If None, creates a default local client.
- `verbose` (bool, optional): Enable verbose output for debugging.

**Returns:**
- pandas.DataFrame with columns `['TF', 'target', 'importance']` representing inferred regulatory links. Each row represents a regulatory relationship with an importance score.

**Algorithm Details:**
- Uses stochastic gradient boosting with early-stopping regularization
- Much faster than GENIE3, especially for large datasets (tens of thousands of observations)
- Extracts important features from trained regression models to identify regulatory relationships
- Recommended as the default choice for most use cases

**Example:**
```python
from arboreto.algo import grnboost2
import pandas as pd

# Load expression data
expression_matrix = pd.read_csv('expression_data.tsv', sep='\t')
tf_list = ['TF1', 'TF2', 'TF3']  # Optional: specify TFs

# Run inference
network = grnboost2(
    expression_data=expression_matrix,
    tf_names=tf_list,
    seed=42  # For reproducibility
)

# Save results
network.to_csv('output_network.tsv', sep='\t', index=False)
```

### GENIE3

Classical gene regulatory network inference using Random Forest regression.

**Function:** `arboreto.algo.genie3()`

**Parameters:**
Same as GRNBoost2 (see above).

**Returns:**
Same format as GRNBoost2 (see above).

**Algorithm Details:**
- Uses Random Forest or ExtraTrees regression models
- Blueprint for multiple regression GRN inference strategy
- More computationally expensive than GRNBoost2
- Better suited for smaller datasets or when maximum accuracy is needed

**When to Use GENIE3 vs GRNBoost2:**
- **Use GRNBoost2:** For large datasets, faster results, or when computational resources are limited
- **Use GENIE3:** For smaller datasets, when following established protocols, or for comparison with published results

## Module Structure

### arboreto.algo

Primary module for typical users. Contains high-level inference functions.

**Main Functions:**
- `grnboost2()` - Fast GRN inference using gradient boosting
- `genie3()` - Classical GRN inference using Random Forest

### arboreto.core

Advanced module for power users. Contains low-level framework components for custom implementations.

**Use cases:**
- Custom inference pipelines
- Algorithm modifications
- Performance tuning

### arboreto.utils

Utility functions for common data processing tasks.

**Key Functions:**
- `load_tf_names(filename)` - Load transcription factor names from file
  - Reads a text file with one TF name per line
  - Returns a list of TF names
  - Example: `tf_names = load_tf_names('transcription_factors.txt')`

## Data Format Requirements

### Input Format

**Expression Matrix:**
- **Format:** pandas DataFrame or numpy ndarray
- **Orientation:** Rows = observations (cells/samples), Columns = genes
- **Convention:** Follows scikit-learn format
- **Gene Names:** Column names (DataFrame) or separate `gene_names` parameter
- **Data Type:** Numeric (float or int)

**Common Mistake:** If data is transposed (genes as rows), use pandas to transpose:
```python
expression_df = pd.read_csv('data.tsv', sep='\t', index_col=0).T
```

**Transcription Factor List:**
- **Format:** Python list of strings or text file (one TF per line)
- **Optional:** If not provided, all genes are considered potential regulators
- **Example:** `['Sox2', 'Oct4', 'Nanog']`

### Output Format

**Network DataFrame:**
- **Columns:**
  - `TF` (str): Transcription factor (regulator) gene name
  - `target` (str): Target gene name
  - `importance` (float): Importance score of the regulatory relationship
- **Interpretation:** Higher importance scores indicate stronger predicted regulatory relationships
- **Sorting:** Typically sorted by importance (descending) for prioritization

**Example Output:**
```
TF      target    importance
Sox2    Gene1     15.234
Oct4    Gene1     12.456
Sox2    Gene2     8.901
```

## Distributed Computing with Dask

### Local Execution (Default)

Arboreto automatically creates a local Dask client if none is provided:

```python
network = grnboost2(expression_data=expr_matrix, tf_names=tf_list)
```

### Custom Local Cluster

For better control over resources or multiple inferences:

```python
from dask.distributed import Client, LocalCluster

# Configure cluster
cluster = LocalCluster(
    n_workers=4,
    threads_per_worker=2,
    memory_limit='4GB'
)
client = Client(cluster)

# Run inference
network = grnboost2(
    expression_data=expr_matrix,
    tf_names=tf_list,
    client_or_address=client
)

# Clean up
client.close()
cluster.close()
```

### Distributed Cluster

For multi-node computation:

**On scheduler node:**
```bash
dask-scheduler --no-bokeh  # Use --no-bokeh to avoid Bokeh errors
```

**On worker nodes:**
```bash
dask-worker scheduler-address:8786 --local-dir /tmp
```

**In Python script:**
```python
from dask.distributed import Client

client = Client('scheduler-address:8786')
network = grnboost2(
    expression_data=expr_matrix,
    tf_names=tf_list,
    client_or_address=client
)
```

### Dask Dashboard

Monitor computation progress via the Dask dashboard:

```python
from dask.distributed import Client, LocalCluster

cluster = LocalCluster(diagnostics_port=8787)
client = Client(cluster)

# Dashboard available at: http://localhost:8787
```

## Reproducibility

To ensure reproducible results across runs:

```python
network = grnboost2(
    expression_data=expr_matrix,
    tf_names=tf_list,
    seed=42  # Fixed seed ensures identical results
)
```

**Note:** Without a seed parameter, results may vary slightly between runs due to randomness in the algorithms.

## Performance Considerations

### Memory Management

- Expression matrices should fit in memory (RAM)
- For very large datasets, consider:
  - Using a machine with more RAM
  - Distributing across multiple nodes
  - Preprocessing to reduce dimensionality

### Worker Configuration

- **Local execution:** Number of workers = number of CPU cores (default)
- **Custom cluster:** Balance workers and threads based on available resources
- **Distributed execution:** Ensure adequate `local_dir` space on worker nodes

### Algorithm Choice

- **GRNBoost2:** ~10-100x faster than GENIE3 for large datasets
- **GENIE3:** More established but slower, better for small datasets (<10k observations)

## Integration with pySCENIC

Arboreto is a core component of the pySCENIC pipeline for single-cell RNA sequencing analysis:

1. **GRN Inference (Arboreto):** Infer regulatory networks using GRNBoost2
2. **Regulon Prediction:** Prune network and identify regulons
3. **Cell Type Identification:** Score regulons across cells

For pySCENIC workflows, arboreto is typically used in the first step to generate the initial regulatory network.

## Common Issues and Solutions

See the main SKILL.md for troubleshooting guidance.
