---
name: anndata
description: "Manipulate AnnData objects for single-cell genomics. Load/save .h5ad files, manage obs/var metadata, layers, embeddings (PCA/UMAP), concatenate datasets, for scRNA-seq workflows."
---

# AnnData

## Overview

AnnData (Annotated Data) is Python's standard for storing and manipulating annotated data matrices, particularly in single-cell genomics. Work with AnnData objects for data creation, manipulation, file I/O, concatenation, and memory-efficient workflows.

## Core Capabilities

### 1. Creating and Structuring AnnData Objects

Create AnnData objects from various data sources and organize multi-dimensional annotations.

**Basic creation:**
```python
import anndata as ad
import numpy as np
from scipy.sparse import csr_matrix

# From dense or sparse arrays
counts = np.random.poisson(1, size=(100, 2000))
adata = ad.AnnData(counts)

# With sparse matrix (memory-efficient)
counts = csr_matrix(np.random.poisson(1, size=(100, 2000)), dtype=np.float32)
adata = ad.AnnData(counts)
```

**With metadata:**
```python
import pandas as pd

obs_meta = pd.DataFrame({
    'cell_type': pd.Categorical(['B', 'T', 'Monocyte'] * 33 + ['B']),
    'batch': ['batch1'] * 50 + ['batch2'] * 50
})
var_meta = pd.DataFrame({
    'gene_name': [f'Gene_{i}' for i in range(2000)],
    'highly_variable': np.random.choice([True, False], 2000)
})

adata = ad.AnnData(counts, obs=obs_meta, var=var_meta)
```

**Understanding the structure:**
- **X**: Primary data matrix (observations × variables)
- **obs**: Row (observation) annotations as DataFrame
- **var**: Column (variable) annotations as DataFrame
- **obsm**: Multi-dimensional observation annotations (e.g., PCA, UMAP coordinates)
- **varm**: Multi-dimensional variable annotations (e.g., gene loadings)
- **layers**: Alternative data matrices with same dimensions as X
- **uns**: Unstructured metadata dictionary
- **obsp/varp**: Pairwise relationship matrices (graphs)

### 2. Adding Annotations and Layers

Organize different data representations and metadata within a single object.

**Cell-level metadata (obs):**
```python
adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)
adata.obs['total_counts'] = adata.X.sum(axis=1)
adata.obs['condition'] = pd.Categorical(['control', 'treated'] * 50)
```

**Gene-level metadata (var):**
```python
adata.var['highly_variable'] = gene_variance > threshold
adata.var['chromosome'] = pd.Categorical(['chr1', 'chr2', ...])
```

**Embeddings (obsm/varm):**
```python
# Dimensionality reduction results
adata.obsm['X_pca'] = pca_coordinates  # Shape: (n_obs, n_components)
adata.obsm['X_umap'] = umap_coordinates  # Shape: (n_obs, 2)
adata.obsm['X_tsne'] = tsne_coordinates

# Gene loadings
adata.varm['PCs'] = principal_components  # Shape: (n_vars, n_components)
```

**Alternative data representations (layers):**
```python
# Store multiple versions
adata.layers['counts'] = raw_counts
adata.layers['log1p'] = np.log1p(adata.X)
adata.layers['scaled'] = (adata.X - mean) / std
```

**Unstructured metadata (uns):**
```python
# Analysis parameters
adata.uns['preprocessing'] = {
    'normalization': 'TPM',
    'min_genes': 200,
    'date': '2024-01-15'
}

# Results
adata.uns['pca'] = {'variance_ratio': variance_explained}
```

### 3. Subsetting and Views

Efficiently subset data while managing memory through views and copies.

**Subsetting operations:**
```python
# By observation/variable names
subset = adata[['Cell_1', 'Cell_10'], ['Gene_5', 'Gene_1900']]

# By boolean masks
b_cells = adata[adata.obs.cell_type == 'B']
high_quality = adata[adata.obs.n_genes > 200]

# By position
first_cells = adata[:100, :]
top_genes = adata[:, :500]

# Combined conditions
filtered = adata[
    (adata.obs.batch == 'batch1') & (adata.obs.n_genes > 200),
    adata.var.highly_variable
]
```

**Understanding views:**
- Subsetting returns **views** by default (memory-efficient, shares data with original)
- Modifying a view affects the original object
- Check with `adata.is_view`
- Convert to independent copy with `.copy()`

```python
# View (memory-efficient)
subset = adata[adata.obs.condition == 'treated']
print(subset.is_view)  # True

# Independent copy
subset_copy = adata[adata.obs.condition == 'treated'].copy()
print(subset_copy.is_view)  # False
```

### 4. File I/O and Backed Mode

Read and write data efficiently, with options for memory-limited environments.

**Writing data:**
```python
# Standard format with compression
adata.write('results.h5ad', compression='gzip')

# Alternative formats
adata.write_zarr('results.zarr')  # For cloud storage
adata.write_loom('results.loom')  # For compatibility
adata.write_csvs('results/')      # As CSV files
```

**Reading data:**
```python
# Load into memory
adata = ad.read_h5ad('results.h5ad')

# Backed mode (disk-backed, memory-efficient)
adata = ad.read_h5ad('large_file.h5ad', backed='r')
print(adata.isbacked)  # True
print(adata.filename)  # Path to file

# Close file connection when done
adata.file.close()
```

**Reading from other formats:**
```python
# 10X format
adata = ad.read_mtx('matrix.mtx')

# CSV
adata = ad.read_csv('data.csv')

# Loom
adata = ad.read_loom('data.loom')
```

**Working with backed mode:**
```python
# Read in backed mode for large files
adata = ad.read_h5ad('large_dataset.h5ad', backed='r')

# Process in chunks
for chunk in adata.chunk_X(chunk_size=1000):
    result = process_chunk(chunk)

# Load to memory if needed
adata_memory = adata.to_memory()
```

### 5. Concatenating Multiple Datasets

Combine multiple AnnData objects with control over how data is merged.

**Basic concatenation:**
```python
# Concatenate observations (most common)
combined = ad.concat([adata1, adata2, adata3], axis=0)

# Concatenate variables (rare)
combined = ad.concat([adata1, adata2], axis=1)
```

**Join strategies:**
```python
# Inner join: only shared variables (no missing data)
combined = ad.concat([adata1, adata2], join='inner')

# Outer join: all variables (fills missing with 0)
combined = ad.concat([adata1, adata2], join='outer')
```

**Tracking data sources:**
```python
# Add source labels
combined = ad.concat(
    [adata1, adata2, adata3],
    label='dataset',
    keys=['exp1', 'exp2', 'exp3']
)
# Creates combined.obs['dataset'] with values 'exp1', 'exp2', 'exp3'

# Make duplicate indices unique
combined = ad.concat(
    [adata1, adata2],
    keys=['batch1', 'batch2'],
    index_unique='-'
)
# Cell names become: Cell_0-batch1, Cell_0-batch2, etc.
```

**Merge strategies for metadata:**
```python
# merge=None: exclude variable annotations (default)
combined = ad.concat([adata1, adata2], merge=None)

# merge='same': keep only identical annotations
combined = ad.concat([adata1, adata2], merge='same')

# merge='first': use first occurrence
combined = ad.concat([adata1, adata2], merge='first')

# merge='unique': keep annotations with single value
combined = ad.concat([adata1, adata2], merge='unique')
```

**Complete example:**
```python
# Load batches
batch1 = ad.read_h5ad('batch1.h5ad')
batch2 = ad.read_h5ad('batch2.h5ad')
batch3 = ad.read_h5ad('batch3.h5ad')

# Concatenate with full tracking
combined = ad.concat(
    [batch1, batch2, batch3],
    axis=0,
    join='outer',              # Keep all genes
    merge='first',             # Use first batch's annotations
    label='batch_id',          # Track source
    keys=['b1', 'b2', 'b3'],  # Custom labels
    index_unique='-'           # Make cell names unique
)
```

### 6. Data Conversion and Extraction

Convert between AnnData and other formats for interoperability.

**To DataFrame:**
```python
# Convert X to DataFrame
df = adata.to_df()

# Convert specific layer
df = adata.to_df(layer='log1p')
```

**Extract vectors:**
```python
# Get 1D arrays from data or annotations
gene_expression = adata.obs_vector('Gene_100')
cell_metadata = adata.obs_vector('n_genes')
```

**Transpose:**
```python
# Swap observations and variables
transposed = adata.T
```

### 7. Memory Optimization

Strategies for working with large datasets efficiently.

**Use sparse matrices:**
```python
from scipy.sparse import csr_matrix

# Check sparsity
density = (adata.X != 0).sum() / adata.X.size
if density < 0.3:  # Less than 30% non-zero
    adata.X = csr_matrix(adata.X)
```

**Convert strings to categoricals:**
```python
# Automatic conversion
adata.strings_to_categoricals()

# Manual conversion (more control)
adata.obs['cell_type'] = pd.Categorical(adata.obs['cell_type'])
```

**Use backed mode:**
```python
# Read without loading into memory
adata = ad.read_h5ad('large_file.h5ad', backed='r')

# Work with subsets
subset = adata[:1000, :500].copy()  # Only this subset in memory
```

**Chunked processing:**
```python
# Process data in chunks
results = []
for chunk in adata.chunk_X(chunk_size=1000):
    result = expensive_computation(chunk)
    results.append(result)
```

## Common Workflows

### Single-Cell RNA-seq Analysis

Complete workflow from loading to analysis:

```python
import anndata as ad
import numpy as np
import pandas as pd

# 1. Load data
adata = ad.read_mtx('matrix.mtx')
adata.obs_names = pd.read_csv('barcodes.tsv', header=None)[0]
adata.var_names = pd.read_csv('genes.tsv', header=None)[0]

# 2. Quality control
adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)
adata.obs['total_counts'] = adata.X.sum(axis=1)
adata = adata[adata.obs.n_genes > 200]
adata = adata[adata.obs.total_counts < 10000]

# 3. Filter genes
min_cells = 3
adata = adata[:, (adata.X > 0).sum(axis=0) >= min_cells]

# 4. Store raw counts
adata.layers['counts'] = adata.X.copy()

# 5. Normalize
adata.X = adata.X / adata.obs.total_counts.values[:, None] * 1e4
adata.X = np.log1p(adata.X)

# 6. Feature selection
gene_var = adata.X.var(axis=0)
adata.var['highly_variable'] = gene_var > np.percentile(gene_var, 90)

# 7. Dimensionality reduction (example with external tools)
# adata.obsm['X_pca'] = compute_pca(adata.X)
# adata.obsm['X_umap'] = compute_umap(adata.obsm['X_pca'])

# 8. Save results
adata.write('analyzed.h5ad', compression='gzip')
```

### Batch Integration

Combining multiple experimental batches:

```python
# Load batches
batches = [ad.read_h5ad(f'batch_{i}.h5ad') for i in range(3)]

# Concatenate with tracking
combined = ad.concat(
    batches,
    axis=0,
    join='outer',
    label='batch',
    keys=['batch_0', 'batch_1', 'batch_2'],
    index_unique='-'
)

# Add batch as numeric for correction algorithms
combined.obs['batch_numeric'] = combined.obs['batch'].cat.codes

# Perform batch correction (with external tools)
# corrected_pca = run_harmony(combined.obsm['X_pca'], combined.obs['batch'])
# combined.obsm['X_pca_corrected'] = corrected_pca

# Save integrated data
combined.write('integrated.h5ad', compression='gzip')
```

### Memory-Efficient Large Dataset Processing

Working with datasets too large for memory:

```python
# Read in backed mode
adata = ad.read_h5ad('huge_dataset.h5ad', backed='r')

# Compute statistics in chunks
total = 0
for chunk in adata.chunk_X(chunk_size=1000):
    total += chunk.sum()

mean_expression = total / (adata.n_obs * adata.n_vars)

# Work with subset
high_quality_cells = adata.obs.n_genes > 1000
subset = adata[high_quality_cells, :].copy()

# Close file
adata.file.close()
```

## Best Practices

### Data Organization

1. **Use layers for different representations**: Store raw counts, normalized, log-transformed, and scaled data in separate layers
2. **Use obsm/varm for multi-dimensional data**: Embeddings, loadings, and other matrix-like annotations
3. **Use uns for metadata**: Analysis parameters, dates, version information
4. **Use categoricals for efficiency**: Convert repeated strings to categorical types

### Subsetting

1. **Understand views vs copies**: Subsetting returns views by default; use `.copy()` when you need independence
2. **Chain conditions efficiently**: Combine boolean masks in a single subsetting operation
3. **Validate after subsetting**: Check dimensions and data integrity

### File I/O

1. **Use compression**: Always use `compression='gzip'` when writing h5ad files
2. **Choose the right format**: H5AD for general use, Zarr for cloud storage, Loom for compatibility
3. **Close backed files**: Always close file connections when done
4. **Use backed mode for large files**: Don't load everything into memory if not needed

### Concatenation

1. **Choose appropriate join**: Inner join for complete cases, outer join to preserve all features
2. **Track sources**: Use `label` and `keys` to track data origin
3. **Handle duplicates**: Use `index_unique` to make observation names unique
4. **Select merge strategy**: Choose appropriate merge strategy for variable annotations

### Memory Management

1. **Use sparse matrices**: For data with <30% non-zero values
2. **Convert to categoricals**: For repeated string values
3. **Process in chunks**: For operations on very large matrices
4. **Use backed mode**: Read large files with `backed='r'`

### Naming Conventions

Follow these conventions for consistency:

- **Embeddings**: `X_pca`, `X_umap`, `X_tsne`
- **Layers**: Descriptive names like `counts`, `log1p`, `scaled`
- **Observations**: Use snake_case like `cell_type`, `n_genes`, `total_counts`
- **Variables**: Use snake_case like `highly_variable`, `gene_name`

## Reference Documentation

For detailed API information, usage patterns, and troubleshooting, refer to the comprehensive reference files in the `references/` directory:

1. **api_reference.md**: Complete API documentation including all classes, methods, and functions with usage examples. Use `grep -r "pattern" references/api_reference.md` to search for specific functions or parameters.

2. **workflows_best_practices.md**: Detailed workflows for common tasks (single-cell analysis, batch integration, large datasets), best practices for memory management, data organization, and common pitfalls to avoid. Use `grep -r "pattern" references/workflows_best_practices.md` to search for specific workflows.

3. **concatenation_guide.md**: Comprehensive guide to concatenation strategies, join types, merge strategies, source tracking, and troubleshooting concatenation issues. Use `grep -r "pattern" references/concatenation_guide.md` to search for concatenation patterns.

## When to Load References

Load reference files into context when:
- Implementing complex concatenation with specific merge strategies
- Troubleshooting errors or unexpected behavior
- Optimizing memory usage for large datasets
- Implementing complete analysis workflows
- Understanding nuances of specific API methods

To search within references without loading them:
```python
# Example: Search for information about backed mode
grep -r "backed mode" references/
```

## Common Error Patterns

### Memory Errors
**Problem**: "MemoryError: Unable to allocate array"
**Solution**: Use backed mode, sparse matrices, or process in chunks

### Dimension Mismatch
**Problem**: "ValueError: operands could not be broadcast together"
**Solution**: Use outer join in concatenation or align indices before operations

### View Modification
**Problem**: "ValueError: assignment destination is read-only"
**Solution**: Convert view to copy with `.copy()` before modification

### File Already Open
**Problem**: "OSError: Unable to open file (file is already open)"
**Solution**: Close previous file connection with `adata.file.close()`
