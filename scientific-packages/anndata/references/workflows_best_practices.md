# AnnData Workflows and Best Practices

## Common Workflows

### 1. Single-Cell RNA-seq Analysis Workflow

#### Loading Data
```python
import anndata as ad
import numpy as np
import pandas as pd

# Load from 10X format
adata = ad.read_mtx('matrix.mtx')
adata.var_names = pd.read_csv('genes.tsv', sep='\t', header=None)[0]
adata.obs_names = pd.read_csv('barcodes.tsv', sep='\t', header=None)[0]

# Or load from pre-processed h5ad
adata = ad.read_h5ad('preprocessed_data.h5ad')
```

#### Quality Control
```python
# Calculate QC metrics
adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)
adata.obs['total_counts'] = adata.X.sum(axis=1)

# Filter cells
adata = adata[adata.obs.n_genes > 200]
adata = adata[adata.obs.total_counts < 10000]

# Filter genes
min_cells = 3
adata = adata[:, (adata.X > 0).sum(axis=0) >= min_cells]
```

#### Normalization and Preprocessing
```python
# Store raw counts
adata.layers['counts'] = adata.X.copy()

# Normalize
adata.X = adata.X / adata.obs.total_counts.values[:, None] * 1e4

# Log transform
adata.layers['log1p'] = np.log1p(adata.X)
adata.X = adata.layers['log1p']

# Identify highly variable genes
gene_variance = adata.X.var(axis=0)
adata.var['highly_variable'] = gene_variance > np.percentile(gene_variance, 90)
```

#### Dimensionality Reduction
```python
# PCA
from sklearn.decomposition import PCA
pca = PCA(n_components=50)
adata.obsm['X_pca'] = pca.fit_transform(adata.X)

# Store PCA variance
adata.uns['pca'] = {'variance_ratio': pca.explained_variance_ratio_}

# UMAP
from umap import UMAP
umap = UMAP(n_components=2)
adata.obsm['X_umap'] = umap.fit_transform(adata.obsm['X_pca'])
```

#### Clustering
```python
# Store cluster assignments
adata.obs['clusters'] = pd.Categorical(['cluster_0', 'cluster_1', ...])

# Store cluster centroids
centroids = np.array([...])
adata.varm['cluster_centroids'] = centroids
```

#### Save Results
```python
# Save complete analysis
adata.write('analyzed_data.h5ad', compression='gzip')
```

### 2. Batch Integration Workflow

```python
import anndata as ad

# Load multiple batches
batch1 = ad.read_h5ad('batch1.h5ad')
batch2 = ad.read_h5ad('batch2.h5ad')
batch3 = ad.read_h5ad('batch3.h5ad')

# Concatenate with batch labels
adata = ad.concat(
    [batch1, batch2, batch3],
    axis=0,
    label='batch',
    keys=['batch1', 'batch2', 'batch3'],
    index_unique='-'
)

# Batch effect correction would go here
# (using external tools like Harmony, Scanorama, etc.)

# Store corrected embeddings
adata.obsm['X_pca_corrected'] = corrected_pca
adata.obsm['X_umap_corrected'] = corrected_umap
```

### 3. Memory-Efficient Large Dataset Workflow

```python
import anndata as ad

# Read in backed mode
adata = ad.read_h5ad('large_dataset.h5ad', backed='r')

# Check backing status
print(f"Is backed: {adata.isbacked}")
print(f"File: {adata.filename}")

# Work with chunks
for chunk in adata.chunk_X(chunk_size=1000):
    # Process chunk
    result = process_chunk(chunk)

# Close file when done
adata.file.close()
```

### 4. Multi-Dataset Comparison Workflow

```python
import anndata as ad

# Load datasets
datasets = {
    'study1': ad.read_h5ad('study1.h5ad'),
    'study2': ad.read_h5ad('study2.h5ad'),
    'study3': ad.read_h5ad('study3.h5ad')
}

# Outer join to keep all genes
combined = ad.concat(
    list(datasets.values()),
    axis=0,
    join='outer',
    label='study',
    keys=list(datasets.keys()),
    merge='first'
)

# Handle missing data
combined.X[np.isnan(combined.X)] = 0

# Add dataset-specific metadata
combined.uns['datasets'] = {
    'study1': {'date': '2023-01', 'n_samples': datasets['study1'].n_obs},
    'study2': {'date': '2023-06', 'n_samples': datasets['study2'].n_obs},
    'study3': {'date': '2024-01', 'n_samples': datasets['study3'].n_obs}
}
```

## Best Practices

### Memory Management

#### Use Sparse Matrices
```python
from scipy.sparse import csr_matrix

# Convert to sparse if data is sparse
if density < 0.3:  # Less than 30% non-zero
    adata.X = csr_matrix(adata.X)
```

#### Use Backed Mode for Large Files
```python
# Read with backing
adata = ad.read_h5ad('large_file.h5ad', backed='r')

# Only load what you need
subset = adata[:1000, :500].copy()  # Now in memory
```

#### Convert Strings to Categoricals
```python
# Efficient storage for repeated strings
adata.strings_to_categoricals()

# Or manually
adata.obs['cell_type'] = pd.Categorical(adata.obs['cell_type'])
```

### Data Organization

#### Use Layers for Different Representations
```python
# Store multiple versions of the data
adata.layers['counts'] = raw_counts
adata.layers['normalized'] = normalized_data
adata.layers['log1p'] = log_transformed_data
adata.layers['scaled'] = scaled_data
```

#### Use obsm/varm for Multi-Dimensional Annotations
```python
# Embeddings
adata.obsm['X_pca'] = pca_coordinates
adata.obsm['X_umap'] = umap_coordinates
adata.obsm['X_tsne'] = tsne_coordinates

# Gene loadings
adata.varm['PCs'] = principal_components
```

#### Use uns for Analysis Metadata
```python
# Store parameters
adata.uns['preprocessing'] = {
    'normalization': 'TPM',
    'min_genes': 200,
    'min_cells': 3,
    'date': '2024-01-15'
}

# Store analysis results
adata.uns['differential_expression'] = {
    'method': 't-test',
    'p_value_threshold': 0.05
}
```

### Subsetting and Views

#### Understand View vs Copy
```python
# Subsetting returns a view
subset = adata[adata.obs.cell_type == 'B cell']  # View
print(subset.is_view)  # True

# Views are memory efficient but modifications affect original
subset.obs['new_column'] = value  # Modifies original adata

# Create independent copy when needed
subset_copy = adata[adata.obs.cell_type == 'B cell'].copy()
```

#### Chain Operations Efficiently
```python
# Bad - creates multiple intermediate views
temp1 = adata[adata.obs.batch == 'batch1']
temp2 = temp1[temp1.obs.n_genes > 200]
result = temp2[:, temp2.var.highly_variable].copy()

# Good - chain operations
result = adata[
    (adata.obs.batch == 'batch1') & (adata.obs.n_genes > 200),
    adata.var.highly_variable
].copy()
```

### File I/O

#### Use Compression
```python
# Save with compression
adata.write('data.h5ad', compression='gzip')
```

#### Choose the Right Format
```python
# H5AD for general use (good compression, fast)
adata.write_h5ad('data.h5ad')

# Zarr for cloud storage and parallel access
adata.write_zarr('data.zarr')

# Loom for compatibility with other tools
adata.write_loom('data.loom')
```

#### Close File Connections
```python
# Use context manager pattern
adata = ad.read_h5ad('file.h5ad', backed='r')
try:
    # Work with data
    process(adata)
finally:
    adata.file.close()
```

### Concatenation

#### Choose Appropriate Join Strategy
```python
# Inner join - only common features (safe, may lose data)
combined = ad.concat([adata1, adata2], join='inner')

# Outer join - all features (keeps all data, may introduce zeros)
combined = ad.concat([adata1, adata2], join='outer')
```

#### Track Data Sources
```python
# Add source labels
combined = ad.concat(
    [adata1, adata2, adata3],
    label='dataset',
    keys=['exp1', 'exp2', 'exp3']
)

# Make indices unique
combined = ad.concat(
    [adata1, adata2, adata3],
    index_unique='-'
)
```

#### Handle Variable-Specific Metadata
```python
# Use merge strategy for var annotations
combined = ad.concat(
    [adata1, adata2],
    merge='same',  # Keep only identical annotations
    join='outer'
)
```

### Naming Conventions

#### Use Consistent Naming
```python
# Embeddings: X_<method>
adata.obsm['X_pca']
adata.obsm['X_umap']
adata.obsm['X_tsne']

# Layers: descriptive names
adata.layers['counts']
adata.layers['log1p']
adata.layers['scaled']

# Observations: snake_case
adata.obs['cell_type']
adata.obs['n_genes']
adata.obs['total_counts']
```

#### Make Indices Unique
```python
# Ensure unique names
adata.obs_names_make_unique()
adata.var_names_make_unique()
```

### Error Handling

#### Validate Data Structure
```python
# Check dimensions
assert adata.n_obs > 0, "No observations in data"
assert adata.n_vars > 0, "No variables in data"

# Check for NaN values
if np.isnan(adata.X).any():
    print("Warning: NaN values detected")

# Check for negative values in count data
if (adata.X < 0).any():
    print("Warning: Negative values in count data")
```

#### Handle Missing Data
```python
# Check for missing annotations
if adata.obs['cell_type'].isna().any():
    print("Warning: Missing cell type annotations")
    # Fill or remove
    adata = adata[~adata.obs['cell_type'].isna()]
```

## Common Pitfalls

### 1. Forgetting to Copy Views
```python
# BAD - modifies original
subset = adata[adata.obs.condition == 'treated']
subset.X = transformed_data  # Changes original adata!

# GOOD
subset = adata[adata.obs.condition == 'treated'].copy()
subset.X = transformed_data  # Only changes subset
```

### 2. Mixing Backed and In-Memory Operations
```python
# BAD - trying to modify backed data
adata = ad.read_h5ad('file.h5ad', backed='r')
adata.X[0, 0] = 100  # Error: can't modify backed data

# GOOD - load to memory first
adata = ad.read_h5ad('file.h5ad', backed='r')
adata = adata.to_memory()
adata.X[0, 0] = 100  # Works
```

### 3. Not Using Categoricals for Metadata
```python
# BAD - stores as strings (memory inefficient)
adata.obs['cell_type'] = ['B cell', 'T cell', ...] * 1000

# GOOD - use categorical
adata.obs['cell_type'] = pd.Categorical(['B cell', 'T cell', ...] * 1000)
```

### 4. Incorrect Concatenation Axis
```python
# Concatenating observations (cells)
combined = ad.concat([adata1, adata2], axis=0)  # Correct

# Concatenating variables (genes) - rare
combined = ad.concat([adata1, adata2], axis=1)  # Less common
```

### 5. Not Preserving Raw Data
```python
# BAD - loses original data
adata.X = normalized_data

# GOOD - preserve original
adata.layers['counts'] = adata.X.copy()
adata.X = normalized_data
```
