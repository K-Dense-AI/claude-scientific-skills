# AnnData API Reference

## Core AnnData Class

The `AnnData` class is the central data structure for storing and manipulating annotated datasets in single-cell genomics and other domains.

### Core Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| **X** | array-like | Primary data matrix (#observations × #variables). Supports NumPy arrays, sparse matrices (CSR/CSC), HDF5 datasets, Zarr arrays, and Dask arrays |
| **obs** | DataFrame | One-dimensional annotation of observations (rows). Length equals observation count |
| **var** | DataFrame | One-dimensional annotation of variables/features (columns). Length equals variable count |
| **uns** | OrderedDict | Unstructured annotation for miscellaneous metadata |
| **obsm** | dict-like | Multi-dimensional observation annotations (structured arrays aligned to observation axis) |
| **varm** | dict-like | Multi-dimensional variable annotations (structured arrays aligned to variable axis) |
| **obsp** | dict-like | Pairwise observation annotations (square matrices representing graphs) |
| **varp** | dict-like | Pairwise variable annotations (graphs between features) |
| **layers** | dict-like | Additional data matrices matching X's dimensions |
| **raw** | AnnData | Stores original versions of X and var before transformations |

### Dimensional Properties

- **n_obs**: Number of observations (sample count)
- **n_vars**: Number of variables/features
- **shape**: Tuple returning (n_obs, n_vars)
- **T**: Transposed view of the entire object

### State Properties

- **isbacked**: Boolean indicating disk-backed storage status
- **is_view**: Boolean identifying whether object is a view of another AnnData
- **filename**: Path to backing .h5ad file; setting this enables disk-backed mode

### Key Methods

#### Construction and Copying
- **`AnnData(X=None, obs=None, var=None, ...)`**: Create new AnnData object
- **`copy(filename=None)`**: Create full copy, optionally stored on disk

#### Subsetting and Views
- **`adata[obs_subset, var_subset]`**: Subset observations and variables (returns view by default)
- **`.copy()`**: Convert view to independent object

#### Data Access
- **`to_df(layer=None)`**: Generate pandas DataFrame representation
- **`obs_vector(k, layer=None)`**: Extract 1D array from X, layers, or annotations
- **`var_vector(k, layer=None)`**: Extract 1D array for a variable
- **`chunk_X(chunk_size)`**: Iterate over data matrix in chunks
- **`chunked_X(chunk_size)`**: Context manager for chunked iteration

#### Transformation
- **`transpose()`**: Return transposed object
- **`concatenate(*adatas, ...)`**: Combine multiple AnnData objects along observation axis
- **`to_memory(copy=False)`**: Load all backed arrays into RAM

#### File I/O
- **`write_h5ad(filename, compression='gzip')`**: Save as .h5ad HDF5 format
- **`write_zarr(store, ...)`**: Export hierarchical Zarr store
- **`write_loom(filename, ...)`**: Output .loom format file
- **`write_csvs(dirname, ...)`**: Write annotations as separate CSV files

#### Data Management
- **`strings_to_categoricals()`**: Convert string annotations to categorical types
- **`rename_categories(key, categories)`**: Update category labels in annotations
- **`obs_names_make_unique(sep='-')`**: Append numeric suffixes to duplicate observation names
- **`var_names_make_unique(sep='-')`**: Append numeric suffixes to duplicate variable names

## Module-Level Functions

### Reading Functions

#### Native Formats
- **`read_h5ad(filename, backed=None, as_sparse=None)`**: Load HDF5-based .h5ad files
- **`read_zarr(store)`**: Access hierarchical Zarr array stores

#### Alternative Formats
- **`read_csv(filename, ...)`**: Import from CSV files
- **`read_excel(filename, ...)`**: Import from Excel files
- **`read_hdf(filename, key)`**: Read from HDF5 files
- **`read_loom(filename, ...)`**: Import from .loom files
- **`read_mtx(filename, ...)`**: Import from Matrix Market format
- **`read_text(filename, ...)`**: Import from text files
- **`read_umi_tools(filename, ...)`**: Import from UMI-tools format

#### Element-Level Access
- **`read_elem(elem)`**: Retrieve specific components from storage
- **`sparse_dataset(group)`**: Generate backed sparse matrix classes

### Combining Operations
- **`concat(adatas, axis=0, join='inner', merge=None, ...)`**: Merge multiple AnnData objects
  - **axis**: 0 (observations) or 1 (variables)
  - **join**: 'inner' (intersection) or 'outer' (union)
  - **merge**: Strategy for non-concatenation axis ('same', 'unique', 'first', 'only', or None)
  - **label**: Column name for source tracking
  - **keys**: Dataset identifiers for source annotation
  - **index_unique**: Separator for making duplicate indices unique

### Writing Functions
- **`write_h5ad(filename, adata, compression='gzip')`**: Export to HDF5 format
- **`write_zarr(store, adata, ...)`**: Save as Zarr hierarchical arrays
- **`write_elem(elem, ...)`**: Write individual components

### Experimental Features
- **`AnnCollection`**: Batch processing for large collections
- **`AnnLoader`**: PyTorch DataLoader integration
- **`concat_on_disk(*adatas, filename, ...)`**: Memory-efficient out-of-core concatenation
- **`read_lazy(filename)`**: Lazy loading with deferred computation
- **`read_dispatched(filename, ...)`**: Custom I/O with callbacks
- **`write_dispatched(filename, ...)`**: Custom writing with callbacks

### Configuration
- **`settings`**: Package-wide configuration object
- **`settings.override(**kwargs)`**: Context manager for temporary settings changes

## Common Usage Patterns

### Creating AnnData Objects

```python
import anndata as ad
import numpy as np
from scipy.sparse import csr_matrix

# From dense array
counts = np.random.poisson(1, size=(100, 2000))
adata = ad.AnnData(counts)

# From sparse matrix
counts = csr_matrix(np.random.poisson(1, size=(100, 2000)), dtype=np.float32)
adata = ad.AnnData(counts)

# With metadata
import pandas as pd
obs_meta = pd.DataFrame({'cell_type': ['B', 'T', 'Monocyte'] * 33 + ['B']})
var_meta = pd.DataFrame({'gene_name': [f'Gene_{i}' for i in range(2000)]})
adata = ad.AnnData(counts, obs=obs_meta, var=var_meta)
```

### Subsetting

```python
# By names
subset = adata[['Cell_1', 'Cell_10'], ['Gene_5', 'Gene_1900']]

# By boolean mask
b_cells = adata[adata.obs.cell_type == 'B']

# By position
first_five = adata[:5, :100]

# Convert view to copy
adata_copy = adata[:5].copy()
```

### Adding Annotations

```python
# Cell-level metadata
adata.obs['batch'] = pd.Categorical(['batch1', 'batch2'] * 50)

# Gene-level metadata
adata.var['highly_variable'] = np.random.choice([True, False], size=adata.n_vars)

# Embeddings
adata.obsm['X_pca'] = np.random.normal(size=(adata.n_obs, 50))
adata.obsm['X_umap'] = np.random.normal(size=(adata.n_obs, 2))

# Alternative data representations
adata.layers['log_transformed'] = np.log1p(adata.X)
adata.layers['scaled'] = (adata.X - adata.X.mean(axis=0)) / adata.X.std(axis=0)

# Unstructured metadata
adata.uns['experiment_date'] = '2024-01-15'
adata.uns['parameters'] = {'min_genes': 200, 'min_cells': 3}
```

### File I/O

```python
# Write to disk
adata.write('my_results.h5ad', compression='gzip')

# Read into memory
adata = ad.read_h5ad('my_results.h5ad')

# Read in backed mode (memory-efficient)
adata = ad.read_h5ad('my_results.h5ad', backed='r')

# Close file connection
adata.file.close()
```

### Concatenation

```python
# Combine multiple datasets
adata1 = ad.AnnData(np.random.poisson(1, size=(100, 2000)))
adata2 = ad.AnnData(np.random.poisson(1, size=(150, 2000)))
adata3 = ad.AnnData(np.random.poisson(1, size=(80, 2000)))

# Simple concatenation
combined = ad.concat([adata1, adata2, adata3], axis=0)

# With source labels
combined = ad.concat(
    [adata1, adata2, adata3],
    axis=0,
    label='dataset',
    keys=['exp1', 'exp2', 'exp3']
)

# Inner join (only shared variables)
combined = ad.concat([adata1, adata2, adata3], axis=0, join='inner')

# Outer join (all variables, pad with zeros)
combined = ad.concat([adata1, adata2, adata3], axis=0, join='outer')
```
