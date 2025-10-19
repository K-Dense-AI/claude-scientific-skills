# AnnData Concatenation Guide

## Overview

The `concat()` function combines multiple AnnData objects through two fundamental operations:
1. **Concatenation**: Stacking sub-elements in order
2. **Merging**: Combining collections into one result

## Basic Concatenation

### Syntax
```python
import anndata as ad

combined = ad.concat(
    adatas,              # List of AnnData objects
    axis=0,              # 0=observations, 1=variables
    join='inner',        # 'inner' or 'outer'
    merge=None,          # Merge strategy for non-concat axis
    label=None,          # Column name for source tracking
    keys=None,           # Dataset identifiers
    index_unique=None,   # Separator for unique indices
    fill_value=None,     # Fill value for missing data
    pairwise=False       # Include pairwise matrices
)
```

### Concatenating Observations (Cells)
```python
# Most common: combining multiple samples/batches
adata1 = ad.AnnData(np.random.rand(100, 2000))
adata2 = ad.AnnData(np.random.rand(150, 2000))
adata3 = ad.AnnData(np.random.rand(80, 2000))

combined = ad.concat([adata1, adata2, adata3], axis=0)
# Result: (330 observations, 2000 variables)
```

### Concatenating Variables (Genes)
```python
# Less common: combining different feature sets
adata1 = ad.AnnData(np.random.rand(100, 1000))
adata2 = ad.AnnData(np.random.rand(100, 500))

combined = ad.concat([adata1, adata2], axis=1)
# Result: (100 observations, 1500 variables)
```

## Join Strategies

### Inner Join (Intersection)

Keeps only shared features across all objects.

```python
# Datasets with different genes
adata1 = ad.AnnData(
    np.random.rand(100, 2000),
    var=pd.DataFrame(index=[f'Gene_{i}' for i in range(2000)])
)
adata2 = ad.AnnData(
    np.random.rand(150, 1800),
    var=pd.DataFrame(index=[f'Gene_{i}' for i in range(200, 2000)])
)

# Inner join: only genes present in both
combined = ad.concat([adata1, adata2], join='inner')
# Result: (250 observations, 1800 variables)
# Only Gene_200 through Gene_1999
```

**Use when:**
- You want to analyze only features measured in all datasets
- Missing features would compromise analysis
- You need a complete case analysis

**Trade-offs:**
- May lose many features
- Ensures no missing data
- Smaller result size

### Outer Join (Union)

Keeps all features from all objects, padding with fill values (default 0).

```python
# Outer join: all genes from both datasets
combined = ad.concat([adata1, adata2], join='outer')
# Result: (250 observations, 2000 variables)
# Missing values filled with 0

# Custom fill value
combined = ad.concat([adata1, adata2], join='outer', fill_value=np.nan)
```

**Use when:**
- You want to preserve all features
- Sparse data is acceptable
- Features are independent

**Trade-offs:**
- Introduces zeros/missing values
- Larger result size
- May need imputation

## Merge Strategies

Merge strategies control how elements on the non-concatenation axis are combined.

### merge=None (Default)

Excludes all non-concatenation axis elements.

```python
# Both datasets have var annotations
adata1.var['gene_type'] = ['protein_coding'] * 2000
adata2.var['gene_type'] = ['protein_coding'] * 1800

# merge=None: var annotations excluded
combined = ad.concat([adata1, adata2], merge=None)
assert 'gene_type' not in combined.var.columns
```

**Use when:**
- Annotations are dataset-specific
- You'll add new annotations after merging

### merge='same'

Keeps only annotations with identical values across datasets.

```python
# Same annotation values
adata1.var['chromosome'] = ['chr1'] * 1000 + ['chr2'] * 1000
adata2.var['chromosome'] = ['chr1'] * 900 + ['chr2'] * 900

# merge='same': keeps chromosome annotation
combined = ad.concat([adata1, adata2], merge='same')
assert 'chromosome' in combined.var.columns
```

**Use when:**
- Annotations should be consistent
- You want to validate consistency
- Shared metadata is important

**Note:** Comparison occurs after index alignment - only shared indices need to match.

### merge='unique'

Includes annotations with a single possible value.

```python
# Unique values per gene
adata1.var['ensembl_id'] = [f'ENSG{i:08d}' for i in range(2000)]
adata2.var['ensembl_id'] = [f'ENSG{i:08d}' for i in range(2000)]

# merge='unique': keeps ensembl_id
combined = ad.concat([adata1, adata2], merge='unique')
```

**Use when:**
- Each feature has a unique identifier
- Annotations are feature-specific

### merge='first'

Takes the first occurrence of each annotation.

```python
# Different annotation versions
adata1.var['description'] = ['desc1'] * 2000
adata2.var['description'] = ['desc2'] * 2000

# merge='first': uses adata1's descriptions
combined = ad.concat([adata1, adata2], merge='first')
# Uses descriptions from adata1
```

**Use when:**
- One dataset has authoritative annotations
- Order matters
- You need a simple resolution strategy

### merge='only'

Retains annotations appearing in exactly one object.

```python
# Dataset-specific annotations
adata1.var['dataset1_specific'] = ['value'] * 2000
adata2.var['dataset2_specific'] = ['value'] * 2000

# merge='only': keeps both (no conflicts)
combined = ad.concat([adata1, adata2], merge='only')
```

**Use when:**
- Datasets have non-overlapping annotations
- You want to preserve all unique metadata

## Source Tracking

### Using label

Add a categorical column to track data origin.

```python
combined = ad.concat(
    [adata1, adata2, adata3],
    label='batch'
)

# Creates obs['batch'] with values 0, 1, 2
print(combined.obs['batch'].cat.categories)  # ['0', '1', '2']
```

### Using keys

Provide custom names for source tracking.

```python
combined = ad.concat(
    [adata1, adata2, adata3],
    label='study',
    keys=['control', 'treatment_a', 'treatment_b']
)

# Creates obs['study'] with custom names
print(combined.obs['study'].unique())  # ['control', 'treatment_a', 'treatment_b']
```

### Making Indices Unique

Append source identifiers to duplicate observation names.

```python
# Both datasets have cells named "Cell_0", "Cell_1", etc.
adata1.obs_names = [f'Cell_{i}' for i in range(100)]
adata2.obs_names = [f'Cell_{i}' for i in range(150)]

# index_unique adds suffix
combined = ad.concat(
    [adata1, adata2],
    keys=['batch1', 'batch2'],
    index_unique='-'
)

# Results in: Cell_0-batch1, Cell_0-batch2, etc.
print(combined.obs_names[:5])
```

## Handling Different Attributes

### X Matrix and Layers

Follows join strategy. Missing values filled according to `fill_value`.

```python
# Both have layers
adata1.layers['counts'] = adata1.X.copy()
adata2.layers['counts'] = adata2.X.copy()

# Concatenates both X and layers
combined = ad.concat([adata1, adata2])
assert 'counts' in combined.layers
```

### obs and var DataFrames

- **obs**: Concatenated along concatenation axis
- **var**: Handled by merge strategy

```python
adata1.obs['cell_type'] = ['B cell'] * 100
adata2.obs['cell_type'] = ['T cell'] * 150

combined = ad.concat([adata1, adata2])
# obs['cell_type'] preserved for all cells
```

### obsm and varm

Multi-dimensional annotations follow same rules as layers.

```python
adata1.obsm['X_pca'] = np.random.rand(100, 50)
adata2.obsm['X_pca'] = np.random.rand(150, 50)

combined = ad.concat([adata1, adata2])
# obsm['X_pca'] concatenated: shape (250, 50)
```

### obsp and varp

Pairwise matrices excluded by default. Enable with `pairwise=True`.

```python
# Distance matrices
adata1.obsp['distances'] = np.random.rand(100, 100)
adata2.obsp['distances'] = np.random.rand(150, 150)

# Excluded by default
combined = ad.concat([adata1, adata2])
assert 'distances' not in combined.obsp

# Include if needed
combined = ad.concat([adata1, adata2], pairwise=True)
# Results in padded block diagonal matrix
```

### uns Dictionary

Merged recursively, applying merge strategy at any nesting depth.

```python
adata1.uns['experiment'] = {'date': '2024-01', 'lab': 'A'}
adata2.uns['experiment'] = {'date': '2024-02', 'lab': 'A'}

# merge='same' keeps 'lab', excludes 'date'
combined = ad.concat([adata1, adata2], merge='same')
# combined.uns['experiment'] = {'lab': 'A'}
```

## Advanced Patterns

### Batch Integration Pipeline

```python
import anndata as ad

# Load batches
batches = [
    ad.read_h5ad(f'batch_{i}.h5ad')
    for i in range(5)
]

# Concatenate with tracking
combined = ad.concat(
    batches,
    axis=0,
    join='outer',
    merge='first',
    label='batch_id',
    keys=[f'batch_{i}' for i in range(5)],
    index_unique='-'
)

# Add batch effects
combined.obs['batch_numeric'] = combined.obs['batch_id'].cat.codes
```

### Multi-Study Meta-Analysis

```python
# Different studies with varying gene coverage
studies = {
    'study_a': ad.read_h5ad('study_a.h5ad'),
    'study_b': ad.read_h5ad('study_b.h5ad'),
    'study_c': ad.read_h5ad('study_c.h5ad')
}

# Outer join to keep all genes
combined = ad.concat(
    list(studies.values()),
    axis=0,
    join='outer',
    label='study',
    keys=list(studies.keys()),
    merge='unique',
    fill_value=0
)

# Track coverage
for study in studies:
    n_genes = studies[study].n_vars
    combined.uns[f'{study}_n_genes'] = n_genes
```

### Incremental Concatenation

```python
# For many datasets, concatenate in batches
chunk_size = 10
all_files = [f'dataset_{i}.h5ad' for i in range(100)]

# Process in chunks
result = None
for i in range(0, len(all_files), chunk_size):
    chunk_files = all_files[i:i+chunk_size]
    chunk_adatas = [ad.read_h5ad(f) for f in chunk_files]
    chunk_combined = ad.concat(chunk_adatas)

    if result is None:
        result = chunk_combined
    else:
        result = ad.concat([result, chunk_combined])
```

### Memory-Efficient On-Disk Concatenation

```python
# Experimental feature for large datasets
from anndata.experimental import concat_on_disk

files = ['dataset1.h5ad', 'dataset2.h5ad', 'dataset3.h5ad']
concat_on_disk(
    files,
    'combined.h5ad',
    join='outer'
)

# Read result in backed mode
combined = ad.read_h5ad('combined.h5ad', backed='r')
```

## Troubleshooting

### Issue: Dimension Mismatch

```python
# Error: shapes don't match
adata1 = ad.AnnData(np.random.rand(100, 2000))
adata2 = ad.AnnData(np.random.rand(150, 1500))

# Solution: use outer join
combined = ad.concat([adata1, adata2], join='outer')
```

### Issue: Memory Error

```python
# Problem: too many large objects in memory
large_adatas = [ad.read_h5ad(f) for f in many_files]

# Solution: read and concatenate incrementally
result = None
for file in many_files:
    adata = ad.read_h5ad(file)
    if result is None:
        result = adata
    else:
        result = ad.concat([result, adata])
        del adata  # Free memory
```

### Issue: Duplicate Indices

```python
# Problem: same cell names in different batches
# Solution: use index_unique
combined = ad.concat(
    [adata1, adata2],
    keys=['batch1', 'batch2'],
    index_unique='-'
)
```

### Issue: Lost Annotations

```python
# Problem: annotations disappear
adata1.var['important'] = values1
adata2.var['important'] = values2

combined = ad.concat([adata1, adata2])  # merge=None by default
# Solution: use appropriate merge strategy
combined = ad.concat([adata1, adata2], merge='first')
```

## Performance Tips

1. **Pre-align indices**: Ensure consistent naming before concatenation
2. **Use sparse matrices**: Convert to sparse before concatenating
3. **Batch operations**: Concatenate in groups for many datasets
4. **Choose inner join**: When possible, to reduce result size
5. **Use categoricals**: Convert string annotations before concatenating
6. **Consider on-disk**: For very large datasets, use `concat_on_disk`
