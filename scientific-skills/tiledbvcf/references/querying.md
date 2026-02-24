# TileDB-VCF Querying Guide

Complete guide to efficiently querying variant data from TileDB-VCF datasets with optimal performance and memory usage.

## Basic Querying

### Opening a Dataset
```python
import tiledbvcf
import pandas as pd

# Open dataset for reading
ds = tiledbvcf.Dataset(uri="my_dataset", mode="r")

# With custom configuration
config = tiledbvcf.ReadConfig(
    memory_budget=2048,
    tiledb_config={"sm.tile_cache_size": "1000000000"}
)
ds = tiledbvcf.Dataset(uri="my_dataset", mode="r", cfg=config)
```

### Simple Region Query
```python
# Query single region
df = ds.read(
    attrs=["sample_name", "pos_start", "pos_end", "alleles"],
    regions=["chr1:1000000-2000000"]
)
print(df.head())
```

### Multi-Region Query
```python
# Query multiple regions efficiently
df = ds.read(
    attrs=["sample_name", "pos_start", "alleles", "fmt_GT"],
    regions=[
        "chr1:1000000-2000000",
        "chr2:500000-1500000",
        "chr22:20000000-21000000"
    ]
)
```

## Sample Filtering

### Query Specific Samples
```python
# Query subset of samples
df = ds.read(
    attrs=["sample_name", "pos_start", "alleles", "fmt_GT"],
    regions=["chr1:1000000-2000000"],
    samples=["SAMPLE_001", "SAMPLE_002", "SAMPLE_003"]
)
```

### Sample Pattern Matching
```python
# Get all sample names first
sample_names = ds.sample_names()
print(f"Total samples: {len(sample_names)}")

# Filter samples by pattern
case_samples = [s for s in sample_names if s.startswith("CASE_")]
control_samples = [s for s in sample_names if s.startswith("CTRL_")]

# Query case samples
case_df = ds.read(
    attrs=["sample_name", "pos_start", "alleles", "fmt_GT"],
    regions=["chr1:1000000-2000000"],
    samples=case_samples
)
```

## Attribute Selection

### Available Attributes
```python
# List all available attributes
attrs = ds.attributes()
print("Available attributes:")
for attr in attrs:
    print(f"  {attr}")
```

### Core Variant Attributes
```python
# Essential variant information
df = ds.read(
    attrs=[
        "sample_name",     # Sample identifier
        "contig",          # Chromosome
        "pos_start",       # Start position (1-based)
        "pos_end",         # End position (1-based)
        "alleles",         # REF and ALT alleles
        "id",              # Variant ID (rs numbers, etc.)
        "qual",            # Quality score
        "filters"          # FILTER field values
    ],
    regions=["chr1:1000000-2000000"]
)
```

### INFO Field Access
```python
# Access INFO fields (prefix with 'info_')
df = ds.read(
    attrs=[
        "sample_name", "pos_start", "alleles",
        "info_AF",         # Allele frequency
        "info_AC",         # Allele count
        "info_AN",         # Allele number
        "info_DP",         # Total depth
        "info_ExcessHet"   # Excess heterozygosity
    ],
    regions=["chr1:1000000-2000000"]
)
```

### FORMAT Field Access
```python
# Access FORMAT fields (prefix with 'fmt_')
df = ds.read(
    attrs=[
        "sample_name", "pos_start", "alleles",
        "fmt_GT",          # Genotype
        "fmt_DP",          # Read depth
        "fmt_GQ",          # Genotype quality
        "fmt_AD",          # Allelic depths
        "fmt_PL"           # Phred-scaled likelihoods
    ],
    regions=["chr1:1000000-2000000"],
    samples=["SAMPLE_001", "SAMPLE_002"]
)
```

## Advanced Filtering

### Quality-Based Filtering
```python
# Query with quality information
df = ds.read(
    attrs=["sample_name", "pos_start", "alleles", "qual", "fmt_GQ", "fmt_DP"],
    regions=["chr1:1000000-2000000"]
)

# Filter high-quality variants in pandas
high_qual_df = df[
    (df["qual"] >= 30) &
    (df["fmt_GQ"] >= 20) &
    (df["fmt_DP"] >= 10)
]
```

### Genotype Filtering
```python
# Query genotype data
df = ds.read(
    attrs=["sample_name", "pos_start", "alleles", "fmt_GT"],
    regions=["chr1:1000000-2000000"]
)

# Filter for heterozygous variants
# fmt_GT format: [allele1, allele2] where -1 = missing
het_variants = df[
    df["fmt_GT"].apply(lambda gt:
        isinstance(gt, list) and len(gt) == 2 and
        gt[0] != gt[1] and -1 not in gt
    )
]
```

## Performance Optimization

### Memory-Efficient Streaming
```python
# For large result sets, use iterator pattern
def stream_variants(dataset, regions, chunk_size=100000):
    """Stream variants in chunks to manage memory"""
    total_variants = 0

    for region in regions:
        print(f"Processing region: {region}")

        # Query region
        df = dataset.read(
            attrs=["sample_name", "pos_start", "alleles", "fmt_GT"],
            regions=[region]
        )

        # Process in chunks
        for i in range(0, len(df), chunk_size):
            chunk = df.iloc[i:i+chunk_size]
            yield chunk
            total_variants += len(chunk)

            if i + chunk_size < len(df):
                print(f"  Processed {i + chunk_size:,} variants...")

    print(f"Total variants processed: {total_variants:,}")

# Usage
regions = ["chr1", "chr2", "chr3"]
for chunk in stream_variants(ds, regions):
    # Process chunk
    pass
```

### Parallel Region Processing
```python
import multiprocessing as mp
from functools import partial

def query_region(dataset_uri, region, samples=None):
    """Query single region - can be parallelized"""
    ds = tiledbvcf.Dataset(uri=dataset_uri, mode="r")

    df = ds.read(
        attrs=["sample_name", "pos_start", "alleles", "fmt_GT"],
        regions=[region],
        samples=samples
    )

    return region, df

def parallel_query(dataset_uri, regions, samples=None, n_processes=4):
    """Query multiple regions in parallel"""
    query_func = partial(query_region, dataset_uri, samples=samples)

    with mp.Pool(n_processes) as pool:
        results = pool.map(query_func, regions)

    # Combine results
    combined_df = pd.concat([df for _, df in results], ignore_index=True)
    return combined_df

# Usage
regions = [f"chr{i}:1000000-2000000" for i in range(1, 23)]
df = parallel_query("my_dataset", regions, n_processes=8)
```

### Query Optimization Strategies
```python
# Optimize tile cache for repeated region access
config = tiledbvcf.ReadConfig(
    memory_budget=4096,
    tiledb_config={
        "sm.tile_cache_size": "2000000000",      # 2GB cache
        "sm.mem.reader.sparse_global_order.ratio_tile_ranges": "0.5",
        "sm.mem.reader.sparse_unordered.ratio_coords": "0.5"
    }
)

ds = tiledbvcf.Dataset(uri="my_dataset", mode="r", cfg=config)

# Combine nearby regions to reduce query overhead
def optimize_regions(regions, max_gap=1000000):
    """Merge nearby regions to reduce query count"""
    if not regions:
        return regions

    # Parse regions (simplified - assumes chr:start-end format)
    parsed = []
    for region in regions:
        chrom, pos_range = region.split(":")
        start, end = map(int, pos_range.split("-"))
        parsed.append((chrom, start, end))

    # Sort by chromosome and position
    parsed.sort()

    merged = []
    current_chrom, current_start, current_end = parsed[0]

    for chrom, start, end in parsed[1:]:
        if chrom == current_chrom and start - current_end <= max_gap:
            # Merge regions
            current_end = max(current_end, end)
        else:
            # Add current region and start new one
            merged.append(f"{current_chrom}:{current_start}-{current_end}")
            current_chrom, current_start, current_end = chrom, start, end

    # Add final region
    merged.append(f"{current_chrom}:{current_start}-{current_end}")

    return merged
```

## Complex Queries

### Population-Specific Queries
```python
# Query specific populations
def query_population(ds, regions, population_file):
    """Query variants for specific population"""
    # Read population assignments
    pop_df = pd.read_csv(population_file)  # sample_id, population

    populations = {}
    for _, row in pop_df.iterrows():
        pop = row['population']
        if pop not in populations:
            populations[pop] = []
        populations[pop].append(row['sample_id'])

    results = {}
    for pop_name, samples in populations.items():
        print(f"Querying {pop_name}: {len(samples)} samples")

        df = ds.read(
            attrs=["sample_name", "pos_start", "alleles", "fmt_GT"],
            regions=regions,
            samples=samples
        )

        results[pop_name] = df

    return results

# Usage
regions = ["chr1:1000000-2000000"]
pop_results = query_population(ds, regions, "population_assignments.csv")
```

### Allele Frequency Calculation
```python
def calculate_allele_frequencies(df):
    """Calculate allele frequencies from genotype data"""
    af_results = []

    # Group by variant position
    for (contig, pos), group in df.groupby(['contig', 'pos_start']):
        alleles = group['alleles'].iloc[0]  # REF and ALT alleles
        genotypes = group['fmt_GT'].values

        # Count alleles
        allele_counts = {}
        total_alleles = 0

        for gt in genotypes:
            if isinstance(gt, list) and -1 not in gt:  # Valid genotype
                for allele_idx in gt:
                    allele_counts[allele_idx] = allele_counts.get(allele_idx, 0) + 1
                    total_alleles += 1

        # Calculate frequencies
        frequencies = {}
        for allele_idx, count in allele_counts.items():
            if allele_idx < len(alleles):
                allele = alleles[allele_idx]
                freq = count / total_alleles if total_alleles > 0 else 0
                frequencies[allele] = freq

        af_results.append({
            'contig': contig,
            'pos': pos,
            'alleles': alleles,
            'frequencies': frequencies,
            'sample_count': len(group)
        })

    return af_results

# Usage
df = ds.read(
    attrs=["contig", "pos_start", "alleles", "fmt_GT"],
    regions=["chr1:1000000-1010000"],
    samples=sample_list[:100]  # Subset for faster calculation
)

af_results = calculate_allele_frequencies(df)
```

### Rare Variant Analysis
```python
def find_rare_variants(ds, regions, samples, max_af=0.01):
    """Find rare variants (MAF < threshold)"""
    # Query with allele frequency info
    df = ds.read(
        attrs=["sample_name", "pos_start", "alleles", "info_AF", "fmt_GT"],
        regions=regions,
        samples=samples
    )

    # Filter for rare variants
    rare_df = df[df["info_AF"] < max_af]

    # Group by variant and count carriers
    rare_summary = []
    for (pos, alleles), group in rare_df.groupby(['pos_start', 'alleles']):
        # Count carriers (samples with non-ref genotypes)
        carriers = []
        for _, row in group.iterrows():
            gt = row['fmt_GT']
            if isinstance(gt, list) and any(allele > 0 for allele in gt if allele != -1):
                carriers.append(row['sample_name'])

        rare_summary.append({
            'pos': pos,
            'alleles': alleles,
            'af': group['info_AF'].iloc[0],
            'carrier_count': len(carriers),
            'carriers': carriers
        })

    return pd.DataFrame(rare_summary)

# Usage
rare_variants = find_rare_variants(
    ds,
    regions=["chr1:1000000-2000000"],
    samples=sample_list,
    max_af=0.005  # 0.5% frequency threshold
)
```

## Cloud Querying

### S3 Dataset Queries
```python
# Configure for S3 access
s3_config = tiledbvcf.ReadConfig(
    memory_budget=2048,
    tiledb_config={
        "vfs.s3.region": "us-east-1",
        "vfs.s3.use_virtual_addressing": "true",
        "sm.tile_cache_size": "1000000000"
    }
)

# Query S3-hosted dataset
ds = tiledbvcf.Dataset(
    uri="s3://my-bucket/vcf-dataset",
    mode="r",
    cfg=s3_config
)

df = ds.read(
    attrs=["sample_name", "pos_start", "alleles", "fmt_GT"],
    regions=["chr1:1000000-2000000"]
)
```

### Federated Queries
```python
def federated_query(dataset_uris, regions, samples):
    """Query multiple datasets and combine results"""
    combined_results = []

    for dataset_uri in dataset_uris:
        print(f"Querying dataset: {dataset_uri}")

        ds = tiledbvcf.Dataset(uri=dataset_uri, mode="r")

        # Check which samples are available in this dataset
        available_samples = set(ds.sample_names())
        query_samples = [s for s in samples if s in available_samples]

        if query_samples:
            df = ds.read(
                attrs=["sample_name", "pos_start", "alleles", "fmt_GT"],
                regions=regions,
                samples=query_samples
            )
            df['dataset'] = dataset_uri
            combined_results.append(df)

    if combined_results:
        return pd.concat(combined_results, ignore_index=True)
    else:
        return pd.DataFrame()

# Usage
datasets = [
    "s3://data-bucket/study1-dataset",
    "s3://data-bucket/study2-dataset",
    "local_dataset"
]
results = federated_query(datasets, ["chr1:1000000-2000000"], ["SAMPLE_001"])
```

## Troubleshooting

### Common Query Issues
```python
# Handle empty results gracefully
def safe_query(ds, **kwargs):
    """Query with error handling"""
    try:
        df = ds.read(**kwargs)
        if df.empty:
            print("Query returned no results")
            return pd.DataFrame()
        return df

    except Exception as e:
        print(f"Query failed: {e}")
        return pd.DataFrame()

# Check dataset contents before querying
def inspect_dataset(ds):
    """Inspect dataset properties"""
    print(f"Sample count: {len(ds.sample_names())}")
    print(f"Sample names (first 10): {ds.sample_names()[:10]}")
    print(f"Available attributes: {ds.attributes()}")

    # Try small test query
    try:
        test_df = ds.read(
            attrs=["sample_name", "pos_start"],
            regions=["chr1:1-1000"]
        )
        print(f"Test query returned {len(test_df)} variants")
    except Exception as e:
        print(f"Test query failed: {e}")

# Usage
inspect_dataset(ds)
```

This comprehensive querying guide covers all aspects of efficiently retrieving and filtering variant data from TileDB-VCF datasets, from basic queries to complex population genomics analyses.