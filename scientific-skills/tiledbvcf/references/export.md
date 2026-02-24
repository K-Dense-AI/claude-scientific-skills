# TileDB-VCF Export Guide

Complete guide to exporting data from TileDB-VCF datasets in various formats for downstream analysis and integration with other genomics tools.

## VCF/BCF Export

### Basic VCF Export
```python
import tiledbvcf

# Open dataset for reading
ds = tiledbvcf.Dataset(uri="my_dataset", mode="r")

# Export specific regions as VCF
ds.export_vcf(
    uri="output.vcf.gz",
    regions=["chr1:1000000-2000000"],
    samples=["SAMPLE_001", "SAMPLE_002", "SAMPLE_003"]
)
```

### BCF Export (Binary VCF)
```python
# Export as compressed BCF for faster processing
ds.export_bcf(
    uri="output.bcf",
    regions=["chr1:1000000-2000000", "chr2:500000-1500000"],
    samples=["SAMPLE_001", "SAMPLE_002"]
)
```

### Large-Scale Export
```python
# Export entire chromosomes efficiently
def export_chromosome(ds, chrom, output_dir, samples=None):
    """Export full chromosome data"""
    output_path = f"{output_dir}/chr{chrom}.bcf"

    print(f"Exporting chromosome {chrom}")
    ds.export_bcf(
        uri=output_path,
        regions=[f"chr{chrom}"],
        samples=samples
    )
    print(f"Exported to {output_path}")

# Export all autosomes
for chrom in range(1, 23):
    export_chromosome(ds, chrom, "exported_data")
```

## TSV Export

### Basic TSV Export
```python
# Export as tab-separated values
ds.export_tsv(
    uri="variants.tsv",
    regions=["chr1:1000000-2000000"],
    samples=["SAMPLE_001", "SAMPLE_002"],
    tsv_fields=["CHR", "POS", "REF", "ALT", "S:GT", "S:DP"]
)
```


## Pandas DataFrame Export

### Query to DataFrame
```python
# Export query results as pandas DataFrame for analysis
df = ds.read(
    attrs=["sample_name", "contig", "pos_start", "alleles", "fmt_GT", "info_AF"],
    regions=["chr1:1000000-2000000"],
    samples=["SAMPLE_001", "SAMPLE_002", "SAMPLE_003"]
)

# Save DataFrame to various formats
df.to_csv("variants.csv", index=False)
df.to_parquet("variants.parquet")
df.to_pickle("variants.pkl")
```

### Processed Data Export
```python
def export_processed_variants(ds, regions, samples, output_file):
    """Export processed variant data with calculated metrics"""
    # Query raw data
    df = ds.read(
        attrs=["sample_name", "contig", "pos_start", "pos_end",
               "alleles", "fmt_GT", "fmt_DP", "fmt_GQ", "info_AF"],
        regions=regions,
        samples=samples
    )

    # Add calculated columns
    df['variant_id'] = df['contig'] + ':' + df['pos_start'].astype(str)

    # Parse genotypes
    def parse_genotype(gt):
        if isinstance(gt, list) and len(gt) == 2:
            if -1 in gt:
                return "missing"
            elif gt[0] == gt[1]:
                return "homozygous"
            else:
                return "heterozygous"
        return "unknown"

    df['genotype_type'] = df['fmt_GT'].apply(parse_genotype)

    # Filter high-quality variants
    high_qual = df[
        (df['fmt_DP'] >= 10) &
        (df['fmt_GQ'] >= 20) &
        (df['genotype_type'] != 'missing')
    ]

    # Export processed data
    high_qual.to_csv(output_file, index=False)
    print(f"Exported {len(high_qual)} high-quality variants to {output_file}")

    return high_qual

# Usage
processed_df = export_processed_variants(
    ds,
    regions=["chr1:1000000-2000000"],
    samples=ds.sample_names()[:50],  # First 50 samples
    output_file="high_quality_variants.csv"
)
```

## Streaming Export for Large Datasets

### Chunked Export
```python
def streaming_export(ds, regions, samples, output_file, chunk_size=100000):
    """Export large datasets in chunks to manage memory"""
    import csv

    total_variants = 0

    with open(output_file, 'w', newline='') as f:
        writer = None
        header_written = False

        for region in regions:
            print(f"Processing region: {region}")

            # Query region
            df = ds.read(
                attrs=["sample_name", "contig", "pos_start", "alleles", "fmt_GT"],
                regions=[region],
                samples=samples
            )

            if df.empty:
                continue

            # Process in chunks
            for i in range(0, len(df), chunk_size):
                chunk = df.iloc[i:i+chunk_size]

                # Write header on first chunk
                if not header_written:
                    writer = csv.writer(f)
                    writer.writerow(chunk.columns)
                    header_written = True

                # Write chunk data
                for _, row in chunk.iterrows():
                    writer.writerow(row.values)

                total_variants += len(chunk)

                if i + chunk_size < len(df):
                    print(f"  Processed {i + chunk_size:,} variants...")

    print(f"Exported {total_variants:,} variants to {output_file}")

# Usage
regions = [f"chr{i}" for i in range(1, 23)]  # All autosomes
streaming_export(ds, regions, ds.sample_names(), "genome_wide_variants.csv")
```

### Parallel Export
```python
import multiprocessing as mp
import os

def export_region_chunk(args):
    """Export single region - for parallel processing"""
    dataset_uri, region, samples, output_dir = args

    # Create separate dataset instance for each process
    ds = tiledbvcf.Dataset(uri=dataset_uri, mode="r")

    # Generate output filename
    region_safe = region.replace(":", "_").replace("-", "_")
    output_file = os.path.join(output_dir, f"variants_{region_safe}.tsv")

    # Export region
    ds.export_tsv(
        uri=output_file,
        regions=[region],
        samples=samples,
        tsv_fields=["CHR", "POS", "REF", "ALT", "S:GT", "S:DP"]
    )

    return region, output_file

def parallel_export(dataset_uri, regions, samples, output_dir, n_processes=4):
    """Export multiple regions in parallel"""
    os.makedirs(output_dir, exist_ok=True)

    # Prepare arguments for parallel processing
    args = [(dataset_uri, region, samples, output_dir) for region in regions]

    # Export in parallel
    with mp.Pool(n_processes) as pool:
        results = pool.map(export_region_chunk, args)

    # Combine results if needed
    output_files = [output_file for _, output_file in results]
    print(f"Exported {len(output_files)} region files to {output_dir}")

    return output_files

# Usage
regions = [f"chr{i}:1-50000000" for i in range(1, 23)]  # First half of each chromosome
output_files = parallel_export(
    dataset_uri="my_dataset",
    regions=regions,
    samples=ds.sample_names()[:100],
    output_dir="parallel_export",
    n_processes=8
)
```

## Integration with Analysis Tools

### PLINK Format Export
```python
def export_for_plink(ds, regions, samples, output_prefix):
    """Export data in format suitable for PLINK analysis"""
    # Query variant data
    df = ds.read(
        attrs=["sample_name", "contig", "pos_start", "id", "alleles", "fmt_GT"],
        regions=regions,
        samples=samples
    )

    # Prepare PLINK-compatible data
    plink_data = []
    for _, row in df.iterrows():
        gt = row['fmt_GT']
        if isinstance(gt, list) and len(gt) == 2 and -1 not in gt:
            # Convert genotype to PLINK format (0/1/2)
            alleles = row['alleles']
            if len(alleles) >= 2:
                ref_allele = alleles[0]
                alt_allele = alleles[1]

                # Count alternative alleles
                alt_count = sum(1 for allele in gt if allele == 1)

                plink_data.append({
                    'sample': row['sample_name'],
                    'chr': row['contig'],
                    'pos': row['pos_start'],
                    'id': row['id'] if row['id'] else f"{row['contig']}_{row['pos_start']}",
                    'ref': ref_allele,
                    'alt': alt_allele,
                    'genotype': alt_count
                })

    # Save as PLINK-compatible format
    plink_df = pd.DataFrame(plink_data)

    # Pivot for PLINK .raw format
    plink_matrix = plink_df.pivot_table(
        index='sample',
        columns=['chr', 'pos', 'id'],
        values='genotype',
        fill_value=-9  # Missing data code
    )

    # Save files
    plink_matrix.to_csv(f"{output_prefix}.raw", sep='\t')

    # Create map file
    map_data = plink_df[['chr', 'id', 'pos']].drop_duplicates()
    map_data['genetic_distance'] = 0  # Placeholder
    map_data = map_data[['chr', 'id', 'genetic_distance', 'pos']]
    map_data.to_csv(f"{output_prefix}.map", sep='\t', header=False, index=False)

    print(f"Exported PLINK files: {output_prefix}.raw, {output_prefix}.map")

# Usage
export_for_plink(
    ds,
    regions=["chr22"],  # Start with smaller chromosome
    samples=ds.sample_names()[:100],
    output_prefix="plink_data"
)
```

### VEP Annotation Preparation
```python
def export_for_vep(ds, regions, output_file):
    """Export variants for VEP (Variant Effect Predictor) annotation"""
    # Query essential variant information
    df = ds.read(
        attrs=["contig", "pos_start", "pos_end", "alleles", "id"],
        regions=regions
    )

    # Prepare VEP input format
    vep_data = []
    for _, row in df.iterrows():
        alleles = row['alleles']
        if len(alleles) >= 2:
            ref = alleles[0]
            for alt in alleles[1:]:  # Can have multiple ALT alleles
                vep_data.append({
                    'chr': row['contig'],
                    'start': row['pos_start'],
                    'end': row['pos_end'],
                    'allele': f"{ref}/{alt}",
                    'strand': '+',
                    'id': row['id'] if row['id'] else '.'
                })

    vep_df = pd.DataFrame(vep_data)

    # Save VEP input format
    vep_df.to_csv(
        output_file,
        sep='\t',
        header=False,
        index=False,
        columns=['chr', 'start', 'end', 'allele', 'strand', 'id']
    )

    print(f"Exported {len(vep_df)} variants for VEP annotation to {output_file}")

# Usage
export_for_vep(ds, ["chr1:1000000-2000000"], "variants_for_vep.txt")
```

## Cloud Export

### S3 Export
```python
def export_to_s3(ds, regions, samples, s3_bucket, s3_prefix):
    """Export data directly to S3"""
    import boto3

    # Configure for S3
    config = tiledbvcf.ReadConfig(
        tiledb_config={
            "vfs.s3.region": "us-east-1",
            "vfs.s3.multipart_part_size": "50MB"
        }
    )

    # Export to S3 paths
    for i, region in enumerate(regions):
        region_safe = region.replace(":", "_").replace("-", "_")
        s3_uri = f"s3://{s3_bucket}/{s3_prefix}/region_{region_safe}.bcf"

        print(f"Exporting region {i+1}/{len(regions)}: {region}")

        ds.export_bcf(
            uri=s3_uri,
            regions=[region],
            samples=samples
        )

        print(f"Exported to {s3_uri}")

# Usage
export_to_s3(
    ds,
    regions=["chr1:1000000-2000000", "chr2:500000-1500000"],
    samples=ds.sample_names()[:50],
    s3_bucket="my-genomics-bucket",
    s3_prefix="exported_variants"
)
```

## Export Validation

### Data Integrity Checks
```python
def validate_export(original_ds, export_file, regions, samples):
    """Validate exported data against original dataset"""
    import pysam

    # Count variants in original dataset
    original_df = original_ds.read(
        attrs=["sample_name", "pos_start"],
        regions=regions,
        samples=samples
    )
    original_count = len(original_df)

    # Count variants in exported file
    try:
        if export_file.endswith('.vcf.gz') or export_file.endswith('.bcf'):
            vcf = pysam.VariantFile(export_file)
            export_count = sum(1 for _ in vcf)
            vcf.close()
        elif export_file.endswith('.tsv') or export_file.endswith('.csv'):
            export_df = pd.read_csv(export_file, sep='\t' if export_file.endswith('.tsv') else ',')
            export_count = len(export_df)
        else:
            print(f"Unknown file format: {export_file}")
            return False

        # Compare counts
        if original_count == export_count:
            print(f"✓ Export validation passed: {export_count} variants")
            return True
        else:
            print(f"✗ Export validation failed: {original_count} original vs {export_count} exported")
            return False

    except Exception as e:
        print(f"✗ Export validation error: {e}")
        return False

# Usage
success = validate_export(
    ds,
    "output.bcf",
    regions=["chr1:1000000-2000000"],
    samples=["SAMPLE_001", "SAMPLE_002"]
)
```

## Best Practices

### Efficient Export Strategies
```python
# 1. Optimize for intended use case
def choose_export_format(use_case, file_size_mb):
    """Choose optimal export format based on use case"""
    if use_case == "downstream_analysis":
        if file_size_mb > 1000:
            return "BCF"  # Compressed binary
        else:
            return "VCF"  # Text format

    elif use_case == "data_sharing":
        return "VCF.gz"  # Standard compressed format

    elif use_case == "statistical_analysis":
        return "TSV"  # Easy to process

    elif use_case == "database_import":
        return "CSV"  # Universal format

    else:
        return "VCF"  # Default

# 2. Batch processing for large exports
def batch_export_by_size(ds, regions, samples, max_variants_per_file=1000000):
    """Export data in batches based on variant count"""
    current_batch = []
    current_count = 0
    batch_num = 1

    for region in regions:
        # Estimate variant count (approximate)
        test_df = ds.read(
            attrs=["pos_start"],
            regions=[region],
            samples=samples[:10]  # Small sample for estimation
        )
        estimated_variants = len(test_df) * len(samples) // 10

        if current_count + estimated_variants > max_variants_per_file and current_batch:
            # Export current batch
            export_batch(ds, current_batch, samples, f"batch_{batch_num}.bcf")
            batch_num += 1
            current_batch = [region]
            current_count = estimated_variants
        else:
            current_batch.append(region)
            current_count += estimated_variants

    # Export final batch
    if current_batch:
        export_batch(ds, current_batch, samples, f"batch_{batch_num}.bcf")

def export_batch(ds, regions, samples, output_file):
    """Export a batch of regions"""
    print(f"Exporting batch to {output_file}")
    ds.export_bcf(uri=output_file, regions=regions, samples=samples)
```

This comprehensive export guide covers all aspects of getting data out of TileDB-VCF in various formats optimized for different downstream analysis workflows.