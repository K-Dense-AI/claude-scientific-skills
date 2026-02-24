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




## Best Practices


This comprehensive export guide covers all aspects of getting data out of TileDB-VCF in various formats optimized for different downstream analysis workflows.