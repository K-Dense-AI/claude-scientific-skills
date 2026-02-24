# TileDB-VCF Ingestion Guide

Complete guide to creating TileDB-VCF datasets and ingesting VCF/BCF files with optimal performance and reliability.

## Dataset Creation

### Basic Dataset Creation
```python
import tiledbvcf

# Create a new dataset
ds = tiledbvcf.Dataset(uri="my_dataset", mode="w")
```

### Advanced Configuration
```python
# Custom configuration for large datasets
config = tiledbvcf.ReadConfig(
    memory_budget=4096,  # MB
    tiledb_config={
        "sm.tile_cache_size": "2000000000",  # 2GB tile cache
        "sm.mem.total_budget": "4000000000",  # 4GB total memory
        "vfs.file.posix_file_permissions": "644"
    }
)

ds = tiledbvcf.Dataset(
    uri="large_dataset",
    mode="w",
    cfg=config
)
```

### Cloud Dataset Creation
```python
# S3 dataset with credentials
config = tiledbvcf.ReadConfig(
    tiledb_config={
        "vfs.s3.aws_access_key_id": "YOUR_KEY",
        "vfs.s3.aws_secret_access_key": "YOUR_SECRET",
        "vfs.s3.region": "us-east-1"
    }
)

ds = tiledbvcf.Dataset(
    uri="s3://my-bucket/vcf-dataset",
    mode="w",
    cfg=config
)
```

## Single Sample Ingestion

### Basic Ingestion
```python
# Ingest a single VCF file
ds.ingest_samples(["sample1.vcf.gz"])

# Multiple files at once
ds.ingest_samples([
    "sample1.vcf.gz",
    "sample2.vcf.gz",
    "sample3.vcf.gz"
])
```

### Custom Sample Names
```python
# Override sample names from VCF headers
ds.ingest_samples(
    ["data/unknown_sample.vcf.gz"],
    sample_names=["SAMPLE_001"]
)
```

### Ingestion with Validation
```python
# Enable additional validation during ingestion
try:
    ds.ingest_samples(
        ["sample1.vcf.gz"],
        contig_fragment_merging=True,  # Merge fragments on same contig
        resume=False  # Start fresh (don't resume)
    )
except Exception as e:
    print(f"Ingestion failed: {e}")
```

## Parallel Ingestion

### Multi-threaded Ingestion
```python
# Configure for parallel ingestion
config = tiledbvcf.ReadConfig(
    tiledb_config={
        "sm.num_async_threads": "8",
        "sm.num_reader_threads": "4",
        "sm.num_writer_threads": "4"
    }
)

ds = tiledbvcf.Dataset(uri="dataset", mode="w", cfg=config)

# Ingest multiple files in parallel
file_list = [f"sample_{i}.vcf.gz" for i in range(1, 101)]
ds.ingest_samples(file_list)
```

### Batched Processing
```python
# Process files in batches to manage memory
import glob

vcf_files = glob.glob("*.vcf.gz")
batch_size = 10

for i in range(0, len(vcf_files), batch_size):
    batch = vcf_files[i:i+batch_size]
    print(f"Processing batch {i//batch_size + 1}: {len(batch)} files")

    try:
        ds.ingest_samples(batch)
        print(f"Successfully ingested batch {i//batch_size + 1}")
    except Exception as e:
        print(f"Error in batch {i//batch_size + 1}: {e}")
        # Continue with next batch
        continue
```

## Incremental Addition

### Adding New Samples
```python
# Open existing dataset and add new samples
ds = tiledbvcf.Dataset(uri="existing_dataset", mode="a")  # append mode

# Add new samples without affecting existing data
ds.ingest_samples(["new_sample1.vcf.gz", "new_sample2.vcf.gz"])
```

### Resuming Interrupted Ingestion
```python
# Resume a previously interrupted ingestion
ds.ingest_samples(
    ["large_sample.vcf.gz"],
    resume=True  # Continue from where it left off
)
```

## Memory Optimization

### Memory Budget Configuration
```python
# Configure memory usage based on system resources
import psutil

# Use 75% of available memory
available_memory = psutil.virtual_memory().available
memory_budget_mb = int((available_memory * 0.75) / (1024 * 1024))

config = tiledbvcf.ReadConfig(
    memory_budget=memory_budget_mb,
    tiledb_config={
        "sm.mem.total_budget": str(int(available_memory * 0.75))
    }
)
```

### Large File Handling
```python
# For very large VCF files (>10GB), use streaming ingestion
config = tiledbvcf.ReadConfig(
    memory_budget=2048,  # Conservative memory usage
    tiledb_config={
        "sm.tile_cache_size": "500000000",  # 500MB cache
        "sm.consolidation.buffer_size": "100000000"  # 100MB buffer
    }
)

# Process large files one at a time
large_files = ["huge_sample1.vcf.gz", "huge_sample2.vcf.gz"]
for vcf_file in large_files:
    print(f"Processing {vcf_file}")
    ds.ingest_samples([vcf_file])
    print(f"Completed {vcf_file}")
```

## Error Handling and Validation

### Comprehensive Error Handling
```python
import logging

logging.basicConfig(level=logging.INFO)

def robust_ingestion(dataset_uri, vcf_files):
    config = tiledbvcf.ReadConfig(memory_budget=2048)

    with tiledbvcf.Dataset(uri=dataset_uri, mode="w", cfg=config) as ds:
        failed_files = []

        for vcf_file in vcf_files:
            try:
                # Validate file exists and is readable
                if not os.path.exists(vcf_file):
                    logging.error(f"File not found: {vcf_file}")
                    failed_files.append(vcf_file)
                    continue

                logging.info(f"Ingesting {vcf_file}")
                ds.ingest_samples([vcf_file])
                logging.info(f"Successfully ingested {vcf_file}")

            except Exception as e:
                logging.error(f"Failed to ingest {vcf_file}: {e}")
                failed_files.append(vcf_file)
                continue

        if failed_files:
            logging.warning(f"Failed to ingest {len(failed_files)} files: {failed_files}")

        return failed_files
```

### Pre-ingestion Validation
```python
import pysam

def validate_vcf_files(vcf_files):
    """Validate VCF files before ingestion"""
    valid_files = []
    invalid_files = []

    for vcf_file in vcf_files:
        try:
            # Basic validation using pysam
            vcf = pysam.VariantFile(vcf_file)

            # Check if file has variants
            try:
                next(iter(vcf))
                valid_files.append(vcf_file)
                print(f"✓ {vcf_file}: Valid")
            except StopIteration:
                print(f"⚠ {vcf_file}: No variants found")
                valid_files.append(vcf_file)  # Empty files are valid

            vcf.close()

        except Exception as e:
            print(f"✗ {vcf_file}: Invalid - {e}")
            invalid_files.append(vcf_file)

    return valid_files, invalid_files

# Use validation before ingestion
vcf_files = ["sample1.vcf.gz", "sample2.vcf.gz"]
valid_files, invalid_files = validate_vcf_files(vcf_files)

if valid_files:
    ds.ingest_samples(valid_files)
```

## Performance Optimization

### I/O Optimization
```python
# Optimize for different storage types
def get_optimized_config(storage_type="local"):
    base_config = {
        "sm.mem.total_budget": "4000000000",
        "sm.tile_cache_size": "1000000000"
    }

    if storage_type == "local":
        # Optimize for local SSD storage
        base_config.update({
            "sm.num_async_threads": "8",
            "vfs.file.enable_filelocks": "true"
        })
    elif storage_type == "s3":
        # Optimize for S3 storage
        base_config.update({
            "vfs.s3.multipart_part_size": "50MB",
            "vfs.s3.max_parallel_ops": "8",
            "vfs.s3.use_multipart_upload": "true"
        })
    elif storage_type == "azure":
        # Optimize for Azure Blob Storage
        base_config.update({
            "vfs.azure.max_parallel_ops": "8",
            "vfs.azure.block_list_block_size": "50MB"
        })

    return tiledbvcf.ReadConfig(
        memory_budget=4096,
        tiledb_config=base_config
    )
```

### Monitoring Ingestion Progress
```python
import time
from pathlib import Path

def ingest_with_progress(dataset, vcf_files):
    """Ingest files with progress monitoring"""
    start_time = time.time()
    total_files = len(vcf_files)

    for i, vcf_file in enumerate(vcf_files, 1):
        file_start = time.time()
        file_size = Path(vcf_file).stat().st_size / (1024*1024)  # MB

        print(f"[{i}/{total_files}] Processing {vcf_file} ({file_size:.1f} MB)")

        try:
            dataset.ingest_samples([vcf_file])
            file_duration = time.time() - file_start

            print(f"  ✓ Completed in {file_duration:.1f}s "
                  f"({file_size/file_duration:.1f} MB/s)")

        except Exception as e:
            print(f"  ✗ Failed: {e}")

    total_duration = time.time() - start_time
    print(f"\nIngestion complete: {total_duration:.1f}s total")
```

## Cloud Storage Patterns

### S3 Ingestion Pipeline
```python
import boto3

def ingest_from_s3_bucket(dataset_uri, bucket, prefix):
    """Ingest VCF files from S3 bucket"""
    s3 = boto3.client('s3')

    # List VCF files in bucket
    response = s3.list_objects_v2(
        Bucket=bucket,
        Prefix=prefix,
        MaxKeys=1000
    )

    vcf_files = [
        f"s3://{bucket}/{obj['Key']}"
        for obj in response.get('Contents', [])
        if obj['Key'].endswith(('.vcf.gz', '.vcf'))
    ]

    print(f"Found {len(vcf_files)} VCF files in s3://{bucket}/{prefix}")

    # Configure for S3
    config = get_optimized_config("s3")

    with tiledbvcf.Dataset(uri=dataset_uri, mode="w", cfg=config) as ds:
        ds.ingest_samples(vcf_files)

# Usage
ingest_from_s3_bucket(
    dataset_uri="s3://my-output-bucket/vcf-dataset",
    bucket="my-input-bucket",
    prefix="vcf_files/"
)
```

## Best Practices

### Dataset Organization
```python
# Organize datasets by study or cohort
study_datasets = {
    "ukb": "s3://genomics-data/ukb-dataset",
    "1kgp": "s3://genomics-data/1kgp-dataset",
    "gnomad": "s3://genomics-data/gnomad-dataset"
}

def create_study_dataset(study_name, vcf_files):
    """Create a dataset for a specific study"""
    dataset_uri = study_datasets[study_name]

    config = tiledbvcf.ReadConfig(
        memory_budget=4096,
        tiledb_config={
            "sm.consolidation.mode": "fragments",
            "sm.consolidation.buffer_size": "200000000"
        }
    )

    with tiledbvcf.Dataset(uri=dataset_uri, mode="w", cfg=config) as ds:
        ds.ingest_samples(vcf_files)
```

### Maintenance and Consolidation
```python
# Consolidate dataset after ingestion for optimal query performance
def consolidate_dataset(dataset_uri):
    """Consolidate dataset fragments for better query performance"""
    config = tiledbvcf.ReadConfig(
        tiledb_config={
            "sm.consolidation.mode": "fragments",
            "sm.consolidation.buffer_size": "1000000000"  # 1GB buffer
        }
    )

    # Note: Consolidation API varies by TileDB-VCF version
    # This is a conceptual example
    print(f"Consolidating dataset: {dataset_uri}")
    # Implementation depends on specific TileDB-VCF version
```

This comprehensive guide covers all aspects of TileDB-VCF ingestion from basic single-file ingestion to complex cloud-based parallel processing workflows. Use the patterns that best fit your data scale and infrastructure requirements.