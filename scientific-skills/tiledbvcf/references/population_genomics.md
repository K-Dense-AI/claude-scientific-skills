# TileDB-VCF Population Genomics Guide

Comprehensive guide to using TileDB-VCF for large-scale population genomics analyses including GWAS, rare variant studies, and population structure analysis.

## GWAS Data Preparation

### Quality Control Pipeline
```python
import tiledbvcf
import pandas as pd
import numpy as np

def gwas_qc_pipeline(ds, regions, samples, output_prefix):
    """Complete quality control pipeline for GWAS data"""
    print("Starting GWAS QC pipeline...")

    # Step 1: Query all relevant data
    print("1. Querying variant data...")
    df = ds.read(
        attrs=[
            "sample_name", "contig", "pos_start", "alleles",
            "fmt_GT", "fmt_DP", "fmt_GQ",
            "info_AF", "info_AC", "info_AN",
            "qual", "filters"
        ],
        regions=regions,
        samples=samples
    )

    print(f"   Initial variants: {len(df):,}")

    # Step 2: Sample-level QC
    print("2. Performing sample-level QC...")
    sample_qc = perform_sample_qc(df)
    good_samples = sample_qc[
        (sample_qc['call_rate'] >= 0.95) &
        (sample_qc['mean_depth'] >= 8) &
        (sample_qc['mean_depth'] <= 50) &
        (sample_qc['het_ratio'] >= 0.8) &
        (sample_qc['het_ratio'] <= 1.2)
    ]['sample'].tolist()

    print(f"   Samples passing QC: {len(good_samples)}/{len(samples)}")

    # Step 3: Variant-level QC
    print("3. Performing variant-level QC...")
    variant_qc = perform_variant_qc(df[df['sample_name'].isin(good_samples)])

    # Step 4: Export clean data
    print("4. Exporting QC'd data...")
    clean_variants = variant_qc[
        (variant_qc['call_rate'] >= 0.95) &
        (variant_qc['hwe_p'] >= 1e-6) &
        (variant_qc['maf'] >= 0.01) &
        (variant_qc['mean_qual'] >= 30)
    ]

    # Save QC results
    sample_qc.to_csv(f"{output_prefix}_sample_qc.csv", index=False)
    variant_qc.to_csv(f"{output_prefix}_variant_qc.csv", index=False)
    clean_variants.to_csv(f"{output_prefix}_clean_variants.csv", index=False)

    print(f"   Final variants for GWAS: {len(clean_variants):,}")
    return clean_variants, good_samples

def perform_sample_qc(df):
    """Calculate sample-level QC metrics"""
    sample_metrics = []

    for sample in df['sample_name'].unique():
        sample_data = df[df['sample_name'] == sample]

        # Calculate metrics
        total_calls = len(sample_data)
        missing_calls = sum(1 for gt in sample_data['fmt_GT']
                          if not isinstance(gt, list) or -1 in gt)
        call_rate = 1 - (missing_calls / total_calls)

        depths = [d for d in sample_data['fmt_DP'] if pd.notna(d) and d > 0]
        mean_depth = np.mean(depths) if depths else 0

        # Heterozygote ratio
        het_count = sum(1 for gt in sample_data['fmt_GT']
                       if isinstance(gt, list) and len(gt) == 2 and
                       gt[0] != gt[1] and -1 not in gt)
        hom_count = sum(1 for gt in sample_data['fmt_GT']
                       if isinstance(gt, list) and len(gt) == 2 and
                       gt[0] == gt[1] and -1 not in gt and gt[0] > 0)
        het_ratio = het_count / hom_count if hom_count > 0 else 0

        sample_metrics.append({
            'sample': sample,
            'total_variants': total_calls,
            'call_rate': call_rate,
            'mean_depth': mean_depth,
            'het_count': het_count,
            'hom_count': hom_count,
            'het_ratio': het_ratio
        })

    return pd.DataFrame(sample_metrics)

def perform_variant_qc(df):
    """Calculate variant-level QC metrics"""
    variant_metrics = []

    for (chrom, pos, alleles), group in df.groupby(['contig', 'pos_start', 'alleles']):
        # Calculate call rate
        total_samples = len(group)
        missing_calls = sum(1 for gt in group['fmt_GT']
                          if not isinstance(gt, list) or -1 in gt)
        call_rate = 1 - (missing_calls / total_samples)

        # Calculate MAF
        allele_counts = {0: 0, 1: 0}  # REF and ALT
        total_alleles = 0

        for gt in group['fmt_GT']:
            if isinstance(gt, list) and -1 not in gt:
                for allele in gt:
                    if allele in allele_counts:
                        allele_counts[allele] += 1
                        total_alleles += 1

        if total_alleles > 0:
            alt_freq = allele_counts[1] / total_alleles
            maf = min(alt_freq, 1 - alt_freq)
        else:
            maf = 0

        # Hardy-Weinberg Equilibrium test (simplified)
        hwe_p = calculate_hwe_p(group['fmt_GT'], alt_freq)

        # Mean quality
        quals = [q for q in group['qual'] if pd.notna(q)]
        mean_qual = np.mean(quals) if quals else 0

        variant_metrics.append({
            'contig': chrom,
            'pos': pos,
            'alleles': alleles,
            'call_rate': call_rate,
            'maf': maf,
            'alt_freq': alt_freq,
            'hwe_p': hwe_p,
            'mean_qual': mean_qual,
            'sample_count': total_samples
        })

    return pd.DataFrame(variant_metrics)

def calculate_hwe_p(genotypes, alt_freq):
    """Simple Hardy-Weinberg equilibrium test"""
    if alt_freq == 0 or alt_freq == 1:
        return 1.0

    # Count observed genotypes
    hom_ref = sum(1 for gt in genotypes
                  if isinstance(gt, list) and gt == [0, 0])
    het = sum(1 for gt in genotypes
              if isinstance(gt, list) and len(gt) == 2 and
              gt[0] != gt[1] and -1 not in gt)
    hom_alt = sum(1 for gt in genotypes
                  if isinstance(gt, list) and gt == [1, 1])

    total = hom_ref + het + hom_alt
    if total == 0:
        return 1.0

    # Expected frequencies under HWE
    p = 1 - alt_freq  # REF frequency
    q = alt_freq      # ALT frequency

    exp_hom_ref = total * p * p
    exp_het = total * 2 * p * q
    exp_hom_alt = total * q * q

    # Chi-square test
    chi2 = 0
    if exp_hom_ref > 0:
        chi2 += (hom_ref - exp_hom_ref) ** 2 / exp_hom_ref
    if exp_het > 0:
        chi2 += (het - exp_het) ** 2 / exp_het
    if exp_hom_alt > 0:
        chi2 += (hom_alt - exp_hom_alt) ** 2 / exp_hom_alt

    # Convert to p-value (simplified - use scipy.stats.chi2 for exact)
    return max(1e-10, 1 - min(0.999, chi2 / 10))  # Rough approximation

# Usage
# clean_vars, good_samples = gwas_qc_pipeline(ds, ["chr22"], sample_list, "gwas_qc")
```

### GWAS Data Export
```python
def export_gwas_data(ds, regions, samples, phenotype_file, output_prefix):
    """Export data in GWAS-ready formats"""
    # Query clean variant data
    df = ds.read(
        attrs=["sample_name", "contig", "pos_start", "id", "alleles", "fmt_GT"],
        regions=regions,
        samples=samples
    )

    # Load phenotype data
    pheno_df = pd.read_csv(phenotype_file)  # sample_id, phenotype, covariates

    # Convert to PLINK format
    plink_data = convert_to_plink_format(df, pheno_df)

    # Save PLINK files
    save_plink_files(plink_data, f"{output_prefix}_plink")

    # Export for other GWAS tools
    export_for_saige(df, pheno_df, f"{output_prefix}_saige")
    export_for_regenie(df, pheno_df, f"{output_prefix}_regenie")

    print(f"GWAS data exported with prefix: {output_prefix}")

def convert_to_plink_format(variant_df, phenotype_df):
    """Convert variant data to PLINK format"""
    # Implementation depends on specific PLINK format requirements
    # This is a simplified version
    pass

def save_plink_files(data, prefix):
    """Save PLINK .bed/.bim/.fam files"""
    # Implementation for PLINK binary format
    pass
```

## Population Structure Analysis

### Principal Component Analysis
```python
def population_pca(ds, regions, samples, output_prefix, n_components=10):
    """Perform PCA for population structure analysis"""
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    print("Performing population structure PCA...")

    # Query variant data for PCA (use common variants)
    df = ds.read(
        attrs=["sample_name", "contig", "pos_start", "alleles", "fmt_GT", "info_AF"],
        regions=regions,
        samples=samples
    )

    # Filter for common variants (MAF > 5%)
    common_variants = df[df['info_AF'].between(0.05, 0.95)]
    print(f"Using {len(common_variants):,} common variants for PCA")

    # Create genotype matrix
    genotype_matrix = create_genotype_matrix(common_variants)

    # Perform PCA
    scaler = StandardScaler()
    scaled_genotypes = scaler.fit_transform(genotype_matrix)

    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(scaled_genotypes)

    # Create results DataFrame
    pca_df = pd.DataFrame(
        pca_result,
        columns=[f'PC{i+1}' for i in range(n_components)],
        index=samples
    )
    pca_df['sample'] = samples

    # Save results
    pca_df.to_csv(f"{output_prefix}_pca.csv", index=False)

    # Plot variance explained
    variance_explained = pd.DataFrame({
        'PC': range(1, n_components + 1),
        'variance_explained': pca.explained_variance_ratio_,
        'cumulative_variance': np.cumsum(pca.explained_variance_ratio_)
    })
    variance_explained.to_csv(f"{output_prefix}_variance_explained.csv", index=False)

    print(f"PCA complete. First 3 PCs explain {variance_explained['cumulative_variance'].iloc[2]:.1%} of variance")

    return pca_df, variance_explained

def create_genotype_matrix(variant_df):
    """Create numerical genotype matrix for PCA"""
    # Pivot to create sample x variant matrix
    genotype_data = []

    for sample in variant_df['sample_name'].unique():
        sample_data = variant_df[variant_df['sample_name'] == sample]
        genotype_row = []

        for _, variant in sample_data.iterrows():
            gt = variant['fmt_GT']
            if isinstance(gt, list) and len(gt) == 2 and -1 not in gt:
                # Count alternative alleles (0, 1, or 2)
                alt_count = sum(1 for allele in gt if allele > 0)
                genotype_row.append(alt_count)
            else:
                # Missing data - use population mean
                genotype_row.append(-1)  # Mark for imputation

        genotype_data.append(genotype_row)

    # Convert to array and handle missing data
    genotype_matrix = np.array(genotype_data, dtype=float)

    # Simple mean imputation for missing values
    for col in range(genotype_matrix.shape[1]):
        col_data = genotype_matrix[:, col]
        valid_values = col_data[col_data >= 0]
        if len(valid_values) > 0:
            mean_val = np.mean(valid_values)
            genotype_matrix[col_data < 0, col] = mean_val

    return genotype_matrix
```

### Population Assignment
```python
def assign_populations(ds, regions, reference_populations, query_samples, output_file):
    """Assign samples to populations based on genetic similarity"""
    print("Performing population assignment...")

    # Load reference population data
    ref_df = pd.read_csv(reference_populations)  # sample, population, PC1, PC2, etc.

    # Query data for assignment
    query_df = ds.read(
        attrs=["sample_name", "contig", "pos_start", "alleles", "fmt_GT"],
        regions=regions,
        samples=query_samples
    )

    # Perform PCA including reference and query samples
    combined_samples = list(ref_df['sample']) + query_samples
    pca_result, _ = population_pca(ds, regions, combined_samples, "temp_pca")

    # Assign populations using k-nearest neighbors
    from sklearn.neighbors import KNeighborsClassifier

    # Prepare training data (reference populations)
    ref_pcs = pca_result[pca_result['sample'].isin(ref_df['sample'])]
    ref_pops = ref_df.set_index('sample').loc[ref_pcs['sample'], 'population']

    # Train classifier
    knn = KNeighborsClassifier(n_neighbors=5)
    knn.fit(ref_pcs[['PC1', 'PC2', 'PC3']], ref_pops)

    # Predict populations for query samples
    query_pcs = pca_result[pca_result['sample'].isin(query_samples)]
    predicted_pops = knn.predict(query_pcs[['PC1', 'PC2', 'PC3']])
    prediction_probs = knn.predict_proba(query_pcs[['PC1', 'PC2', 'PC3']])

    # Create results
    assignment_results = pd.DataFrame({
        'sample': query_samples,
        'predicted_population': predicted_pops,
        'confidence': np.max(prediction_probs, axis=1)
    })

    assignment_results.to_csv(output_file, index=False)
    print(f"Population assignments saved to {output_file}")

    return assignment_results
```

## Rare Variant Analysis

### Burden Testing
```python
def rare_variant_burden_test(ds, regions, cases, controls, gene_annotations, output_file):
    """Perform rare variant burden testing"""
    print("Performing rare variant burden testing...")

    # Load gene annotations
    genes_df = pd.read_csv(gene_annotations)  # gene, chr, start, end

    burden_results = []

    for _, gene in genes_df.iterrows():
        gene_region = f"{gene['chr']}:{gene['start']}-{gene['end']}"

        print(f"Testing gene: {gene['gene']}")

        # Query rare variants in gene region
        df = ds.read(
            attrs=["sample_name", "pos_start", "alleles", "fmt_GT", "info_AF"],
            regions=[gene_region],
            samples=cases + controls
        )

        # Filter for rare variants (MAF < 1%)
        rare_variants = df[df['info_AF'] < 0.01]

        if len(rare_variants) == 0:
            continue

        # Calculate burden scores
        case_burden = calculate_burden_score(rare_variants, cases)
        control_burden = calculate_burden_score(rare_variants, controls)

        # Perform statistical test
        from scipy.stats import mannwhitneyu
        statistic, p_value = mannwhitneyu(case_burden, control_burden, alternative='two-sided')

        burden_results.append({
            'gene': gene['gene'],
            'chr': gene['chr'],
            'start': gene['start'],
            'end': gene['end'],
            'rare_variant_count': len(rare_variants),
            'case_mean_burden': np.mean(case_burden),
            'control_mean_burden': np.mean(control_burden),
            'p_value': p_value,
            'statistic': statistic
        })

    # Adjust for multiple testing
    results_df = pd.DataFrame(burden_results)
    if len(results_df) > 0:
        from statsmodels.stats.multitest import multipletests
        _, corrected_p, _, _ = multipletests(results_df['p_value'], method='bonferroni')
        results_df['corrected_p_value'] = corrected_p

        results_df.to_csv(output_file, index=False)
        print(f"Burden test results saved to {output_file}")

    return results_df

def calculate_burden_score(variant_df, samples):
    """Calculate burden score for each sample"""
    burden_scores = []

    for sample in samples:
        sample_variants = variant_df[variant_df['sample_name'] == sample]
        score = 0

        for _, variant in sample_variants.iterrows():
            gt = variant['fmt_GT']
            if isinstance(gt, list) and len(gt) == 2 and -1 not in gt:
                # Count alternative alleles
                alt_count = sum(1 for allele in gt if allele > 0)
                score += alt_count

        burden_scores.append(score)

    return burden_scores
```

### Variant Annotation Integration
```python
def annotate_rare_variants(ds, regions, annotation_file, output_file):
    """Annotate rare variants with functional predictions"""
    print("Annotating rare variants...")

    # Query rare variants
    df = ds.read(
        attrs=["contig", "pos_start", "pos_end", "alleles", "id", "info_AF"],
        regions=regions
    )

    rare_variants = df[df['info_AF'] < 0.05]  # 5% threshold

    # Load functional annotations
    annotations = pd.read_csv(annotation_file)  # chr, pos, consequence, sift, polyphen

    # Merge with annotations
    annotated = rare_variants.merge(
        annotations,
        left_on=['contig', 'pos_start'],
        right_on=['chr', 'pos'],
        how='left'
    )

    # Classify by functional impact
    def classify_impact(row):
        if pd.isna(row['consequence']):
            return 'unknown'
        elif 'stop_gained' in row['consequence'] or 'frameshift' in row['consequence']:
            return 'high_impact'
        elif 'missense' in row['consequence']:
            return 'moderate_impact'
        else:
            return 'low_impact'

    annotated['impact'] = annotated.apply(classify_impact, axis=1)

    # Save annotated variants
    annotated.to_csv(output_file, index=False)

    # Summary statistics
    impact_summary = annotated['impact'].value_counts()
    print("Functional impact summary:")
    for impact, count in impact_summary.items():
        print(f"  {impact}: {count:,} variants")

    return annotated
```

## Allele Frequency Analysis

### Population-Specific Allele Frequencies
```python
def calculate_population_allele_frequencies(ds, regions, population_file, output_prefix):
    """Calculate allele frequencies for each population"""
    print("Calculating population-specific allele frequencies...")

    # Load population assignments
    pop_df = pd.read_csv(population_file)  # sample, population

    populations = pop_df['population'].unique()
    results = {}

    for population in populations:
        print(f"Processing population: {population}")

        pop_samples = pop_df[pop_df['population'] == population]['sample'].tolist()

        # Query variant data
        df = ds.read(
            attrs=["contig", "pos_start", "alleles", "fmt_GT"],
            regions=regions,
            samples=pop_samples
        )

        # Calculate allele frequencies
        af_results = []

        for (chrom, pos, alleles), group in df.groupby(['contig', 'pos_start', 'alleles']):
            allele_counts = {i: 0 for i in range(len(alleles))}
            total_alleles = 0

            for gt in group['fmt_GT']:
                if isinstance(gt, list) and -1 not in gt:
                    for allele in gt:
                        if allele in allele_counts:
                            allele_counts[allele] += 1
                            total_alleles += 1

            if total_alleles > 0:
                frequencies = {alleles[i]: count/total_alleles
                             for i, count in allele_counts.items()
                             if i < len(alleles)}

                af_results.append({
                    'chr': chrom,
                    'pos': pos,
                    'ref': alleles[0],
                    'alt': ','.join(alleles[1:]) if len(alleles) > 1 else '',
                    'ref_freq': frequencies.get(alleles[0], 0),
                    'alt_freq': sum(frequencies.get(alleles[i], 0) for i in range(1, len(alleles))),
                    'sample_count': len(group)
                })

        af_df = pd.DataFrame(af_results)
        af_df.to_csv(f"{output_prefix}_{population}_frequencies.csv", index=False)

        results[population] = af_df
        print(f"  Calculated frequencies for {len(af_df):,} variants")

    return results

# Compare allele frequencies between populations
def compare_population_frequencies(pop_frequencies, output_file):
    """Compare allele frequencies between populations"""
    print("Comparing population allele frequencies...")

    populations = list(pop_frequencies.keys())
    if len(populations) < 2:
        print("Need at least 2 populations for comparison")
        return

    # Merge frequency data
    merged = pop_frequencies[populations[0]].copy()
    merged = merged.rename(columns={'alt_freq': f'{populations[0]}_freq'})

    for pop in populations[1:]:
        pop_df = pop_frequencies[pop][['chr', 'pos', 'alt_freq']].copy()
        pop_df = pop_df.rename(columns={'alt_freq': f'{pop}_freq'})

        merged = merged.merge(pop_df, on=['chr', 'pos'], how='inner')

    # Calculate Fst and frequency differences
    for i, pop1 in enumerate(populations):
        for pop2 in populations[i+1:]:
            freq1_col = f'{pop1}_freq'
            freq2_col = f'{pop2}_freq'

            # Frequency difference
            merged[f'freq_diff_{pop1}_{pop2}'] = abs(merged[freq1_col] - merged[freq2_col])

            # Simplified Fst calculation
            merged[f'fst_{pop1}_{pop2}'] = calculate_fst(merged[freq1_col], merged[freq2_col])

    merged.to_csv(output_file, index=False)
    print(f"Population frequency comparisons saved to {output_file}")

    return merged

def calculate_fst(freq1, freq2):
    """Calculate simplified Fst between two populations"""
    # Simplified Fst = (Ht - Hs) / Ht where Ht is total het, Hs is within-pop het
    p_total = (freq1 + freq2) / 2
    ht = 2 * p_total * (1 - p_total)
    hs = (2 * freq1 * (1 - freq1) + 2 * freq2 * (1 - freq2)) / 2

    fst = np.where(ht > 0, (ht - hs) / ht, 0)
    return np.clip(fst, 0, 1)  # Fst should be between 0 and 1
```

## Integration with Analysis Tools

### SAIGE Integration
```python
def prepare_saige_input(ds, regions, samples, phenotype_file, output_prefix):
    """Prepare input files for SAIGE analysis"""
    # Export genotype data in PLINK format
    export_plink_format(ds, regions, samples, f"{output_prefix}_geno")

    # Prepare phenotype file
    pheno_df = pd.read_csv(phenotype_file)
    saige_pheno = pheno_df[['sample_id', 'phenotype']].copy()
    saige_pheno.columns = ['IID', 'y_binary']
    saige_pheno['FID'] = saige_pheno['IID']  # Family ID same as individual ID
    saige_pheno = saige_pheno[['FID', 'IID', 'y_binary']]

    saige_pheno.to_csv(f"{output_prefix}_phenotype.txt", sep='\t', index=False)

    print(f"SAIGE input files prepared with prefix: {output_prefix}")

def export_plink_format(ds, regions, samples, output_prefix):
    """Export data in PLINK format for SAIGE"""
    # Implementation specific to PLINK .bed/.bim/.fam format
    # This is a placeholder - actual implementation would depend on
    # specific PLINK format requirements
    pass
```

### Integration with Cloud Computing
```python
def run_cloud_gwas_pipeline(ds, regions, samples, phenotype_file, cloud_config):
    """Run GWAS pipeline on cloud computing platform"""
    import subprocess

    # Export data to cloud storage
    cloud_prefix = cloud_config['storage_prefix']

    print("Exporting data to cloud...")
    ds.export_bcf(
        uri=f"{cloud_prefix}/variants.bcf",
        regions=regions,
        samples=samples
    )

    # Upload phenotype file
    subprocess.run([
        'aws', 's3', 'cp', phenotype_file,
        f"{cloud_prefix}/phenotypes.txt"
    ])

    # Submit GWAS job
    job_config = {
        'job_name': 'gwas_analysis',
        'input_data': f"{cloud_prefix}/variants.bcf",
        'phenotypes': f"{cloud_prefix}/phenotypes.txt",
        'output_path': f"{cloud_prefix}/results/",
        'instance_type': 'r5.xlarge',
        'max_runtime': '2h'
    }

    print("Submitting GWAS job to cloud...")
    # Implementation would depend on specific cloud platform
    # (AWS Batch, Google Cloud, etc.)

    return job_config
```

## Best Practices for Population Genomics

### Workflow Optimization
```python
def optimized_population_workflow():
    """Best practices for population genomics with TileDB-VCF"""
    best_practices = {
        'data_organization': [
            "Organize by population/cohort in separate datasets",
            "Use consistent sample naming conventions",
            "Maintain metadata in separate files"
        ],

        'quality_control': [
            "Perform sample QC before variant QC",
            "Use population-appropriate QC thresholds",
            "Document all filtering steps"
        ],

        'analysis_strategy': [
            "Use common variants for PCA/population structure",
            "Filter rare variants by functional annotation",
            "Validate results in independent cohorts"
        ],

        'computational_efficiency': [
            "Process chromosomes in parallel",
            "Use appropriate memory budgets for large cohorts",
            "Cache intermediate results for iterative analysis"
        ],

        'reproducibility': [
            "Version control all analysis scripts",
            "Document software versions and parameters",
            "Use containerized environments for complex analyses"
        ]
    }

    return best_practices

def monitoring_pipeline_progress(total_steps):
    """Monitor progress of long-running population genomics pipelines"""
    import time
    from datetime import datetime

    def log_step(step_num, description, start_time=None):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        if start_time:
            duration = time.time() - start_time
            print(f"[{timestamp}] Step {step_num}/{total_steps}: {description} - Completed in {duration:.1f}s")
        else:
            print(f"[{timestamp}] Step {step_num}/{total_steps}: {description} - Starting...")

        return time.time()

    return log_step

# Usage example
# log_step = monitoring_pipeline_progress(5)
# start_time = log_step(1, "Loading data")
# # ... do work ...
# log_step(1, "Loading data", start_time)
```

This comprehensive population genomics guide demonstrates how to leverage TileDB-VCF for complex analyses from GWAS preparation through rare variant studies and population structure analysis.