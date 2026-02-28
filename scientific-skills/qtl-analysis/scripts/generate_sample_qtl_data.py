# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

"""
Generate sample QTL and genomic data for testing visualizations
"""
import pandas as pd
import numpy as np
import os

def generate_phenotype_data(n_samples=200, seed=42):
    """Generate sample phenotype data with genotypes and environments"""
    np.random.seed(seed)
    
    data = pd.DataFrame({
        'Sample_ID': [f'Sample_{i:03d}' for i in range(1, n_samples + 1)],
        'Genotype': np.random.choice(['AA', 'AB', 'BB'], n_samples),
        'Environment': np.random.choice(['Control', 'Stress', 'Optimal'], n_samples),
        'Plant_Height': np.concatenate([
            np.random.normal(100, 15, int(n_samples * 0.35)),
            np.random.normal(120, 18, int(n_samples * 0.35)),
            np.random.normal(95, 12, n_samples - int(n_samples * 0.7))
        ]),
        'Yield': np.concatenate([
            np.random.normal(50, 10, int(n_samples * 0.35)),
            np.random.normal(65, 12, int(n_samples * 0.35)),
            np.random.normal(45, 8, n_samples - int(n_samples * 0.7))
        ]),
        'Quality': np.random.normal(85, 8, n_samples),
        'Disease_Resistance': np.random.normal(60, 20, n_samples),
        'Flowering_Time': np.random.normal(45, 5, n_samples)
    })
    
    return data

def generate_gwas_results(n_snps=10000, n_chrom=5, seed=42):
    """Generate simulated GWAS results"""
    np.random.seed(seed)
    
    snp_data = []
    snp_counter = 1
    
    for chrom in range(1, n_chrom + 1):
        n_markers = np.random.randint(1500, 2500)
        positions = np.sort(np.random.randint(0, 100000000, n_markers))
        
        # Generate p-values with some significant hits
        if chrom == 3:
            p_values = np.random.uniform(0.0001, 0.05, n_markers)
            # Add peak
            peak_idx = n_markers // 2
            p_values[peak_idx-10:peak_idx+10] = np.random.uniform(1e-10, 1e-6, 20)
        elif chrom == 5:
            p_values = np.random.uniform(0.001, 0.1, n_markers)
            peak_idx = int(n_markers * 0.6)
            p_values[peak_idx-15:peak_idx+15] = np.random.uniform(1e-8, 1e-5, 30)
        else:
            p_values = np.random.uniform(0.1, 1.0, n_markers)
        
        for i in range(n_markers):
            snp_data.append({
                'SNP': f'rs{snp_counter + i}',
                'CHR': chrom,
                'BP': positions[i],
                'P': p_values[i]
            })
        
        snp_counter += n_markers
    
    return pd.DataFrame(snp_data)

def generate_qtl_scan(n_chrom=5, seed=42):
    """Generate QTL scan results with LOD scores"""
    np.random.seed(seed)
    
    scan_data = []
    
    for chrom in range(1, n_chrom + 1):
        chrom_length = np.random.randint(80, 150)
        n_markers = np.random.randint(50, 100)
        positions = np.linspace(0, chrom_length, n_markers)
        
        # Generate LOD scores
        if chrom == 3:
            lod = np.random.exponential(2, n_markers)
            peak_pos = 45
            lod += 15 * np.exp(-((positions - peak_pos) / 10) ** 2)
        elif chrom == 5:
            lod = np.random.exponential(1.5, n_markers)
            peak_pos = 60
            lod += 8 * np.exp(-((positions - peak_pos) / 8) ** 2)
        else:
            lod = np.random.exponential(1, n_markers)
        
        for i in range(n_markers):
            scan_data.append({
                'Chromosome': chrom,
                'Position': positions[i],
                'LOD': lod[i],
                'Marker': f'Marker_{chrom}_{i+1}'
            })
    
    return pd.DataFrame(scan_data)

def generate_genotype_matrix(n_samples=50, n_markers=20, seed=42):
    """Generate genotype matrix (0, 1, 2 coding)"""
    np.random.seed(seed)
    
    genotype_matrix = np.random.choice([0, 1, 2], size=(n_samples, n_markers))
    
    df = pd.DataFrame(
        genotype_matrix,
        index=[f'Sample_{i:03d}' for i in range(1, n_samples + 1)],
        columns=[f'Marker_{i+1}' for i in range(n_markers)]
    )
    
    return df

def generate_gxe_data(n_genotypes=5, n_envs=4, seed=42):
    """Generate genotype x environment data"""
    np.random.seed(seed)
    
    genotypes = [f'G{i}' for i in range(1, n_genotypes + 1)]
    environments = [f'Env{i}' for i in range(1, n_envs + 1)]
    
    gxe_data = []
    for geno in genotypes:
        for env in environments:
            base = {'G1': 50, 'G2': 55, 'G3': 45, 'G4': 60, 'G5': 52}[geno]
            env_eff = {'Env1': 0, 'Env2': 5, 'Env3': -3, 'Env4': 8}[env]
            
            # Add GxE interaction
            if geno == 'G4' and env == 'Env4':
                interaction = 10
            elif geno == 'G3' and env == 'Env3':
                interaction = -8
            else:
                interaction = np.random.normal(0, 3)
            
            yield_val = base + env_eff + interaction + np.random.normal(0, 4)
            
            gxe_data.append({
                'Genotype': geno,
                'Environment': env,
                'Yield': yield_val,
                'Biomass': yield_val * 0.8 + np.random.normal(0, 5),
                'Quality': np.random.normal(80, 10)
            })
    
    return pd.DataFrame(gxe_data)

def generate_effect_sizes(n_qtl=8, n_traits=4, seed=42):
    """Generate QTL effect size data"""
    np.random.seed(seed)
    
    traits = ['Height', 'Yield', 'Quality', 'Disease_Res']
    effects = []
    
    for i in range(1, n_qtl + 1):
        for trait in traits:
            effects.append({
                'QTL': f'QTL_{i}',
                'Trait': trait,
                'Effect': np.random.normal(0, 2),
                'SE': np.random.uniform(0.5, 1.5),
                'Chromosome': np.random.randint(1, 6),
                'Position': np.random.uniform(0, 100)
            })
    
    df = pd.DataFrame(effects)
    df['Lower'] = df['Effect'] - 1.96 * df['SE']
    df['Upper'] = df['Effect'] + 1.96 * df['SE']
    
    return df

def main():
    """Generate all sample datasets"""
    output_dir = 'sample_qtl_data'
    os.makedirs(output_dir, exist_ok=True)
    
    print("Generating sample QTL datasets...")
    
    # 1. Phenotype data
    print("  - phenotype_data.csv")
    phenotype_data = generate_phenotype_data(n_samples=200)
    phenotype_data.to_csv(f'{output_dir}/phenotype_data.csv', index=False)
    
    # 2. GWAS results
    print("  - gwas_results.csv")
    gwas_results = generate_gwas_results(n_snps=10000)
    gwas_results.to_csv(f'{output_dir}/gwas_results.csv', index=False)
    
    # 3. QTL scan
    print("  - qtl_scan.csv")
    qtl_scan = generate_qtl_scan()
    qtl_scan.to_csv(f'{output_dir}/qtl_scan.csv', index=False)
    
    # 4. Genotype matrix
    print("  - genotype_matrix.csv")
    genotype_matrix = generate_genotype_matrix()
    genotype_matrix.to_csv(f'{output_dir}/genotype_matrix.csv')
    
    # 5. GxE data
    print("  - gxe_data.csv")
    gxe_data = generate_gxe_data()
    gxe_data.to_csv(f'{output_dir}/gxe_data.csv', index=False)
    
    # 6. Effect sizes
    print("  - effect_sizes.csv")
    effect_sizes = generate_effect_sizes()
    effect_sizes.to_csv(f'{output_dir}/effect_sizes.csv', index=False)
    
    print(f"\nAll datasets saved to '{output_dir}/' directory")
    print("\nDataset summaries:")
    print(f"  Phenotype data: {len(phenotype_data)} samples, {len(phenotype_data.columns)} traits")
    print(f"  GWAS results: {len(gwas_results)} SNPs across {gwas_results['CHR'].nunique()} chromosomes")
    print(f"  QTL scan: {len(qtl_scan)} markers")
    print(f"  Genotype matrix: {genotype_matrix.shape[0]} samples × {genotype_matrix.shape[1]} markers")
    print(f"  GxE data: {len(gxe_data)} genotype-environment combinations")
    print(f"  Effect sizes: {len(effect_sizes)} QTL-trait combinations")

if __name__ == '__main__':
    main()
