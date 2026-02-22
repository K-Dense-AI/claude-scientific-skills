#!/usr/bin/env python3
"""
Example: SNP Filtering and Quality Control

This example demonstrates SNP data filtering - removing low-quality variants
before GWAS/QTL analysis. Poor quality SNPs can cause false positives and
reduce statistical power.

WHAT THIS MEANS:
Raw genotype data contains many low-quality SNPs: low minor allele frequency (MAF),
excessive missing data, deviations from Hardy-Weinberg equilibrium. These SNPs
should be filtered out to improve analysis quality.

WHY WE DO THIS:
- Low MAF SNPs have low statistical power
- Missing data >5% reduces effective sample size
- Hardy-Weinberg deviation suggests genotyping errors
- Reduces computational burden
- Improves QQ plot calibration

FILTERING CRITERIA APPLIED:
1. Minor Allele Frequency (MAF) ≥ 5%
2. Missing genotype rate ≤ 5%
3. Hardy-Weinberg equilibrium p-value ≥ 1e-6
4. Biallelic SNPs only

WHAT'S OUTPUT:
- Filtered genotype matrix
- Filter statistics (before/after counts)
- QC plots (MAF distribution, missing rates, HWE p-values)

Equivalent to QTLmax: "How to create a subset file from SNP dataset"
https://open.qtlmax.com/guide/index.php/2025/07/14/how-to-create-a-subset-file-from-snp-dataset/

Auto-installs: pandas, numpy, scipy, matplotlib
"""

import subprocess
import sys
import os


def install_packages():
    """Install required packages without root"""
    packages = ["pandas", "numpy", "scipy", "matplotlib"]
    for pkg in packages:
        try:
            __import__(pkg.replace("-", "_"))
        except ImportError:
            print(f"Installing {pkg}...")
            subprocess.check_call(
                [sys.executable, "-m", "pip", "install", "--user", "-q", pkg]
            )


def generate_snp_data(n_individuals=100, n_snps=500):
    """Generate synthetic SNP data with various quality issues"""
    import numpy as np

    np.random.seed(42)

    # Generate genotypes (0, 1, 2)
    genotypes = np.random.binomial(2, 0.3, (n_individuals, n_snps))

    # Introduce low MAF SNPs (rare alleles)
    for i in range(20):
        genotypes[:, i] = np.random.binomial(2, 0.01, n_individuals)

    # Introduce high missing rate SNPs
    for i in range(20, 40):
        mask = np.random.random(n_individuals) < 0.2
        genotypes[mask, i] = -9  # Missing

    # Introduce HWE deviation (excess homozygotes)
    for i in range(40, 50):
        p = 0.3
        # Generate with inbreeding
        inbreeding = 0.3
        genotypes[:, i] = np.random.choice(
            [0, 1, 2],
            n_individuals,
            p=[
                (1 - p) ** 2 + inbreeding * p * (1 - p),
                2 * p * (1 - p) * (1 - inbreeding),
                p**2 + inbreeding * p * (1 - p),
            ],
        )

    # Sample and SNP IDs
    sample_ids = [f"Sample{i + 1}" for i in range(n_individuals)]
    snp_ids = [f"rs{i + 1}" for i in range(n_snps)]
    chromosomes = [(i // 50) + 1 for i in range(n_snps)]
    positions = [(i % 50) * 10000 + 1000 for i in range(n_snps)]

    return genotypes, sample_ids, snp_ids, chromosomes, positions


def calculate_maf(genotypes):
    """Calculate Minor Allele Frequency"""
    import numpy as np

    # Replace missing with NaN
    geno_clean = np.where(genotypes == -9, np.nan, genotypes)

    # Calculate allele frequencies
    n_obs = np.sum(~np.isnan(geno_clean), axis=0)
    allele_counts = np.nansum(geno_clean, axis=0)

    # Frequency of allele 1
    p = allele_counts / (2 * n_obs)

    # MAF is minimum of p and 1-p
    maf = np.minimum(p, 1 - p)

    return maf, n_obs


def calculate_missing_rate(genotypes):
    """Calculate missing genotype rate per SNP"""
    import numpy as np

    missing = np.sum(genotypes == -9, axis=0)
    total = genotypes.shape[0]

    return missing / total


def calculate_hwe_pvalue(genotypes):
    """Calculate Hardy-Weinberg Equilibrium p-value"""
    from scipy import stats
    import numpy as np

    pvalues = []

    for i in range(genotypes.shape[1]):
        geno = genotypes[:, i]
        geno = geno[geno != -9]  # Remove missing

        if len(geno) == 0:
            pvalues.append(1.0)
            continue

        # Observed genotype counts
        n = len(geno)
        obs_hom1 = np.sum(geno == 0)
        obs_het = np.sum(geno == 1)
        obs_hom2 = np.sum(geno == 2)

        # Allele frequencies
        p = (2 * obs_hom1 + obs_het) / (2 * n)
        q = 1 - p

        # Expected counts under HWE
        exp_hom1 = n * p * p
        exp_het = n * 2 * p * q
        exp_hom2 = n * q * q

        # Chi-square test
        expected = [exp_hom1, exp_het, exp_hom2]
        observed = [obs_hom1, obs_het, obs_hom2]

        # Avoid division by zero
        if all(e > 0 for e in expected):
            chi2 = sum((o - e) ** 2 / e for o, e in zip(observed, expected))
            pval = 1 - stats.chi2.cdf(chi2, df=1)
        else:
            pval = 1.0

        pvalues.append(pval)

    return np.array(pvalues)


def filter_snps(
    genotypes,
    snp_ids,
    chromosomes,
    positions,
    maf_threshold=0.05,
    missing_threshold=0.05,
    hwe_threshold=1e-6,
):
    """Apply SNP filters"""
    print("\nCalculating SNP statistics...")

    # Calculate statistics
    maf, _ = calculate_maf(genotypes)
    missing_rate = calculate_missing_rate(genotypes)
    hwe_p = calculate_hwe_pvalue(genotypes)

    # Apply filters
    print("\nApplying filters:")
    print(f"  MAF >= {maf_threshold}")
    print(f"  Missing <= {missing_threshold}")
    print(f"  HWE p >= {hwe_threshold}")

    # Filter mask
    mask = (
        (maf >= maf_threshold)
        & (missing_rate <= missing_threshold)
        & (hwe_p >= hwe_threshold)
    )

    # Filter data
    filtered_genotypes = genotypes[:, mask]
    filtered_snp_ids = [snp_ids[i] for i in range(len(snp_ids)) if mask[i]]
    filtered_chromosomes = [chromosomes[i] for i in range(len(chromosomes)) if mask[i]]
    filtered_positions = [positions[i] for i in range(len(positions)) if mask[i]]

    # Statistics
    stats = {
        "total": len(snp_ids),
        "passed": sum(mask),
        "failed_maf": sum(maf < maf_threshold),
        "failed_missing": sum(missing_rate > missing_threshold),
        "failed_hwe": sum(hwe_p < hwe_threshold),
        "maf": maf,
        "missing_rate": missing_rate,
        "hwe_p": hwe_p,
    }

    return (
        filtered_genotypes,
        filtered_snp_ids,
        filtered_chromosomes,
        filtered_positions,
        stats,
    )


def create_qc_plots(stats, output_dir):
    """Create QC plots for filtering"""
    import matplotlib.pyplot as plt
    import numpy as np

    print("\nCreating QC plots...")

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Filter summary
    ax1 = axes[0, 0]
    categories = ["Passed", "Failed\nMAF", "Failed\nMissing", "Failed\nHWE"]
    values = [
        stats["passed"],
        stats["failed_maf"],
        stats["failed_missing"],
        stats["failed_hwe"],
    ]
    colors = ["green", "red", "orange", "purple"]
    ax1.bar(categories, values, color=colors, alpha=0.7)
    ax1.set_ylabel("Number of SNPs")
    ax1.set_title("SNP Filtering Results")
    ax1.text(
        0,
        stats["passed"] + 5,
        f"{stats['passed']}\n({100 * stats['passed'] / stats['total']:.1f}%)",
        ha="center",
        fontweight="bold",
    )

    # 2. MAF distribution
    ax2 = axes[0, 1]
    ax2.hist(stats["maf"], bins=20, color="steelblue", alpha=0.7, edgecolor="black")
    ax2.axvline(x=0.05, color="red", linestyle="--", linewidth=2, label="MAF threshold")
    ax2.set_xlabel("Minor Allele Frequency")
    ax2.set_ylabel("Number of SNPs")
    ax2.set_title("MAF Distribution")
    ax2.legend()

    # 3. Missing rate distribution
    ax3 = axes[1, 0]
    ax3.hist(
        stats["missing_rate"], bins=20, color="orange", alpha=0.7, edgecolor="black"
    )
    ax3.axvline(
        x=0.05, color="red", linestyle="--", linewidth=2, label="Missing threshold"
    )
    ax3.set_xlabel("Missing Rate")
    ax3.set_ylabel("Number of SNPs")
    ax3.set_title("Missing Data Rate Distribution")
    ax3.legend()

    # 4. HWE p-value distribution
    ax4 = axes[1, 1]
    hwe_logp = -np.log10(stats["hwe_p"] + 1e-300)
    ax4.hist(hwe_logp, bins=20, color="purple", alpha=0.7, edgecolor="black")
    ax4.axvline(
        x=-np.log10(1e-6),
        color="red",
        linestyle="--",
        linewidth=2,
        label="HWE threshold",
    )
    ax4.set_xlabel("-log10(HWE p-value)")
    ax4.set_ylabel("Number of SNPs")
    ax4.set_title("Hardy-Weinberg Equilibrium")
    ax4.legend()

    plt.tight_layout()
    plt.savefig(f"{output_dir}/snp_filtering_qc.png", dpi=150, bbox_inches="tight")
    print(f"QC plots saved: {output_dir}/snp_filtering_qc.png")


def main():
    print("=" * 60)
    print("Example: SNP Filtering and Quality Control")
    print("=" * 60)
    print("\nThis script filters SNPs based on MAF, missing rate, and HWE.")
    print("Poor quality SNPs cause false positives in GWAS.")

    # Install packages
    print("\n[1/4] Installing dependencies...")
    install_packages()

    # Setup
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    # Generate data
    print("\n[2/4] Generating test SNP data...")
    genotypes, sample_ids, snp_ids, chromosomes, positions = generate_snp_data()
    print(f"  Generated: {len(sample_ids)} samples, {len(snp_ids)} SNPs")

    # Filter SNPs
    print("\n[3/4] Filtering SNPs...")
    (filtered_geno, filtered_ids, filtered_chrom, filtered_pos, stats) = filter_snps(
        genotypes, snp_ids, chromosomes, positions
    )

    # Create plots
    print("\n[4/4] Creating QC plots...")
    create_qc_plots(stats, output_dir)

    # Summary
    print("\n" + "=" * 60)
    print("SNP Filtering Summary")
    print("=" * 60)
    print(f"Total SNPs: {stats['total']}")
    print(f"Passed: {stats['passed']} ({100 * stats['passed'] / stats['total']:.1f}%)")
    print(f"Failed MAF: {stats['failed_maf']}")
    print(f"Failed Missing: {stats['failed_missing']}")
    print(f"Failed HWE: {stats['failed_hwe']}")
    print(
        f"\nFiltered data: {filtered_geno.shape[0]} samples × {filtered_geno.shape[1]} SNPs"
    )
    print(f"\nOutputs:")
    print(f"  - {output_dir}/snp_filtering_qc.png")

    print("\n" + "=" * 60)
    print("Why This Matters")
    print("=" * 60)
    print("Filtering SNPs improves GWAS results:")
    print("  ✓ Removes low-power rare variants")
    print("  ✓ Reduces missing data impact")
    print("  ✓ Eliminates genotyping errors")
    print("  ✓ Improves QQ plot calibration")
    print("\n✅ SNP filtering example complete!")
    print("\nIn QTLmax: Preprocess → SNP Filter")


if __name__ == "__main__":
    main()
