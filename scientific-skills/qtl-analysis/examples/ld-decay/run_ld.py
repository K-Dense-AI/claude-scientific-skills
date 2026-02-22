#!/usr/bin/env python3
"""
Example: LD Decay Analysis

This example demonstrates Linkage Disequilibrium (LD) analysis including:
- Pairwise LD calculation (r²)
- LD decay curve
- LD heatmap for genomic window

Equivalent to QTLmax: "How to draw a linkage disequilibrium decay curve"
https://open.qtlmax.com/guide/index.php/2026/02/09/linkage-disequilibrium-decay-curve/

Auto-installs: pandas, numpy, matplotlib, scipy
"""

import subprocess
import sys
import os


def install_packages():
    """Install required packages without root"""
    packages = ["pandas", "numpy", "matplotlib", "scipy"]
    for pkg in packages:
        try:
            __import__(pkg.replace("-", "_"))
        except ImportError:
            print(f"Installing {pkg}...")
            subprocess.check_call(
                [sys.executable, "-m", "pip", "install", "--user", "-q", pkg]
            )


def generate_linkage_data(n_individuals=200, n_snps=500):
    """Generate synthetic genotype data for LD analysis"""
    import numpy as np

    np.random.seed(42)

    # Simulate chromosome with realistic LD structure
    # Use autoregressive model for LD decay
    genotypes = np.zeros((n_individuals, n_snps))

    # First SNP
    genotypes[:, 0] = np.random.binomial(2, 0.3, n_individuals)

    # Subsequent SNPs with decaying correlation
    for i in range(1, n_snps):
        # Decay factor based on distance
        decay = 0.98**i  # Slow decay for demonstration
        # Generate correlated genotype
        parent = genotypes[:, max(0, i - 1)]
        switch_prob = 1 - decay
        new_geno = parent.copy()
        mask = np.random.random(n_individuals) < switch_prob
        new_geno[mask] = np.random.binomial(2, 0.3, mask.sum())
        genotypes[:, i] = new_geno

    # Clip to valid range
    genotypes = np.clip(genotypes, 0, 2).astype(int)

    # Positions (in kb)
    positions = np.arange(n_snps) * 10  # 10kb intervals

    snp_ids = [f"rs{i}" for i in range(n_snps)]

    return genotypes, positions, snp_ids


def calculate_pairwise_ld(genotypes, max_lag=50):
    """Calculate pairwise LD (r²) between nearby SNPs"""
    import numpy as np

    print("Calculating pairwise LD...")

    n_snps = genotypes.shape[1]
    ld_values = []
    distances = []

    for lag in range(1, min(max_lag + 1, n_snps)):
        r2_values = []

        for i in range(n_snps - lag):
            # Calculate r² between SNP i and SNP i+lag
            geno1 = genotypes[:, i]
            geno2 = genotypes[:, i + lag]

            # Skip if either is monomorphic
            if geno1.var() == 0 or geno2.var() == 0:
                continue

            # Calculate correlation (r)
            r = np.corrcoef(geno1, geno2)[0, 1]

            # r²
            if not np.isnan(r):
                r2_values.append(r**2)

        if r2_values:
            ld_values.append(np.mean(r2_values))
            distances.append(lag * 10)  # 10kb intervals

    return np.array(distances), np.array(ld_values)


def calculate_ld_matrix(genotypes, window_size=20):
    """Calculate LD matrix for a genomic window"""
    import numpy as np

    print(f"Calculating LD matrix ({window_size}x{window_size})...")

    n_snps = min(window_size, genotypes.shape[1])
    ld_matrix = np.ones((n_snps, n_snps))

    for i in range(n_snps):
        for j in range(i + 1, n_snps):
            geno1 = genotypes[:, i]
            geno2 = genotypes[:, j]

            if geno1.var() == 0 or geno2.var() == 0:
                continue

            r = np.corrcoef(geno1, geno2)[0, 1]
            if not np.isnan(r):
                ld_matrix[i, j] = r**2
                ld_matrix[j, i] = r**2

    return ld_matrix


def create_ld_plots(distances, ld_values, ld_matrix, output_dir):
    """Create LD visualization plots"""
    import numpy as np
    import matplotlib.pyplot as plt

    print("Creating LD visualizations...")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # 1. LD Decay Curve
    ax1 = axes[0]
    ax1.scatter(distances, ld_values, s=20, alpha=0.6)

    # Fit exponential decay
    from scipy.optimize import curve_fit

    def exp_decay(x, a, b):
        return a * np.exp(-b * x)

    try:
        popt, _ = curve_fit(exp_decay, distances, ld_values, p0=[1, 0.01], maxfev=5000)
        x_fit = np.linspace(0, max(distances), 100)
        y_fit = exp_decay(x_fit, *popt)
        ax1.plot(
            x_fit,
            y_fit,
            "r-",
            linewidth=2,
            label=f"Fit: r² = {popt[0]:.2f} × exp(-{popt[1]:.4f} × d)",
        )
        ax1.legend()
    except:
        pass

    ax1.set_xlabel("Distance (kb)")
    ax1.set_ylabel("Mean r²")
    ax1.set_title("LD Decay Curve")
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, 1)

    # 2. LD Heatmap
    ax2 = axes[1]
    im = ax2.imshow(ld_matrix, cmap="YlOrRd", aspect="auto", vmin=0, vmax=1)
    ax2.set_xlabel("SNP Index")
    ax2.set_ylabel("SNP Index")
    ax2.set_title("LD Heatmap (20 SNP window)")
    plt.colorbar(im, ax=ax2, label="r²")

    plt.tight_layout()
    plt.savefig(f"{output_dir}/ld_analysis.png", dpi=150, bbox_inches="tight")
    print(f"LD plots saved: {output_dir}/ld_analysis.png")


def main():
    print("=" * 60)
    print("Example: Linkage Disequilibrium (LD) Analysis")
    print("=" * 60)

    # Install packages
    print("\n[1/5] Installing dependencies...")
    install_packages()

    # Setup
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    # Generate data
    print("\n[2/5] Generating genotype data...")
    genotypes, positions, snp_ids = generate_linkage_data()
    print(f"  Generated: {genotypes.shape[0]} individuals, {genotypes.shape[1]} SNPs")

    # Calculate LD decay
    print("\n[3/5] Calculating LD decay...")
    distances, ld_values = calculate_pairwise_ld(genotypes, max_lag=50)
    print(f"  Calculated {len(ld_values)} pairwise comparisons")

    # Calculate LD matrix
    print("\n[4/5] Calculating LD matrix...")
    ld_matrix = calculate_ld_matrix(genotypes, window_size=20)

    # Create plots
    print("\n[5/5] Creating visualizations...")
    create_ld_plots(distances, ld_values, ld_matrix, output_dir)

    # Summary
    print("\n" + "=" * 40)
    print("Summary")
    print("=" * 40)
    print(f"Individuals: {genotypes.shape[0]}")
    print(f"SNPs: {genotypes.shape[1]}")
    print(f"Distance range: {distances.min()}-{distances.max()} kb")
    print(f"Initial r²: {ld_values[0]:.3f}")
    print(f"r² at 200kb: {ld_values[min(20, len(ld_values) - 1)]:.3f}")
    print(f"\nOutputs:")
    print(f"  - {output_dir}/ld_analysis.png")
    print("\n✅ LD analysis example complete!")
    print("\nIn QTLmax: Linkage Disequilibrium → Decay curve")


if __name__ == "__main__":
    main()
