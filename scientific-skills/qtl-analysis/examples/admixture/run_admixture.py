#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Author: Clayton Young <Clayton@SuperiorByteWorks.com>
# LinkedIn: https://linkedin.com/in/claytoneyoung/
# GitHub: https://github.com/borealBytes

#!/usr/bin/env python3
"""
Example: Admixture Analysis

This example demonstrates population admixture analysis using PCA and clustering.
Identifies ancestral populations and estimates individual ancestry proportions.

Equivalent to QTLmax: "How to calculate admixture and plotting admixture bar chart"
https://open.qtlmax.com/guide/index.php/2025/07/10/calculating-admixture-and-plotting-admixture-bar-chart/

Auto-installs: pandas, numpy, scikit-learn, matplotlib
"""

import subprocess
import sys
import os


def install_packages():
    """Install required packages without root"""
    packages = ["pandas", "numpy", "scikit-learn", "matplotlib"]
    for pkg in packages:
        try:
            __import__(pkg.replace("-", "_"))
        except ImportError:
            print(f"Installing {pkg}...")
            subprocess.check_call(
                [sys.executable, "-m", "pip", "install", "--user", "-q", pkg]
            )


def generate_admixed_population(n_individuals=150, n_snps=1000, n_pops=3):
    """Generate synthetic admixed population data"""
    import numpy as np
    import pandas as pd

    np.random.seed(42)

    # Generate K ancestral populations with different allele frequencies
    ancestral_freqs = []
    for k in range(n_pops):
        freqs = np.random.uniform(0.1, 0.9, n_snps)
        ancestral_freqs.append(freqs)

    # Generate individuals with mixed ancestry
    genotypes = np.zeros((n_individuals, n_snps))
    ancestry = np.zeros((n_individuals, n_pops))

    for i in range(n_individuals):
        # Random ancestry proportions (Dirichlet distribution)
        alpha = np.ones(n_pops)
        props = np.random.dirichlet(alpha)
        ancestry[i] = props

        # Generate genotype based on ancestry
        for k in range(n_pops):
            # Binomial based on population frequency
            freq = ancestral_freqs[k]
            geno_k = np.random.binomial(2, freq, n_snps)
            genotypes[i] += props[k] * geno_k

    genotypes = np.clip(genotypes, 0, 2).astype(int)

    # Sample IDs
    sample_ids = [f"Ind{i + 1}" for i in range(n_individuals)]
    snp_ids = [f"snp{i + 1}" for i in range(n_snps)]

    return genotypes, sample_ids, snp_ids, ancestry


def run_admixture_analysis(genotypes, sample_ids, snp_ids, output_dir, n_pops=3):
    """Run PCA and estimate ancestry"""
    import numpy as np
    import pandas as pd
    from sklearn.decomposition import PCA

    print("Running admixture analysis...")

    # Center and scale
    geno_mean = genotypes.mean(axis=0)
    geno_std = genotypes.std(axis=0) + 1e-8
    geno_scaled = (genotypes - geno_mean) / geno_std

    # PCA
    pca = PCA(n_components=n_pops)
    pc_scores = pca.fit_transform(geno_scaled)

    # Estimate ancestry using PC scores normalized to sum to 1
    # This is a simplified proxy for ancestry proportions
    ancestry = np.abs(pc_scores)
    ancestry = ancestry / (ancestry.sum(axis=1, keepdims=True) + 1e-8)

    # Save results
    results = pd.DataFrame(
        {
            "sample": sample_ids,
            "PC1": pc_scores[:, 0],
            "PC2": pc_scores[:, 1],
            "PC3": pc_scores[:, 2] if n_pops > 2 else 0,
            "Pop1": ancestry[:, 0],
            "Pop2": ancestry[:, 1],
            "Pop3": ancestry[:, 2],
        }
    )
    results.to_csv(f"{output_dir}/admixture_results.csv", index=False)

    return results, pca


def create_admixture_plot(results, pca, output_dir):
    """Create admixture bar chart and PCA plot"""
    import numpy as np
    import matplotlib.pyplot as plt

    print("Creating admixture visualizations...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1. PCA plot
    ax1 = axes[0, 0]
    scatter = ax1.scatter(
        results["PC1"],
        results["PC2"],
        c=results[["Pop1", "Pop2", "Pop3"]].values.argmax(axis=1),
        cmap="Set1",
        s=50,
        alpha=0.7,
    )
    ax1.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
    ax1.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)")
    ax1.set_title("PCA - Population Structure")
    ax1.grid(True, alpha=0.3)

    # 2. Stacked bar chart (admixture)
    ax2 = axes[0, 1]
    n_ind = len(results)
    x = np.arange(n_ind)
    width = 0.8

    colors = ["#e41a1c", "#377eb8", "#4daf4a"]
    bottom = np.zeros(n_ind)

    for k in range(3):
        ax2.bar(
            x,
            results[f"Pop{k + 1}"],
            width,
            bottom=bottom,
            color=colors[k],
            label=f"Population {k + 1}",
        )
        bottom += results[f"Pop{k + 1}"].values

    ax2.set_xlabel("Individual")
    ax2.set_ylabel("Ancestry Proportion")
    ax2.set_title("Admixture Plot (K=3)")
    ax2.legend(loc="upper right")

    # 3. Sorted admixture
    ax3 = axes[1, 0]
    sorted_idx = results["Pop1"].argsort()[::-1]
    bottom = np.zeros(n_ind)

    for k in range(3):
        ax3.bar(
            x,
            results.iloc[sorted_idx][f"Pop{k + 1}"],
            width,
            bottom=bottom,
            color=colors[k],
        )
        bottom += results.iloc[sorted_idx][f"Pop{k + 1}"].values

    ax3.set_xlabel("Individual (sorted by Pop1)")
    ax3.set_ylabel("Ancestry Proportion")
    ax3.set_title("Admixture Plot (Sorted)")

    # 4. Variance explained
    ax4 = axes[1, 1]
    var_explained = pca.explained_variance_ratio_[:5] * 100
    ax4.bar(range(1, len(var_explained) + 1), var_explained, color="steelblue")
    ax4.set_xlabel("Principal Component")
    ax4.set_ylabel("Variance Explained (%)")
    ax4.set_title("PCA Variance Explained")
    ax4.set_xticks(range(1, len(var_explained) + 1))

    plt.tight_layout()
    plt.savefig(f"{output_dir}/admixture_plot.png", dpi=150, bbox_inches="tight")
    print(f"Admixture plot saved: {output_dir}/admixture_plot.png")


def main():
    print("=" * 60)
    print("Example: Admixture Analysis")
    print("=" * 60)

    # Install packages
    print("\n[1/4] Installing dependencies...")
    install_packages()

    # Setup
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    # Generate data
    print("\n[2/4] Generating admixed population...")
    genotypes, sample_ids, snp_ids, true_ancestry = generate_admixed_population(
        n_pops=3
    )
    print(f"  Generated: {len(sample_ids)} individuals, {len(snp_ids)} SNPs")

    # Run analysis
    print("\n[3/4] Running admixture analysis...")
    results, pca = run_admixture_analysis(
        genotypes, sample_ids, snp_ids, output_dir, n_pops=3
    )

    # Create plots
    print("\n[4/4] Creating visualizations...")
    create_admixture_plot(results, pca, output_dir)

    # Summary
    print("\n" + "=" * 40)
    print("Summary")
    print("=" * 40)
    print(f"Individuals: {len(sample_ids)}")
    print(f"SNPs: {len(snp_ids)}")
    print(f"Ancestral populations (K): 3")
    print(f"PC1 variance: {pca.explained_variance_ratio_[0] * 100:.1f}%")
    print(f"PC2 variance: {pca.explained_variance_ratio_[1] * 100:.1f}%")
    print(f"\nOutputs:")
    print(f"  - {output_dir}/admixture_results.csv")
    print(f"  - {output_dir}/admixture_plot.png")
    print("\n✅ Admixture example complete!")
    print("\nIn QTLmax: Admixture → Calculate → Bar chart")


if __name__ == "__main__":
    main()
