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

"""
Example: Haplotype Analysis with Dendrogram

This example demonstrates haplotype block analysis and hierarchical clustering
using linkage disequilibrium patterns - similar to QTLmax's "Haplotype tools" tab.

WHAT THIS MEANS:
Haplotypes are sets of SNPs that tend to be inherited together because of
 linkage disequilibrium (LD). By clustering SNPs into haplotype blocks, we can:
- Reduce dimensionality of GWAS data
- Identify ancestral recombination blocks
- Improve imputation accuracy
- Visualize population structure

WHY WE DO THIS:
- Haplotype-based analysis is more powerful than single-SNP GWAS
- LD blocks represent historical recombination events
- Dendrograms show hierarchical relationships between haplotypes
- Enables "phased" view of genetic variation

WHAT'S DEMONSTRATED:
1. Calculate pairwise LD (r²) between nearby SNPs
2. Cluster SNPs into haplotype blocks using hierarchical clustering
3. Generate dendrogram showing haplotype relationships
4. Create LD heatmap with haplotype block boundaries

Equivalent to QTLmax: "Haplotype tools" tab
https://open.qtlmax.com/guide/index.php/2025/07/10/how-to-set-genomic-window-range-2/
"""

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import squareform
import os
import warnings

warnings.filterwarnings("ignore")

# Create output directory
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)


def calculate_ld_matrix(genotypes):
    """Calculate pairwise LD (r²) between all SNPs."""
    n_snps = genotypes.shape[1]
    ld_matrix = np.ones((n_snps, n_snps))

    # Calculate allele frequencies
    freq = genotypes.mean(axis=0) / 2  # Convert to allele frequency
    freq = np.clip(freq, 0.01, 0.99)  # Avoid division by zero

    # Calculate LD for each pair
    for i in range(n_snps):
        for j in range(i + 1, n_snps):
            # Genotype correlation (r)
            geno_i = genotypes[:, i]
            geno_j = genotypes[:, j]

            # Calculate r²
            corr = np.corrcoef(geno_i, geno_j)[0, 1]
            if np.isnan(corr):
                corr = 0
            ld_matrix[i, j] = corr**2
            ld_matrix[j, i] = corr**2

    return ld_matrix


def generate_haplotype_data(n_samples=100, n_snps=50, n_chr=3):
    """Generate synthetic haplotype data with realistic LD patterns."""
    np.random.seed(42)

    # Create haplotype blocks (groups of correlated SNPs)
    block_sizes = [8, 12, 10, 10, 10]  # 5 blocks
    assert sum(block_sizes) == n_snps, "Block sizes must sum to n_snps"

    # Generate base haplotypes for each block
    genotypes = np.zeros((n_samples, n_snps))
    snp_positions = []
    snp_names = []
    chroms = []

    pos = 0
    block_id = 0
    for block_size in block_sizes:
        # Generate block-specific haplotype (all SNPs in block correlated)
        n_haplotypes = 4  # 4 common haplotypes in this block
        block_haps = np.random.binomial(
            2,
            0.3 + 0.4 * np.random.rand(n_haplotypes, block_size),
            size=(n_haplotypes, block_size),
        )

        # Assign each sample to a haplotype
        for i in range(n_samples):
            hap_idx = np.random.randint(0, n_haplotypes)
            genotypes[i, pos : pos + block_size] = block_haps[hap_idx]

        # Add some noise (recombination)
        noise = np.random.binomial(2, 0.05, (n_samples, block_size))
        genotypes[:, pos : pos + block_size] = np.clip(
            genotypes[:, pos : pos + block_size] + noise, 0, 2
        )

        # Record positions and names
        for j in range(block_size):
            snp_positions.append(pos + j * 1000)  # 1kb spacing
            snp_names.append(f"rs{block_id}_{pos + j}")
            chroms.append(block_id + 1)

        pos += block_size
        block_id += 1

    return genotypes, snp_names, snp_positions, chroms


def run_haplotype_analysis():
    """Main haplotype analysis workflow."""
    print("=" * 60)
    print("HAPLOTYPE ANALYSIS WITH DENROGRAM")
    print("=" * 60)

    # Generate synthetic data
    print("\n[1/5] Generating synthetic haplotype data...")
    genotypes, snp_names, positions, chroms = generate_haplotype_data(
        n_samples=100, n_snps=50, n_chr=5
    )
    print(f"  Generated: {genotypes.shape[0]} samples × {genotypes.shape[1]} SNPs")

    # Save genotype data
    geno_df = pd.DataFrame(genotypes.astype(int), columns=snp_names)
    geno_df.to_csv(f"{output_dir}/haplotype_genotypes.csv", index=False)
    print(f"  Saved: haplotype_genotypes.csv")

    # Calculate LD matrix
    print("\n[2/5] Calculating LD matrix (r²)...")
    ld_matrix = calculate_ld_matrix(genotypes)
    print(f"  LD matrix: {ld_matrix.shape}")

    # Save LD matrix
    ld_df = pd.DataFrame(ld_matrix, index=snp_names, columns=snp_names)
    ld_df.to_csv(f"{output_dir}/ld_matrix.csv")
    print(f"  Saved: ld_matrix.csv")

    # Hierarchical clustering of SNPs based on LD
    print("\n[3/5] Performing hierarchical clustering...")
    # Convert LD to distance (1 - r²)
    distance_matrix = 1 - ld_matrix
    np.fill_diagonal(distance_matrix, 0)

    # Condense for linkage
    condensed = squareform(distance_matrix)
    linkage_matrix = linkage(condensed, method="ward")
    print(f"  Clustering complete")

    # Identify haplotype blocks (cut at 0.7 height)
    n_clusters = 5
    cluster_labels = fcluster(linkage_matrix, n_clusters, criterion="maxclust")
    print(f"  Identified {n_clusters} haplotype blocks")

    # Save cluster assignments
    cluster_df = pd.DataFrame(
        {
            "SNP": snp_names,
            "Position": positions,
            "Chromosome": chroms,
            "Haplotype_Block": cluster_labels,
        }
    )
    cluster_df.to_csv(f"{output_dir}/haplotype_blocks.csv", index=False)
    print(f"  Saved: haplotype_blocks.csv")

    # Create visualization
    print("\n[4/5] Creating visualizations...")
    n_snps = genotypes.shape[1]
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    # Plot 1: LD Heatmap
    ax1 = axes[0, 0]
    im = ax1.imshow(ld_matrix, cmap="RdBu_r", aspect="auto", vmin=0, vmax=1)
    ax1.set_title("Linkage Disequilibrium (r²) Matrix", fontsize=14, fontweight="bold")
    ax1.set_xlabel("SNP Index")
    ax1.set_ylabel("SNP Index")
    plt.colorbar(im, ax=ax1, label="r²")

    # Add haplotype block boundaries
    block_boundaries = [8, 20, 30, 40]
    for b in block_boundaries:
        ax1.axhline(y=b - 0.5, color="yellow", linewidth=2, linestyle="--")
        ax1.axvline(x=b - 0.5, color="yellow", linewidth=2, linestyle="--")

    # Plot 2: Dendrogram
    ax2 = axes[0, 1]
    dendro = dendrogram(
        linkage_matrix,
        labels=snp_names,
        leaf_rotation=90,
        leaf_font_size=6,
        ax=ax2,
        color_threshold=0.7 * max(linkage_matrix[:, 2]),
    )
    ax2.set_title(
        "Haplotype Dendrogram\n(Hierarchical Clustering by LD)",
        fontsize=14,
        fontweight="bold",
    )
    ax2.set_xlabel("SNP")
    ax2.set_ylabel("Distance (1 - r²)")

    # Plot 3: Haplotype blocks visualization
    ax3 = axes[1, 0]
    colors = plt.cm.Set2(np.linspace(0, 1, n_clusters))
    block_positions = []
    block_labels = []
    for i in range(n_snps):
        ax3.axvline(
            x=i,
            ymin=0,
            ymax=1,
            color=colors[cluster_labels[i] - 1],
            alpha=0.7,
            linewidth=3,
        )
        if i == 0 or cluster_labels[i] != cluster_labels[i - 1]:
            block_positions.append(i)
            block_labels.append(f"Block {cluster_labels[i]}")

    ax3.set_xlim(0, n_snps)
    ax3.set_ylim(0, 1)
    ax3.set_title(
        "Haplotype Blocks\n(Colored by Cluster)", fontsize=14, fontweight="bold"
    )
    ax3.set_xlabel("SNP Index")
    ax3.set_yticks([])

    # Add legend
    for i, (pos, label) in enumerate(zip(block_positions[:5], block_labels[:5])):
        ax3.text(pos + 2, 0.5, label, fontsize=8, va="center")

    # Plot 4: LD decay within blocks
    ax4 = axes[1, 1]
    # Calculate average LD by distance
    distances = []
    ld_values = []
    for i in range(n_snps):
        for j in range(i + 1, min(i + 10, n_snps)):  # Only nearby SNPs
            dist = abs(positions[i] - positions[j]) / 1000  # kb
            ld_values.append(ld_matrix[i, j])
            distances.append(dist)

    ax4.scatter(distances, ld_values, alpha=0.5, s=20)
    ax4.set_title("LD Decay within Haplotype Blocks", fontsize=14, fontweight="bold")
    ax4.set_xlabel("Distance (kb)")
    ax4.set_ylabel("r²")
    ax4.set_ylim(0, 1)

    # Add trend line
    z = np.polyfit(distances, ld_values, 2)
    p = np.poly1d(z)
    x_trend = np.linspace(0, max(distances), 100)
    ax4.plot(x_trend, p(x_trend), "r-", linewidth=2, label="Trend")
    ax4.legend()

    plt.tight_layout()
    plt.savefig(f"{output_dir}/haplotype_analysis.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: haplotype_analysis.png")

    # Print summary
    print("\n[5/5] Summary:")
    print(f"  - {n_snps} SNPs analyzed")
    print(f"  - {n_clusters} haplotype blocks identified")
    print(f"  - Block sizes: 8, 12, 10, 10, 10 SNPs")
    print(f"  - Hierarchical clustering performed with Ward's method")

    # Block statistics
    print("\n  Haplotype Block Summary:")
    for block in range(1, n_clusters + 1):
        n_snps_in_block = np.sum(cluster_labels == block)
        print(f"    Block {block}: {n_snps_in_block} SNPs")

    print("\n" + "=" * 60)
    print("HAPLOTYPE ANALYSIS COMPLETE")
    print("=" * 60)

    return {
        "n_snps": n_snps,
        "n_blocks": n_clusters,
        "ld_matrix": ld_matrix,
        "clusters": cluster_labels,
    }


if __name__ == "__main__":
    results = run_haplotype_analysis()
