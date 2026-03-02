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
Example: Qmapper - Physical Mapping with SNPs (Ideogram)

This example demonstrates chromosome ideogram visualization with SNP positioning -
similar to QTLmax's Qmapper tool in the "Genome browse" tab.

WHAT THIS MEANS:
Physical maps (ideograms) are graphical representations of chromosomes showing:
- Chromosome length and centromere position
- SNP locations as tick marks or bands
- QTL regions highlighted
- Linkage group alignment with Manhattan plots

WHY WE DO THIS:
- Visualize spatial distribution of SNPs on chromosomes
- Identify physical positions of QTL
- Compare genetic and physical maps
- Publication-quality chromosome figures

WHAT'S DEMONSTRATED:
1. Create chromosome ideogram templates
2. Position SNPs on physical map
3. Add QTL regions with confidence intervals
4. Create publication-quality figure with multiple chromosomes

Equivalent to QTLmax: "Qmapper" tool
https://open.qtlmax.com/guide/index.php/2026/02/09/physical-mapping-of-snps/
"""

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import warnings

warnings.filterwarnings("ignore")

# Create output directory
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)


def generate_chromosome_data():
    """Generate synthetic chromosome and SNP data."""
    np.random.seed(42)

    # Define chromosomes with lengths (in Mb - like rice genome)
    chromosomes = {
        "Chr1": 43.3,
        "Chr2": 35.9,
        "Chr3": 36.6,
        "Chr4": 35.2,
        "Chr5": 30.4,
    }

    # Generate SNPs on each chromosome
    snp_data = []
    snp_id = 0

    for chrom, length in chromosomes.items():
        # Distribute SNPs randomly but with some clustering
        n_snps = np.random.randint(15, 30)

        # Create clusters of SNPs
        n_clusters = np.random.randint(2, 5)
        cluster_centers = np.sort(np.random.uniform(0, length, n_clusters))

        for cluster in cluster_centers:
            # SNPs around cluster center
            cluster_snps = max(3, n_snps // n_clusters)
            positions = np.random.normal(cluster, length * 0.05, cluster_snps)
            positions = np.clip(positions, 0.1, length - 0.1)

            for pos in positions:
                # Determine if SNP is significant (QTL region)
                is_qtl = abs(pos - cluster) < length * 0.02

                snp_data.append(
                    {
                        "SNP": f"rs{snp_id}",
                        "Chromosome": chrom,
                        "Position_Mb": round(pos, 2),
                        "P_value": 10 ** (-np.random.uniform(2, 10))
                        if is_qtl
                        else np.random.uniform(0.01, 1.0),
                        "Effect": np.random.uniform(-0.5, 0.5) if is_qtl else 0,
                        "In_QTL": is_qtl,
                    }
                )
                snp_id += 1

    return pd.DataFrame(snp_data), chromosomes


def draw_ideogram(ax, chrom_name, length, snps, show_labels=True):
    """Draw a single chromosome ideogram with SNP positions."""
    # Chromosome body (rounded rectangle)
    chrom_height = 0.4
    chrom_width = 0.8

    # Draw chromosome
    ax.add_patch(
        plt.Rectangle(
            (0.1, 0.3),
            chrom_width,
            chrom_height,
            facecolor="#E8E8E8",
            edgecolor="#333333",
            linewidth=2,
            zorder=1,
        )
    )

    # Draw centromere (constriction)
    centromere_pos = length * 0.45  # Approximate centromere position
    centromere_width = 0.15
    ax.add_patch(
        plt.Rectangle(
            (0.1 + (centromere_pos / length) * chrom_width - centromere_width / 2, 0.3),
            centromere_width,
            chrom_height,
            facecolor="#888888",
            edgecolor="#333333",
            linewidth=1,
            zorder=2,
        )
    )

    # Get SNPs for this chromosome
    chrom_snps = snps[snps["Chromosome"] == chrom_name].copy()

    if len(chrom_snps) > 0:
        # Scale positions to chromosome width
        x_positions = 0.1 + (chrom_snps["Position_Mb"].values / length) * chrom_width

        # Color by significance
        colors = []
        for _, row in chrom_snps.iterrows():
            if row["P_value"] < 5e-8:
                colors.append("#FF0000")  # Red for genome-wide significant
            elif row["P_value"] < 1e-5:
                colors.append("#FF6600")  # Orange for suggestive
            elif row["In_QTL"]:
                colors.append("#0066FF")  # Blue for QTL region
            else:
                colors.append("#666666")  # Gray for other

        # Draw SNP markers
        y_base = 0.5
        for x, color in zip(x_positions, colors):
            ax.plot(
                [x, x],
                [y_base - 0.15, y_base + 0.15],
                color=color,
                linewidth=2,
                zorder=3,
            )

        # Add QTL region highlights
        qtl_snps = chrom_snps[chrom_snps["In_QTL"]]
        if len(qtl_snps) > 0:
            qtl_start = qtl_snps["Position_Mb"].min()
            qtl_end = qtl_snps["Position_Mb"].max()
            qtl_x_start = 0.1 + (qtl_start / length) * chrom_width
            qtl_x_end = 0.1 + (qtl_end / length) * chrom_width

            # Semi-transparent QTL region
            ax.add_patch(
                plt.Rectangle(
                    (qtl_x_start, 0.25),
                    qtl_x_end - qtl_x_start,
                    0.5,
                    facecolor="#FF0000",
                    alpha=0.15,
                    edgecolor="#FF0000",
                    linewidth=1,
                    linestyle="--",
                    zorder=0,
                )
            )

    # Add chromosome label
    ax.text(
        0.5,
        0.05,
        f"{chrom_name} ({length:.1f} Mb)",
        ha="center",
        va="bottom",
        fontsize=10,
        fontweight="bold",
    )

    # Set limits and remove axes
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.axis("off")


def run_qmapper():
    """Main Qmapper workflow."""
    print("=" * 60)
    print("QAPPER: PHYSICAL MAPPING WITH SNPs")
    print("=" * 60)

    # Generate data
    print("\n[1/4] Generating chromosome and SNP data...")
    snps_df, chromosomes = generate_chromosome_data()
    print(f"  Generated: {len(snps_df)} SNPs on {len(chromosomes)} chromosomes")

    # Save SNP data
    snps_df.to_csv(f"{output_dir}/snp_positions.csv", index=False)
    print(f"  Saved: snp_positions.csv")

    # Create chromosome map
    print("\n[2/4] Creating chromosome ideogram...")
    n_chr = len(chromosomes)
    fig, axes = plt.subplots(n_chr, 1, figsize=(14, n_chr * 1.5))

    if n_chr == 1:
        axes = [axes]

    for idx, (chrom, length) in enumerate(chromosomes.items()):
        draw_ideogram(axes[idx], chrom, length, snps_df)
        axes[idx].set_title(f"", fontsize=1)  # Title reserved for labels

    plt.suptitle(
        "Qmapper: Physical Map of SNP Markers", fontsize=16, fontweight="bold", y=0.98
    )
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f"{output_dir}/ideogram.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: ideogram.png")

    # Create combined figure with Manhattan-style overlay
    print("\n[3/4] Creating combined visualization...")
    fig, axes = plt.subplots(
        2, 1, figsize=(14, 8), gridspec_kw={"height_ratios": [1, 3]}
    )

    # Top: Ideogram
    ax_ideo = axes[0]
    chrom_list = list(chromosomes.items())
    for idx, (chrom, length) in enumerate(chrom_list):
        # Normalize to 0-1 range
        x_offset = idx * 1.2

        # Draw chromosome
        ax_ideo.add_patch(
            plt.Rectangle(
                (x_offset, 0.4),
                0.8,
                0.2,
                facecolor="#E8E8E8",
                edgecolor="#333333",
                linewidth=2,
            )
        )

        # Add SNPs
        chrom_snps = snps_df[snps_df["Chromosome"] == chrom]
        if len(chrom_snps) > 0:
            x_positions = (
                x_offset + 0.1 + (chrom_snps["Position_Mb"].values / length) * 0.6
            )
            colors = [
                "#FF0000" if p < 5e-8 else "#666666" for p in chrom_snps["P_value"]
            ]
            for x, c in zip(x_positions, colors):
                ax_ideo.plot([x, x], [0.55, 0.65], color=c, linewidth=2)

        ax_ideo.text(x_offset + 0.5, 0.25, chrom, ha="center", fontsize=8)

    ax_ideo.set_xlim(-0.2, len(chromosomes) * 1.2 - 0.2)
    ax_ideo.set_ylim(0, 1)
    ax_ideo.set_title(
        "Chromosome Ideogram with SNP Positions", fontsize=12, fontweight="bold"
    )
    ax_ideo.set_yticks([])
    ax_ideo.set_xticks([])
    ax_ideo.spines["top"].set_visible(False)
    ax_ideo.spines["right"].set_visible(False)
    ax_ideo.spines["bottom"].set_visible(False)
    ax_ideo.spines["left"].set_visible(False)

    # Bottom: Manhattan plot aligned with ideogram
    ax_man = axes[1]

    # Prepare Manhattan data
    all_chroms = []
    all_pos = []
    all_pvals = []
    chrom_offsets = {}

    offset = 0
    for chrom, length in chrom_list:
        chrom_data = snps_df[snps_df["Chromosome"] == chrom].copy()
        if len(chrom_data) > 0:
            chrom_data["PlotPos"] = chrom_data["Position_Mb"] + offset
            all_chroms.append(chrom_data["PlotPos"].values)
            all_pos.extend(chrom_data["Position_Mb"].values)
            all_pvals.extend(-np.log10(chrom_data["P_value"].values))
        chrom_offsets[chrom] = offset
        offset += length + 5  # Add gap between chromosomes

    # Flatten
    plot_positions = []
    p_values = []
    for c in all_chroms:
        plot_positions.extend(c)
    plot_positions = np.array(plot_positions)
    p_values = np.array(all_pvals)

    # Color by chromosome
    colors = ["#1f77b4" if i % 2 == 0 else "#aec7e8" for i in range(len(chrom_list))]

    for i, (chrom, length) in enumerate(chrom_list):
        mask = (plot_positions >= chrom_offsets[chrom]) & (
            plot_positions < chrom_offsets[chrom] + length + 5
        )
        ax_man.scatter(
            plot_positions[mask],
            p_values[mask],
            c=colors[i % len(colors)],
            s=20,
            alpha=0.7,
        )

    # Add significance line
    ax_man.axhline(
        y=-np.log10(5e-8),
        color="red",
        linestyle="--",
        linewidth=1.5,
        label="Genome-wide (5e-8)",
    )
    ax_man.axhline(
        y=-np.log10(1e-5),
        color="blue",
        linestyle="--",
        linewidth=1,
        label="Suggestive (1e-5)",
    )

    ax_man.set_xlabel("Chromosome", fontsize=12)
    ax_man.set_ylabel("-log10(p)", fontsize=12)
    ax_man.set_title(
        "Manhattan Plot (Aligned with Physical Map)", fontsize=12, fontweight="bold"
    )
    ax_man.legend(loc="upper right")

    # Add chromosome labels
    for chrom, length in chrom_list:
        mid = chrom_offsets[chrom] + length / 2
        ax_man.text(mid, 0, chrom, ha="center", fontsize=9)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/qmapper_combined.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: qmapper_combined.png")

    # Summary
    print("\n[4/4] Summary:")
    print(f"  - {len(snps_df)} SNPs mapped to {len(chromosomes)} chromosomes")
    print(f"  - {sum(snps_df['P_value'] < 5e-8)} genome-wide significant SNPs")
    print(f"  - {sum(snps_df['In_QTL'])} SNPs in QTL regions")

    print("\n" + "=" * 60)
    print("QAPPER ANALYSIS COMPLETE")
    print("=" * 60)

    return {"n_snps": len(snps_df), "chromosomes": len(chromosomes)}


if __name__ == "__main__":
    results = run_qmapper()
