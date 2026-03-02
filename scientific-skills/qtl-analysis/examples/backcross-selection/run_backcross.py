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
Example: Backcross Selection for Breeding

This example demonstrates backcross breeding with marker-assisted selection -
similar to QTLmax's "Backcross" tab.

WHAT THIS MEANS:
Backcrossing is a breeding technique where offspring are repeatedly crossed
back to one parent (the "recurrent" parent) to recover that parent's genome
while retaining a specific trait from the other parent (the "donor").

WHY WE DO THIS:
- Introduce specific traits (e.g., disease resistance) from donor to elite variety
- Each backcross adds ~50% of recurrent parent genome
- After 7-8 generations, have ~99% recurrent parent but keep target trait
- Marker-assisted selection tracks target gene

WHAT'S DEMONSTRATED:
1. Generate pedigree with F1 and multiple BC generations
2. Calculate similarity to recurrent parent
3. Select best offspring for next generation
4. Visualize progress across generations

Equivalent to QTLmax: "Backcross" tab
https://open.qtlmax.com/guide/index.php/2025/07/11/backcross-selection/
"""

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
import os
import warnings

warnings.filterwarnings("ignore")

# Create output directory
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)


def generate_backcross_data(n_snps=100, n_per_gen=30):
    """Generate synthetic backcross breeding data."""
    np.random.seed(42)

    # Define SNPs where parents differ
    recurrent_parent = np.random.choice([0, 2], n_snps)  # Reference alleles
    donor_parent = np.where(recurrent_parent == 0, 2, 0)  # Alternative alleles

    # Mark target trait region (SNPs 40-50)
    target_region = slice(40, 50)

    data = []

    # Generation 0: Parents
    for i in range(100):
        data.append(
            {
                "Sample": f"Recurrent_{i}",
                "Generation": 0,
                "Type": "Recurrent Parent",
                "Similarity": 100.0,
                **{f"SNP_{j}": recurrent_parent[j] for j in range(n_snps)},
            }
        )

    for i in range(100):
        data.append(
            {
                "Sample": f"Donor_{i}",
                "Generation": 0,
                "Type": "Donor Parent",
                "Similarity": 0.0,
                **{f"SNP_{j}": donor_parent[j] for j in range(n_snps)},
            }
        )

    # Generation 1: F1 (50% from each parent)
    for i in range(n_per_gen):
        genotype = np.where(
            np.random.random(n_snps) < 0.5, recurrent_parent, donor_parent
        )
        similarity = 100 * np.sum(genotype == recurrent_parent) / n_snps

        record = {
            "Sample": f"F1_{i:03d}",
            "Generation": 1,
            "Type": "F1",
            "Similarity": similarity,
        }
        record.update({f"SNP_{j}": int(genotype[j]) for j in range(n_snps)})
        data.append(record)

    # Generations 2-7: BC1-BC6
    for gen in range(2, 8):  # BC1 to BC6
        for i in range(n_per_gen):
            # Random parent from previous generation
            parent_idx = np.random.randint(0, n_per_gen)
            parent_genotype = np.array(
                [data[-n_per_gen + parent_idx][f"SNP_{j}"] for j in range(n_snps)]
            )

            # Cross with recurrent parent
            genotype = np.where(
                np.random.random(n_snps) < 0.5, recurrent_parent, parent_genotype
            )

            # Ensure target trait is retained (heterozygous in target region)
            genotype[target_region] = 1

            similarity = 100 * np.sum(genotype == recurrent_parent) / n_snps

            record = {
                "Sample": f"BC{gen - 1}_{i:03d}",
                "Generation": gen - 1,
                "Type": f"BC{gen - 1}",
                "Similarity": similarity,
            }
            record.update({f"SNP_{j}": int(genotype[j]) for j in range(n_snps)})
            data.append(record)

    df = pd.DataFrame(data)
    return df, recurrent_parent, donor_parent


def run_backcross_selection():
    """Main backcross selection workflow."""
    print("=" * 60)
    print("BACKCROSS SELECTION")
    print("=" * 60)

    # Generate data
    print("\n[1/5] Generating backcross breeding data...")
    df, recurrent, donor = generate_backcross_data(n_snps=100, n_per_gen=30)
    print(f"  Generated: {len(df)} individuals across 7 generations")

    # Save full data
    df.to_csv(f"{output_dir}/backcross_data.csv", index=False)
    print(f"  Saved: backcross_data.csv")

    # Calculate statistics by generation
    print("\n[2/5] Calculating generation statistics...")
    gen_stats = df.groupby("Generation")["Similarity"].agg(
        ["mean", "std", "min", "max"]
    )
    print(gen_stats)

    # Expected values: F1=50%, BC1=75%, BC2=87.5%, etc.
    expected = [100, 50, 75, 87.5, 93.75, 96.875, 98.4375, 99.21875]

    # Save statistics
    gen_stats.to_csv(f"{output_dir}/generation_statistics.csv")
    print(f"  Saved: generation_statistics.csv")

    # Select best candidates from each generation
    print("\n[3/5] Selecting best candidates...")
    selected = []
    for gen in range(1, 8):
        gen_df = df[df["Generation"] == gen]
        # Top 20% by similarity
        threshold = gen_df["Similarity"].quantile(0.8)
        best = gen_df[gen_df["Similarity"] >= threshold]
        selected.append(best)
        print(
            f"  Generation {gen}: Selected {len(best)}/{len(gen_df)} candidates (threshold: {threshold:.1f}%)"
        )

    selected_df = pd.concat(selected)
    selected_df.to_csv(f"{output_dir}/selected_candidates.csv", index=False)
    print(f"  Saved: selected_candidates.csv")

    # Create visualization
    print("\n[4/5] Creating visualizations...")
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Plot 1: Pedigree diagram
    ax1 = axes[0, 0]
    ax1.set_xlim(0, 10)
    ax1.set_ylim(0, 10)

    # Draw pedigree flow
    generations = ["P", "F1", "BC1", "BC2", "BC3", "BC4", "BC5", "BC6"]
    y_pos = [8, 7, 6, 5, 4, 3, 2, 1]
    x_pos = [2, 5, 5, 5, 5, 5, 5, 5]

    # Parents
    ax1.add_patch(
        mpatches.FancyBboxPatch(
            (1, 8.5),
            1.5,
            0.8,
            boxstyle="round,pad=0.1",
            facecolor="lightgreen",
            edgecolor="green",
            linewidth=2,
        )
    )
    ax1.text(
        2.75,
        8.9,
        "Recurrent Parent",
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
    )

    ax1.add_patch(
        mpatches.FancyBboxPatch(
            (6, 8.5),
            1.5,
            0.8,
            boxstyle="round,pad=0.1",
            facecolor="lightcoral",
            edgecolor="red",
            linewidth=2,
        )
    )
    ax1.text(
        6.75,
        8.9,
        "Donor Parent",
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
    )

    # Arrows
    ax1.arrow(
        2.75, 8.5, 1.5, -0.5, head_width=0.15, head_length=0.1, fc="black", ec="black"
    )
    ax1.arrow(
        6.75, 8.5, -2.5, -0.5, head_width=0.15, head_length=0.1, fc="black", ec="black"
    )

    # Generations
    for i, gen in enumerate(generations[1:]):
        ax1.add_patch(
            mpatches.FancyBboxPatch(
                (4.5, 8.5 - (i + 1) * 0.9),
                1.5,
                0.6,
                boxstyle="round,pad=0.1",
                facecolor="lightblue",
                edgecolor="blue",
                linewidth=1.5,
            )
        )
        ax1.text(
            5.25,
            8.2 - (i + 1) * 0.9,
            gen,
            ha="center",
            va="center",
            fontsize=9,
            fontweight="bold",
        )
        if i > 0:
            ax1.arrow(
                5.25,
                8.5 - i * 0.9,
                0,
                -0.2,
                head_width=0.1,
                head_length=0.08,
                fc="black",
                ec="black",
            )

    # Cross arrows
    for i in range(1, 7):
        y = 8.2 - (i + 1) * 0.9
        ax1.arrow(
            2.75,
            8.5 - i * 0.9,
            2.3,
            y - (8.5 - i * 0.9) + 0.2,
            head_width=0.1,
            head_length=0.08,
            fc="gray",
            ec="gray",
            linestyle="--",
            alpha=0.7,
        )

    ax1.set_title("Backcross Pedigree", fontsize=14, fontweight="bold")
    ax1.axis("off")

    # Plot 2: Similarity heatmap (subset of SNPs)
    ax2 = axes[0, 1]

    # Sample individuals from each generation
    sample_ids = []
    for gen in range(1, 8):
        gen_df = df[df["Generation"] == gen]
        sample_ids.extend(
            gen_df["Sample"].iloc[::3].head(5).tolist()
        )  # Sample every 3rd

    # Get genotypes for sampled individuals
    snp_cols = [f"SNP_{i}" for i in range(30)]  # Show first 30 SNPs
    heatmap_data = df[df["Sample"].isin(sample_ids)][snp_cols].values

    # Color by similarity to recurrent
    heatmap_colored = np.zeros((*heatmap_data.shape, 3))
    for i in range(heatmap_data.shape[0]):
        for j in range(heatmap_data.shape[1]):
            val = heatmap_data[i, j]
            if val == recurrent[j]:  # Matches recurrent
                heatmap_colored[i, j] = [0.7, 1.0, 0.7]  # Green
            elif val == donor[j]:  # Matches donor
                heatmap_colored[i, j] = [1.0, 0.7, 0.7]  # Red
            else:  # Heterozygous
                heatmap_colored[i, j] = [1.0, 1.0, 0.7]  # Yellow

    ax2.imshow(heatmap_colored, aspect="auto", interpolation="nearest")
    ax2.set_title(
        "Genotype Heatmap\n(Green=Recurrent, Red=Donor, Yellow=Het)",
        fontsize=12,
        fontweight="bold",
    )
    ax2.set_xlabel("SNP")
    ax2.set_ylabel("Individual")
    ax2.set_yticks([])

    # Add generation labels
    gen_positions = []
    for gen in range(1, 8):
        gen_df = df[df["Generation"] == gen]
        gen_samples = gen_df[gen_df["Sample"].isin(sample_ids)]
        if len(gen_samples) > 0:
            idx = df[df["Sample"].isin(sample_ids)].index.get_indexer(
                gen_samples.index
            )[0]
            gen_positions.append((idx, f"BC{gen - 1}" if gen > 1 else "F1"))

    for pos, label in gen_positions:
        ax2.text(-2, pos, label, ha="right", va="center", fontsize=8)

    # Plot 3: Similarity progress
    ax3 = axes[1, 0]

    generations_list = df["Generation"].unique()
    gen_means = [
        df[df["Generation"] == g]["Similarity"].mean() for g in generations_list
    ]
    gen_stds = [df[df["Generation"] == g]["Similarity"].std() for g in generations_list]

    ax3.plot(
        generations_list,
        expected[: len(generations_list)],
        "r--",
        linewidth=2,
        label="Expected",
        marker="o",
        markersize=6,
    )
    ax3.errorbar(
        generations_list,
        gen_means,
        yerr=gen_stds,
        fmt="b-o",
        capsize=5,
        capthick=2,
        linewidth=2,
        markersize=6,
        label="Observed (mean ± SD)",
    )

    # Add selection threshold
    ax3.axhline(
        y=93,
        color="orange",
        linestyle="--",
        linewidth=1.5,
        label="Selection threshold (93%)",
    )
    ax3.axhline(
        y=99, color="green", linestyle="--", linewidth=1.5, label="Target (99%)"
    )

    ax3.set_xlabel("Generation")
    ax3.set_ylabel("% Similarity to Recurrent Parent")
    ax3.set_title("Backcross Progress", fontsize=12, fontweight="bold")
    ax3.legend(loc="lower right")
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim(0, 105)

    # Add generation labels
    ax3.set_xticks(generations_list)
    ax3.set_xticklabels(
        ["P", "F1", "BC1", "BC2", "BC3", "BC4", "BC5", "BC6"][: len(generations_list)]
    )

    # Plot 4: Selection candidates
    ax4 = axes[1, 1]

    # Distribution of similarity by generation
    for gen in [1, 3, 5, 7]:  # F1, BC2, BC4, BC6
        if gen in generations_list:
            gen_data = df[df["Generation"] == gen]["Similarity"]
            ax4.hist(
                gen_data,
                bins=10,
                alpha=0.6,
                label=f"{'F1' if gen == 1 else f'BC{gen - 1}'} (n={len(gen_data)})",
                edgecolor="black",
                linewidth=0.5,
            )

    ax4.axvline(
        x=93, color="orange", linestyle="--", linewidth=2, label="Selection threshold"
    )
    ax4.set_xlabel("% Similarity to Recurrent Parent")
    ax4.set_ylabel("Number of Individuals")
    ax4.set_title("Selection Distribution", fontsize=12, fontweight="bold")
    ax4.legend(loc="upper left")
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/backcross_selection.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: backcross_selection.png")

    # Summary
    print("\n[5/5] Summary:")
    print(f"  - Total individuals: {len(df)}")
    print(f"  - SNPs analyzed: 100")
    print(f"  - Generations: Parents + F1 + BC1-BC6 (7 total)")
    print(f"  - Selected candidates: {len(selected_df)} (top 20% per generation)")
    print(f"  - Target similarity: 93-99%")

    print("\n  Generation Progress:")
    for gen in range(1, 8):
        gen_df = df[df["Generation"] == gen]
        mean_sim = gen_df["Similarity"].mean()
        expected_sim = expected[gen]
        print(
            f"    {'F1' if gen == 1 else f'BC{gen - 1}'}: {mean_sim:.1f}% (expected: {expected_sim:.1f}%)"
        )

    print("\n" + "=" * 60)
    print("BACKCROSS SELECTION COMPLETE")
    print("=" * 60)

    return {"generations": 7, "selected": len(selected_df)}


if __name__ == "__main__":
    results = run_backcross_selection()
