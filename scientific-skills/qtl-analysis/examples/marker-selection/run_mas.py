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
Example: Marker-Assisted Selection (MAS)

This example demonstrates MAS - using significant markers to select individuals
for breeding based on their genotype at target loci.

Equivalent to QTLmax: "Marker assisted selection (MAS)"
https://open.qtlmax.com/guide/index.php/2025/07/10/marker-assisted-selection-mas/

Auto-installs: pandas, numpy, matplotlib
"""

import subprocess
import sys
import os


def install_packages():
    """Install required packages without root"""
    packages = ["pandas", "numpy", "matplotlib"]
    for pkg in packages:
        try:
            __import__(pkg.replace("-", "_"))
        except ImportError:
            print(f"Installing {pkg}...")
            subprocess.check_call(
                [sys.executable, "-m", "pip", "install", "--user", "-q", pkg]
            )


def generate_mas_data(n_candidates=100, n_snps=50):
    """Generate candidate genotypes for MAS"""
    import numpy as np

    np.random.seed(42)

    # Define target QTL markers with known effects
    target_snps = [5, 15, 25, 35]  # Indices of causal SNPs
    target_effects = [0.5, -0.3, 0.8, 0.4]  # Effect sizes
    target_alleles = [1, 1, 0, 1]  # Desired allele (1 = alternate)

    # Generate genotypes
    genotypes = np.random.binomial(2, 0.3, (n_candidates, n_snps))

    # Generate phenotypes from target markers + noise
    breeding_values = np.zeros(n_candidates)
    for snp_idx, effect, desired_allele in zip(
        target_snps, target_effects, target_alleles
    ):
        # Calculate allele dosage (0, 1, 2)
        allele_match = (genotypes[:, snp_idx] == desired_allele).astype(float)
        breeding_values += allele_match * effect

    phenotype = breeding_values + np.random.normal(0, 0.5, n_candidates)

    # Sample IDs
    sample_ids = [f"Calf{i + 1}" for i in range(n_candidates)]
    snp_ids = [f"rs{i}" for i in range(n_snps)]

    return (
        genotypes,
        sample_ids,
        snp_ids,
        target_snps,
        target_effects,
        target_alleles,
        phenotype,
    )


def run_mas(genotypes, sample_ids, snp_ids, target_snps, target_alleles, output_dir):
    """Run Marker-Assisted Selection"""
    import pandas as pd

    print("Running Marker-Assisted Selection...")

    # Create selection criteria
    # For each target SNP, select individuals with desired genotype
    results = []

    for i, sample_id in enumerate(sample_ids):
        row = {"sample": sample_id, "phenotype": genotypes[i]}

        # Calculate breeding value from markers
        marker_score = 0
        selected = True

        for snp_idx, desired_allele in zip(target_snps, target_alleles):
            geno = genotypes[i, snp_idx]
            row[f"{snp_ids[snp_idx]}_geno"] = geno

            # Score: +1 if has desired allele, 0 otherwise
            allele_match = 1 if geno == desired_allele else 0
            marker_score += allele_match

            # If any target SNP has wrong genotype, not selected
            if geno != desired_allele:
                selected = False

        row["marker_score"] = marker_score
        row["selected"] = selected
        results.append(row)

    results_df = pd.DataFrame(results)
    results_df.to_csv(f"{output_dir}/mas_results.csv", index=False)

    return results_df


def create_mas_plot(results_df, snp_ids, target_snps, output_dir):
    """Create MAS visualization"""
    import matplotlib.pyplot as plt
    import numpy as np

    print("Creating MAS visualizations...")

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Selection distribution
    ax1 = axes[0, 0]
    selected = results_df[results_df["selected"] == True]
    not_selected = results_df[results_df["selected"] == False]

    ax1.bar(
        ["Selected", "Not Selected"],
        [len(selected), len(not_selected)],
        color=["green", "red"],
        alpha=0.7,
    )
    ax1.set_ylabel("Number of Candidates")
    ax1.set_title("MAS Selection Results")
    for i, v in enumerate([len(selected), len(not_selected)]):
        ax1.text(i, v + 1, str(v), ha="center", fontweight="bold")

    # 2. Marker score distribution
    ax2 = axes[0, 1]
    scores = results_df["marker_score"]
    ax2.hist(
        scores,
        bins=range(0, 5),
        align="left",
        color="steelblue",
        alpha=0.7,
        edgecolor="black",
    )
    ax2.set_xlabel("Marker Score (target alleles present)")
    ax2.set_ylabel("Frequency")
    ax2.set_title("Distribution of Marker Scores")
    ax2.axvline(x=len(target_snps), color="red", linestyle="--", label="Maximum score")
    ax2.legend()

    # 3. Target SNP genotypes heatmap
    ax3 = axes[1, 0]
    target_data = results_df[[f"{snp_ids[snp]}_geno" for snp in target_snps[:4]]]
    im = ax3.imshow(target_data.values, cmap="RdYlGn", aspect="auto", vmin=0, vmax=2)
    ax3.set_yticks([])
    ax3.set_xticks(range(len(target_snps[:4])))
    ax3.set_xticklabels([snp_ids[s] for s in target_snps[:4]], rotation=45)
    ax3.set_xlabel("Target SNPs")
    ax3.set_ylabel("Candidates")
    ax3.set_title("Target SNP Genotypes\n(0=AA, 1=Aa, 2=aa)")
    plt.colorbar(im, ax=ax3)

    # 4. Marker score vs phenotype (if available)
    ax4 = axes[1, 1]
    pheno_values = results_df["phenotype"].apply(
        lambda x: x.mean() if isinstance(x, np.ndarray) else 0
    )
    marker_scores = results_df["marker_score"]

    colors = ["green" if s else "red" for s in results_df["selected"]]
    ax4.scatter(marker_scores, pheno_values, c=colors, alpha=0.6, s=50)
    ax4.set_xlabel("Marker Score")
    ax4.set_ylabel("Phenotype Value")
    ax4.set_title("Marker Score vs Phenotype\n(Green=Selected, Red=Not Selected)")
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/mas_selection.png", dpi=150, bbox_inches="tight")
    print(f"MAS plot saved: {output_dir}/mas_selection.png")


def main():
    print("=" * 60)
    print("Example: Marker-Assisted Selection (MAS)")
    print("=" * 60)

    # Install packages
    print("\n[1/4] Installing dependencies...")
    install_packages()

    # Setup
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    # Generate data
    print("\n[2/4] Generating MAS candidate data...")
    (
        genotypes,
        sample_ids,
        snp_ids,
        target_snps,
        target_effects,
        target_alleles,
        phenotype,
    ) = generate_mas_data()
    print(f"  Generated: {len(sample_ids)} candidates, {len(snp_ids)} SNPs")
    print(f"  Target markers: {[snp_ids[s] for s in target_snps]}")

    # Run MAS
    print("\n[3/4] Running Marker-Assisted Selection...")
    results_df = run_mas(
        genotypes, sample_ids, snp_ids, target_snps, target_alleles, output_dir
    )

    selected_count = results_df["selected"].sum()
    print(f"  Selected: {selected_count}/{len(sample_ids)}")

    # Create plots
    print("\n[4/4] Creating visualizations...")
    create_mas_plot(results_df, snp_ids, target_snps, output_dir)

    # Summary
    print("\n" + "=" * 40)
    print("Summary")
    print("=" * 40)
    print(f"Candidates: {len(sample_ids)}")
    print(f"Target SNPs: {len(target_snps)}")
    print(f"Selected: {selected_count} ({100 * selected_count / len(sample_ids):.1f}%)")
    print(f"\nOutputs:")
    print(f"  - {output_dir}/mas_results.csv")
    print(f"  - {output_dir}/mas_selection.png")
    print("\n✅ MAS example complete!")
    print("\nIn QTLmax: Marker-Assisted Selection → MAS")


if __name__ == "__main__":
    main()
