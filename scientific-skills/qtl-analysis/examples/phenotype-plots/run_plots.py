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
Example: Phenotype Visualization Plots

This example demonstrates comprehensive phenotype data visualization:
- Box plots (show distribution, outliers, quartiles)
- Density plots (show probability distribution)
- Violin plots (show distribution shape)
- Heatmaps (show patterns in multivariate data)
- Scatter plot matrices (show correlations)

WHAT THIS MEANS:
Phenotype data exploration is crucial before QTL analysis. These plots help:
1. Identify outliers and data quality issues
2. Understand trait distributions
3. Detect population structure effects
4. Find correlations between traits
5. Visualize group differences

WHY WE DO THIS:
- Box plots: Show median, quartiles, outliers
- Density plots: Reveal distribution shape (normal, skewed, bimodal)
- Violin plots: Combine box plot with density
- Heatmaps: Show patterns across samples and traits
- Scatter matrices: Find trait correlations

Equivalent to QTLmax features:
- "How to draw a box plot"
- "How to draw a density plot"
- "How to draw a violin plot"
- "How to draw a heatmap"
- "How to draw a scatter plot matrix"

Auto-installs: pandas, numpy, matplotlib, seaborn
"""

import subprocess
import sys
import os


def install_packages():
    """Install required packages without root"""
    packages = ["pandas", "numpy", "matplotlib", "seaborn"]
    for pkg in packages:
        try:
            __import__(pkg.replace("-", "_"))
        except ImportError:
            print(f"Installing {pkg}...")
            subprocess.check_call(
                [sys.executable, "-m", "pip", "install", "--user", "-q", pkg]
            )


def generate_phenotype_data(n_samples=150, n_traits=5):
    """Generate synthetic phenotype data"""
    import numpy as np
    import pandas as pd

    np.random.seed(42)

    # Define trait names
    trait_names = ["Height", "Weight", "Yield", "Protein", "Oil"]

    # Generate correlated trait data
    # Create correlation structure
    mean = [100, 50, 200, 15, 5]
    std = [10, 8, 30, 3, 1]

    # Correlation matrix
    corr = np.array(
        [
            [1.0, 0.3, 0.5, 0.2, -0.1],
            [0.3, 1.0, 0.4, 0.1, 0.0],
            [0.5, 0.4, 1.0, 0.3, -0.2],
            [0.2, 0.1, 0.3, 1.0, 0.4],
            [-0.1, 0.0, -0.2, 0.4, 1.0],
        ]
    )

    # Convert to covariance
    cov = np.outer(std, std) * corr

    # Generate data
    data = np.random.multivariate_normal(mean, cov, n_samples)

    # Add some outliers
    outlier_idx = np.random.choice(n_samples, 5, replace=False)
    data[outlier_idx, 0] += 30  # Tall outliers

    # Create DataFrame
    df = pd.DataFrame(data, columns=trait_names)

    # Add group labels (e.g., different populations or treatments)
    df["Population"] = np.random.choice(["PopA", "PopB", "PopC"], n_samples)
    df["Sample"] = [f"Sample{i + 1}" for i in range(n_samples)]

    return df


def create_phenotype_plots(df, output_dir):
    """Create comprehensive phenotype plots"""
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np

    print("\nCreating phenotype visualizations...")

    # Set style
    sns.set_style("whitegrid")

    # Create figure with subplots
    fig = plt.figure(figsize=(16, 12))

    # 1. Box plots by population
    ax1 = plt.subplot(3, 3, 1)
    trait_cols = ["Height", "Weight", "Yield", "Protein", "Oil"]
    df_melted = df.melt(id_vars=["Population"], value_vars=trait_cols)
    sns.boxplot(data=df_melted, x="variable", y="value", hue="Population", ax=ax1)
    ax1.set_title("Box Plots: Traits by Population")
    ax1.set_xlabel("Trait")
    ax1.set_ylabel("Value")
    ax1.tick_params(axis="x", rotation=45)

    # 2. Density plot (Height)
    ax2 = plt.subplot(3, 3, 2)
    for pop in ["PopA", "PopB", "PopC"]:
        subset = df[df["Population"] == pop]["Height"]
        sns.kdeplot(data=subset, label=pop, ax=ax2, fill=True, alpha=0.3)
    ax2.set_title("Density Plot: Height Distribution")
    ax2.set_xlabel("Height")
    ax2.set_ylabel("Density")
    ax2.legend()

    # 3. Violin plot (Yield)
    ax3 = plt.subplot(3, 3, 3)
    sns.violinplot(data=df, x="Population", y="Yield", ax=ax3)
    ax3.set_title("Violin Plot: Yield by Population")

    # 4. Heatmap of trait correlations
    ax4 = plt.subplot(3, 3, 4)
    corr_matrix = df[trait_cols].corr()
    sns.heatmap(
        corr_matrix, annot=True, cmap="RdBu_r", center=0, square=True, ax=ax4, fmt=".2f"
    )
    ax4.set_title("Trait Correlation Heatmap")

    # 5. Sample heatmap (first 20 samples)
    ax5 = plt.subplot(3, 3, 5)
    sample_data = df[trait_cols].head(20).T
    sns.heatmap(sample_data, cmap="YlOrRd", ax=ax5, cbar_kws={"label": "Value"})
    ax5.set_title("Sample-Trait Heatmap (first 20)")
    ax5.set_xlabel("Sample")

    # 6. Scatter plot: Height vs Weight
    ax6 = plt.subplot(3, 3, 6)
    for pop, color in zip(["PopA", "PopB", "PopC"], ["red", "green", "blue"]):
        subset = df[df["Population"] == pop]
        ax6.scatter(
            subset["Height"], subset["Weight"], label=pop, alpha=0.6, s=50, c=color
        )
    ax6.set_xlabel("Height")
    ax6.set_ylabel("Weight")
    ax6.set_title("Scatter: Height vs Weight")
    ax6.legend()
    ax6.grid(True, alpha=0.3)

    # 7. Distribution of all traits
    ax7 = plt.subplot(3, 3, 7)
    df[trait_cols].hist(bins=15, ax=plt.subplot(3, 3, 7))
    # This will be handled differently

    # Actually create separate small multiples for distributions
    ax7.clear()
    positions = [1, 2, 3, 4, 5]
    bp = ax7.boxplot(
        [df[t].values for t in trait_cols],
        positions=positions,
        labels=trait_cols,
        patch_artist=True,
    )
    for patch in bp["boxes"]:
        patch.set_facecolor("lightblue")
    ax7.set_title("All Traits: Box Plot Summary")
    ax7.set_ylabel("Standardized Values")

    # 8. Population comparison (bar chart)
    ax8 = plt.subplot(3, 3, 8)
    pop_stats = df.groupby("Population")[trait_cols].mean()
    pop_stats.plot(kind="bar", ax=ax8)
    ax8.set_title("Mean Traits by Population")
    ax8.set_xlabel("Population")
    ax8.set_ylabel("Mean Value")
    ax8.tick_params(axis="x", rotation=0)
    ax8.legend(loc="upper right", fontsize=8)

    # 9. Outlier detection
    ax9 = plt.subplot(3, 3, 9)
    z_scores = np.abs((df[trait_cols] - df[trait_cols].mean()) / df[trait_cols].std())
    outliers = (z_scores > 2.5).sum(axis=1)
    ax9.hist(outliers, bins=range(0, 7), color="orange", alpha=0.7, edgecolor="black")
    ax9.set_xlabel("Number of Outlier Traits")
    ax9.set_ylabel("Number of Samples")
    ax9.set_title("Outlier Detection")
    ax9.axvline(x=2, color="red", linestyle="--", label="Outlier threshold")
    ax9.legend()

    plt.tight_layout()
    plt.savefig(f"{output_dir}/phenotype_plots.png", dpi=150, bbox_inches="tight")
    print(f"Phenotype plots saved: {output_dir}/phenotype_plots.png")


def create_scatter_matrix(df, output_dir):
    """Create scatter plot matrix"""
    import matplotlib.pyplot as plt
    import seaborn as sns

    print("Creating scatter plot matrix...")

    trait_cols = ["Height", "Weight", "Yield", "Protein", "Oil"]

    # Create pairplot
    g = sns.pairplot(
        df, vars=trait_cols, hue="Population", diag_kind="kde", plot_kws={"alpha": 0.6}
    )
    g.fig.suptitle("Scatter Plot Matrix: Trait Correlations", y=1.02)

    plt.savefig(f"{output_dir}/scatter_matrix.png", dpi=150, bbox_inches="tight")
    print(f"Scatter matrix saved: {output_dir}/scatter_matrix.png")


def main():
    print("=" * 60)
    print("Example: Phenotype Visualization Plots")
    print("=" * 60)
    print("\nThis script creates comprehensive phenotype visualizations.")
    print("Explore your data before running QTL analysis.")

    # Install packages
    print("\n[1/3] Installing dependencies...")
    install_packages()

    # Setup
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    # Generate data
    print("\n[2/3] Generating phenotype data...")
    df = generate_phenotype_data()
    print(f"  Generated: {len(df)} samples, 5 traits")
    print(f"\n  Trait summary:")
    print(df[["Height", "Weight", "Yield", "Protein", "Oil"]].describe().round(2))

    # Create plots
    print("\n[3/3] Creating visualizations...")
    create_phenotype_plots(df, output_dir)
    create_scatter_matrix(df, output_dir)

    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print("Created plots:")
    print("  ✓ Box plots - show quartiles and outliers")
    print("  ✓ Density plots - show distribution shapes")
    print("  ✓ Violin plots - distribution + summary statistics")
    print("  ✓ Heatmaps - trait correlations and sample patterns")
    print("  ✓ Scatter matrix - pairwise correlations")
    print(f"\nOutputs:")
    print(f"  - {output_dir}/phenotype_plots.png")
    print(f"  - {output_dir}/scatter_matrix.png")
    print("\nKey Insights:")
    print("  • Check for outliers before GWAS")
    print("  • Verify trait distributions are suitable")
    print("  • Look for population structure effects")
    print("  • Identify highly correlated traits")
    print("\n✅ Phenotype plots example complete!")
    print("\nIn QTLmax: Visualization → Phenotype plots")


if __name__ == "__main__":
    main()
