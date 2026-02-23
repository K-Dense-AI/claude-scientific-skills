#!/usr/bin/env python3
"""
Example: K-means Clustering for Population Structure

This example demonstrates K-means clustering to identify subpopulations
within a sample set based on genotype data. Useful when you don't know
the number of populations beforehand.

WHAT THIS MEANS:
K-means is an unsupervised learning algorithm that partitions samples into
K clusters based on genetic similarity. Unlike PCA, which shows continuous
variation, K-means assigns discrete cluster labels.

WHY WE DO THIS:
- Identify cryptic population structure
- Assign individuals to ancestry groups
- Quality control (detect outliers, mislabeled samples)
- Pre-GWAS stratification correction
- When admixture proportions aren't known

K-MEANS ALGORITHM:
1. Randomly initialize K cluster centroids
2. Assign each sample to nearest centroid
3. Update centroids as mean of assigned samples
4. Repeat until convergence

WHAT'S OUTPUT:
- Cluster assignments for each individual
- Cluster centroids (mean genotype profiles)
- Within-cluster sum of squares (WCSS)
- Visualization of clusters in PCA space

Equivalent to QTLmax: "Clustering subpopulations using K-means"
https://open.qtlmax.com/guide/index.php/2025/07/24/clustering-subpopulations-using-the-k-means-algorithm/

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


def generate_population_data(n_individuals=150, n_snps=1000):
    """Generate data from 3 distinct populations"""
    import numpy as np

    np.random.seed(42)

    # Generate 3 populations with distinct allele frequencies
    pop_freqs = [
        np.random.uniform(0.1, 0.3, n_snps),  # Population 1
        np.random.uniform(0.4, 0.6, n_snps),  # Population 2
        np.random.uniform(0.7, 0.9, n_snps),  # Population 3
    ]

    genotypes = []
    true_labels = []

    for pop_idx, freq in enumerate(pop_freqs):
        n_pop = n_individuals // 3
        for _ in range(n_pop):
            geno = np.random.binomial(2, freq)
            genotypes.append(geno)
            true_labels.append(f"Pop{pop_idx + 1}")

    return np.array(genotypes), true_labels


def run_kmeans(genotypes, k=3):
    """Run K-means clustering"""
    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA

    print(f"\nRunning K-means with K={k}...")

    # Standardize
    geno_std = (genotypes - genotypes.mean(axis=0)) / (genotypes.std(axis=0) + 1e-8)

    # K-means
    kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
    cluster_labels = kmeans.fit_predict(geno_std)

    # PCA for visualization
    pca = PCA(n_components=2)
    pc_scores = pca.fit_transform(geno_std)

    # Calculate WCSS for different K (elbow method)
    wcss = []
    for k_test in range(1, 10):
        km = KMeans(n_clusters=k_test, random_state=42, n_init=10)
        km.fit(geno_std)
        wcss.append(km.inertia_)

    return cluster_labels, kmeans, pc_scores, wcss


def create_cluster_plots(
    genotypes, true_labels, cluster_labels, pc_scores, wcss, output_dir
):
    """Create clustering visualization"""
    import matplotlib.pyplot as plt
    import numpy as np

    print("\nCreating cluster plots...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1. True labels
    ax1 = axes[0, 0]
    unique_true = list(set(true_labels))
    colors_true = ["red", "green", "blue"]
    for label, color in zip(unique_true, colors_true):
        mask = np.array(true_labels) == label
        ax1.scatter(
            pc_scores[mask, 0],
            pc_scores[mask, 1],
            c=color,
            label=label,
            alpha=0.6,
            s=50,
        )
    ax1.set_xlabel("PC1")
    ax1.set_ylabel("PC2")
    ax1.set_title("True Population Structure")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # 2. K-means clusters
    ax2 = axes[0, 1]
    colors_pred = ["red", "green", "blue"]
    for k in range(3):
        mask = cluster_labels == k
        ax2.scatter(
            pc_scores[mask, 0],
            pc_scores[mask, 1],
            c=colors_pred[k],
            label=f"Cluster {k + 1}",
            alpha=0.6,
            s=50,
        )
    ax2.set_xlabel("PC1")
    ax2.set_ylabel("PC2")
    ax2.set_title("K-means Clustering Results")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # 3. Elbow plot
    ax3 = axes[1, 0]
    ax3.plot(range(1, 10), wcss, "bo-", linewidth=2, markersize=8)
    ax3.set_xlabel("Number of Clusters (K)")
    ax3.set_ylabel("Within-Cluster Sum of Squares (WCSS)")
    ax3.set_title("Elbow Method for Optimal K")
    ax3.axvline(x=3, color="red", linestyle="--", label="Optimal K=3")
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # 4. Cluster comparison
    ax4 = axes[1, 1]

    # Create confusion-like matrix
    from sklearn.metrics import confusion_matrix

    cm = confusion_matrix(
        true_labels,
        [f"Cluster{i + 1}" for i in cluster_labels],
        labels=unique_true + [f"Cluster{i + 1}" for i in range(3)],
    )

    # Simplified: just show cluster sizes
    cluster_sizes = [np.sum(cluster_labels == i) for i in range(3)]
    ax4.bar(
        ["Cluster 1", "Cluster 2", "Cluster 3"],
        cluster_sizes,
        color=["red", "green", "blue"],
        alpha=0.7,
    )
    ax4.set_xlabel("Cluster")
    ax4.set_ylabel("Number of Individuals")
    ax4.set_title("Cluster Sizes")

    plt.tight_layout()
    plt.savefig(f"{output_dir}/kmeans_clustering.png", dpi=150, bbox_inches="tight")
    print(f"Clustering plots saved: {output_dir}/kmeans_clustering.png")


def main():
    print("=" * 60)
    print("Example: K-means Clustering")
    print("=" * 60)
    print("\nThis script demonstrates K-means for population clustering.")

    # Install packages
    print("\n[1/4] Installing dependencies...")
    install_packages()

    # Setup
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    # Generate data
    print("\n[2/4] Generating population data...")
    genotypes, true_labels = generate_population_data()
    print(f"  Generated: {len(genotypes)} individuals, {genotypes.shape[1]} SNPs")
    print(f"  True populations: {len(set(true_labels))}")

    # Run K-means
    print("\n[3/4] Running K-means clustering...")
    cluster_labels, kmeans, pc_scores, wcss = run_kmeans(genotypes, k=3)

    # Calculate accuracy
    from sklearn.metrics import adjusted_rand_score

    ari = adjusted_rand_score(true_labels, cluster_labels)
    print(f"  Adjusted Rand Index: {ari:.3f} (1.0 = perfect agreement)")

    # Create plots
    print("\n[4/4] Creating visualizations...")
    create_cluster_plots(
        genotypes, true_labels, cluster_labels, pc_scores, wcss, output_dir
    )

    # Summary
    print("\n" + "=" * 60)
    print("K-means Summary")
    print("=" * 60)
    print(f"Individuals: {len(genotypes)}")
    print(f"SNPs: {genotypes.shape[1]}")
    print(f"Clusters (K): 3")
    print(f"Agreement with truth: {ari:.1%}")

    for k in range(3):
        n = np.sum(cluster_labels == k)
        print(f"  Cluster {k + 1}: {n} individuals")

    print(f"\nOutputs:")
    print(f"  - {output_dir}/kmeans_clustering.png")
    print("\n✅ K-means example complete!")
    print("\nIn QTLmax: Population Structure → K-means")


if __name__ == "__main__":
    main()
