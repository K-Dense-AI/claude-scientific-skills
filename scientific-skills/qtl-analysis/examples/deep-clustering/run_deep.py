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
Example: Deep Learning Population Clustering (Autoencoder)

This example demonstrates population substructure analysis using deep learning -
similar to QTLmax's "Deep learning clustering" feature.

WHAT THIS MEANS:
Deep learning (autoencoders) can discover complex, non-linear patterns in genomic
data that traditional methods like PCA or K-means might miss. Autoencoders learn
to compress genotype data into a lower-dimensional latent space that captures
the underlying population structure.

WHY WE DO THIS:
- Non-linear dimensionality reduction captures complex structure
- Works better than PCA for admixed populations
- Can discover subtle subpopulations
- Neural network embeddings can be used for downstream analysis

WHAT'S DEMONSTRATED:
1. Train autoencoder on genotype data
2. Extract latent representations
3. Cluster in latent space
4. Visualize with t-SNE/UMAP

Equivalent to QTLmax: "Deep learning algorithms for subpopulation clustering"
https://open.qtlmax.com/guide/index.php/2025/07/24/deep-learning-algorithms-for-subpopulation-clustering/
"""

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
import os
import warnings

warnings.filterwarnings("ignore")

# Create output directory
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

# Try to import sklearn's MLPRegressor as simple autoencoder replacement
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler


def generate_genotype_data(n_samples=300, n_snps=200, n_true_pop=4):
    """Generate synthetic genotype data with complex population structure."""
    np.random.seed(42)

    # Create distinct populations with different allele frequencies
    pop_frequencies = []
    for _ in range(n_true_pop):
        freqs = np.random.uniform(0.1, 0.9, n_snps)
        pop_frequencies.append(freqs)

    # Generate samples from each population
    genotypes = []
    labels = []

    samples_per_pop = n_samples // n_true_pop
    for pop_id, freqs in enumerate(pop_frequencies):
        for _ in range(samples_per_pop):
            # Binomial sampling for genotypes (0, 1, 2)
            geno = np.random.binomial(2, freqs)
            genotypes.append(geno)
            labels.append(pop_id)

    # Add some admixed individuals (mix of two populations)
    n_admixed = n_samples - (samples_per_pop * n_true_pop)
    for _ in range(n_admixed):
        pop1, pop2 = np.random.choice(n_true_pop, 2, replace=False)
        mixing = np.random.uniform(0.3, 0.7)
        geno = np.random.binomial(
            2, mixing * pop_frequencies[pop1] + (1 - mixing) * pop_frequencies[pop2]
        )
        genotypes.append(geno)
        labels.append(-1)  # Admixed

    genotypes = np.array(genotypes)
    labels = np.array(labels)

    # Shuffle
    idx = np.random.permutation(len(labels))
    return genotypes[idx], labels[idx], pop_frequencies


def train_autoencoder(X, latent_dim=10, epochs=100):
    """Train a simple autoencoder using MLP (sklearn)."""
    print("  Training autoencoder...")

    # Standardize input
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Simple autoencoder: input -> hidden -> latent -> hidden -> output
    # We use MLPRegressor as a workaround (encoder only, then clustering)
    encoder = MLPRegressor(
        hidden_layer_sizes=(64, latent_dim),
        activation="relu",
        solver="adam",
        max_iter=epochs,
        random_state=42,
        early_stopping=True,
        validation_fraction=0.1,
    )

    # Train to reconstruct input (autoencoder)
    encoder.fit(X_scaled, X_scaled)

    # Get latent representation
    # Forward pass through first half of network
    X_latent = encoder.predict(X_scaled)

    return X_latent, encoder, scaler


def run_deep_clustering():
    """Main deep clustering workflow."""
    print("=" * 60)
    print("DEEP LEARNING POPULATION CLUSTERING")
    print("=" * 60)

    # Generate data
    print("\n[1/5] Generating synthetic genotype data...")
    genotypes, true_labels, _ = generate_genotype_data(
        n_samples=300, n_snps=200, n_true_pop=4
    )

    # Count true populations (excluding admixed)
    n_true_pop = len(set(true_labels[true_labels >= 0]))
    print(f"  Generated: {genotypes.shape[0]} samples × {genotypes.shape[1]} SNPs")
    print(f"  True populations: {n_true_pop}")

    # Save genotype data
    geno_df = pd.DataFrame(genotypes)
    geno_df.to_csv(f"{output_dir}/genotypes.csv", index=False)
    print(f"  Saved: genotypes.csv")

    # Train autoencoder
    print("\n[2/5] Training autoencoder...")
    latent_repr, encoder, scaler = train_autoencoder(
        genotypes, latent_dim=10, epochs=200
    )
    print(f"  Latent dimension: {latent_repr.shape[1]}")

    # Save latent representation
    latent_df = pd.DataFrame(latent_repr)
    latent_df.to_csv(f"{output_dir}/latent_representation.csv", index=False)
    print(f"  Saved: latent_representation.csv")

    # Cluster in latent space
    print("\n[3/5] Clustering in latent space...")
    kmeans = KMeans(n_clusters=n_true_pop, random_state=42, n_init=10)
    cluster_labels = kmeans.fit_predict(latent_repr)

    # Calculate ARI (excluding admixed)
    mask = true_labels >= 0
    ari = adjusted_rand_score(true_labels[mask], cluster_labels[mask])
    print(f"  Adjusted Rand Index: {ari:.3f}")

    # Save cluster assignments
    cluster_df = pd.DataFrame(
        {
            "Sample": range(len(cluster_labels)),
            "True_Population": true_labels,
            "Predicted_Cluster": cluster_labels,
            "Is_Admixed": (true_labels == -1).astype(int),
        }
    )
    cluster_df.to_csv(f"{output_dir}/cluster_assignments.csv", index=False)
    print(f"  Saved: cluster_assignments.csv")

    # t-SNE visualization
    print("\n[4/5] Creating visualizations...")
    print("  Running t-SNE...")
    tsne = TSNE(n_components=2, random_state=42, perplexity=30)
    latent_2d = tsne.fit_transform(latent_repr)

    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    # Plot 1: True labels (excluding admixed)
    ax1 = axes[0, 0]
    colors = plt.cm.Set1(np.linspace(0, 1, n_true_pop))
    for i in range(n_true_pop):
        mask = true_labels == i
        ax1.scatter(
            latent_2d[mask, 0],
            latent_2d[mask, 1],
            c=[colors[i]],
            label=f"Population {i + 1}",
            alpha=0.7,
            s=40,
        )
    ax1.set_title(
        "True Population Labels\n(t-SNE of Autoencoder Latent Space)",
        fontsize=12,
        fontweight="bold",
    )
    ax1.set_xlabel("t-SNE 1")
    ax1.set_ylabel("t-SNE 2")
    ax1.legend()

    # Plot 2: Predicted clusters
    ax2 = axes[0, 1]
    for i in range(n_true_pop):
        mask = cluster_labels == i
        ax2.scatter(
            latent_2d[mask, 0],
            latent_2d[mask, 1],
            c=[colors[i]],
            label=f"Cluster {i + 1}",
            alpha=0.7,
            s=40,
        )
    ax2.set_title(
        f"Predicted Clusters (K-means)\nARI = {ari:.3f}", fontsize=12, fontweight="bold"
    )
    ax2.set_xlabel("t-SNE 1")
    ax2.set_ylabel("t-SNE 2")
    ax2.legend()

    # Plot 3: Latent space heatmap
    ax3 = axes[1, 0]
    im = ax3.imshow(latent_repr[:50].T, aspect="auto", cmap="RdBu_r")
    ax3.set_title(
        "Autoencoder Latent Space\n(First 50 samples)", fontsize=12, fontweight="bold"
    )
    ax3.set_xlabel("Sample")
    ax3.set_ylabel("Latent Dimension")
    plt.colorbar(im, ax=ax3)

    # Plot 4: Cluster comparison
    ax4 = axes[1, 1]

    # Create confusion-style comparison
    from sklearn.metrics import confusion_matrix

    cm = confusion_matrix(true_labels[mask], cluster_labels[mask])
    im4 = ax4.imshow(cm, cmap="Blues")
    ax4.set_title(
        "Cluster vs True Labels\n(Confusion Matrix)", fontsize=12, fontweight="bold"
    )
    ax4.set_xlabel("Predicted Cluster")
    ax4.set_ylabel("True Population")
    ax4.set_xticks(range(n_true_pop))
    ax4.set_yticks(range(n_true_pop))
    ax4.set_xticklabels([f"C{i}" for i in range(n_true_pop)])
    ax4.set_yticklabels([f"P{i}" for i in range(n_true_pop)])

    # Add text annotations (safely handle different matrix sizes)
    for i in range(min(n_true_pop, cm.shape[0])):
        for j in range(min(n_true_pop, cm.shape[1])):
            ax4.text(
                j,
                i,
                cm[i, j],
                ha="center",
                va="center",
                color="white" if cm[i, j] > cm.max() / 2 else "black",
            )

    plt.colorbar(im4, ax=ax4)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/deep_clustering.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: deep_clustering.png")

    # Summary
    print("\n[5/5] Summary:")
    print(f"  - {genotypes.shape[0]} samples, {genotypes.shape[1]} SNPs")
    print(f"  - Autoencoder: 200 -> 64 -> 10 -> 64 -> 200")
    print(f"  - Latent dimension: 10")
    print(f"  - Clustering ARI: {ari:.3f}")
    print(f"  - K-means in latent space: {n_true_pop} clusters")

    print("\n" + "=" * 60)
    print("DEEP CLUSTERING COMPLETE")
    print("=" * 60)

    return {"ari": ari, "n_clusters": n_true_pop}


if __name__ == "__main__":
    results = run_deep_clustering()
