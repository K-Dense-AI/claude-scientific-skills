#!/usr/bin/env python3
"""
Clustering analysis script with multiple algorithms and evaluation.
Demonstrates k-means, DBSCAN, and hierarchical clustering with visualization.
"""

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns


def scale_data(X):
    """
    Scale features using StandardScaler.
    ALWAYS scale data before clustering!

    Args:
        X: Feature matrix

    Returns:
        Scaled feature matrix and fitted scaler
    """
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    return X_scaled, scaler


def find_optimal_k(X_scaled, k_range=range(2, 11)):
    """
    Find optimal number of clusters using elbow method and silhouette score.

    Args:
        X_scaled: Scaled feature matrix
        k_range: Range of k values to try

    Returns:
        Dictionary with inertias and silhouette scores
    """
    inertias = []
    silhouette_scores = []

    for k in k_range:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = kmeans.fit_predict(X_scaled)
        inertias.append(kmeans.inertia_)
        silhouette_scores.append(silhouette_score(X_scaled, labels))

    return {
        'k_values': list(k_range),
        'inertias': inertias,
        'silhouette_scores': silhouette_scores
    }


def plot_elbow_silhouette(results):
    """
    Plot elbow method and silhouette scores.

    Args:
        results: Dictionary from find_optimal_k
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Elbow plot
    ax1.plot(results['k_values'], results['inertias'], 'bo-')
    ax1.set_xlabel('Number of clusters (k)')
    ax1.set_ylabel('Inertia')
    ax1.set_title('Elbow Method')
    ax1.grid(True, alpha=0.3)

    # Silhouette plot
    ax2.plot(results['k_values'], results['silhouette_scores'], 'ro-')
    ax2.set_xlabel('Number of clusters (k)')
    ax2.set_ylabel('Silhouette Score')
    ax2.set_title('Silhouette Score vs k')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('elbow_silhouette.png', dpi=300, bbox_inches='tight')
    print("Saved elbow and silhouette plots to 'elbow_silhouette.png'")
    plt.close()


def evaluate_clustering(X_scaled, labels, algorithm_name):
    """
    Evaluate clustering using multiple metrics.

    Args:
        X_scaled: Scaled feature matrix
        labels: Cluster labels
        algorithm_name: Name of clustering algorithm

    Returns:
        Dictionary with evaluation metrics
    """
    # Filter out noise points for DBSCAN (-1 labels)
    mask = labels != -1
    X_filtered = X_scaled[mask]
    labels_filtered = labels[mask]

    n_clusters = len(set(labels_filtered))
    n_noise = list(labels).count(-1)

    results = {
        'algorithm': algorithm_name,
        'n_clusters': n_clusters,
        'n_noise': n_noise
    }

    # Calculate metrics if we have valid clusters
    if n_clusters > 1:
        results['silhouette'] = silhouette_score(X_filtered, labels_filtered)
        results['davies_bouldin'] = davies_bouldin_score(X_filtered, labels_filtered)
        results['calinski_harabasz'] = calinski_harabasz_score(X_filtered, labels_filtered)
    else:
        results['silhouette'] = None
        results['davies_bouldin'] = None
        results['calinski_harabasz'] = None

    return results


def perform_kmeans(X_scaled, n_clusters=3):
    """
    Perform k-means clustering.

    Args:
        X_scaled: Scaled feature matrix
        n_clusters: Number of clusters

    Returns:
        Fitted KMeans model and labels
    """
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    labels = kmeans.fit_predict(X_scaled)
    return kmeans, labels


def perform_dbscan(X_scaled, eps=0.5, min_samples=5):
    """
    Perform DBSCAN clustering.

    Args:
        X_scaled: Scaled feature matrix
        eps: Maximum distance between neighbors
        min_samples: Minimum points to form dense region

    Returns:
        Fitted DBSCAN model and labels
    """
    dbscan = DBSCAN(eps=eps, min_samples=min_samples)
    labels = dbscan.fit_predict(X_scaled)
    return dbscan, labels


def perform_hierarchical(X_scaled, n_clusters=3, linkage='ward'):
    """
    Perform hierarchical clustering.

    Args:
        X_scaled: Scaled feature matrix
        n_clusters: Number of clusters
        linkage: Linkage criterion ('ward', 'complete', 'average', 'single')

    Returns:
        Fitted AgglomerativeClustering model and labels
    """
    hierarchical = AgglomerativeClustering(n_clusters=n_clusters, linkage=linkage)
    labels = hierarchical.fit_predict(X_scaled)
    return hierarchical, labels


def visualize_clusters_2d(X_scaled, labels, algorithm_name, method='pca'):
    """
    Visualize clusters in 2D using PCA or t-SNE.

    Args:
        X_scaled: Scaled feature matrix
        labels: Cluster labels
        algorithm_name: Name of algorithm for title
        method: 'pca' or 'tsne'
    """
    # Reduce to 2D
    if method == 'pca':
        pca = PCA(n_components=2, random_state=42)
        X_2d = pca.fit_transform(X_scaled)
        variance = pca.explained_variance_ratio_
        xlabel = f'PC1 ({variance[0]:.1%} variance)'
        ylabel = f'PC2 ({variance[1]:.1%} variance)'
    else:
        from sklearn.manifold import TSNE
        # Use PCA first to speed up t-SNE
        pca = PCA(n_components=min(50, X_scaled.shape[1]), random_state=42)
        X_pca = pca.fit_transform(X_scaled)
        tsne = TSNE(n_components=2, random_state=42, perplexity=30)
        X_2d = tsne.fit_transform(X_pca)
        xlabel = 't-SNE 1'
        ylabel = 't-SNE 2'

    # Plot
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(X_2d[:, 0], X_2d[:, 1], c=labels, cmap='viridis', alpha=0.6, s=50)
    plt.colorbar(scatter, label='Cluster')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(f'{algorithm_name} Clustering ({method.upper()})')
    plt.grid(True, alpha=0.3)

    filename = f'{algorithm_name.lower().replace(" ", "_")}_{method}.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Saved visualization to '{filename}'")
    plt.close()


def main():
    """
    Example clustering analysis workflow.
    """
    # Load your data here
    # X = load_data()

    # Example with synthetic data
    from sklearn.datasets import make_blobs
    X, y_true = make_blobs(
        n_samples=500,
        n_features=10,
        centers=4,
        cluster_std=1.0,
        random_state=42
    )

    print(f"Dataset shape: {X.shape}")

    # Scale data (ALWAYS scale for clustering!)
    print("\nScaling data...")
    X_scaled, scaler = scale_data(X)

    # Find optimal k
    print("\nFinding optimal number of clusters...")
    results = find_optimal_k(X_scaled)
    plot_elbow_silhouette(results)

    # Based on elbow/silhouette, choose optimal k
    optimal_k = 4  # Adjust based on plots

    # Perform k-means
    print(f"\nPerforming k-means with k={optimal_k}...")
    kmeans, kmeans_labels = perform_kmeans(X_scaled, n_clusters=optimal_k)
    kmeans_results = evaluate_clustering(X_scaled, kmeans_labels, 'K-Means')

    # Perform DBSCAN
    print("\nPerforming DBSCAN...")
    dbscan, dbscan_labels = perform_dbscan(X_scaled, eps=0.5, min_samples=5)
    dbscan_results = evaluate_clustering(X_scaled, dbscan_labels, 'DBSCAN')

    # Perform hierarchical clustering
    print("\nPerforming hierarchical clustering...")
    hierarchical, hier_labels = perform_hierarchical(X_scaled, n_clusters=optimal_k)
    hier_results = evaluate_clustering(X_scaled, hier_labels, 'Hierarchical')

    # Print results
    print("\n" + "="*60)
    print("CLUSTERING RESULTS")
    print("="*60)

    for results in [kmeans_results, dbscan_results, hier_results]:
        print(f"\n{results['algorithm']}:")
        print(f"  Clusters: {results['n_clusters']}")
        if results['n_noise'] > 0:
            print(f"  Noise points: {results['n_noise']}")
        if results['silhouette']:
            print(f"  Silhouette Score: {results['silhouette']:.3f}")
            print(f"  Davies-Bouldin Index: {results['davies_bouldin']:.3f} (lower is better)")
            print(f"  Calinski-Harabasz Index: {results['calinski_harabasz']:.1f} (higher is better)")

    # Visualize clusters
    print("\nCreating visualizations...")
    visualize_clusters_2d(X_scaled, kmeans_labels, 'K-Means', method='pca')
    visualize_clusters_2d(X_scaled, dbscan_labels, 'DBSCAN', method='pca')
    visualize_clusters_2d(X_scaled, hier_labels, 'Hierarchical', method='pca')

    print("\nClustering analysis complete!")


if __name__ == "__main__":
    main()
