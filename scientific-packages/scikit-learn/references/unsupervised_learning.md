# Unsupervised Learning in scikit-learn

## Overview
Unsupervised learning discovers patterns in data without labeled targets. Main tasks include clustering (grouping similar samples), dimensionality reduction (reducing feature count), and anomaly detection (finding outliers).

## Clustering Algorithms

### K-Means

Groups data into k clusters by minimizing within-cluster variance.

**Algorithm**:
1. Initialize k centroids (k-means++ initialization recommended)
2. Assign each point to nearest centroid
3. Update centroids to mean of assigned points
4. Repeat until convergence

```python
from sklearn.cluster import KMeans

kmeans = KMeans(
    n_clusters=3,
    init='k-means++',  # Smart initialization
    n_init=10,         # Number of times to run with different seeds
    max_iter=300,
    random_state=42
)
labels = kmeans.fit_predict(X)
centroids = kmeans.cluster_centers_
```

**Use cases**:
- Customer segmentation
- Image compression
- Data preprocessing (clustering as features)

**Strengths**:
- Fast and scalable
- Simple to understand
- Works well with spherical clusters

**Limitations**:
- Assumes spherical clusters of similar size
- Sensitive to initialization (mitigated by k-means++)
- Must specify k beforehand
- Sensitive to outliers

**Choosing k**: Use elbow method, silhouette score, or domain knowledge

**Variants**:
- **MiniBatchKMeans**: Faster for large datasets, uses mini-batches
- **KMeans with n_init='auto'**: Adaptive number of initializations

### DBSCAN

Density-Based Spatial Clustering of Applications with Noise. Identifies clusters as dense regions separated by sparse areas.

```python
from sklearn.cluster import DBSCAN

dbscan = DBSCAN(
    eps=0.5,           # Maximum distance between neighbors
    min_samples=5,     # Minimum points to form dense region
    metric='euclidean'
)
labels = dbscan.fit_predict(X)
# -1 indicates noise/outliers
```

**Use cases**:
- Arbitrary cluster shapes
- Outlier detection
- When cluster count is unknown
- Geographic/spatial data

**Strengths**:
- Discovers arbitrary-shaped clusters
- Automatically detects outliers
- Doesn't require specifying number of clusters
- Robust to outliers

**Limitations**:
- Struggles with varying densities
- Sensitive to eps and min_samples parameters
- Not deterministic (border points may vary)

**Parameter tuning**:
- `eps`: Plot k-distance graph, look for elbow
- `min_samples`: Rule of thumb: 2 * dimensions

### HDBSCAN

Hierarchical DBSCAN that handles variable cluster densities.

```python
from sklearn.cluster import HDBSCAN

hdbscan = HDBSCAN(
    min_cluster_size=5,
    min_samples=None,  # Defaults to min_cluster_size
    metric='euclidean'
)
labels = hdbscan.fit_predict(X)
```

**Advantages over DBSCAN**:
- Handles variable density clusters
- More robust parameter selection
- Provides cluster membership probabilities
- Hierarchical structure

**Use cases**: When DBSCAN struggles with varying densities

### Hierarchical Clustering

Builds nested cluster hierarchies using agglomerative (bottom-up) approach.

```python
from sklearn.cluster import AgglomerativeClustering

agg_clust = AgglomerativeClustering(
    n_clusters=3,
    linkage='ward',  # 'ward', 'complete', 'average', 'single'
    metric='euclidean'
)
labels = agg_clust.fit_predict(X)

# Visualize with dendrogram
from scipy.cluster.hierarchy import dendrogram, linkage as scipy_linkage
import matplotlib.pyplot as plt

linkage_matrix = scipy_linkage(X, method='ward')
dendrogram(linkage_matrix)
plt.show()
```

**Linkage methods**:
- `ward`: Minimizes variance (only with Euclidean) - **most common**
- `complete`: Maximum distance between clusters
- `average`: Average distance between clusters
- `single`: Minimum distance between clusters

**Use cases**:
- When hierarchical structure is meaningful
- Taxonomy/phylogenetic trees
- When visualization is important (dendrograms)

**Strengths**:
- No need to specify k initially (cut dendrogram at desired level)
- Produces hierarchy of clusters
- Deterministic

**Limitations**:
- Computationally expensive (O(n²) to O(n³))
- Not suitable for large datasets
- Cannot undo previous merges

### Spectral Clustering

Performs dimensionality reduction using affinity matrix before clustering.

```python
from sklearn.cluster import SpectralClustering

spectral = SpectralClustering(
    n_clusters=3,
    affinity='rbf',  # 'rbf', 'nearest_neighbors', 'precomputed'
    gamma=1.0,
    n_neighbors=10,
    random_state=42
)
labels = spectral.fit_predict(X)
```

**Use cases**:
- Non-convex clusters
- Image segmentation
- Graph clustering
- When similarity matrix is available

**Strengths**:
- Handles non-convex clusters
- Works with similarity matrices
- Often better than k-means for complex shapes

**Limitations**:
- Computationally expensive
- Requires specifying number of clusters
- Memory intensive

### Mean Shift

Discovers clusters through iterative centroid updates based on density.

```python
from sklearn.cluster import MeanShift, estimate_bandwidth

# Estimate bandwidth
bandwidth = estimate_bandwidth(X, quantile=0.2, n_samples=500)

mean_shift = MeanShift(bandwidth=bandwidth)
labels = mean_shift.fit_predict(X)
cluster_centers = mean_shift.cluster_centers_
```

**Use cases**:
- When cluster count is unknown
- Computer vision applications
- Object tracking

**Strengths**:
- Automatically determines number of clusters
- Handles arbitrary shapes
- No assumptions about cluster shape

**Limitations**:
- Computationally expensive
- Very sensitive to bandwidth parameter
- Doesn't scale well

### Affinity Propagation

Uses message-passing between samples to identify exemplars.

```python
from sklearn.cluster import AffinityPropagation

affinity_prop = AffinityPropagation(
    damping=0.5,       # Damping factor (0.5-1.0)
    preference=None,   # Self-preference (controls number of clusters)
    random_state=42
)
labels = affinity_prop.fit_predict(X)
exemplars = affinity_prop.cluster_centers_indices_
```

**Use cases**:
- When number of clusters is unknown
- When exemplars (representative samples) are needed

**Strengths**:
- Automatically determines number of clusters
- Identifies exemplar samples
- No initialization required

**Limitations**:
- Very slow: O(n²t) where t is iterations
- Not suitable for large datasets
- Memory intensive

### Gaussian Mixture Models (GMM)

Probabilistic model assuming data comes from mixture of Gaussian distributions.

```python
from sklearn.mixture import GaussianMixture

gmm = GaussianMixture(
    n_components=3,
    covariance_type='full',  # 'full', 'tied', 'diag', 'spherical'
    random_state=42
)
labels = gmm.fit_predict(X)
probabilities = gmm.predict_proba(X)  # Soft clustering
```

**Covariance types**:
- `full`: Each component has its own covariance matrix
- `tied`: All components share same covariance
- `diag`: Diagonal covariance (independent features)
- `spherical`: Spherical covariance (isotropic)

**Use cases**:
- When soft clustering is needed (probabilities)
- When clusters have different shapes/sizes
- Generative modeling
- Density estimation

**Strengths**:
- Provides probabilities (soft clustering)
- Can handle elliptical clusters
- Generative model (can sample new data)
- Model selection with BIC/AIC

**Limitations**:
- Assumes Gaussian distributions
- Sensitive to initialization
- Can converge to local optima

**Model selection**:
```python
from sklearn.mixture import GaussianMixture
import numpy as np

n_components_range = range(2, 10)
bic_scores = []

for n in n_components_range:
    gmm = GaussianMixture(n_components=n, random_state=42)
    gmm.fit(X)
    bic_scores.append(gmm.bic(X))

optimal_n = n_components_range[np.argmin(bic_scores)]
```

### BIRCH

Builds Clustering Feature Tree for memory-efficient processing of large datasets.

```python
from sklearn.cluster import Birch

birch = Birch(
    n_clusters=3,
    threshold=0.5,
    branching_factor=50
)
labels = birch.fit_predict(X)
```

**Use cases**:
- Very large datasets
- Streaming data
- Memory constraints

**Strengths**:
- Memory efficient
- Single pass over data
- Incremental learning

## Dimensionality Reduction

### Principal Component Analysis (PCA)

Finds orthogonal components that explain maximum variance.

```python
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# Specify number of components
pca = PCA(n_components=2, random_state=42)
X_transformed = pca.fit_transform(X)

print("Explained variance ratio:", pca.explained_variance_ratio_)
print("Total variance explained:", pca.explained_variance_ratio_.sum())

# Or specify variance to retain
pca = PCA(n_components=0.95)  # Keep 95% of variance
X_transformed = pca.fit_transform(X)
print(f"Components needed: {pca.n_components_}")

# Visualize explained variance
plt.plot(np.cumsum(pca.explained_variance_ratio_))
plt.xlabel('Number of components')
plt.ylabel('Cumulative explained variance')
plt.show()
```

**Use cases**:
- Visualization (reduce to 2-3 dimensions)
- Remove multicollinearity
- Noise reduction
- Speed up training
- Feature extraction

**Strengths**:
- Fast and efficient
- Reduces multicollinearity
- Works well for linear relationships
- Interpretable components

**Limitations**:
- Only linear transformations
- Sensitive to scaling (always standardize first!)
- Components may be hard to interpret

**Variants**:
- **IncrementalPCA**: For datasets that don't fit in memory
- **KernelPCA**: Non-linear dimensionality reduction
- **SparsePCA**: Sparse loadings for interpretability

### t-SNE

t-Distributed Stochastic Neighbor Embedding for visualization.

```python
from sklearn.manifold import TSNE

tsne = TSNE(
    n_components=2,
    perplexity=30,      # Balance local vs global structure (5-50)
    learning_rate='auto',
    n_iter=1000,
    random_state=42
)
X_embedded = tsne.fit_transform(X)

# Visualize
plt.scatter(X_embedded[:, 0], X_embedded[:, 1], c=y)
plt.show()
```

**Use cases**:
- Visualization only (do not use for preprocessing!)
- Exploring high-dimensional data
- Finding clusters visually

**Important notes**:
- **Only for visualization**, not for preprocessing
- Each run produces different results (use random_state for reproducibility)
- Slow for large datasets
- Cannot transform new data (no transform() method)

**Parameter tuning**:
- `perplexity`: 5-50, larger for larger datasets
- Lower perplexity = focus on local structure
- Higher perplexity = focus on global structure

### UMAP

Uniform Manifold Approximation and Projection (requires umap-learn package).

**Advantages over t-SNE**:
- Preserves global structure better
- Faster
- Can transform new data
- Can be used for preprocessing (not just visualization)

### Truncated SVD (LSA)

Similar to PCA but works with sparse matrices (e.g., TF-IDF).

```python
from sklearn.decomposition import TruncatedSVD

svd = TruncatedSVD(n_components=100, random_state=42)
X_reduced = svd.fit_transform(X_sparse)
```

**Use cases**:
- Text data (after TF-IDF)
- Sparse matrices
- Latent Semantic Analysis (LSA)

### Non-negative Matrix Factorization (NMF)

Factorizes data into non-negative components.

```python
from sklearn.decomposition import NMF

nmf = NMF(n_components=10, init='nndsvd', random_state=42)
W = nmf.fit_transform(X)  # Document-topic matrix
H = nmf.components_        # Topic-word matrix
```

**Use cases**:
- Topic modeling
- Audio source separation
- Image processing
- When non-negativity is important (e.g., counts)

**Strengths**:
- Interpretable components (additive, non-negative)
- Sparse representations

### Independent Component Analysis (ICA)

Separates multivariate signal into independent components.

```python
from sklearn.decomposition import FastICA

ica = FastICA(n_components=10, random_state=42)
X_independent = ica.fit_transform(X)
```

**Use cases**:
- Blind source separation
- Signal processing
- Feature extraction when independence is expected

### Factor Analysis

Models observed variables as linear combinations of latent factors plus noise.

```python
from sklearn.decomposition import FactorAnalysis

fa = FactorAnalysis(n_components=5, random_state=42)
X_factors = fa.fit_transform(X)
```

**Use cases**:
- When noise is heteroscedastic
- Latent variable modeling
- Psychology/social science research

**Difference from PCA**: Models noise explicitly, assumes features have independent noise

## Anomaly Detection

### One-Class SVM

Learns boundary around normal data.

```python
from sklearn.svm import OneClassSVM

oc_svm = OneClassSVM(
    nu=0.1,           # Proportion of outliers expected
    kernel='rbf',
    gamma='auto'
)
oc_svm.fit(X_train)
predictions = oc_svm.predict(X_test)  # 1 for inliers, -1 for outliers
```

**Use cases**:
- Novelty detection
- When only normal data is available for training

### Isolation Forest

Isolates outliers using random forests.

```python
from sklearn.ensemble import IsolationForest

iso_forest = IsolationForest(
    contamination=0.1,  # Expected proportion of outliers
    random_state=42
)
predictions = iso_forest.fit_predict(X)  # 1 for inliers, -1 for outliers
scores = iso_forest.score_samples(X)     # Anomaly scores
```

**Use cases**:
- General anomaly detection
- Works well with high-dimensional data
- Fast and scalable

**Strengths**:
- Fast
- Effective in high dimensions
- Low memory requirements

### Local Outlier Factor (LOF)

Detects outliers based on local density deviation.

```python
from sklearn.neighbors import LocalOutlierFactor

lof = LocalOutlierFactor(
    n_neighbors=20,
    contamination=0.1
)
predictions = lof.fit_predict(X)  # 1 for inliers, -1 for outliers
scores = lof.negative_outlier_factor_  # Anomaly scores (negative)
```

**Use cases**:
- Finding local outliers
- When global methods fail

## Clustering Evaluation

### With Ground Truth Labels

When true labels are available (for validation):

**Adjusted Rand Index (ARI)**:
```python
from sklearn.metrics import adjusted_rand_score
ari = adjusted_rand_score(y_true, y_pred)
# Range: [-1, 1], 1 = perfect, 0 = random
```

**Normalized Mutual Information (NMI)**:
```python
from sklearn.metrics import normalized_mutual_info_score
nmi = normalized_mutual_info_score(y_true, y_pred)
# Range: [0, 1], 1 = perfect
```

**V-Measure**:
```python
from sklearn.metrics import v_measure_score
v = v_measure_score(y_true, y_pred)
# Range: [0, 1], harmonic mean of homogeneity and completeness
```

### Without Ground Truth Labels

When true labels are unavailable (unsupervised evaluation):

**Silhouette Score**:
Measures how similar objects are to their own cluster vs other clusters.

```python
from sklearn.metrics import silhouette_score, silhouette_samples
import matplotlib.pyplot as plt

score = silhouette_score(X, labels)
# Range: [-1, 1], higher is better
# >0.7: Strong structure
# 0.5-0.7: Reasonable structure
# 0.25-0.5: Weak structure
# <0.25: No substantial structure

# Per-sample scores for detailed analysis
sample_scores = silhouette_samples(X, labels)

# Visualize silhouette plot
for i in range(n_clusters):
    cluster_scores = sample_scores[labels == i]
    cluster_scores.sort()
    plt.barh(range(len(cluster_scores)), cluster_scores)
plt.axvline(x=score, color='red', linestyle='--')
plt.show()
```

**Davies-Bouldin Index**:
```python
from sklearn.metrics import davies_bouldin_score
db = davies_bouldin_score(X, labels)
# Lower is better, 0 = perfect
```

**Calinski-Harabasz Index** (Variance Ratio Criterion):
```python
from sklearn.metrics import calinski_harabasz_score
ch = calinski_harabasz_score(X, labels)
# Higher is better
```

**Inertia** (K-Means specific):
```python
inertia = kmeans.inertia_
# Sum of squared distances to nearest cluster center
# Use for elbow method
```

### Elbow Method (K-Means)

```python
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

inertias = []
K_range = range(2, 11)

for k in K_range:
    kmeans = KMeans(n_clusters=k, random_state=42)
    kmeans.fit(X)
    inertias.append(kmeans.inertia_)

plt.plot(K_range, inertias, 'bo-')
plt.xlabel('Number of clusters (k)')
plt.ylabel('Inertia')
plt.title('Elbow Method')
plt.show()
# Look for "elbow" where inertia starts decreasing more slowly
```

## Best Practices

### Clustering Algorithm Selection

**Use K-Means when**:
- Clusters are spherical and similar size
- Speed is important
- Data is not too high-dimensional

**Use DBSCAN when**:
- Arbitrary cluster shapes
- Number of clusters unknown
- Outlier detection needed

**Use Hierarchical when**:
- Hierarchy is meaningful
- Small to medium datasets
- Visualization is important

**Use GMM when**:
- Soft clustering needed
- Clusters have different shapes/sizes
- Probabilistic interpretation needed

**Use Spectral Clustering when**:
- Non-convex clusters
- Have similarity matrix
- Moderate dataset size

### Preprocessing for Clustering

1. **Always scale features**: Use StandardScaler or MinMaxScaler
2. **Handle outliers**: Remove or use robust algorithms (DBSCAN, HDBSCAN)
3. **Reduce dimensionality if needed**: PCA for speed, careful with interpretation
4. **Check for categorical variables**: Encode appropriately or use specialized algorithms

### Dimensionality Reduction Guidelines

**For preprocessing/feature extraction**:
- PCA (linear relationships)
- TruncatedSVD (sparse data)
- NMF (non-negative data)

**For visualization only**:
- t-SNE (preserves local structure)
- UMAP (preserves both local and global structure)

**Always**:
- Standardize features before PCA
- Use appropriate n_components (elbow plot, explained variance)
- Don't use t-SNE for anything except visualization

### Common Pitfalls

1. **Not scaling data**: Most algorithms sensitive to scale
2. **Using t-SNE for preprocessing**: Only for visualization!
3. **Overfitting cluster count**: Too many clusters = overfitting noise
4. **Ignoring outliers**: Can severely affect centroid-based methods
5. **Wrong metric**: Euclidean assumes all features equally important
6. **Not validating results**: Always check with multiple metrics and domain knowledge
7. **PCA without standardization**: Components dominated by high-variance features
