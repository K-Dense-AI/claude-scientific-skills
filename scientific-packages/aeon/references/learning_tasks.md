# Learning Tasks: Classification, Regression, Clustering, and Similarity Search

This reference provides comprehensive details on supervised and unsupervised learning tasks for time series collections.

## Time Series Classification

Time series classification (TSC) assigns labels to entire sequences. Aeon provides diverse algorithm families with unique strengths.

### Algorithm Categories

#### 1. Convolution-Based Classifiers
Transform time series using random convolutional kernels:

**ROCKET (RAndom Convolutional KErnel Transform)**
- Ultra-fast feature extraction via random kernels
- 10,000+ kernels generate discriminative features
- Linear classifier on extracted features

```python
from aeon.classification.convolution_based import RocketClassifier

clf = RocketClassifier(num_kernels=10000, n_jobs=-1)
clf.fit(X_train, y_train)
predictions = clf.predict(X_test)
probabilities = clf.predict_proba(X_test)
```

**Variants:**
- `MiniRocketClassifier`: Faster, streamlined version
- `MultiRocketClassifier`: Multivariate extensions
- `Arsenal`: Ensemble of ROCKET transformers
- `Hydra`: Dictionary-based convolution variant

#### 2. Deep Learning Classifiers
Neural networks specialized for time series:

**InceptionTime**
- Ensemble of Inception modules
- Captures patterns at multiple scales
- State-of-the-art on UCR benchmarks

```python
from aeon.classification.deep_learning import InceptionTimeClassifier

clf = InceptionTimeClassifier(n_epochs=200, batch_size=64)
clf.fit(X_train, y_train)
```

**Other architectures:**
- `ResNetClassifier`: Residual connections
- `FCNClassifier`: Fully Convolutional Networks
- `CNNClassifier`: Standard convolutional architecture
- `LITEClassifier`: Lightweight networks
- `MLPClassifier`: Multi-layer perceptrons
- `TapNetClassifier`: Attentional prototype networks

#### 3. Dictionary-Based Classifiers
Symbolic representations and bag-of-words approaches:

**BOSS (Bag of SFA Symbols)**
- Converts series to symbolic words
- Histogram-based classification
- Effective for shape patterns

```python
from aeon.classification.dictionary_based import BOSSEnsemble

clf = BOSSEnsemble(max_ensemble_size=500)
clf.fit(X_train, y_train)
```

**Other dictionary methods:**
- `TemporalDictionaryEnsemble (TDE)`: Enhanced BOSS with temporal info
- `WEASEL`: Word ExtrAction for time SEries cLassification
- `MUSE`: MUltivariate Symbolic Extension
- `MrSEQL`: Multiple Representations SEQuence Learner

#### 4. Distance-Based Classifiers
Leverage time series-specific distance metrics:

**K-Nearest Neighbors with DTW**
- Dynamic Time Warping handles temporal shifts
- Effective for shape-based similarity

```python
from aeon.classification.distance_based import KNeighborsTimeSeriesClassifier

clf = KNeighborsTimeSeriesClassifier(
    distance="dtw",
    n_neighbors=5
)
clf.fit(X_train, y_train)
```

**Other distance methods:**
- `ElasticEnsemble`: Ensemble of elastic distances
- `ProximityForest`: Tree-based with elastic measures
- `ProximityTree`: Single tree variant
- `ShapeDTW`: DTW with shape descriptors

#### 5. Feature-Based Classifiers
Extract statistical and domain-specific features:

**Catch22**
- 22 time series features
- Canonical Time-series CHaracteristics
- Fast and interpretable

```python
from aeon.classification.feature_based import Catch22Classifier

clf = Catch22Classifier(estimator=RandomForestClassifier())
clf.fit(X_train, y_train)
```

**Other feature methods:**
- `FreshPRINCEClassifier`: Fresh Pipelines with Random Interval and Catch22 Features
- `SignatureClassifier`: Path signature features
- `TSFreshClassifier`: Comprehensive feature extraction (slower, more features)
- `SummaryClassifier`: Simple summary statistics

#### 6. Interval-Based Classifiers
Analyze discriminative time intervals:

**Time Series Forest (TSF)**
- Random intervals + summary statistics
- Random forest on extracted features

```python
from aeon.classification.interval_based import TimeSeriesForestClassifier

clf = TimeSeriesForestClassifier(n_estimators=500)
clf.fit(X_train, y_train)
```

**Other interval methods:**
- `CanonicalIntervalForest (CIF)`: Canonical Interval Forest
- `DrCIF`: Diverse Representation CIF
- `RISE`: Random Interval Spectral Ensemble
- `RandomIntervalClassifier`: Basic random interval approach
- `STSF`: Shapelet Transform Interval Forest

#### 7. Shapelet-Based Classifiers
Discover discriminative subsequences:

**Shapelets**: Small subsequences that best distinguish classes

```python
from aeon.classification.shapelet_based import ShapeletTransformClassifier

clf = ShapeletTransformClassifier(
    n_shapelet_samples=10000,
    max_shapelets=20
)
clf.fit(X_train, y_train)
```

**Other shapelet methods:**
- `LearningShapeletClassifier`: Gradient-based learning
- `SASTClassifier`: Shapelet-Attention Subsequence Transform

#### 8. Hybrid Ensembles
Combine multiple algorithm families:

**HIVE-COTE (Hierarchical Vote Collective of Transformation-based Ensembles)**
- State-of-the-art accuracy
- Combines shapelets, intervals, dictionaries, and spectral features
- V2 uses ROCKET and improved components

```python
from aeon.classification.hybrid import HIVECOTEV2

clf = HIVECOTEV2(n_jobs=-1)  # Slow but highly accurate
clf.fit(X_train, y_train)
```

### Algorithm Selection Guide

**Fast and accurate (default choice):**
- `RocketClassifier` or `MiniRocketClassifier`

**Maximum accuracy (slow):**
- `HIVECOTEV2` or `InceptionTimeClassifier`

**Interpretable:**
- `Catch22Classifier` or `ShapeletTransformClassifier`

**Multivariate focus:**
- `MultiRocketClassifier` or `MUSE`

**Small datasets:**
- `KNeighborsTimeSeriesClassifier` with DTW

### Classification Workflow

```python
from aeon.classification.convolution_based import RocketClassifier
from aeon.datasets import load_arrow_head
from sklearn.metrics import accuracy_score, classification_report

# Load data
X_train, y_train = load_arrow_head(split="train")
X_test, y_test = load_arrow_head(split="test")

# Train classifier
clf = RocketClassifier(n_jobs=-1)
clf.fit(X_train, y_train)

# Evaluate
y_pred = clf.predict(X_test)
accuracy = accuracy_score(y_test, y_pred)
print(f"Accuracy: {accuracy:.3f}")
print(classification_report(y_test, y_pred))
```

## Time Series Regression

Time series regression predicts continuous values from sequences. Most classification algorithms have regression equivalents.

### Regression Algorithms

Available regressors mirror classification structure:
- `RocketRegressor`, `MiniRocketRegressor`, `MultiRocketRegressor`
- `InceptionTimeRegressor`, `ResNetRegressor`, `FCNRegressor`
- `KNeighborsTimeSeriesRegressor`
- `Catch22Regressor`, `FreshPRINCERegressor`
- `TimeSeriesForestRegressor`, `DrCIFRegressor`

### Regression Workflow

```python
from aeon.regression.convolution_based import RocketRegressor
from sklearn.metrics import mean_squared_error, r2_score

# Train regressor
reg = RocketRegressor(num_kernels=10000)
reg.fit(X_train, y_train_continuous)

# Predict and evaluate
y_pred = reg.predict(X_test)
mse = mean_squared_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)
print(f"MSE: {mse:.3f}, R²: {r2:.3f}")
```

## Time Series Clustering

Clustering groups similar time series without labels.

### Clustering Algorithms

**TimeSeriesKMeans**
- K-means with time series distances
- Supports DTW, Euclidean, and other metrics

```python
from aeon.clustering import TimeSeriesKMeans

clusterer = TimeSeriesKMeans(
    n_clusters=3,
    distance="dtw",
    n_init=10
)
clusterer.fit(X_collection)
labels = clusterer.labels_
```

**TimeSeriesKMedoids**
- Uses actual series as cluster centers
- More robust to outliers

```python
from aeon.clustering import TimeSeriesKMedoids

clusterer = TimeSeriesKMedoids(
    n_clusters=3,
    distance="euclidean"
)
clusterer.fit(X_collection)
```

**Other clustering methods:**
- `TimeSeriesKernelKMeans`: Kernel-based clustering
- `ElasticSOM`: Self-organizing maps with elastic distances

### Clustering Workflow

```python
from aeon.clustering import TimeSeriesKMeans
from aeon.distances import dtw_distance
import numpy as np

# Cluster time series
clusterer = TimeSeriesKMeans(n_clusters=4, distance="dtw")
clusterer.fit(X_train)

# Get cluster labels
labels = clusterer.predict(X_test)

# Compute cluster centers
centers = clusterer.cluster_centers_

# Evaluate clustering quality (if ground truth available)
from sklearn.metrics import adjusted_rand_score
ari = adjusted_rand_score(y_true, labels)
```

## Similarity Search

Similarity search finds motifs, nearest neighbors, and repeated patterns.

### Key Concepts

**Motifs**: Frequently repeated subsequences within a time series
**Matrix Profile**: Data structure encoding nearest neighbor distances for all subsequences

### Similarity Search Methods

**Matrix Profile**
- Efficient motif discovery
- Change point detection
- Anomaly detection

```python
from aeon.similarity_search import MatrixProfile

mp = MatrixProfile(window_size=50)
profile = mp.fit_transform(X_series)

# Find top motif
motif_idx = np.argmin(profile)
```

**Query Search**
- Find nearest neighbors to a query subsequence
- Useful for template matching

```python
from aeon.similarity_search import QuerySearch

searcher = QuerySearch(distance="euclidean")
distances, indices = searcher.search(X_series, query_subsequence)
```

### Similarity Search Workflow

```python
from aeon.similarity_search import MatrixProfile
import numpy as np

# Compute matrix profile
mp = MatrixProfile(window_size=100)
profile, profile_index = mp.fit_transform(X_series)

# Find top-k motifs (lowest profile values)
k = 3
motif_indices = np.argsort(profile)[:k]

# Find anomalies (highest profile values)
anomaly_indices = np.argsort(profile)[-k:]
```

## Ensemble and Composition Tools

### Voting Ensembles
```python
from aeon.classification.ensemble import WeightedEnsembleClassifier
from aeon.classification.convolution_based import RocketClassifier
from aeon.classification.distance_based import KNeighborsTimeSeriesClassifier

ensemble = WeightedEnsembleClassifier(
    estimators=[
        ('rocket', RocketClassifier()),
        ('knn', KNeighborsTimeSeriesClassifier())
    ]
)
ensemble.fit(X_train, y_train)
```

### Pipelines
```python
from sklearn.pipeline import Pipeline
from aeon.transformations.collection import Catch22
from sklearn.ensemble import RandomForestClassifier

pipeline = Pipeline([
    ('features', Catch22()),
    ('classifier', RandomForestClassifier())
])
pipeline.fit(X_train, y_train)
```

## Model Selection and Validation

### Cross-Validation
```python
from sklearn.model_selection import cross_val_score
from aeon.classification.convolution_based import RocketClassifier

clf = RocketClassifier()
scores = cross_val_score(clf, X_train, y_train, cv=5)
print(f"CV Accuracy: {scores.mean():.3f} (+/- {scores.std():.3f})")
```

### Grid Search
```python
from sklearn.model_selection import GridSearchCV
from aeon.classification.distance_based import KNeighborsTimeSeriesClassifier

param_grid = {
    'n_neighbors': [1, 3, 5, 7],
    'distance': ['dtw', 'euclidean', 'erp']
}

clf = KNeighborsTimeSeriesClassifier()
grid_search = GridSearchCV(clf, param_grid, cv=5)
grid_search.fit(X_train, y_train)
print(f"Best params: {grid_search.best_params_}")
```

## Discovery Functions

Find available estimators programmatically:

```python
from aeon.utils.discovery import all_estimators

# Get all classifiers
classifiers = all_estimators(type_filter="classifier")

# Get all regressors
regressors = all_estimators(type_filter="regressor")

# Get all clusterers
clusterers = all_estimators(type_filter="clusterer")

# Filter by tag (e.g., multivariate capable)
mv_classifiers = all_estimators(
    type_filter="classifier",
    filter_tags={"capability:multivariate": True}
)
```
