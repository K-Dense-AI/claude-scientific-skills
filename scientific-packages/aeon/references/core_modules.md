# Core Modules: Transformations, Distances, Networks, Datasets, and Benchmarking

This reference provides comprehensive details on foundational modules that support aeon's learning tasks.

## Transformations

Transformations convert time series into alternative representations for feature extraction, preprocessing, or visualization.

### Two Types of Transformers

**Collection Transformers**: Process entire collections of time series
- Input: `(n_cases, n_channels, n_timepoints)`
- Output: Features, transformed collections, or tabular data

**Series Transformers**: Work on individual time series
- Input: Single time series
- Output: Transformed single series

### Collection-Level Transformations

#### ROCKET (RAndom Convolutional KErnel Transform)

Fast feature extraction via random convolutional kernels:

```python
from aeon.transformations.collection.convolution_based import Rocket

rocket = Rocket(num_kernels=10000, n_jobs=-1)
X_transformed = rocket.fit_transform(X_train)
# Output shape: (n_cases, 2 * num_kernels)
```

**Variants:**
```python
from aeon.transformations.collection.convolution_based import (
    MiniRocket,
    MultiRocket,
    Hydra
)

# MiniRocket: Faster, streamlined version
minirocket = MiniRocket(num_kernels=10000)
X_features = minirocket.fit_transform(X_train)

# MultiRocket: Multivariate extensions
multirocket = MultiRocket(num_kernels=10000)
X_features = multirocket.fit_transform(X_train)

# Hydra: Dictionary-based convolution
hydra = Hydra(n_kernels=8)
X_features = hydra.fit_transform(X_train)
```

#### Catch22

22 canonical time series features:

```python
from aeon.transformations.collection.feature_based import Catch22

catch22 = Catch22(n_jobs=-1)
X_features = catch22.fit_transform(X_train)
# Output shape: (n_cases, 22)
```

**Feature categories:**
- Distribution (mean, variance, skewness)
- Autocorrelation properties
- Entropy measures
- Nonlinear dynamics
- Spectral properties

#### TSFresh

Comprehensive feature extraction (779 features):

```python
from aeon.transformations.collection.feature_based import TSFresh

tsfresh = TSFresh(
    default_fc_parameters="comprehensive",
    n_jobs=-1
)
X_features = tsfresh.fit_transform(X_train)
```

**Warning**: Slow on large datasets; use Catch22 for faster alternative

#### FreshPRINCE

Fresh Pipelines with Random Interval and Catch22 Features:

```python
from aeon.transformations.collection.feature_based import FreshPRINCE

freshprince = FreshPRINCE(n_intervals=50, n_jobs=-1)
X_features = freshprince.fit_transform(X_train)
```

#### Shapelet Transform

Extract discriminative subsequences:

```python
from aeon.transformations.collection.shapelet_based import ShapeletTransform

shapelet = ShapeletTransform(
    n_shapelet_samples=10000,
    max_shapelets=20,
    n_jobs=-1
)
X_features = shapelet.fit_transform(X_train, y_train)
# Requires labels for supervised shapelet discovery
```

**Random Shapelet Transform**:
```python
from aeon.transformations.collection.shapelet_based import RandomShapeletTransform

rst = RandomShapeletTransform(n_shapelets=1000)
X_features = rst.fit_transform(X_train)
```

#### SAST (Shapelet-Attention Subsequence Transform)

Attention-based shapelet discovery:

```python
from aeon.transformations.collection.shapelet_based import SAST

sast = SAST(window_size=0.1, n_shapelets=100)
X_features = sast.fit_transform(X_train, y_train)
```

#### Symbolic Representations

**SAX (Symbolic Aggregate approXimation)**:
```python
from aeon.transformations.collection.dictionary_based import SAX

sax = SAX(n_segments=8, alphabet_size=4)
X_symbolic = sax.fit_transform(X_train)
```

**PAA (Piecewise Aggregate Approximation)**:
```python
from aeon.transformations.collection.dictionary_based import PAA

paa = PAA(n_segments=10)
X_approximated = paa.fit_transform(X_train)
```

**SFA (Symbolic Fourier Approximation)**:
```python
from aeon.transformations.collection.dictionary_based import SFA

sfa = SFA(word_length=8, alphabet_size=4)
X_symbolic = sfa.fit_transform(X_train)
```

#### Channel Selection and Operations

**Channel Selection**:
```python
from aeon.transformations.collection.channel_selection import ChannelSelection

selector = ChannelSelection(channels=[0, 2, 5])
X_selected = selector.fit_transform(X_train)
```

**Channel Scoring**:
```python
from aeon.transformations.collection.channel_selection import ChannelScorer

scorer = ChannelScorer()
scores = scorer.fit_transform(X_train, y_train)
```

#### Data Balancing

**SMOTE (Synthetic Minority Over-sampling)**:
```python
from aeon.transformations.collection.smote import SMOTE

smote = SMOTE(k_neighbors=5)
X_resampled, y_resampled = smote.fit_resample(X_train, y_train)
```

**ADASYN**:
```python
from aeon.transformations.collection.smote import ADASYN

adasyn = ADASYN(n_neighbors=5)
X_resampled, y_resampled = adasyn.fit_resample(X_train, y_train)
```

### Series-Level Transformations

#### Smoothing Filters

**Moving Average**:
```python
from aeon.transformations.series.moving_average import MovingAverage

ma = MovingAverage(window_size=5)
X_smoothed = ma.fit_transform(X_series)
```

**Exponential Smoothing**:
```python
from aeon.transformations.series.exponent import ExponentTransformer

exp_smooth = ExponentTransformer(power=0.5)
X_smoothed = exp_smooth.fit_transform(X_series)
```

**Savitzky-Golay Filter**:
```python
from aeon.transformations.series.savgol import SavitzkyGolay

savgol = SavitzkyGolay(window_length=11, polyorder=3)
X_smoothed = savgol.fit_transform(X_series)
```

**Gaussian Filter**:
```python
from aeon.transformations.series.gaussian import GaussianFilter

gaussian = GaussianFilter(sigma=2.0)
X_smoothed = gaussian.fit_transform(X_series)
```

#### Statistical Transforms

**Box-Cox Transformation**:
```python
from aeon.transformations.series.boxcox import BoxCoxTransformer

boxcox = BoxCoxTransformer()
X_transformed = boxcox.fit_transform(X_series)
```

**AutoCorrelation**:
```python
from aeon.transformations.series.acf import AutoCorrelationTransformer

acf = AutoCorrelationTransformer(n_lags=40)
X_acf = acf.fit_transform(X_series)
```

**PCA (Principal Component Analysis)**:
```python
from aeon.transformations.series.pca import PCATransformer

pca = PCATransformer(n_components=3)
X_reduced = pca.fit_transform(X_series)
```

#### Approximation Methods

**Discrete Fourier Transform (DFT)**:
```python
from aeon.transformations.series.fourier import FourierTransform

dft = FourierTransform()
X_freq = dft.fit_transform(X_series)
```

**Piecewise Linear Approximation (PLA)**:
```python
from aeon.transformations.series.pla import PLA

pla = PLA(n_segments=10)
X_approx = pla.fit_transform(X_series)
```

#### Anomaly Detection Transform

**DOBIN (Distance-based Outlier BasIs using Neighbors)**:
```python
from aeon.transformations.series.dobin import DOBIN

dobin = DOBIN()
X_transformed = dobin.fit_transform(X_series)
```

### Transformation Pipelines

Chain transformers together:

```python
from sklearn.pipeline import Pipeline
from aeon.transformations.collection import Catch22, PCA

pipeline = Pipeline([
    ('features', Catch22()),
    ('reduce', PCA(n_components=10))
])
X_transformed = pipeline.fit_transform(X_train)
```

## Distance Metrics

Specialized distance functions for time series similarity measurement.

### Distance Categories

#### Warping-Based Distances

**DTW (Dynamic Time Warping)**:
```python
from aeon.distances import dtw_distance, dtw_pairwise_distance

# Compute distance between two series
dist = dtw_distance(series1, series2, window=0.2)

# Pairwise distances for a collection
dist_matrix = dtw_pairwise_distance(X_collection)

# Get alignment path
from aeon.distances import dtw_alignment_path
path = dtw_alignment_path(series1, series2)

# Get cost matrix
from aeon.distances import dtw_cost_matrix
cost = dtw_cost_matrix(series1, series2)
```

**DTW Variants**:
```python
from aeon.distances import (
    wdtw_distance,    # Weighted DTW
    ddtw_distance,    # Derivative DTW
    wddtw_distance,   # Weighted Derivative DTW
    adtw_distance,    # Amerced DTW
    shape_dtw_distance # Shape DTW
)

# Weighted DTW (penalize warping)
dist = wdtw_distance(series1, series2, g=0.05)

# Derivative DTW (compare shapes)
dist = ddtw_distance(series1, series2)

# Shape DTW (with shape descriptors)
dist = shape_dtw_distance(series1, series2)
```

**DTW Parameters**:
- `window`: Sakoe-Chiba band constraint (0.0-1.0)
- `g`: Penalty weight for warping distances

#### Edit Distances

**ERP (Edit distance with Real Penalty)**:
```python
from aeon.distances import erp_distance

dist = erp_distance(series1, series2, g=0.0, window=None)
```

**EDR (Edit Distance on Real sequences)**:
```python
from aeon.distances import edr_distance

dist = edr_distance(series1, series2, epsilon=0.1, window=None)
```

**LCSS (Longest Common SubSequence)**:
```python
from aeon.distances import lcss_distance

dist = lcss_distance(series1, series2, epsilon=1.0, window=None)
```

**TWE (Time Warp Edit)**:
```python
from aeon.distances import twe_distance

dist = twe_distance(series1, series2, penalty=0.1, stiffness=0.001)
```

#### Standard Metrics

```python
from aeon.distances import (
    euclidean_distance,
    manhattan_distance,
    minkowski_distance,
    squared_distance
)

# Euclidean distance
dist = euclidean_distance(series1, series2)

# Manhattan (L1) distance
dist = manhattan_distance(series1, series2)

# Minkowski distance
dist = minkowski_distance(series1, series2, p=3)

# Squared Euclidean
dist = squared_distance(series1, series2)
```

#### Specialized Distances

**MSM (Move-Split-Merge)**:
```python
from aeon.distances import msm_distance

dist = msm_distance(series1, series2, c=1.0)
```

**SBD (Shape-Based Distance)**:
```python
from aeon.distances import sbd_distance

dist = sbd_distance(series1, series2)
```

### Unified Distance Interface

```python
from aeon.distances import distance, pairwise_distance

# Compute any distance by name
dist = distance(series1, series2, metric="dtw", window=0.1)

# Pairwise distance matrix
dist_matrix = pairwise_distance(X_collection, metric="euclidean")

# Get available distance names
from aeon.distances import get_distance_function_names
available_distances = get_distance_function_names()
```

### Distance Selection Guide

**Fast and accurate**:
- Euclidean for aligned series
- Squared for even faster computation

**Handle temporal shifts**:
- DTW for general warping
- WDTW to penalize excessive warping

**Shape-based similarity**:
- DDTW or Shape DTW
- SBD for normalized shape comparison

**Robust to noise**:
- ERP, EDR, or LCSS

**Multivariate**:
- DTW supports multivariate via independent/dependent alignment

## Deep Learning Networks

Neural network architectures specialized for time series.

### Network Architectures

#### InceptionTime
Ensemble of Inception modules capturing multi-scale patterns:

```python
from aeon.networks import InceptionNetwork
from aeon.classification.deep_learning import InceptionTimeClassifier

# Use via classifier
clf = InceptionTimeClassifier(
    n_epochs=200,
    batch_size=64,
    n_ensemble=5
)

# Or use network directly
network = InceptionNetwork(
    n_classes=3,
    n_channels=1,
    n_timepoints=100
)
```

#### ResNet
Residual networks with skip connections:

```python
from aeon.networks import ResNetNetwork
from aeon.classification.deep_learning import ResNetClassifier

clf = ResNetClassifier(
    n_epochs=200,
    batch_size=64,
    n_res_blocks=3
)
```

#### FCN (Fully Convolutional Network)
```python
from aeon.networks import FCNNetwork
from aeon.classification.deep_learning import FCNClassifier

clf = FCNClassifier(
    n_epochs=200,
    batch_size=64,
    n_conv_layers=3
)
```

#### CNN
Standard convolutional architecture:

```python
from aeon.classification.deep_learning import CNNClassifier

clf = CNNClassifier(
    n_epochs=100,
    batch_size=32,
    kernel_size=7,
    n_filters=32
)
```

#### TapNet
Attentional prototype networks:

```python
from aeon.classification.deep_learning import TapNetClassifier

clf = TapNetClassifier(
    n_epochs=200,
    batch_size=64
)
```

#### MLP (Multi-Layer Perceptron)
```python
from aeon.classification.deep_learning import MLPClassifier

clf = MLPClassifier(
    n_epochs=100,
    batch_size=32,
    hidden_layer_sizes=[500]
)
```

#### LITE (Light Inception with boosTing tEchnique)
Lightweight ensemble network:

```python
from aeon.classification.deep_learning import LITEClassifier

clf = LITEClassifier(
    n_epochs=100,
    batch_size=64
)
```

### Training Configuration

```python
from aeon.classification.deep_learning import InceptionTimeClassifier

clf = InceptionTimeClassifier(
    n_epochs=200,
    batch_size=64,
    learning_rate=0.001,
    use_bias=True,
    verbose=1
)
clf.fit(X_train, y_train)
```

**Common parameters:**
- `n_epochs`: Training iterations
- `batch_size`: Samples per gradient update
- `learning_rate`: Optimizer learning rate
- `verbose`: Training output verbosity
- `callbacks`: Keras callbacks (early stopping, etc.)

## Datasets

Load built-in datasets and access UCR/UEA archives.

### Built-in Datasets

```python
from aeon.datasets import (
    load_arrow_head,
    load_airline,
    load_gunpoint,
    load_italy_power_demand,
    load_basic_motions,
    load_japanese_vowels
)

# Classification dataset
X_train, y_train = load_arrow_head(split="train")
X_test, y_test = load_arrow_head(split="test")

# Forecasting dataset (univariate series)
y = load_airline()

# Multivariate classification
X_train, y_train = load_basic_motions(split="train")
print(X_train.shape)  # (n_cases, n_channels, n_timepoints)
```

### UCR/UEA Archives

Access 100+ benchmark datasets:

```python
from aeon.datasets import load_from_tsfile, load_classification

# Load UCR/UEA dataset by name
X_train, y_train = load_classification("GunPoint", split="train")
X_test, y_test = load_classification("GunPoint", split="test")

# Load from local .ts file
X, y = load_from_tsfile("data/my_dataset_TRAIN.ts")
```

### Dataset Information

```python
from aeon.datasets import get_dataset_meta_data

# Get metadata about a dataset
info = get_dataset_meta_data("GunPoint")
print(info)
# {'n_cases': 150, 'n_timepoints': 150, 'n_classes': 2, ...}
```

### Custom Dataset Format

Save/load custom datasets in aeon format:

```python
from aeon.datasets import write_to_tsfile, load_from_tsfile

# Save
write_to_tsfile(
    X_train,
    "my_dataset_TRAIN.ts",
    y=y_train,
    problem_name="MyDataset"
)

# Load
X, y = load_from_tsfile("my_dataset_TRAIN.ts")
```

## Benchmarking

Tools for reproducible evaluation and comparison.

### Benchmarking Utilities

```python
from aeon.benchmarking import benchmark_estimator

# Benchmark a classifier on multiple datasets
results = benchmark_estimator(
    estimator=RocketClassifier(),
    datasets=["GunPoint", "ArrowHead", "ItalyPowerDemand"],
    n_resamples=10
)
```

### Result Storage and Comparison

```python
from aeon.benchmarking import (
    write_results_to_csv,
    read_results_from_csv,
    compare_results
)

# Save results
write_results_to_csv(results, "results.csv")

# Load and compare
results_rocket = read_results_from_csv("results_rocket.csv")
results_inception = read_results_from_csv("results_inception.csv")

comparison = compare_results(
    [results_rocket, results_inception],
    estimator_names=["ROCKET", "InceptionTime"]
)
```

### Critical Difference Diagrams

Visualize statistical significance of differences:

```python
from aeon.benchmarking.results_plotting import plot_critical_difference_diagram

plot_critical_difference_diagram(
    results_dict={
        'ROCKET': results_rocket,
        'InceptionTime': results_inception,
        'BOSS': results_boss
    },
    dataset_names=["GunPoint", "ArrowHead", "ItalyPowerDemand"]
)
```

## Discovery and Tags

### Finding Estimators

```python
from aeon.utils.discovery import all_estimators

# Get all classifiers
classifiers = all_estimators(type_filter="classifier")

# Get all transformers
transformers = all_estimators(type_filter="transformer")

# Filter by capability tags
multivariate_classifiers = all_estimators(
    type_filter="classifier",
    filter_tags={"capability:multivariate": True}
)
```

### Checking Estimator Tags

```python
from aeon.utils.tags import all_tags_for_estimator
from aeon.classification.convolution_based import RocketClassifier

tags = all_tags_for_estimator(RocketClassifier)
print(tags)
# {'capability:multivariate': True, 'X_inner_type': ['numpy3D'], ...}
```

### Common Tags

- `capability:multivariate`: Handles multivariate series
- `capability:unequal_length`: Handles variable-length series
- `capability:missing_values`: Handles missing data
- `algorithm_type`: Algorithm family (e.g., "convolution", "distance")
- `python_dependencies`: Required packages
