---
name: aeon
description: Time series machine learning toolkit for classification, regression, clustering, forecasting, anomaly detection, segmentation, and similarity search. Use this skill when working with temporal data, performing time series analysis, building predictive models on sequential data, or implementing workflows that involve distance metrics (DTW), transformations (ROCKET, Catch22), or deep learning for time series. Applicable for tasks like ECG classification, stock price forecasting, sensor anomaly detection, or activity recognition from wearable devices.
---

# Aeon

## Overview

Aeon is a comprehensive Python toolkit for time series machine learning, providing state-of-the-art algorithms and classical techniques for analyzing temporal data. Use this skill when working with sequential/temporal data across seven primary learning tasks: classification, regression, clustering, forecasting, anomaly detection, segmentation, and similarity search.

## When to Use This Skill

Apply this skill when:
- Classifying or predicting from time series data (e.g., ECG classification, activity recognition)
- Forecasting future values in temporal sequences (e.g., stock prices, energy demand)
- Detecting anomalies in sensor streams or operational data
- Clustering temporal patterns or discovering motifs
- Segmenting time series into meaningful regions (change point detection)
- Computing distances between time series using specialized metrics (DTW, MSM, ERP)
- Extracting features from temporal data using ROCKET, Catch22, TSFresh, or shapelets
- Building deep learning models for time series with specialized architectures

## Core Capabilities

### 1. Time Series Classification
Classify labeled time series using diverse algorithm families:
- **Convolution-based**: ROCKET, MiniRocket, MultiRocket, Arsenal, Hydra
- **Deep learning**: InceptionTime, ResNet, FCN, TimeCNN, LITE
- **Dictionary-based**: BOSS, TDE, WEASEL, MrSEQL (symbolic representations)
- **Distance-based**: KNN with elastic distances, Elastic Ensemble, Proximity Forest
- **Feature-based**: Catch22, FreshPRINCE, Signature classifiers
- **Interval-based**: CIF, DrCIF, RISE, Random Interval variants
- **Shapelet-based**: Learning Shapelet, SAST
- **Hybrid ensembles**: HIVE-COTE V1/V2

Example:
```python
from aeon.classification.convolution_based import RocketClassifier
from aeon.datasets import load_arrow_head

X_train, y_train = load_arrow_head(split="train")
X_test, y_test = load_arrow_head(split="test")

clf = RocketClassifier()
clf.fit(X_train, y_train)
accuracy = clf.score(X_test, y_test)
```

### 2. Time Series Regression
Predict continuous values from time series using adapted classification algorithms:
```python
from aeon.regression.convolution_based import RocketRegressor

reg = RocketRegressor()
reg.fit(X_train, y_train_continuous)
predictions = reg.predict(X_test)
```

### 3. Forecasting
Predict future values using statistical and deep learning models:
- Statistical: ARIMA, ETS, Theta, TAR, AutoTAR, TVP
- Naive baselines: NaiveForecaster with seasonal strategies
- Deep learning: TCN (Temporal Convolutional Networks)
- Regression-based: RegressionForecaster with sliding windows

Example:
```python
from aeon.forecasting.naive import NaiveForecaster

forecaster = NaiveForecaster(strategy="last")
forecaster.fit(y_train)
y_pred = forecaster.predict(fh=[1, 2, 3])  # forecast 3 steps ahead
```

### 4. Anomaly Detection
Identify outliers in time series data:
- **Distance-based**: KMeansAD, CBLOF, LOF, STOMP, LeftSTAMPi, MERLIN, ROCKAD
- **Distribution-based**: COPOD, DWT_MLEAD
- **Outlier detection**: IsolationForest, OneClassSVM, STRAY
- **Collection adapters**: ClassificationAdapter, OutlierDetectionAdapter

Example:
```python
from aeon.anomaly_detection import STOMP

detector = STOMP(window_size=50)
anomaly_scores = detector.fit_predict(X_series)
```

### 5. Clustering
Group similar time series without labels:
```python
from aeon.clustering import TimeSeriesKMeans

clusterer = TimeSeriesKMeans(n_clusters=3, distance="dtw")
clusterer.fit(X_collection)
labels = clusterer.predict(X_new)
```

### 6. Segmentation
Divide time series into distinct regions or identify change points:
```python
from aeon.segmentation import ClaSPSegmenter

segmenter = ClaSPSegmenter()
change_points = segmenter.fit_predict(X_series)
```

### 7. Similarity Search
Find motifs and nearest neighbors in time series collections using specialized distance metrics and matrix profile techniques.

### 8. Transformations
Preprocess and extract features from time series:
- **Collection transformers**: ROCKET, Catch22, TSFresh, Shapelet, SAX, PAA, SFA
- **Series transformers**: Moving Average, Box-Cox, PCA, Fourier, Savitzky-Golay
- **Channel operations**: Selection, scoring, balancing
- **Data balancing**: SMOTE, ADASYN

Example:
```python
from aeon.transformations.collection.convolution_based import Rocket

rocket = Rocket(num_kernels=10000)
X_transformed = rocket.fit_transform(X_train)
```

### 9. Distance Metrics
Compute specialized time series distances:
- **Warping**: DTW, WDTW, DDTW, WDDTW, Shape DTW, ADTW
- **Edit distances**: ERP, EDR, LCSS, TWE
- **Standard**: Euclidean, Manhattan, Minkowski, Squared
- **Specialized**: MSM, SBD

Example:
```python
from aeon.distances import dtw_distance, pairwise_distance

dist = dtw_distance(series1, series2)
dist_matrix = pairwise_distance(X_collection, metric="dtw")
```

## Installation

Install aeon using pip:
```bash
# Core dependencies only
pip install -U aeon

# All optional dependencies
pip install -U "aeon[all_extras]"
```

Or using conda:
```bash
conda create -n aeon-env -c conda-forge aeon
conda activate aeon-env
```

**Requirements**: Python 3.9, 3.10, 3.11, or 3.12

## Data Format

Aeon uses standardized data shapes:
- **Collections**: `(n_cases, n_channels, n_timepoints)` as NumPy arrays or pandas DataFrames
- **Single series**: NumPy arrays or pandas Series
- **Variable-length**: Supported with padding or specialized handling

Load example datasets:
```python
from aeon.datasets import load_arrow_head, load_airline

# Classification dataset
X_train, y_train = load_arrow_head(split="train")

# Forecasting dataset
y = load_airline()
```

## Workflow Patterns

### Pipeline Construction
Combine transformers and estimators using scikit-learn pipelines:
```python
from sklearn.pipeline import Pipeline
from aeon.transformations.collection import Catch22
from aeon.classification.distance_based import KNeighborsTimeSeriesClassifier

pipeline = Pipeline([
    ('features', Catch22()),
    ('classifier', KNeighborsTimeSeriesClassifier())
])
pipeline.fit(X_train, y_train)
```

### Discovery and Tags
Find estimators programmatically:
```python
from aeon.utils.discovery import all_estimators

# Find all classifiers
classifiers = all_estimators(type_filter="classifier")

# Find all forecasters
forecasters = all_estimators(type_filter="forecaster")
```

## References

The skill includes modular reference files with comprehensive details:

### references/learning_tasks.md
In-depth coverage of classification, regression, clustering, and similarity search, including algorithm categories, use cases, and code patterns.

### references/temporal_analysis.md
Detailed information on forecasting, anomaly detection, and segmentation tasks with model descriptions and workflows.

### references/core_modules.md
Comprehensive documentation of transformations, distances, networks, datasets, and benchmarking utilities.

### references/workflows.md
Common workflow patterns, pipeline examples, cross-validation strategies, and integration with scikit-learn.

Load these reference files as needed for detailed information on specific modules or workflows.
