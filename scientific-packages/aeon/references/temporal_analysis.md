# Temporal Analysis: Forecasting, Anomaly Detection, and Segmentation

This reference provides comprehensive details on forecasting future values, detecting anomalies, and segmenting time series.

## Forecasting

Forecasting predicts future values in a time series based on historical patterns.

### Forecasting Concepts

**Forecasting horizon (fh)**: Number of steps ahead to predict
- Absolute: `fh=[1, 2, 3]` (predict steps 1, 2, 3)
- Relative: `fh=ForecastingHorizon([1, 2, 3], is_relative=True)`

**Exogenous variables**: External features that influence predictions

### Statistical Forecasters

#### ARIMA (AutoRegressive Integrated Moving Average)
Classical time series model combining AR, differencing, and MA components:

```python
from aeon.forecasting.arima import ARIMA

forecaster = ARIMA(
    order=(1, 1, 1),  # (p, d, q)
    seasonal_order=(1, 1, 1, 12),  # (P, D, Q, s)
    suppress_warnings=True
)
forecaster.fit(y_train)
y_pred = forecaster.predict(fh=[1, 2, 3, 4, 5])
```

**Parameters:**
- `p`: AR order (lags)
- `d`: Differencing order
- `q`: MA order (moving average)
- `P, D, Q, s`: Seasonal components

#### ETS (Exponential Smoothing)
State space model for trend and seasonality:

```python
from aeon.forecasting.ets import ETS

forecaster = ETS(
    error="add",
    trend="add",
    seasonal="add",
    sp=12  # seasonal period
)
forecaster.fit(y_train)
y_pred = forecaster.predict(fh=[1, 2, 3])
```

**Model types:**
- Error: "add" (additive) or "mul" (multiplicative)
- Trend: "add", "mul", or None
- Seasonal: "add", "mul", or None

#### Theta Forecaster
Simple, effective method using exponential smoothing:

```python
from aeon.forecasting.theta import ThetaForecaster

forecaster = ThetaForecaster(deseasonalize=True, sp=12)
forecaster.fit(y_train)
y_pred = forecaster.predict(fh=np.arange(1, 13))
```

#### TAR (Threshold AutoRegressive)
Non-linear autoregressive model with regime switching:

```python
from aeon.forecasting.tar import TAR

forecaster = TAR(
    delay=1,
    threshold=0.0,
    order_below=2,
    order_above=2
)
forecaster.fit(y_train)
y_pred = forecaster.predict(fh=[1, 2, 3])
```

**AutoTAR**: Automatically optimizes threshold:
```python
from aeon.forecasting.tar import AutoTAR

forecaster = AutoTAR(max_order=5)
forecaster.fit(y_train)
```

#### TVP (Time-Varying Parameter)
Kalman filter-based forecaster with dynamic coefficients:

```python
from aeon.forecasting.tvp import TVP

forecaster = TVP(
    order=2,
    use_exog=False
)
forecaster.fit(y_train)
y_pred = forecaster.predict(fh=[1, 2, 3])
```

### Naive Baselines

Simple forecasting strategies for benchmarking:

```python
from aeon.forecasting.naive import NaiveForecaster

# Last value
forecaster = NaiveForecaster(strategy="last")
forecaster.fit(y_train)
y_pred = forecaster.predict(fh=[1, 2, 3])

# Seasonal naive (use value from same season last year)
forecaster = NaiveForecaster(strategy="seasonal_last", sp=12)
forecaster.fit(y_train)

# Mean
forecaster = NaiveForecaster(strategy="mean")
forecaster.fit(y_train)

# Drift (linear trend from first to last)
forecaster = NaiveForecaster(strategy="drift")
forecaster.fit(y_train)
```

**Strategies:**
- `"last"`: Repeat last observed value
- `"mean"`: Use mean of training data
- `"seasonal_last"`: Repeat value from previous season
- `"drift"`: Linear extrapolation

### Deep Learning Forecasters

#### TCN (Temporal Convolutional Network)
Deep learning with dilated causal convolutions:

```python
from aeon.forecasting.deep_learning import TCNForecaster

forecaster = TCNForecaster(
    n_epochs=100,
    batch_size=32,
    kernel_size=3,
    n_filters=64,
    dilation_rate=2
)
forecaster.fit(y_train)
y_pred = forecaster.predict(fh=[1, 2, 3, 4, 5])
```

### Regression-Based Forecasting

Transform forecasting into a supervised learning problem:

```python
from aeon.forecasting.compose import RegressionForecaster
from sklearn.ensemble import RandomForestRegressor

forecaster = RegressionForecaster(
    regressor=RandomForestRegressor(n_estimators=100),
    window_length=10
)
forecaster.fit(y_train)
y_pred = forecaster.predict(fh=[1, 2, 3])
```

### Forecasting Workflow

```python
from aeon.forecasting.arima import ARIMA
from aeon.datasets import load_airline
from sklearn.metrics import mean_absolute_error, mean_squared_error
import numpy as np

# Load data
y = load_airline()

# Split train/test
split_point = int(len(y) * 0.8)
y_train, y_test = y[:split_point], y[split_point:]

# Fit forecaster
forecaster = ARIMA(order=(2, 1, 2), suppress_warnings=True)
forecaster.fit(y_train)

# Predict
fh = np.arange(1, len(y_test) + 1)
y_pred = forecaster.predict(fh=fh)

# Evaluate
mae = mean_absolute_error(y_test, y_pred)
rmse = np.sqrt(mean_squared_error(y_test, y_pred))
print(f"MAE: {mae:.3f}, RMSE: {rmse:.3f}")
```

### Forecasting with Exogenous Variables

```python
from aeon.forecasting.arima import ARIMA

# X contains exogenous features
forecaster = ARIMA(order=(1, 1, 1))
forecaster.fit(y_train, X=X_train)

# Must provide future exogenous values
y_pred = forecaster.predict(fh=[1, 2, 3], X=X_future)
```

### Multi-Step Forecasting Strategies

**Direct**: Train separate model for each horizon
**Recursive**: Use predictions as inputs for next step
**DirRec**: Combine both strategies

```python
from aeon.forecasting.compose import DirectReductionForecaster
from sklearn.linear_model import Ridge

forecaster = DirectReductionForecaster(
    regressor=Ridge(),
    window_length=10
)
forecaster.fit(y_train)
y_pred = forecaster.predict(fh=[1, 2, 3, 4, 5])
```

## Anomaly Detection

Anomaly detection identifies unusual patterns or outliers in time series data.

### Anomaly Detection Types

**Point anomalies**: Single unusual values
**Contextual anomalies**: Values anomalous given context
**Collective anomalies**: Sequences of unusual behavior

### Distance-Based Anomaly Detectors

#### STOMP (Scalable Time series Ordered-search Matrix Profile)
Matrix profile-based anomaly detection:

```python
from aeon.anomaly_detection import STOMP

detector = STOMP(window_size=50)
anomaly_scores = detector.fit_predict(X_series)

# High scores indicate anomalies
threshold = np.percentile(anomaly_scores, 95)
anomalies = anomaly_scores > threshold
```

#### LeftSTAMPi
Incremental matrix profile for streaming data:

```python
from aeon.anomaly_detection import LeftSTAMPi

detector = LeftSTAMPi(window_size=50)
anomaly_scores = detector.fit_predict(X_series)
```

#### MERLIN
Matrix profile with range constraints:

```python
from aeon.anomaly_detection import MERLIN

detector = MERLIN(window_size=50, k=3)
anomaly_scores = detector.fit_predict(X_series)
```

#### KMeansAD
K-means clustering-based anomaly detection:

```python
from aeon.anomaly_detection import KMeansAD

detector = KMeansAD(n_clusters=5, window_size=50)
anomaly_scores = detector.fit_predict(X_series)
```

#### CBLOF (Cluster-Based Local Outlier Factor)
```python
from aeon.anomaly_detection import CBLOF

detector = CBLOF(n_clusters=8, alpha=0.9)
anomaly_scores = detector.fit_predict(X_series)
```

#### LOF (Local Outlier Factor)
Density-based outlier detection:

```python
from aeon.anomaly_detection import LOF

detector = LOF(n_neighbors=20, window_size=50)
anomaly_scores = detector.fit_predict(X_series)
```

#### ROCKAD
ROCKET-based anomaly detection:

```python
from aeon.anomaly_detection import ROCKAD

detector = ROCKAD(num_kernels=1000, window_size=50)
anomaly_scores = detector.fit_predict(X_series)
```

### Distribution-Based Anomaly Detectors

#### COPOD (Copula-Based Outlier Detection)
```python
from aeon.anomaly_detection import COPOD

detector = COPOD(window_size=50)
anomaly_scores = detector.fit_predict(X_series)
```

#### DWT_MLEAD
Discrete Wavelet Transform with Machine Learning:

```python
from aeon.anomaly_detection import DWT_MLEAD

detector = DWT_MLEAD(window_size=50, wavelet='db4')
anomaly_scores = detector.fit_predict(X_series)
```

### Outlier Detection Methods

#### IsolationForest
Ensemble tree-based isolation:

```python
from aeon.anomaly_detection import IsolationForest

detector = IsolationForest(
    n_estimators=100,
    window_size=50,
    contamination=0.1
)
anomaly_scores = detector.fit_predict(X_series)
```

#### OneClassSVM
Support vector machine for novelty detection:

```python
from aeon.anomaly_detection import OneClassSVM

detector = OneClassSVM(
    kernel='rbf',
    nu=0.1,
    window_size=50
)
anomaly_scores = detector.fit_predict(X_series)
```

#### STRAY (Search TRace AnomalY)
```python
from aeon.anomaly_detection import STRAY

detector = STRAY(alpha=0.05)
anomaly_scores = detector.fit_predict(X_series)
```

### Collection Anomaly Detection

Detect anomalous time series within a collection:

```python
from aeon.anomaly_detection import ClassificationAdapter
from aeon.classification.convolution_based import RocketClassifier

detector = ClassificationAdapter(
    classifier=RocketClassifier()
)
detector.fit(X_normal)  # Train on normal data
anomaly_labels = detector.predict(X_test)  # 1 = anomaly, 0 = normal
```

### Anomaly Detection Workflow

```python
from aeon.anomaly_detection import STOMP
import numpy as np
import matplotlib.pyplot as plt

# Detect anomalies
detector = STOMP(window_size=100)
anomaly_scores = detector.fit_predict(X_series)

# Identify anomalies (top 5%)
threshold = np.percentile(anomaly_scores, 95)
anomaly_indices = np.where(anomaly_scores > threshold)[0]

# Visualize
plt.figure(figsize=(12, 6))
plt.subplot(2, 1, 1)
plt.plot(X_series[0, 0, :])
plt.scatter(anomaly_indices, X_series[0, 0, anomaly_indices],
            color='red', label='Anomalies', zorder=5)
plt.legend()
plt.title('Time Series with Detected Anomalies')

plt.subplot(2, 1, 2)
plt.plot(anomaly_scores)
plt.axhline(threshold, color='red', linestyle='--', label='Threshold')
plt.legend()
plt.title('Anomaly Scores')
plt.tight_layout()
plt.show()
```

## Segmentation

Segmentation divides time series into distinct regions or identifies change points.

### Segmentation Concepts

**Change points**: Locations where statistical properties change
**Segments**: Homogeneous regions between change points
**Applications**: Regime detection, event identification, structural breaks

### Segmentation Algorithms

#### ClaSP (Classification Score Profile)
Discover change points using classification performance:

```python
from aeon.segmentation import ClaSPSegmenter

segmenter = ClaSPSegmenter(
    n_segments=3,
    period_length=10
)
change_points = segmenter.fit_predict(X_series)
print(f"Change points at indices: {change_points}")
```

**How it works:**
- Slides a window over the series
- Computes classification score for left vs. right segments
- High scores indicate change points

#### FLUSS (Fast Low-cost Unipotent Semantic Segmentation)
Matrix profile-based segmentation:

```python
from aeon.segmentation import FLUSSSegmenter

segmenter = FLUSSSegmenter(
    n_segments=5,
    window_size=50
)
change_points = segmenter.fit_predict(X_series)
```

#### BinSeg (Binary Segmentation)
Recursive splitting for change point detection:

```python
from aeon.segmentation import BinSegSegmenter

segmenter = BinSegSegmenter(
    n_segments=4,
    model="l2"  # cost function
)
change_points = segmenter.fit_predict(X_series)
```

**Models:**
- `"l2"`: Least squares (continuous data)
- `"l1"`: Absolute deviation (robust to outliers)
- `"rbf"`: Radial basis function
- `"ar"`: Autoregressive model

#### HMM (Hidden Markov Model) Segmentation
Probabilistic state-based segmentation:

```python
from aeon.segmentation import HMMSegmenter

segmenter = HMMSegmenter(
    n_states=3,
    covariance_type="full"
)
segmenter.fit(X_series)
states = segmenter.predict(X_series)
```

### Segmentation Workflow

```python
from aeon.segmentation import ClaSPSegmenter
import matplotlib.pyplot as plt

# Detect change points
segmenter = ClaSPSegmenter(n_segments=4)
change_points = segmenter.fit_predict(X_series)

# Visualize segments
plt.figure(figsize=(12, 4))
plt.plot(X_series[0, 0, :])
for cp in change_points:
    plt.axvline(cp, color='red', linestyle='--', alpha=0.7)
plt.title('Time Series Segmentation')
plt.xlabel('Time')
plt.ylabel('Value')
plt.show()

# Extract segments
segments = []
prev_cp = 0
for cp in np.append(change_points, len(X_series[0, 0, :])):
    segment = X_series[0, 0, prev_cp:cp]
    segments.append(segment)
    prev_cp = cp
```

### Multi-Variate Segmentation

```python
from aeon.segmentation import ClaSPSegmenter

# X_multivariate has shape (1, n_channels, n_timepoints)
segmenter = ClaSPSegmenter(n_segments=3)
change_points = segmenter.fit_predict(X_multivariate)
```

## Combining Forecasting, Anomaly Detection, and Segmentation

### Robust Forecasting with Anomaly Detection

```python
from aeon.forecasting.arima import ARIMA
from aeon.anomaly_detection import IsolationForest

# Detect and remove anomalies
detector = IsolationForest(window_size=50, contamination=0.1)
anomaly_scores = detector.fit_predict(X_series)
normal_mask = anomaly_scores < np.percentile(anomaly_scores, 90)

# Forecast on cleaned data
y_clean = y_train[normal_mask]
forecaster = ARIMA(order=(2, 1, 2))
forecaster.fit(y_clean)
y_pred = forecaster.predict(fh=[1, 2, 3])
```

### Segmentation-Based Forecasting

```python
from aeon.segmentation import ClaSPSegmenter
from aeon.forecasting.arima import ARIMA

# Segment time series
segmenter = ClaSPSegmenter(n_segments=3)
change_points = segmenter.fit_predict(X_series)

# Forecast using most recent segment
last_segment_start = change_points[-1]
y_recent = y_train[last_segment_start:]

forecaster = ARIMA(order=(1, 1, 1))
forecaster.fit(y_recent)
y_pred = forecaster.predict(fh=[1, 2, 3])
```

## Discovery Functions

Find available forecasters, detectors, and segmenters:

```python
from aeon.utils.discovery import all_estimators

# Get all forecasters
forecasters = all_estimators(type_filter="forecaster")

# Get all anomaly detectors
detectors = all_estimators(type_filter="anomaly-detector")

# Get all segmenters
segmenters = all_estimators(type_filter="segmenter")
```
