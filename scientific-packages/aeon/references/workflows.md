# Common Workflows and Integration Patterns

This reference provides end-to-end workflows, best practices, and integration patterns for using aeon effectively.

## Complete Classification Workflow

### Basic Classification Pipeline

```python
# 1. Import required modules
from aeon.classification.convolution_based import RocketClassifier
from aeon.datasets import load_arrow_head
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns

# 2. Load and inspect data
X_train, y_train = load_arrow_head(split="train")
X_test, y_test = load_arrow_head(split="test")

print(f"Training shape: {X_train.shape}")  # (n_cases, n_channels, n_timepoints)
print(f"Unique classes: {np.unique(y_train)}")

# 3. Train classifier
clf = RocketClassifier(num_kernels=10000, n_jobs=-1)
clf.fit(X_train, y_train)

# 4. Make predictions
y_pred = clf.predict(X_test)
y_proba = clf.predict_proba(X_test)

# 5. Evaluate performance
accuracy = accuracy_score(y_test, y_pred)
print(f"Accuracy: {accuracy:.3f}")
print("\nClassification Report:")
print(classification_report(y_test, y_pred))

# 6. Visualize confusion matrix
cm = confusion_matrix(y_test, y_pred)
plt.figure(figsize=(8, 6))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
plt.title('Confusion Matrix')
plt.ylabel('True Label')
plt.xlabel('Predicted Label')
plt.show()
```

### Feature Extraction + Classifier Pipeline

```python
from sklearn.pipeline import Pipeline
from aeon.transformations.collection import Catch22
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score

# Create pipeline
pipeline = Pipeline([
    ('features', Catch22(n_jobs=-1)),
    ('classifier', RandomForestClassifier(n_estimators=500, n_jobs=-1))
])

# Cross-validation
scores = cross_val_score(pipeline, X_train, y_train, cv=5, scoring='accuracy')
print(f"CV Accuracy: {scores.mean():.3f} (+/- {scores.std():.3f})")

# Train on full training set
pipeline.fit(X_train, y_train)

# Evaluate on test set
accuracy = pipeline.score(X_test, y_test)
print(f"Test Accuracy: {accuracy:.3f}")
```

### Multi-Algorithm Comparison

```python
from aeon.classification.convolution_based import RocketClassifier, MiniRocketClassifier
from aeon.classification.distance_based import KNeighborsTimeSeriesClassifier
from aeon.classification.feature_based import Catch22Classifier
from aeon.classification.interval_based import TimeSeriesForestClassifier
import time

classifiers = {
    'ROCKET': RocketClassifier(num_kernels=10000),
    'MiniRocket': MiniRocketClassifier(),
    'KNN-DTW': KNeighborsTimeSeriesClassifier(distance='dtw', n_neighbors=5),
    'Catch22': Catch22Classifier(),
    'TSF': TimeSeriesForestClassifier(n_estimators=200)
}

results = {}
for name, clf in classifiers.items():
    start_time = time.time()
    clf.fit(X_train, y_train)
    train_time = time.time() - start_time

    start_time = time.time()
    accuracy = clf.score(X_test, y_test)
    test_time = time.time() - start_time

    results[name] = {
        'accuracy': accuracy,
        'train_time': train_time,
        'test_time': test_time
    }

# Display results
import pandas as pd
df_results = pd.DataFrame(results).T
df_results = df_results.sort_values('accuracy', ascending=False)
print(df_results)
```

## Complete Forecasting Workflow

### Univariate Forecasting

```python
from aeon.forecasting.arima import ARIMA
from aeon.forecasting.naive import NaiveForecaster
from aeon.datasets import load_airline
from sklearn.metrics import mean_absolute_error, mean_squared_error
import numpy as np
import matplotlib.pyplot as plt

# 1. Load data
y = load_airline()

# 2. Train/test split (temporal)
split_point = int(len(y) * 0.8)
y_train, y_test = y[:split_point], y[split_point:]

# 3. Create baseline (naive forecaster)
baseline = NaiveForecaster(strategy="last")
baseline.fit(y_train)
y_pred_baseline = baseline.predict(fh=np.arange(1, len(y_test) + 1))

# 4. Train ARIMA model
forecaster = ARIMA(order=(2, 1, 2), seasonal_order=(1, 1, 1, 12))
forecaster.fit(y_train)
y_pred = forecaster.predict(fh=np.arange(1, len(y_test) + 1))

# 5. Evaluate
mae_baseline = mean_absolute_error(y_test, y_pred_baseline)
mae_arima = mean_absolute_error(y_test, y_pred)
rmse_baseline = np.sqrt(mean_squared_error(y_test, y_pred_baseline))
rmse_arima = np.sqrt(mean_squared_error(y_test, y_pred))

print(f"Baseline - MAE: {mae_baseline:.2f}, RMSE: {rmse_baseline:.2f}")
print(f"ARIMA    - MAE: {mae_arima:.2f}, RMSE: {rmse_arima:.2f}")

# 6. Visualize
plt.figure(figsize=(12, 6))
plt.plot(y_train.index, y_train, label='Train', alpha=0.7)
plt.plot(y_test.index, y_test, label='Test (Actual)', alpha=0.7)
plt.plot(y_test.index, y_pred, label='ARIMA Forecast', linestyle='--')
plt.plot(y_test.index, y_pred_baseline, label='Baseline', linestyle=':', alpha=0.5)
plt.legend()
plt.title('Forecasting Results')
plt.xlabel('Time')
plt.ylabel('Value')
plt.show()
```

### Forecast with Confidence Intervals

```python
from aeon.forecasting.arima import ARIMA

forecaster = ARIMA(order=(2, 1, 2))
forecaster.fit(y_train)

# Predict with prediction intervals
y_pred = forecaster.predict(fh=np.arange(1, len(y_test) + 1))
pred_interval = forecaster.predict_interval(
    fh=np.arange(1, len(y_test) + 1),
    coverage=0.95
)

# Visualize with confidence bands
plt.figure(figsize=(12, 6))
plt.plot(y_test.index, y_test, label='Actual')
plt.plot(y_test.index, y_pred, label='Forecast')
plt.fill_between(
    y_test.index,
    pred_interval.iloc[:, 0],
    pred_interval.iloc[:, 1],
    alpha=0.3,
    label='95% Confidence'
)
plt.legend()
plt.show()
```

### Multi-Step Ahead Forecasting

```python
from aeon.forecasting.compose import DirectReductionForecaster
from sklearn.ensemble import GradientBoostingRegressor

# Convert to supervised learning problem
forecaster = DirectReductionForecaster(
    regressor=GradientBoostingRegressor(n_estimators=100),
    window_length=12
)
forecaster.fit(y_train)

# Forecast multiple steps
fh = np.arange(1, 13)  # 12 months ahead
y_pred = forecaster.predict(fh=fh)
```

## Complete Anomaly Detection Workflow

```python
from aeon.anomaly_detection import STOMP
from aeon.datasets import load_airline
import numpy as np
import matplotlib.pyplot as plt

# 1. Load data
y = load_airline()
X_series = y.values.reshape(1, 1, -1)  # Convert to aeon format

# 2. Detect anomalies
detector = STOMP(window_size=50)
anomaly_scores = detector.fit_predict(X_series)

# 3. Identify anomalies (top 5%)
threshold = np.percentile(anomaly_scores, 95)
anomaly_indices = np.where(anomaly_scores > threshold)[0]

# 4. Visualize
fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

# Plot time series with anomalies
axes[0].plot(y.values, label='Time Series')
axes[0].scatter(
    anomaly_indices,
    y.values[anomaly_indices],
    color='red',
    s=100,
    label='Anomalies',
    zorder=5
)
axes[0].set_ylabel('Value')
axes[0].legend()
axes[0].set_title('Time Series with Detected Anomalies')

# Plot anomaly scores
axes[1].plot(anomaly_scores, label='Anomaly Score')
axes[1].axhline(threshold, color='red', linestyle='--', label='Threshold')
axes[1].set_xlabel('Time')
axes[1].set_ylabel('Score')
axes[1].legend()
axes[1].set_title('Anomaly Scores')

plt.tight_layout()
plt.show()

# 5. Extract anomalous segments
print(f"Found {len(anomaly_indices)} anomalies")
for idx in anomaly_indices[:5]:  # Show first 5
    print(f"Anomaly at index {idx}, value: {y.values[idx]:.2f}")
```

## Complete Clustering Workflow

```python
from aeon.clustering import TimeSeriesKMeans
from aeon.datasets import load_basic_motions
from sklearn.metrics import silhouette_score, davies_bouldin_score
import matplotlib.pyplot as plt

# 1. Load data
X_train, y_train = load_basic_motions(split="train")

# 2. Determine optimal number of clusters (elbow method)
inertias = []
silhouettes = []
K = range(2, 11)

for k in K:
    clusterer = TimeSeriesKMeans(n_clusters=k, distance="euclidean", n_init=5)
    labels = clusterer.fit_predict(X_train)
    inertias.append(clusterer.inertia_)
    silhouettes.append(silhouette_score(X_train.reshape(len(X_train), -1), labels))

# Plot elbow curve
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
axes[0].plot(K, inertias, 'bo-')
axes[0].set_xlabel('Number of Clusters')
axes[0].set_ylabel('Inertia')
axes[0].set_title('Elbow Method')

axes[1].plot(K, silhouettes, 'ro-')
axes[1].set_xlabel('Number of Clusters')
axes[1].set_ylabel('Silhouette Score')
axes[1].set_title('Silhouette Analysis')
plt.tight_layout()
plt.show()

# 3. Cluster with optimal k
optimal_k = 4
clusterer = TimeSeriesKMeans(n_clusters=optimal_k, distance="dtw", n_init=10)
labels = clusterer.fit_predict(X_train)

# 4. Visualize clusters
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.ravel()

for cluster_id in range(optimal_k):
    cluster_indices = np.where(labels == cluster_id)[0]
    ax = axes[cluster_id]

    # Plot all series in cluster
    for idx in cluster_indices[:20]:  # Plot up to 20 series
        ax.plot(X_train[idx, 0, :], alpha=0.3, color='blue')

    # Plot cluster center
    ax.plot(clusterer.cluster_centers_[cluster_id, 0, :],
            color='red', linewidth=2, label='Center')
    ax.set_title(f'Cluster {cluster_id} (n={len(cluster_indices)})')
    ax.legend()

plt.tight_layout()
plt.show()
```

## Cross-Validation Strategies

### Standard K-Fold Cross-Validation

```python
from sklearn.model_selection import cross_val_score, StratifiedKFold
from aeon.classification.convolution_based import RocketClassifier

clf = RocketClassifier()

# Stratified K-Fold (preserves class distribution)
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
scores = cross_val_score(clf, X_train, y_train, cv=cv, scoring='accuracy')

print(f"Cross-validation scores: {scores}")
print(f"Mean accuracy: {scores.mean():.3f} (+/- {scores.std():.3f})")
```

### Time Series Cross-Validation (for forecasting)

```python
from sklearn.model_selection import TimeSeriesSplit
from aeon.forecasting.arima import ARIMA
from sklearn.metrics import mean_squared_error
import numpy as np

# Time-aware split (no future data leakage)
tscv = TimeSeriesSplit(n_splits=5)
mse_scores = []

for train_idx, test_idx in tscv.split(y):
    y_train_cv, y_test_cv = y.iloc[train_idx], y.iloc[test_idx]

    forecaster = ARIMA(order=(2, 1, 2))
    forecaster.fit(y_train_cv)

    fh = np.arange(1, len(y_test_cv) + 1)
    y_pred = forecaster.predict(fh=fh)

    mse = mean_squared_error(y_test_cv, y_pred)
    mse_scores.append(mse)

print(f"CV MSE: {np.mean(mse_scores):.3f} (+/- {np.std(mse_scores):.3f})")
```

## Hyperparameter Tuning

### Grid Search

```python
from sklearn.model_selection import GridSearchCV
from aeon.classification.distance_based import KNeighborsTimeSeriesClassifier

# Define parameter grid
param_grid = {
    'n_neighbors': [1, 3, 5, 7, 9],
    'distance': ['dtw', 'euclidean', 'erp', 'msm'],
    'distance_params': [{'window': 0.1}, {'window': 0.2}, None]
}

# Grid search with cross-validation
clf = KNeighborsTimeSeriesClassifier()
grid_search = GridSearchCV(
    clf,
    param_grid,
    cv=5,
    scoring='accuracy',
    n_jobs=-1,
    verbose=2
)

grid_search.fit(X_train, y_train)

print(f"Best parameters: {grid_search.best_params_}")
print(f"Best CV score: {grid_search.best_score_:.3f}")
print(f"Test accuracy: {grid_search.score(X_test, y_test):.3f}")
```

### Random Search

```python
from sklearn.model_selection import RandomizedSearchCV
from scipy.stats import randint, uniform

param_distributions = {
    'n_neighbors': randint(1, 20),
    'distance': ['dtw', 'euclidean', 'ddtw'],
    'distance_params': [{'window': w} for w in np.linspace(0.0, 0.5, 10)]
}

clf = KNeighborsTimeSeriesClassifier()
random_search = RandomizedSearchCV(
    clf,
    param_distributions,
    n_iter=50,
    cv=5,
    scoring='accuracy',
    n_jobs=-1,
    random_state=42
)

random_search.fit(X_train, y_train)
print(f"Best parameters: {random_search.best_params_}")
```

## Integration with scikit-learn

### Using aeon in scikit-learn Pipelines

```python
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from aeon.transformations.collection import Catch22
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.ensemble import RandomForestClassifier

pipeline = Pipeline([
    ('features', Catch22()),
    ('scaler', StandardScaler()),
    ('feature_selection', SelectKBest(f_classif, k=15)),
    ('classifier', RandomForestClassifier(n_estimators=500))
])

pipeline.fit(X_train, y_train)
accuracy = pipeline.score(X_test, y_test)
```

### Voting Ensemble with scikit-learn

```python
from sklearn.ensemble import VotingClassifier
from aeon.classification.convolution_based import RocketClassifier
from aeon.classification.distance_based import KNeighborsTimeSeriesClassifier
from aeon.classification.feature_based import Catch22Classifier

ensemble = VotingClassifier(
    estimators=[
        ('rocket', RocketClassifier()),
        ('knn', KNeighborsTimeSeriesClassifier()),
        ('catch22', Catch22Classifier())
    ],
    voting='soft',
    n_jobs=-1
)

ensemble.fit(X_train, y_train)
accuracy = ensemble.score(X_test, y_test)
```

### Stacking with Meta-Learner

```python
from sklearn.ensemble import StackingClassifier
from sklearn.linear_model import LogisticRegression
from aeon.classification.convolution_based import MiniRocketClassifier
from aeon.classification.interval_based import TimeSeriesForestClassifier

stacking = StackingClassifier(
    estimators=[
        ('minirocket', MiniRocketClassifier()),
        ('tsf', TimeSeriesForestClassifier(n_estimators=100))
    ],
    final_estimator=LogisticRegression(),
    cv=5
)

stacking.fit(X_train, y_train)
accuracy = stacking.score(X_test, y_test)
```

## Data Preprocessing

### Handling Variable-Length Series

```python
from aeon.transformations.collection import PaddingTransformer

# Pad series to equal length
padder = PaddingTransformer(pad_length=None, fill_value=0)
X_padded = padder.fit_transform(X_variable_length)
```

### Handling Missing Values

```python
from aeon.transformations.series import Imputer

imputer = Imputer(method='mean')
X_imputed = imputer.fit_transform(X_with_missing)
```

### Normalization

```python
from aeon.transformations.collection import Normalizer

normalizer = Normalizer(method='z-score')
X_normalized = normalizer.fit_transform(X_train)
```

## Model Persistence

### Saving and Loading Models

```python
import pickle
from aeon.classification.convolution_based import RocketClassifier

# Train and save
clf = RocketClassifier()
clf.fit(X_train, y_train)

with open('rocket_model.pkl', 'wb') as f:
    pickle.dump(clf, f)

# Load and predict
with open('rocket_model.pkl', 'rb') as f:
    loaded_clf = pickle.load(f)

predictions = loaded_clf.predict(X_test)
```

### Using joblib (recommended for large models)

```python
import joblib

# Save
joblib.dump(clf, 'rocket_model.joblib')

# Load
loaded_clf = joblib.load('rocket_model.joblib')
```

## Visualization Utilities

### Plotting Time Series

```python
from aeon.visualisation import plot_series
import matplotlib.pyplot as plt

# Plot multiple series
fig, ax = plt.subplots(figsize=(12, 6))
plot_series(X_train[0], X_train[1], X_train[2], labels=['Series 1', 'Series 2', 'Series 3'], ax=ax)
plt.title('Time Series Visualization')
plt.show()
```

### Plotting Distance Matrices

```python
from aeon.distances import pairwise_distance
import seaborn as sns

dist_matrix = pairwise_distance(X_train[:50], metric="dtw")

plt.figure(figsize=(10, 8))
sns.heatmap(dist_matrix, cmap='viridis', square=True)
plt.title('DTW Distance Matrix')
plt.show()
```

## Performance Optimization Tips

1. **Use n_jobs=-1** for parallel processing:
   ```python
   clf = RocketClassifier(num_kernels=10000, n_jobs=-1)
   ```

2. **Use MiniRocket instead of ROCKET** for faster training:
   ```python
   clf = MiniRocketClassifier()  # 75% faster
   ```

3. **Reduce num_kernels** for faster training:
   ```python
   clf = RocketClassifier(num_kernels=2000)  # Default is 10000
   ```

4. **Use Catch22 instead of TSFresh**:
   ```python
   transform = Catch22()  # Much faster, fewer features
   ```

5. **Window constraints for DTW**:
   ```python
   clf = KNeighborsTimeSeriesClassifier(
       distance='dtw',
       distance_params={'window': 0.1}  # Constrain warping
   )
   ```

## Best Practices

1. **Always use train/test split** with time series ordering preserved
2. **Use stratified splits** for classification to maintain class balance
3. **Start with fast algorithms** (ROCKET, MiniRocket) before trying slow ones
4. **Use cross-validation** to estimate generalization performance
5. **Benchmark against naive baselines** to establish minimum performance
6. **Normalize/standardize** when using distance-based methods
7. **Use appropriate distance metrics** for your data characteristics
8. **Save trained models** to avoid retraining
9. **Monitor training time** and computational resources
10. **Visualize results** to understand model behavior
