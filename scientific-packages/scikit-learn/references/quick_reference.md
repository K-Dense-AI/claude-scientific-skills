# Scikit-learn Quick Reference

## Essential Imports

```python
# Core
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV
from sklearn.pipeline import Pipeline, make_pipeline
from sklearn.compose import ColumnTransformer

# Preprocessing
from sklearn.preprocessing import (
    StandardScaler, MinMaxScaler, RobustScaler,
    OneHotEncoder, OrdinalEncoder, LabelEncoder,
    PolynomialFeatures
)
from sklearn.impute import SimpleImputer

# Models - Classification
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import (
    RandomForestClassifier,
    GradientBoostingClassifier,
    HistGradientBoostingClassifier
)
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier

# Models - Regression
from sklearn.linear_model import LinearRegression, Ridge, Lasso
from sklearn.ensemble import (
    RandomForestRegressor,
    GradientBoostingRegressor,
    HistGradientBoostingRegressor
)

# Clustering
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering
from sklearn.mixture import GaussianMixture

# Dimensionality Reduction
from sklearn.decomposition import PCA, NMF, TruncatedSVD
from sklearn.manifold import TSNE

# Metrics
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    confusion_matrix, classification_report,
    mean_squared_error, r2_score, mean_absolute_error
)
```

## Basic Workflow Template

### Classification

```python
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report

# Split data
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y
)

# Scale features
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Train model
model = RandomForestClassifier(n_estimators=100, random_state=42)
model.fit(X_train_scaled, y_train)

# Predict and evaluate
y_pred = model.predict(X_test_scaled)
print(classification_report(y_test, y_pred))
```

### Regression

```python
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score

# Split data
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

# Scale features
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Train model
model = RandomForestRegressor(n_estimators=100, random_state=42)
model.fit(X_train_scaled, y_train)

# Predict and evaluate
y_pred = model.predict(X_test_scaled)
print(f"RMSE: {mean_squared_error(y_test, y_pred, squared=False):.3f}")
print(f"R²: {r2_score(y_test, y_pred):.3f}")
```

### With Pipeline (Recommended)

```python
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, cross_val_score

# Create pipeline
pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('classifier', RandomForestClassifier(n_estimators=100, random_state=42))
])

# Split and train
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)
pipeline.fit(X_train, y_train)

# Evaluate
score = pipeline.score(X_test, y_test)
cv_scores = cross_val_score(pipeline, X_train, y_train, cv=5)
print(f"Test accuracy: {score:.3f}")
print(f"CV accuracy: {cv_scores.mean():.3f} (+/- {cv_scores.std():.3f})")
```

## Common Preprocessing Patterns

### Numeric Data

```python
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline

numeric_transformer = Pipeline([
    ('imputer', SimpleImputer(strategy='median')),
    ('scaler', StandardScaler())
])
```

### Categorical Data

```python
from sklearn.preprocessing import OneHotEncoder
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline

categorical_transformer = Pipeline([
    ('imputer', SimpleImputer(strategy='constant', fill_value='missing')),
    ('onehot', OneHotEncoder(handle_unknown='ignore'))
])
```

### Mixed Data with ColumnTransformer

```python
from sklearn.compose import ColumnTransformer

numeric_features = ['age', 'income', 'credit_score']
categorical_features = ['country', 'occupation']

preprocessor = ColumnTransformer(
    transformers=[
        ('num', numeric_transformer, numeric_features),
        ('cat', categorical_transformer, categorical_features)
    ])

# Complete pipeline
from sklearn.ensemble import RandomForestClassifier
pipeline = Pipeline([
    ('preprocessor', preprocessor),
    ('classifier', RandomForestClassifier())
])
```

## Model Selection Cheat Sheet

### Quick Decision Tree

```
Is it supervised?
├─ Yes
│  ├─ Predicting categories? → Classification
│  │  ├─ Start with: LogisticRegression (baseline)
│  │  ├─ Then try: RandomForestClassifier
│  │  └─ Best performance: HistGradientBoostingClassifier
│  └─ Predicting numbers? → Regression
│     ├─ Start with: LinearRegression/Ridge (baseline)
│     ├─ Then try: RandomForestRegressor
│     └─ Best performance: HistGradientBoostingRegressor
└─ No
   ├─ Grouping similar items? → Clustering
   │  ├─ Know # clusters: KMeans
   │  └─ Unknown # clusters: DBSCAN or HDBSCAN
   ├─ Reducing dimensions?
   │  ├─ For preprocessing: PCA
   │  └─ For visualization: t-SNE or UMAP
   └─ Finding outliers? → IsolationForest or LocalOutlierFactor
```

### Algorithm Selection by Data Size

- **Small (<1K samples)**: Any algorithm
- **Medium (1K-100K)**: Random Forests, Gradient Boosting, Neural Networks
- **Large (>100K)**: SGDClassifier/Regressor, HistGradientBoosting, LinearSVC

### When to Scale Features

**Always scale**:
- SVM, Neural Networks
- K-Nearest Neighbors
- Linear/Logistic Regression (with regularization)
- PCA, LDA
- Any gradient descent algorithm

**Don't need to scale**:
- Tree-based (Decision Trees, Random Forests, Gradient Boosting)
- Naive Bayes

## Hyperparameter Tuning

### GridSearchCV

```python
from sklearn.model_selection import GridSearchCV

param_grid = {
    'n_estimators': [100, 200, 500],
    'max_depth': [10, 20, None],
    'min_samples_split': [2, 5, 10]
}

grid_search = GridSearchCV(
    RandomForestClassifier(random_state=42),
    param_grid,
    cv=5,
    scoring='f1_weighted',
    n_jobs=-1
)

grid_search.fit(X_train, y_train)
best_model = grid_search.best_estimator_
print(f"Best params: {grid_search.best_params_}")
```

### RandomizedSearchCV (Faster)

```python
from sklearn.model_selection import RandomizedSearchCV
from scipy.stats import randint, uniform

param_distributions = {
    'n_estimators': randint(100, 1000),
    'max_depth': randint(5, 50),
    'min_samples_split': randint(2, 20)
}

random_search = RandomizedSearchCV(
    RandomForestClassifier(random_state=42),
    param_distributions,
    n_iter=50,  # Number of combinations to try
    cv=5,
    n_jobs=-1,
    random_state=42
)

random_search.fit(X_train, y_train)
```

### Pipeline with GridSearchCV

```python
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV

pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('svm', SVC())
])

param_grid = {
    'svm__C': [0.1, 1, 10],
    'svm__kernel': ['rbf', 'linear'],
    'svm__gamma': ['scale', 'auto']
}

grid = GridSearchCV(pipeline, param_grid, cv=5)
grid.fit(X_train, y_train)
```

## Cross-Validation

### Basic Cross-Validation

```python
from sklearn.model_selection import cross_val_score

scores = cross_val_score(model, X, y, cv=5, scoring='accuracy')
print(f"Accuracy: {scores.mean():.3f} (+/- {scores.std():.3f})")
```

### Multiple Metrics

```python
from sklearn.model_selection import cross_validate

scoring = ['accuracy', 'precision_weighted', 'recall_weighted', 'f1_weighted']
results = cross_validate(model, X, y, cv=5, scoring=scoring)

for metric in scoring:
    scores = results[f'test_{metric}']
    print(f"{metric}: {scores.mean():.3f} (+/- {scores.std():.3f})")
```

### Custom CV Strategies

```python
from sklearn.model_selection import StratifiedKFold, TimeSeriesSplit

# For imbalanced classification
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# For time series
cv = TimeSeriesSplit(n_splits=5)

scores = cross_val_score(model, X, y, cv=cv)
```

## Common Metrics

### Classification

```python
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score,
    precision_score, recall_score, f1_score,
    confusion_matrix, classification_report,
    roc_auc_score
)

# Basic metrics
accuracy = accuracy_score(y_true, y_pred)
f1 = f1_score(y_true, y_pred, average='weighted')

# Comprehensive report
print(classification_report(y_true, y_pred))

# ROC AUC (requires probabilities)
y_proba = model.predict_proba(X_test)[:, 1]
auc = roc_auc_score(y_true, y_proba)
```

### Regression

```python
from sklearn.metrics import (
    mean_squared_error,
    mean_absolute_error,
    r2_score
)

mse = mean_squared_error(y_true, y_pred)
rmse = mean_squared_error(y_true, y_pred, squared=False)
mae = mean_absolute_error(y_true, y_pred)
r2 = r2_score(y_true, y_pred)

print(f"RMSE: {rmse:.3f}")
print(f"MAE: {mae:.3f}")
print(f"R²: {r2:.3f}")
```

## Feature Engineering

### Polynomial Features

```python
from sklearn.preprocessing import PolynomialFeatures

poly = PolynomialFeatures(degree=2, include_bias=False)
X_poly = poly.fit_transform(X)
# [x1, x2] → [x1, x2, x1², x1·x2, x2²]
```

### Feature Selection

```python
from sklearn.feature_selection import (
    SelectKBest, f_classif,
    RFE,
    SelectFromModel
)

# Univariate selection
selector = SelectKBest(f_classif, k=10)
X_selected = selector.fit_transform(X, y)

# Recursive feature elimination
from sklearn.ensemble import RandomForestClassifier
rfe = RFE(RandomForestClassifier(), n_features_to_select=10)
X_selected = rfe.fit_transform(X, y)

# Model-based selection
selector = SelectFromModel(
    RandomForestClassifier(n_estimators=100),
    threshold='median'
)
X_selected = selector.fit_transform(X, y)
```

### Feature Importance

```python
# Tree-based models
model = RandomForestClassifier()
model.fit(X_train, y_train)
importances = model.feature_importances_

# Visualize
import matplotlib.pyplot as plt
indices = np.argsort(importances)[::-1]
plt.bar(range(X.shape[1]), importances[indices])
plt.xticks(range(X.shape[1]), feature_names[indices], rotation=90)
plt.show()

# Permutation importance (works for any model)
from sklearn.inspection import permutation_importance
result = permutation_importance(model, X_test, y_test, n_repeats=10)
importances = result.importances_mean
```

## Clustering

### K-Means

```python
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

# Always scale for k-means
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Fit k-means
kmeans = KMeans(n_clusters=3, random_state=42)
labels = kmeans.fit_predict(X_scaled)

# Evaluate
from sklearn.metrics import silhouette_score
score = silhouette_score(X_scaled, labels)
print(f"Silhouette score: {score:.3f}")
```

### Elbow Method

```python
inertias = []
K_range = range(2, 11)

for k in K_range:
    kmeans = KMeans(n_clusters=k, random_state=42)
    kmeans.fit(X_scaled)
    inertias.append(kmeans.inertia_)

plt.plot(K_range, inertias, 'bo-')
plt.xlabel('k')
plt.ylabel('Inertia')
plt.show()
```

### DBSCAN

```python
from sklearn.cluster import DBSCAN

dbscan = DBSCAN(eps=0.5, min_samples=5)
labels = dbscan.fit_predict(X_scaled)

# -1 indicates noise/outliers
n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
n_noise = list(labels).count(-1)
print(f"Clusters: {n_clusters}, Noise points: {n_noise}")
```

## Dimensionality Reduction

### PCA

```python
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Always scale before PCA
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Specify n_components
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

# Or specify variance to retain
pca = PCA(n_components=0.95)  # Keep 95% variance
X_pca = pca.fit_transform(X_scaled)

print(f"Explained variance: {pca.explained_variance_ratio_}")
print(f"Components needed: {pca.n_components_}")
```

### t-SNE (Visualization Only)

```python
from sklearn.manifold import TSNE

# Reduce to 50 dimensions with PCA first (recommended)
pca = PCA(n_components=50)
X_pca = pca.fit_transform(X_scaled)

# Apply t-SNE
tsne = TSNE(n_components=2, random_state=42, perplexity=30)
X_tsne = tsne.fit_transform(X_pca)

# Visualize
plt.scatter(X_tsne[:, 0], X_tsne[:, 1], c=y, cmap='viridis')
plt.colorbar()
plt.show()
```

## Saving and Loading Models

```python
import joblib

# Save model
joblib.dump(model, 'model.pkl')

# Save pipeline
joblib.dump(pipeline, 'pipeline.pkl')

# Load
model = joblib.load('model.pkl')
pipeline = joblib.load('pipeline.pkl')

# Use loaded model
y_pred = model.predict(X_new)
```

## Common Pitfalls and Solutions

### Data Leakage
❌ **Wrong**: Fit on all data before split
```python
scaler = StandardScaler().fit(X)
X_train, X_test = train_test_split(scaler.transform(X))
```

✅ **Correct**: Use pipeline or fit only on train
```python
X_train, X_test = train_test_split(X)
pipeline = Pipeline([('scaler', StandardScaler()), ('model', model)])
pipeline.fit(X_train, y_train)
```

### Not Scaling
❌ **Wrong**: Using SVM without scaling
```python
svm = SVC()
svm.fit(X_train, y_train)
```

✅ **Correct**: Scale for SVM
```python
pipeline = Pipeline([('scaler', StandardScaler()), ('svm', SVC())])
pipeline.fit(X_train, y_train)
```

### Wrong Metric for Imbalanced Data
❌ **Wrong**: Using accuracy for 99:1 imbalance
```python
accuracy = accuracy_score(y_true, y_pred)  # Can be misleading
```

✅ **Correct**: Use appropriate metrics
```python
f1 = f1_score(y_true, y_pred, average='weighted')
balanced_acc = balanced_accuracy_score(y_true, y_pred)
```

### Not Using Stratification
❌ **Wrong**: Random split for imbalanced data
```python
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
```

✅ **Correct**: Stratify for imbalanced classes
```python
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, stratify=y
)
```

## Performance Tips

1. **Use n_jobs=-1** for parallel processing (RandomForest, GridSearchCV)
2. **Use HistGradientBoosting** for large datasets (>10K samples)
3. **Use MiniBatchKMeans** for large clustering tasks
4. **Use IncrementalPCA** for data that doesn't fit in memory
5. **Use sparse matrices** for high-dimensional sparse data (text)
6. **Cache transformers** in pipelines during grid search
7. **Use RandomizedSearchCV** instead of GridSearchCV for large parameter spaces
8. **Reduce dimensionality** with PCA before applying expensive algorithms
