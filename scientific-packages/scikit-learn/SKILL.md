---
name: scikit-learn
description: "ML toolkit. Classification, regression, clustering, PCA, preprocessing, pipelines, GridSearch, cross-validation, RandomForest, SVM, for general machine learning workflows."
---

# Scikit-learn: Machine Learning in Python

## Overview

Scikit-learn is Python's premier machine learning library, offering simple and efficient tools for predictive data analysis. Apply this skill for classification, regression, clustering, dimensionality reduction, model selection, preprocessing, and hyperparameter optimization.

## When to Use This Skill

This skill should be used when:
- Building classification models (spam detection, image recognition, medical diagnosis)
- Creating regression models (price prediction, forecasting, trend analysis)
- Performing clustering analysis (customer segmentation, pattern discovery)
- Reducing dimensionality (PCA, t-SNE for visualization)
- Preprocessing data (scaling, encoding, imputation)
- Evaluating model performance (cross-validation, metrics)
- Tuning hyperparameters (grid search, random search)
- Creating machine learning pipelines
- Detecting anomalies or outliers
- Implementing ensemble methods

## Core Machine Learning Workflow

### Standard ML Pipeline

Follow this general workflow for supervised learning tasks:

1. **Data Preparation**
   - Load and explore data
   - Split into train/test sets
   - Handle missing values
   - Encode categorical features
   - Scale/normalize features

2. **Model Selection**
   - Start with baseline model
   - Try more complex models
   - Use domain knowledge to guide selection

3. **Model Training**
   - Fit model on training data
   - Use pipelines to prevent data leakage
   - Apply cross-validation

4. **Model Evaluation**
   - Evaluate on test set
   - Use appropriate metrics
   - Analyze errors

5. **Model Optimization**
   - Tune hyperparameters
   - Feature engineering
   - Ensemble methods

6. **Deployment**
   - Save model using joblib
   - Create prediction pipeline
   - Monitor performance

### Classification Quick Start

```python
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
from sklearn.pipeline import Pipeline

# Create pipeline (prevents data leakage)
pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('classifier', RandomForestClassifier(n_estimators=100, random_state=42))
])

# Split data (use stratify for imbalanced classes)
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y
)

# Train
pipeline.fit(X_train, y_train)

# Evaluate
y_pred = pipeline.predict(X_test)
print(classification_report(y_test, y_pred))

# Cross-validation for robust evaluation
from sklearn.model_selection import cross_val_score
scores = cross_val_score(pipeline, X_train, y_train, cv=5)
print(f"CV Accuracy: {scores.mean():.3f} (+/- {scores.std():.3f})")
```

### Regression Quick Start

```python
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.pipeline import Pipeline

# Create pipeline
pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('regressor', RandomForestRegressor(n_estimators=100, random_state=42))
])

# Split data
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

# Train
pipeline.fit(X_train, y_train)

# Evaluate
y_pred = pipeline.predict(X_test)
rmse = mean_squared_error(y_test, y_pred, squared=False)
r2 = r2_score(y_test, y_pred)
print(f"RMSE: {rmse:.3f}, R²: {r2:.3f}")
```

## Algorithm Selection Guide

### Classification Algorithms

**Start with baseline**: LogisticRegression
- Fast, interpretable, works well for linearly separable data
- Good for high-dimensional data (text classification)

**General-purpose**: RandomForestClassifier
- Handles non-linear relationships
- Robust to outliers
- Provides feature importance
- Good default choice

**Best performance**: HistGradientBoostingClassifier
- State-of-the-art for tabular data
- Fast on large datasets (>10K samples)
- Often wins Kaggle competitions

**Special cases**:
- **Small datasets (<1K)**: SVC with RBF kernel
- **Very large datasets (>100K)**: SGDClassifier or LinearSVC
- **Interpretability critical**: LogisticRegression or DecisionTreeClassifier
- **Probabilistic predictions**: GaussianNB or calibrated models
- **Text classification**: LogisticRegression with TfidfVectorizer

### Regression Algorithms

**Start with baseline**: LinearRegression or Ridge
- Fast, interpretable
- Works well when relationships are linear

**General-purpose**: RandomForestRegressor
- Handles non-linear relationships
- Robust to outliers
- Good default choice

**Best performance**: HistGradientBoostingRegressor
- State-of-the-art for tabular data
- Fast on large datasets

**Special cases**:
- **Regularization needed**: Ridge (L2) or Lasso (L1 + feature selection)
- **Very large datasets**: SGDRegressor
- **Outliers present**: HuberRegressor or RANSAC

### Clustering Algorithms

**Known number of clusters**: KMeans
- Fast and scalable
- Assumes spherical clusters

**Unknown number of clusters**: DBSCAN or HDBSCAN
- Handles arbitrary shapes
- Automatic outlier detection

**Hierarchical relationships**: AgglomerativeClustering
- Creates hierarchy of clusters
- Good for visualization (dendrograms)

**Soft clustering (probabilities)**: GaussianMixture
- Provides cluster probabilities
- Handles elliptical clusters

### Dimensionality Reduction

**Preprocessing/feature extraction**: PCA
- Fast and efficient
- Linear transformation
- ALWAYS standardize first

**Visualization only**: t-SNE or UMAP
- Preserves local structure
- Non-linear
- DO NOT use for preprocessing

**Sparse data (text)**: TruncatedSVD
- Works with sparse matrices
- Latent Semantic Analysis

**Non-negative data**: NMF
- Interpretable components
- Topic modeling

## Working with Different Data Types

### Numeric Features

**Continuous features**:
1. Check distribution
2. Handle outliers (remove, clip, or use RobustScaler)
3. Scale using StandardScaler (most algorithms) or MinMaxScaler (neural networks)

**Count data**:
1. Consider log transformation or sqrt
2. Scale after transformation

**Skewed data**:
1. Use PowerTransformer (Yeo-Johnson or Box-Cox)
2. Or QuantileTransformer for stronger normalization

### Categorical Features

**Low cardinality (<10 categories)**:
```python
from sklearn.preprocessing import OneHotEncoder
encoder = OneHotEncoder(drop='first', sparse_output=True)
```

**High cardinality (>10 categories)**:
```python
from sklearn.preprocessing import TargetEncoder
encoder = TargetEncoder()
# Uses target statistics, prevents leakage with cross-fitting
```

**Ordinal relationships**:
```python
from sklearn.preprocessing import OrdinalEncoder
encoder = OrdinalEncoder(categories=[['small', 'medium', 'large']])
```

### Text Data

```python
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.naive_bayes import MultinomialNB
from sklearn.pipeline import Pipeline

text_pipeline = Pipeline([
    ('tfidf', TfidfVectorizer(max_features=1000, stop_words='english')),
    ('classifier', MultinomialNB())
])

text_pipeline.fit(X_train_text, y_train)
```

### Mixed Data Types

```python
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline

# Define feature types
numeric_features = ['age', 'income', 'credit_score']
categorical_features = ['country', 'occupation']

# Separate preprocessing pipelines
numeric_transformer = Pipeline([
    ('imputer', SimpleImputer(strategy='median')),
    ('scaler', StandardScaler())
])

categorical_transformer = Pipeline([
    ('imputer', SimpleImputer(strategy='constant', fill_value='missing')),
    ('onehot', OneHotEncoder(handle_unknown='ignore', sparse_output=True))
])

# Combine with ColumnTransformer
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

pipeline.fit(X_train, y_train)
```

## Model Evaluation

### Classification Metrics

**Balanced datasets**: Use accuracy or F1-score

**Imbalanced datasets**: Use balanced_accuracy, F1-weighted, or ROC-AUC
```python
from sklearn.metrics import balanced_accuracy_score, f1_score, roc_auc_score

balanced_acc = balanced_accuracy_score(y_true, y_pred)
f1 = f1_score(y_true, y_pred, average='weighted')

# ROC-AUC requires probabilities
y_proba = model.predict_proba(X_test)
auc = roc_auc_score(y_true, y_proba, multi_class='ovr')
```

**Cost-sensitive**: Define custom scorer or adjust decision threshold

**Comprehensive report**:
```python
from sklearn.metrics import classification_report, confusion_matrix

print(classification_report(y_true, y_pred))
print(confusion_matrix(y_true, y_pred))
```

### Regression Metrics

**Standard use**: RMSE and R²
```python
from sklearn.metrics import mean_squared_error, r2_score

rmse = mean_squared_error(y_true, y_pred, squared=False)
r2 = r2_score(y_true, y_pred)
```

**Outliers present**: Use MAE (robust to outliers)
```python
from sklearn.metrics import mean_absolute_error
mae = mean_absolute_error(y_true, y_pred)
```

**Percentage errors matter**: Use MAPE
```python
from sklearn.metrics import mean_absolute_percentage_error
mape = mean_absolute_percentage_error(y_true, y_pred)
```

### Cross-Validation

**Standard approach** (5-10 folds):
```python
from sklearn.model_selection import cross_val_score

scores = cross_val_score(model, X, y, cv=5, scoring='accuracy')
print(f"CV Score: {scores.mean():.3f} (+/- {scores.std():.3f})")
```

**Imbalanced classes** (use stratification):
```python
from sklearn.model_selection import StratifiedKFold

cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
scores = cross_val_score(model, X, y, cv=cv)
```

**Time series** (respect temporal order):
```python
from sklearn.model_selection import TimeSeriesSplit

cv = TimeSeriesSplit(n_splits=5)
scores = cross_val_score(model, X, y, cv=cv)
```

**Multiple metrics**:
```python
from sklearn.model_selection import cross_validate

scoring = ['accuracy', 'precision_weighted', 'recall_weighted', 'f1_weighted']
results = cross_validate(model, X, y, cv=5, scoring=scoring)

for metric in scoring:
    scores = results[f'test_{metric}']
    print(f"{metric}: {scores.mean():.3f}")
```

## Hyperparameter Tuning

### Grid Search (Exhaustive)

```python
from sklearn.model_selection import GridSearchCV

param_grid = {
    'n_estimators': [100, 200, 500],
    'max_depth': [10, 20, 30, None],
    'min_samples_split': [2, 5, 10]
}

grid_search = GridSearchCV(
    RandomForestClassifier(random_state=42),
    param_grid,
    cv=5,
    scoring='f1_weighted',
    n_jobs=-1,  # Use all CPU cores
    verbose=1
)

grid_search.fit(X_train, y_train)

print(f"Best parameters: {grid_search.best_params_}")
print(f"Best CV score: {grid_search.best_score_:.3f}")

# Use best model
best_model = grid_search.best_estimator_
test_score = best_model.score(X_test, y_test)
```

### Random Search (Faster)

```python
from sklearn.model_selection import RandomizedSearchCV
from scipy.stats import randint, uniform

param_distributions = {
    'n_estimators': randint(100, 1000),
    'max_depth': randint(5, 50),
    'min_samples_split': randint(2, 20),
    'max_features': uniform(0.1, 0.9)
}

random_search = RandomizedSearchCV(
    RandomForestClassifier(random_state=42),
    param_distributions,
    n_iter=100,  # Number of combinations to try
    cv=5,
    scoring='f1_weighted',
    n_jobs=-1,
    random_state=42
)

random_search.fit(X_train, y_train)
```

### Pipeline Hyperparameter Tuning

```python
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('svm', SVC())
])

# Use double underscore for nested parameters
param_grid = {
    'svm__C': [0.1, 1, 10, 100],
    'svm__kernel': ['rbf', 'linear'],
    'svm__gamma': ['scale', 'auto', 0.001, 0.01]
}

grid_search = GridSearchCV(pipeline, param_grid, cv=5, n_jobs=-1)
grid_search.fit(X_train, y_train)
```

## Feature Engineering and Selection

### Feature Importance

```python
# Tree-based models have built-in feature importance
from sklearn.ensemble import RandomForestClassifier

model = RandomForestClassifier(n_estimators=100)
model.fit(X_train, y_train)

importances = model.feature_importances_
feature_importance_df = pd.DataFrame({
    'feature': feature_names,
    'importance': importances
}).sort_values('importance', ascending=False)

# Permutation importance (works for any model)
from sklearn.inspection import permutation_importance

result = permutation_importance(
    model, X_test, y_test,
    n_repeats=10,
    random_state=42,
    n_jobs=-1
)

importance_df = pd.DataFrame({
    'feature': feature_names,
    'importance': result.importances_mean,
    'std': result.importances_std
}).sort_values('importance', ascending=False)
```

### Feature Selection Methods

**Univariate selection**:
```python
from sklearn.feature_selection import SelectKBest, f_classif

selector = SelectKBest(f_classif, k=10)
X_selected = selector.fit_transform(X, y)
selected_features = selector.get_support(indices=True)
```

**Recursive Feature Elimination**:
```python
from sklearn.feature_selection import RFECV
from sklearn.ensemble import RandomForestClassifier

selector = RFECV(
    RandomForestClassifier(n_estimators=100),
    step=1,
    cv=5,
    n_jobs=-1
)
X_selected = selector.fit_transform(X, y)
print(f"Optimal features: {selector.n_features_}")
```

**Model-based selection**:
```python
from sklearn.feature_selection import SelectFromModel

selector = SelectFromModel(
    RandomForestClassifier(n_estimators=100),
    threshold='median'  # or '0.5*mean', or specific value
)
X_selected = selector.fit_transform(X, y)
```

### Polynomial Features

```python
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import Ridge
from sklearn.pipeline import Pipeline

pipeline = Pipeline([
    ('poly', PolynomialFeatures(degree=2, include_bias=False)),
    ('scaler', StandardScaler()),
    ('ridge', Ridge())
])

pipeline.fit(X_train, y_train)
```

## Common Patterns and Best Practices

### Always Use Pipelines

Pipelines prevent data leakage and ensure proper workflow:

✅ **Correct**:
```python
pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('model', LogisticRegression())
])
pipeline.fit(X_train, y_train)
y_pred = pipeline.predict(X_test)
```

❌ **Wrong** (data leakage):
```python
scaler = StandardScaler().fit(X)  # Fit on all data!
X_train, X_test = train_test_split(scaler.transform(X))
```

### Stratify for Imbalanced Classes

```python
# Always use stratify for classification with imbalanced classes
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, stratify=y, random_state=42
)
```

### Scale When Necessary

**Scale for**: SVM, Neural Networks, KNN, Linear Models with regularization, PCA, Gradient Descent

**Don't scale for**: Tree-based models (Random Forest, Gradient Boosting), Naive Bayes

### Handle Missing Values

```python
from sklearn.impute import SimpleImputer

# Numeric: use median (robust to outliers)
imputer = SimpleImputer(strategy='median')

# Categorical: use constant value or most_frequent
imputer = SimpleImputer(strategy='constant', fill_value='missing')
```

### Use Appropriate Metrics

- **Balanced classification**: accuracy, F1
- **Imbalanced classification**: balanced_accuracy, F1-weighted, ROC-AUC
- **Regression with outliers**: MAE instead of RMSE
- **Cost-sensitive**: custom scorer

### Set Random States

```python
# For reproducibility
model = RandomForestClassifier(random_state=42)
X_train, X_test, y_train, y_test = train_test_split(
    X, y, random_state=42
)
```

### Use Parallel Processing

```python
# Use all CPU cores
model = RandomForestClassifier(n_jobs=-1)
grid_search = GridSearchCV(model, param_grid, n_jobs=-1)
```

## Unsupervised Learning

### Clustering Workflow

```python
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

# Always scale for clustering
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Elbow method to find optimal k
inertias = []
silhouette_scores = []
K_range = range(2, 11)

for k in K_range:
    kmeans = KMeans(n_clusters=k, random_state=42)
    labels = kmeans.fit_predict(X_scaled)
    inertias.append(kmeans.inertia_)
    silhouette_scores.append(silhouette_score(X_scaled, labels))

# Plot and choose k
import matplotlib.pyplot as plt
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
ax1.plot(K_range, inertias, 'bo-')
ax1.set_xlabel('k')
ax1.set_ylabel('Inertia')
ax2.plot(K_range, silhouette_scores, 'ro-')
ax2.set_xlabel('k')
ax2.set_ylabel('Silhouette Score')
plt.show()

# Fit final model
optimal_k = 5  # Based on elbow/silhouette
kmeans = KMeans(n_clusters=optimal_k, random_state=42)
labels = kmeans.fit_predict(X_scaled)
```

### Dimensionality Reduction

```python
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# ALWAYS scale before PCA
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Specify variance to retain
pca = PCA(n_components=0.95)  # Keep 95% of variance
X_pca = pca.fit_transform(X_scaled)

print(f"Original features: {X.shape[1]}")
print(f"Reduced features: {pca.n_components_}")
print(f"Variance explained: {pca.explained_variance_ratio_.sum():.3f}")

# Visualize explained variance
import matplotlib.pyplot as plt
plt.plot(np.cumsum(pca.explained_variance_ratio_))
plt.xlabel('Number of components')
plt.ylabel('Cumulative explained variance')
plt.show()
```

### Visualization with t-SNE

```python
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

# Reduce to 50 dimensions with PCA first (faster)
pca = PCA(n_components=min(50, X.shape[1]))
X_pca = pca.fit_transform(X_scaled)

# Apply t-SNE (only for visualization!)
tsne = TSNE(n_components=2, random_state=42, perplexity=30)
X_tsne = tsne.fit_transform(X_pca)

# Visualize
plt.scatter(X_tsne[:, 0], X_tsne[:, 1], c=y, cmap='viridis', alpha=0.6)
plt.colorbar()
plt.title('t-SNE Visualization')
plt.show()
```

## Saving and Loading Models

```python
import joblib

# Save model or pipeline
joblib.dump(model, 'model.pkl')
joblib.dump(pipeline, 'pipeline.pkl')

# Load
loaded_model = joblib.load('model.pkl')
loaded_pipeline = joblib.load('pipeline.pkl')

# Use loaded model
predictions = loaded_model.predict(X_new)
```

## Reference Documentation

This skill includes comprehensive reference files:

- **`references/supervised_learning.md`**: Detailed coverage of all classification and regression algorithms, parameters, use cases, and selection guidelines
- **`references/preprocessing.md`**: Complete guide to data preprocessing including scaling, encoding, imputation, transformations, and best practices
- **`references/model_evaluation.md`**: In-depth coverage of cross-validation strategies, metrics, hyperparameter tuning, and validation techniques
- **`references/unsupervised_learning.md`**: Comprehensive guide to clustering, dimensionality reduction, anomaly detection, and evaluation methods
- **`references/pipelines_and_composition.md`**: Complete guide to Pipeline, ColumnTransformer, FeatureUnion, custom transformers, and composition patterns
- **`references/quick_reference.md`**: Quick lookup guide with code snippets, common patterns, and decision trees for algorithm selection

Read these files when:
- Need detailed parameter explanations for specific algorithms
- Comparing multiple algorithms for a task
- Understanding evaluation metrics in depth
- Building complex preprocessing workflows
- Troubleshooting common issues

Example search patterns:
```python
# To find information about specific algorithms
grep -r "GradientBoosting" references/

# To find preprocessing techniques
grep -r "OneHotEncoder" references/preprocessing.md

# To find evaluation metrics
grep -r "f1_score" references/model_evaluation.md
```

## Common Pitfalls to Avoid

1. **Data leakage**: Always use pipelines, fit only on training data
2. **Not scaling**: Scale for distance-based algorithms (SVM, KNN, Neural Networks)
3. **Wrong metrics**: Use appropriate metrics for imbalanced data
4. **Not using cross-validation**: Single train-test split can be misleading
5. **Forgetting stratification**: Stratify for imbalanced classification
6. **Using t-SNE for preprocessing**: t-SNE is for visualization only!
7. **Not setting random_state**: Results won't be reproducible
8. **Ignoring class imbalance**: Use stratification, appropriate metrics, or resampling
9. **PCA without scaling**: Components will be dominated by high-variance features
10. **Testing on training data**: Always evaluate on held-out test set
