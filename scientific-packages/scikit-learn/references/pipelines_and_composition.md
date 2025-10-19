# Pipelines and Composite Estimators in scikit-learn

## Overview
Pipelines chain multiple estimators into a single unit, ensuring proper workflow sequencing and preventing data leakage. As the documentation states: "Pipeline can be used to chain multiple estimators into one. This is useful as there is often a fixed sequence of steps in processing the data, for example feature selection, normalization and classification."

## Pipeline Basics

### Creating Pipelines

```python
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression

# Method 1: List of (name, estimator) tuples
pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('pca', PCA(n_components=10)),
    ('classifier', LogisticRegression())
])

# Method 2: Using make_pipeline (auto-generates names)
from sklearn.pipeline import make_pipeline
pipeline = make_pipeline(
    StandardScaler(),
    PCA(n_components=10),
    LogisticRegression()
)
```

### Using Pipelines

```python
# Fit and predict like any estimator
pipeline.fit(X_train, y_train)
y_pred = pipeline.predict(X_test)
score = pipeline.score(X_test, y_test)

# Access steps
pipeline.named_steps['scaler']
pipeline.steps[0]  # Returns ('scaler', StandardScaler(...))
pipeline[0]        # Returns StandardScaler(...) object
pipeline['scaler'] # Returns StandardScaler(...) object

# Get final estimator
pipeline[-1]  # Returns LogisticRegression(...) object
```

### Pipeline Rules

**All steps except the last must be transformers** (have `fit()` and `transform()` methods).

**The final step** can be:
- Predictor (classifier/regressor) with `fit()` and `predict()`
- Transformer with `fit()` and `transform()`
- Any estimator with at least `fit()`

### Pipeline Benefits

1. **Convenience**: Single `fit()` and `predict()` call
2. **Prevents data leakage**: Ensures proper fit/transform on train/test
3. **Joint parameter selection**: Tune all steps together with GridSearchCV
4. **Reproducibility**: Encapsulates entire workflow

## Accessing and Setting Parameters

### Nested Parameters

Access step parameters using `stepname__parameter` syntax:

```python
from sklearn.model_selection import GridSearchCV

pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('clf', LogisticRegression())
])

# Grid search over pipeline parameters
param_grid = {
    'scaler__with_mean': [True, False],
    'clf__C': [0.1, 1.0, 10.0],
    'clf__penalty': ['l1', 'l2']
}

grid_search = GridSearchCV(pipeline, param_grid, cv=5)
grid_search.fit(X_train, y_train)
```

### Setting Parameters

```python
# Set parameters
pipeline.set_params(clf__C=10.0, scaler__with_std=False)

# Get parameters
params = pipeline.get_params()
```

## Caching Intermediate Results

Cache fitted transformers to avoid recomputation:

```python
from tempfile import mkdtemp
from shutil import rmtree

# Create cache directory
cachedir = mkdtemp()

pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('pca', PCA(n_components=10)),
    ('clf', LogisticRegression())
], memory=cachedir)

# When doing grid search, scaler and PCA only fit once per fold
grid_search = GridSearchCV(pipeline, param_grid, cv=5)
grid_search.fit(X_train, y_train)

# Clean up cache
rmtree(cachedir)

# Or use joblib for persistent caching
from joblib import Memory
memory = Memory(location='./cache', verbose=0)
pipeline = Pipeline([...], memory=memory)
```

**When to use caching**:
- Expensive transformations (PCA, feature selection)
- Grid search over final estimator parameters only
- Multiple experiments with same preprocessing

## ColumnTransformer

Apply different transformations to different columns (essential for heterogeneous data).

### Basic Usage

```python
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler, OneHotEncoder

# Define which transformations for which columns
preprocessor = ColumnTransformer(
    transformers=[
        ('num', StandardScaler(), ['age', 'income', 'credit_score']),
        ('cat', OneHotEncoder(), ['country', 'occupation'])
    ],
    remainder='drop'  # What to do with remaining columns
)

X_transformed = preprocessor.fit_transform(X)
```

### Column Selection Methods

```python
# Method 1: Column names (list of strings)
('num', StandardScaler(), ['age', 'income'])

# Method 2: Column indices (list of integers)
('num', StandardScaler(), [0, 1, 2])

# Method 3: Boolean mask
('num', StandardScaler(), [True, True, False, True, False])

# Method 4: Slice
('num', StandardScaler(), slice(0, 3))

# Method 5: make_column_selector (by dtype or pattern)
from sklearn.compose import make_column_selector as selector

preprocessor = ColumnTransformer([
    ('num', StandardScaler(), selector(dtype_include='number')),
    ('cat', OneHotEncoder(), selector(dtype_include='object'))
])

# Select by pattern
selector(pattern='.*_score$')  # All columns ending with '_score'
```

### Remainder Parameter

Controls what happens to columns not specified:

```python
# Drop remaining columns (default)
remainder='drop'

# Pass through remaining columns unchanged
remainder='passthrough'

# Apply transformer to remaining columns
remainder=StandardScaler()
```

### Full Pipeline with ColumnTransformer

```python
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.ensemble import RandomForestClassifier

# Separate preprocessing for numeric and categorical
numeric_features = ['age', 'income', 'credit_score']
categorical_features = ['country', 'occupation', 'education']

numeric_transformer = Pipeline(steps=[
    ('imputer', SimpleImputer(strategy='median')),
    ('scaler', StandardScaler())
])

categorical_transformer = Pipeline(steps=[
    ('imputer', SimpleImputer(strategy='constant', fill_value='missing')),
    ('onehot', OneHotEncoder(handle_unknown='ignore'))
])

preprocessor = ColumnTransformer(
    transformers=[
        ('num', numeric_transformer, numeric_features),
        ('cat', categorical_transformer, categorical_features)
    ])

# Complete pipeline
clf = Pipeline(steps=[
    ('preprocessor', preprocessor),
    ('classifier', RandomForestClassifier())
])

clf.fit(X_train, y_train)
y_pred = clf.predict(X_test)

# Grid search over preprocessing and model parameters
param_grid = {
    'preprocessor__num__imputer__strategy': ['mean', 'median'],
    'preprocessor__cat__onehot__max_categories': [10, 20, None],
    'classifier__n_estimators': [100, 200],
    'classifier__max_depth': [10, 20, None]
}

grid_search = GridSearchCV(clf, param_grid, cv=5)
grid_search.fit(X_train, y_train)
```

## FeatureUnion

Combine multiple transformer outputs by concatenating features side-by-side.

```python
from sklearn.pipeline import FeatureUnion
from sklearn.decomposition import PCA
from sklearn.feature_selection import SelectKBest

# Combine PCA and feature selection
combined_features = FeatureUnion([
    ('pca', PCA(n_components=10)),
    ('univ_select', SelectKBest(k=5))
])

X_features = combined_features.fit_transform(X, y)
# Result: 15 features (10 from PCA + 5 from SelectKBest)

# In a pipeline
pipeline = Pipeline([
    ('features', combined_features),
    ('classifier', LogisticRegression())
])
```

### FeatureUnion with Transformers on Different Data

```python
from sklearn.pipeline import FeatureUnion
from sklearn.preprocessing import FunctionTransformer
import numpy as np

def get_numeric_data(X):
    return X[:, :3]  # First 3 columns

def get_text_data(X):
    return X[:, 3]   # 4th column (text)

from sklearn.feature_extraction.text import TfidfVectorizer

combined = FeatureUnion([
    ('numeric_features', Pipeline([
        ('selector', FunctionTransformer(get_numeric_data)),
        ('scaler', StandardScaler())
    ])),
    ('text_features', Pipeline([
        ('selector', FunctionTransformer(get_text_data)),
        ('tfidf', TfidfVectorizer())
    ]))
])
```

**Note**: ColumnTransformer is usually more convenient than FeatureUnion for heterogeneous data.

## Common Pipeline Patterns

### Classification Pipeline

```python
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.svm import SVC

pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('feature_selection', SelectKBest(f_classif, k=10)),
    ('classifier', SVC(kernel='rbf'))
])
```

### Regression Pipeline

```python
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.linear_model import Ridge

pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('poly', PolynomialFeatures(degree=2)),
    ('ridge', Ridge(alpha=1.0))
])
```

### Text Classification Pipeline

```python
from sklearn.pipeline import Pipeline
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.naive_bayes import MultinomialNB

pipeline = Pipeline([
    ('tfidf', TfidfVectorizer(max_features=1000)),
    ('classifier', MultinomialNB())
])

# Works directly with text
pipeline.fit(X_train_text, y_train)
y_pred = pipeline.predict(X_test_text)
```

### Image Processing Pipeline

```python
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.neural_network import MLPClassifier

pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('pca', PCA(n_components=100)),
    ('mlp', MLPClassifier(hidden_layer_sizes=(100, 50)))
])
```

### Dimensionality Reduction + Clustering

```python
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('pca', PCA(n_components=10)),
    ('kmeans', KMeans(n_clusters=5))
])

labels = pipeline.fit_predict(X)
```

## Custom Transformers

### Using FunctionTransformer

```python
from sklearn.preprocessing import FunctionTransformer
import numpy as np

# Log transformation
log_transformer = FunctionTransformer(np.log1p)

# Custom function
def custom_transform(X):
    # Your transformation logic
    return X_transformed

custom_transformer = FunctionTransformer(custom_transform)

# In pipeline
pipeline = Pipeline([
    ('log', log_transformer),
    ('scaler', StandardScaler()),
    ('model', LinearRegression())
])
```

### Creating Custom Transformer Class

```python
from sklearn.base import BaseEstimator, TransformerMixin

class CustomTransformer(BaseEstimator, TransformerMixin):
    def __init__(self, parameter=1.0):
        self.parameter = parameter

    def fit(self, X, y=None):
        # Learn parameters from X
        self.learned_param_ = X.mean()  # Example
        return self

    def transform(self, X):
        # Transform X using learned parameters
        return X * self.parameter - self.learned_param_

    # Optional: for pipelines that need inverse transform
    def inverse_transform(self, X):
        return (X + self.learned_param_) / self.parameter

# Use in pipeline
pipeline = Pipeline([
    ('custom', CustomTransformer(parameter=2.0)),
    ('model', LinearRegression())
])
```

**Key requirements**:
- Inherit from `BaseEstimator` and `TransformerMixin`
- Implement `fit()` and `transform()` methods
- `fit()` must return `self`
- Use trailing underscore for learned attributes (`learned_param_`)
- Constructor parameters should be stored as attributes

### Transformer for Pandas DataFrames

```python
from sklearn.base import BaseEstimator, TransformerMixin
import pandas as pd

class DataFrameTransformer(BaseEstimator, TransformerMixin):
    def __init__(self, columns=None):
        self.columns = columns

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        if isinstance(X, pd.DataFrame):
            if self.columns:
                return X[self.columns].values
            return X.values
        return X
```

## Visualization

### Display Pipeline in Jupyter

```python
from sklearn import set_config

# Enable HTML display
set_config(display='diagram')

# Now displaying the pipeline shows interactive diagram
pipeline
```

### Print Pipeline Structure

```python
from sklearn.utils import estimator_html_repr

# Get HTML representation
html = estimator_html_repr(pipeline)

# Or just print
print(pipeline)
```

## Advanced Patterns

### Conditional Transformations

```python
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, FunctionTransformer

def conditional_scale(X, scale=True):
    if scale:
        return StandardScaler().fit_transform(X)
    return X

pipeline = Pipeline([
    ('conditional_scaler', FunctionTransformer(
        conditional_scale,
        kw_args={'scale': True}
    )),
    ('model', LogisticRegression())
])
```

### Multiple Preprocessing Paths

```python
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline

# Different preprocessing for different feature types
preprocessor = ColumnTransformer([
    # Numeric: impute + scale
    ('num_standard', Pipeline([
        ('imputer', SimpleImputer(strategy='mean')),
        ('scaler', StandardScaler())
    ]), ['age', 'income']),

    # Numeric: impute + log + scale
    ('num_skewed', Pipeline([
        ('imputer', SimpleImputer(strategy='median')),
        ('log', FunctionTransformer(np.log1p)),
        ('scaler', StandardScaler())
    ]), ['price', 'revenue']),

    # Categorical: impute + one-hot
    ('cat', Pipeline([
        ('imputer', SimpleImputer(strategy='constant', fill_value='missing')),
        ('onehot', OneHotEncoder(handle_unknown='ignore'))
    ]), ['category', 'region']),

    # Text: TF-IDF
    ('text', TfidfVectorizer(), 'description')
])
```

### Feature Engineering Pipeline

```python
from sklearn.base import BaseEstimator, TransformerMixin

class FeatureEngineer(BaseEstimator, TransformerMixin):
    def fit(self, X, y=None):
        return self

    def transform(self, X):
        X = X.copy()
        # Add engineered features
        X['age_income_ratio'] = X['age'] / (X['income'] + 1)
        X['total_score'] = X['score1'] + X['score2'] + X['score3']
        return X

pipeline = Pipeline([
    ('engineer', FeatureEngineer()),
    ('preprocessor', preprocessor),
    ('model', RandomForestClassifier())
])
```

## Best Practices

### Always Use Pipelines When

1. **Preprocessing is needed**: Scaling, encoding, imputation
2. **Cross-validation**: Ensures proper fit/transform split
3. **Hyperparameter tuning**: Joint optimization of preprocessing and model
4. **Production deployment**: Single object to serialize
5. **Multiple steps**: Any workflow with >1 step

### Pipeline Do's

- ✅ Fit pipeline only on training data
- ✅ Use ColumnTransformer for heterogeneous data
- ✅ Cache expensive transformations during grid search
- ✅ Use make_pipeline for simple cases
- ✅ Set verbose=True to debug issues
- ✅ Use remainder='passthrough' when appropriate

### Pipeline Don'ts

- ❌ Fit preprocessing on full dataset before split (data leakage!)
- ❌ Manually transform test data (use pipeline.predict())
- ❌ Forget to handle missing values before scaling
- ❌ Mix pandas DataFrames and arrays inconsistently
- ❌ Skip using pipelines for "just one preprocessing step"

### Data Leakage Prevention

```python
# ❌ WRONG - Data leakage
scaler = StandardScaler().fit(X)  # Fit on all data
X_train, X_test, y_train, y_test = train_test_split(X, y)
X_train_scaled = scaler.transform(X_train)
X_test_scaled = scaler.transform(X_test)

# ✅ CORRECT - No leakage with pipeline
pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('model', LogisticRegression())
])

X_train, X_test, y_train, y_test = train_test_split(X, y)
pipeline.fit(X_train, y_train)  # Scaler fits only on train
y_pred = pipeline.predict(X_test)  # Scaler transforms only on test

# ✅ CORRECT - No leakage in cross-validation
scores = cross_val_score(pipeline, X, y, cv=5)
# Each fold: scaler fits on train folds, transforms on test fold
```

### Debugging Pipelines

```python
# Examine intermediate outputs
pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('pca', PCA(n_components=10)),
    ('model', LogisticRegression())
])

# Fit pipeline
pipeline.fit(X_train, y_train)

# Get output after scaling
X_scaled = pipeline.named_steps['scaler'].transform(X_train)

# Get output after PCA
X_pca = pipeline[:-1].transform(X_train)  # All steps except last

# Or build partial pipeline
partial_pipeline = Pipeline(pipeline.steps[:-1])
X_transformed = partial_pipeline.transform(X_train)
```

### Saving and Loading Pipelines

```python
import joblib

# Save pipeline
joblib.dump(pipeline, 'model_pipeline.pkl')

# Load pipeline
pipeline = joblib.load('model_pipeline.pkl')

# Use loaded pipeline
y_pred = pipeline.predict(X_new)
```

## Common Errors and Solutions

**Error**: `ValueError: could not convert string to float`
- **Cause**: Categorical features not encoded
- **Solution**: Add OneHotEncoder or OrdinalEncoder to pipeline

**Error**: `All intermediate steps should be transformers`
- **Cause**: Non-transformer in non-final position
- **Solution**: Ensure only last step is predictor

**Error**: `X has different number of features than during fitting`
- **Cause**: Different columns in train and test
- **Solution**: Ensure consistent column handling, use `handle_unknown='ignore'` in OneHotEncoder

**Error**: Different results in cross-validation vs train-test split
- **Cause**: Data leakage (fitting preprocessing on all data)
- **Solution**: Always use Pipeline for preprocessing

**Error**: Pipeline too slow during grid search
- **Solution**: Use caching with `memory` parameter
