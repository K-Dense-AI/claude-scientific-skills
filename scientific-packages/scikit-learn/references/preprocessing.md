# Data Preprocessing in scikit-learn

## Overview
Preprocessing transforms raw data into a format suitable for machine learning algorithms. Many algorithms require standardized or normalized data to perform well.

## Standardization and Scaling

### StandardScaler
Removes mean and scales to unit variance (z-score normalization).

**Formula**: `z = (x - μ) / σ`

**Use cases**:
- Most ML algorithms (especially SVM, neural networks, PCA)
- When features have different units or scales
- When assuming Gaussian-like distribution

**Important**: Fit only on training data, then transform both train and test sets.

```python
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)  # Use same parameters
```

### MinMaxScaler
Scales features to a specified range, typically [0, 1].

**Formula**: `X_scaled = (X - X_min) / (X_max - X_min)`

**Use cases**:
- When bounded range is needed
- Neural networks (often prefer [0, 1] range)
- When distribution is not Gaussian
- Image pixel values

**Parameters**:
- `feature_range`: Tuple (min, max), default (0, 1)

**Warning**: Sensitive to outliers since it uses min/max.

### MaxAbsScaler
Scales to [-1, 1] by dividing by maximum absolute value.

**Use cases**:
- Sparse data (preserves sparsity)
- Data already centered at zero
- When sign of values is meaningful

**Advantage**: Doesn't shift/center the data, preserves zero entries.

### RobustScaler
Uses median and interquartile range (IQR) instead of mean and standard deviation.

**Formula**: `X_scaled = (X - median) / IQR`

**Use cases**:
- When outliers are present
- When StandardScaler produces skewed results
- Robust statistics preferred

**Parameters**:
- `quantile_range`: Tuple (q_min, q_max), default (25.0, 75.0)

## Normalization

### normalize() function and Normalizer
Scales individual samples (rows) to unit norm, not features (columns).

**Use cases**:
- Text classification (TF-IDF vectors)
- When similarity metrics (dot product, cosine) are used
- When each sample should have equal weight

**Norms**:
- `l1`: Manhattan norm (sum of absolutes = 1)
- `l2`: Euclidean norm (sum of squares = 1) - **most common**
- `max`: Maximum absolute value = 1

**Key difference from scalers**: Operates on rows (samples), not columns (features).

```python
from sklearn.preprocessing import Normalizer
normalizer = Normalizer(norm='l2')
X_normalized = normalizer.transform(X)
```

## Encoding Categorical Features

### OrdinalEncoder
Converts categories to integers (0 to n_categories - 1).

**Use cases**:
- Ordinal relationships exist (small < medium < large)
- Preprocessing before other transformations
- Tree-based algorithms (which can handle integers)

**Parameters**:
- `handle_unknown`: 'error' or 'use_encoded_value'
- `unknown_value`: Value for unknown categories
- `encoded_missing_value`: Value for missing data

```python
from sklearn.preprocessing import OrdinalEncoder
encoder = OrdinalEncoder()
X_encoded = encoder.fit_transform(X_categorical)
```

### OneHotEncoder
Creates binary columns for each category.

**Use cases**:
- Nominal categories (no order)
- Linear models, neural networks
- When category relationships shouldn't be assumed

**Parameters**:
- `drop`: 'first', 'if_binary', array-like (prevents multicollinearity)
- `sparse_output`: True (default, memory efficient) or False
- `handle_unknown`: 'error', 'ignore', 'infrequent_if_exist'
- `min_frequency`: Group infrequent categories
- `max_categories`: Limit number of categories

**High cardinality handling**:
```python
encoder = OneHotEncoder(min_frequency=100, handle_unknown='infrequent_if_exist')
# Groups categories appearing < 100 times into 'infrequent' category
```

**Memory tip**: Use `sparse_output=True` (default) for high-cardinality features.

### TargetEncoder
Uses target statistics to encode categories.

**Use cases**:
- High-cardinality categorical features (zip codes, user IDs)
- When linear relationships with target are expected
- Often improves performance over one-hot encoding

**How it works**:
- Replaces category with mean of target for that category
- Uses cross-fitting during fit_transform() to prevent target leakage
- Applies smoothing to handle rare categories

**Parameters**:
- `smooth`: Smoothing parameter for rare categories
- `cv`: Cross-validation strategy

**Warning**: Only for supervised learning. Requires target variable.

```python
from sklearn.preprocessing import TargetEncoder
encoder = TargetEncoder()
X_encoded = encoder.fit_transform(X_categorical, y)
```

### LabelEncoder
Encodes target labels into integers 0 to n_classes - 1.

**Use cases**: Encoding target variable for classification (not features!)

**Important**: Use `LabelEncoder` for targets, not features. For features, use OrdinalEncoder or OneHotEncoder.

### Binarizer
Converts numeric values to binary (0 or 1) based on threshold.

**Use cases**: Creating binary features from continuous values

## Non-linear Transformations

### QuantileTransformer
Maps features to uniform or normal distribution using rank transformation.

**Use cases**:
- Unusual distributions (bimodal, heavy tails)
- Reducing outlier impact
- When normal distribution is desired

**Parameters**:
- `output_distribution`: 'uniform' (default) or 'normal'
- `n_quantiles`: Number of quantiles (default: min(1000, n_samples))

**Effect**: Strong transformation that reduces outlier influence and makes data more Gaussian-like.

### PowerTransformer
Applies parametric monotonic transformation to make data more Gaussian.

**Methods**:
- `yeo-johnson`: Works with positive and negative values (default)
- `box-cox`: Only positive values

**Use cases**:
- Skewed distributions
- When Gaussian assumption is important
- Variance stabilization

**Advantage**: Less radical than QuantileTransformer, preserves more of original relationships.

## Discretization

### KBinsDiscretizer
Bins continuous features into discrete intervals.

**Strategies**:
- `uniform`: Equal-width bins
- `quantile`: Equal-frequency bins
- `kmeans`: K-means clustering to determine bins

**Encoding**:
- `ordinal`: Integer encoding (0 to n_bins - 1)
- `onehot`: One-hot encoding
- `onehot-dense`: Dense one-hot encoding

**Use cases**:
- Making linear models handle non-linear relationships
- Reducing noise in features
- Making features more interpretable

```python
from sklearn.preprocessing import KBinsDiscretizer
disc = KBinsDiscretizer(n_bins=5, encode='onehot', strategy='quantile')
X_binned = disc.fit_transform(X)
```

## Feature Generation

### PolynomialFeatures
Generates polynomial and interaction features.

**Parameters**:
- `degree`: Polynomial degree
- `interaction_only`: Only multiplicative interactions (no x²)
- `include_bias`: Include constant feature

**Use cases**:
- Adding non-linearity to linear models
- Feature engineering
- Polynomial regression

**Warning**: Number of features grows rapidly: (n+d)!/d!n! for degree d.

```python
from sklearn.preprocessing import PolynomialFeatures
poly = PolynomialFeatures(degree=2, include_bias=False)
X_poly = poly.fit_transform(X)
# [x1, x2] → [x1, x2, x1², x1·x2, x2²]
```

### SplineTransformer
Generates B-spline basis functions.

**Use cases**:
- Smooth non-linear transformations
- Alternative to PolynomialFeatures (less oscillation at boundaries)
- Generalized additive models (GAMs)

**Parameters**:
- `n_knots`: Number of knots
- `degree`: Spline degree
- `knots`: Knot positions ('uniform', 'quantile', or array)

## Missing Value Handling

### SimpleImputer
Imputes missing values with various strategies.

**Strategies**:
- `mean`: Mean of column (numeric only)
- `median`: Median of column (numeric only)
- `most_frequent`: Mode (numeric or categorical)
- `constant`: Fill with constant value

**Parameters**:
- `strategy`: Imputation strategy
- `fill_value`: Value when strategy='constant'
- `missing_values`: What represents missing (np.nan, None, specific value)

```python
from sklearn.impute import SimpleImputer
imputer = SimpleImputer(strategy='median')
X_imputed = imputer.fit_transform(X)
```

### KNNImputer
Imputes using k-nearest neighbors.

**Use cases**: When relationships between features should inform imputation

**Parameters**:
- `n_neighbors`: Number of neighbors
- `weights`: 'uniform' or 'distance'

### IterativeImputer
Models each feature with missing values as function of other features.

**Use cases**:
- Complex relationships between features
- When multiple features have missing values
- Higher quality imputation (but slower)

**Parameters**:
- `estimator`: Estimator for regression (default: BayesianRidge)
- `max_iter`: Maximum iterations

## Function Transformers

### FunctionTransformer
Applies custom function to data.

**Use cases**:
- Custom transformations in pipelines
- Log transformation, square root, etc.
- Domain-specific preprocessing

```python
from sklearn.preprocessing import FunctionTransformer
import numpy as np

log_transformer = FunctionTransformer(np.log1p, validate=True)
X_log = log_transformer.transform(X)
```

## Best Practices

### Feature Scaling Guidelines

**Always scale**:
- SVM, neural networks
- K-nearest neighbors
- Linear/Logistic regression with regularization
- PCA, LDA
- Gradient descent-based algorithms

**Don't need to scale**:
- Tree-based algorithms (Decision Trees, Random Forests, Gradient Boosting)
- Naive Bayes

### Pipeline Integration

Always use preprocessing within pipelines to prevent data leakage:

```python
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('classifier', LogisticRegression())
])

pipeline.fit(X_train, y_train)  # Scaler fit only on train data
y_pred = pipeline.predict(X_test)  # Scaler transform only on test data
```

### Common Transformations by Data Type

**Numeric - Continuous**:
- StandardScaler (most common)
- MinMaxScaler (neural networks)
- RobustScaler (outliers present)
- PowerTransformer (skewed data)

**Numeric - Count Data**:
- sqrt or log transformation
- QuantileTransformer
- StandardScaler after transformation

**Categorical - Low Cardinality (<10 categories)**:
- OneHotEncoder

**Categorical - High Cardinality (>10 categories)**:
- TargetEncoder (supervised)
- Frequency encoding
- OneHotEncoder with min_frequency parameter

**Categorical - Ordinal**:
- OrdinalEncoder

**Text**:
- CountVectorizer or TfidfVectorizer
- Normalizer after vectorization

### Data Leakage Prevention

1. **Fit only on training data**: Never include test data when fitting preprocessors
2. **Use pipelines**: Ensures proper fit/transform separation
3. **Cross-validation**: Use Pipeline with cross_val_score() for proper evaluation
4. **Target encoding**: Use cv parameter in TargetEncoder for cross-fitting

```python
# WRONG - data leakage
scaler = StandardScaler().fit(X_full)
X_train_scaled = scaler.transform(X_train)
X_test_scaled = scaler.transform(X_test)

# CORRECT - no leakage
scaler = StandardScaler().fit(X_train)
X_train_scaled = scaler.transform(X_train)
X_test_scaled = scaler.transform(X_test)
```

## Preprocessing Checklist

Before modeling:
1. Handle missing values (imputation or removal)
2. Encode categorical variables appropriately
3. Scale/normalize numeric features (if needed for algorithm)
4. Handle outliers (RobustScaler, clipping, removal)
5. Create additional features if beneficial (PolynomialFeatures, domain knowledge)
6. Check for data leakage in preprocessing steps
7. Wrap everything in a Pipeline
