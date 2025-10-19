# Model Evaluation and Selection in scikit-learn

## Overview
Model evaluation assesses how well models generalize to unseen data. Scikit-learn provides three main APIs for evaluation:
1. **Estimator score methods**: Built-in evaluation (accuracy for classifiers, R² for regressors)
2. **Scoring parameter**: Used in cross-validation and hyperparameter tuning
3. **Metric functions**: Specialized evaluation in `sklearn.metrics`

## Cross-Validation

Cross-validation evaluates model performance by splitting data into multiple train/test sets. This addresses overfitting: "a model that would just repeat the labels of the samples that it has just seen would have a perfect score but would fail to predict anything useful on yet-unseen data."

### Basic Cross-Validation

```python
from sklearn.model_selection import cross_val_score
from sklearn.linear_model import LogisticRegression

model = LogisticRegression()
scores = cross_val_score(model, X, y, cv=5, scoring='accuracy')
print(f"Accuracy: {scores.mean():.3f} (+/- {scores.std():.3f})")
```

### Cross-Validation Strategies

#### For i.i.d. Data

**KFold**: Standard k-fold cross-validation
- Splits data into k equal folds
- Each fold used once as test set
- `n_splits`: Number of folds (typically 5 or 10)

```python
from sklearn.model_selection import KFold
cv = KFold(n_splits=5, shuffle=True, random_state=42)
```

**RepeatedKFold**: Repeats KFold with different randomization
- More robust estimation
- Computationally expensive

**LeaveOneOut (LOO)**: Each sample is a test set
- Maximum training data usage
- Very computationally expensive
- High variance in estimates
- Use only for small datasets (<1000 samples)

**ShuffleSplit**: Random train/test splits
- Flexible train/test sizes
- Can control number of iterations
- Good for quick evaluation

```python
from sklearn.model_selection import ShuffleSplit
cv = ShuffleSplit(n_splits=10, test_size=0.2, random_state=42)
```

#### For Imbalanced Classes

**StratifiedKFold**: Preserves class proportions in each fold
- Essential for imbalanced datasets
- Default for classification in cross_val_score()

```python
from sklearn.model_selection import StratifiedKFold
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
```

**StratifiedShuffleSplit**: Stratified random splits

#### For Grouped Data

Use when samples are not independent (e.g., multiple measurements from same subject).

**GroupKFold**: Groups don't appear in both train and test
```python
from sklearn.model_selection import GroupKFold
cv = GroupKFold(n_splits=5)
scores = cross_val_score(model, X, y, groups=groups, cv=cv)
```

**StratifiedGroupKFold**: Combines stratification with group separation

**LeaveOneGroupOut**: Each group becomes a test set

#### For Time Series

**TimeSeriesSplit**: Expanding window approach
- Successive training sets are supersets
- Respects temporal ordering
- No data leakage from future to past

```python
from sklearn.model_selection import TimeSeriesSplit
cv = TimeSeriesSplit(n_splits=5)
for train_idx, test_idx in cv.split(X):
    # Train on indices 0 to t, test on t+1 to t+k
    pass
```

### Cross-Validation Functions

**cross_val_score**: Returns array of scores
```python
scores = cross_val_score(model, X, y, cv=5, scoring='f1_weighted')
```

**cross_validate**: Returns multiple metrics and timing
```python
results = cross_validate(
    model, X, y, cv=5,
    scoring=['accuracy', 'f1_weighted', 'roc_auc'],
    return_train_score=True,
    return_estimator=True  # Returns fitted estimators
)
print(results['test_accuracy'])
print(results['fit_time'])
```

**cross_val_predict**: Returns predictions for model blending/visualization
```python
from sklearn.model_selection import cross_val_predict
y_pred = cross_val_predict(model, X, y, cv=5)
# Use for confusion matrix, error analysis, etc.
```

## Hyperparameter Tuning

### GridSearchCV
Exhaustively searches all parameter combinations.

```python
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier

param_grid = {
    'n_estimators': [100, 200, 500],
    'max_depth': [10, 20, 30, None],
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 2, 4]
}

grid_search = GridSearchCV(
    RandomForestClassifier(random_state=42),
    param_grid,
    cv=5,
    scoring='f1_weighted',
    n_jobs=-1,  # Use all CPU cores
    verbose=2
)

grid_search.fit(X_train, y_train)
print("Best parameters:", grid_search.best_params_)
print("Best score:", grid_search.best_score_)

# Use best model
best_model = grid_search.best_estimator_
```

**When to use**:
- Small parameter spaces
- When computational resources allow
- When exhaustive search is desired

### RandomizedSearchCV
Samples parameter combinations from distributions.

```python
from sklearn.model_selection import RandomizedSearchCV
from scipy.stats import randint, uniform

param_distributions = {
    'n_estimators': randint(100, 1000),
    'max_depth': randint(5, 50),
    'min_samples_split': randint(2, 20),
    'min_samples_leaf': randint(1, 10),
    'max_features': uniform(0.1, 0.9)
}

random_search = RandomizedSearchCV(
    RandomForestClassifier(random_state=42),
    param_distributions,
    n_iter=100,  # Number of parameter settings sampled
    cv=5,
    scoring='f1_weighted',
    n_jobs=-1,
    random_state=42
)

random_search.fit(X_train, y_train)
```

**When to use**:
- Large parameter spaces
- When budget is limited
- Often finds good parameters faster than GridSearchCV

**Advantage**: "Budget can be chosen independent of the number of parameters and possible values"

### Successive Halving

**HalvingGridSearchCV** and **HalvingRandomSearchCV**: Tournament-style selection

**How it works**:
1. Start with many candidates, minimal resources
2. Eliminate poor performers
3. Increase resources for remaining candidates
4. Repeat until best candidates found

**When to use**:
- Large parameter spaces
- Expensive model training
- When many parameter combinations are clearly inferior

```python
from sklearn.experimental import enable_halving_search_cv
from sklearn.model_selection import HalvingGridSearchCV

halving_search = HalvingGridSearchCV(
    estimator,
    param_grid,
    factor=3,  # Proportion of candidates eliminated each round
    cv=5
)
```

## Classification Metrics

### Accuracy-Based Metrics

**Accuracy**: Proportion of correct predictions
```python
from sklearn.metrics import accuracy_score
accuracy = accuracy_score(y_true, y_pred)
```

**When to use**: Balanced datasets only
**When NOT to use**: Imbalanced datasets (misleading)

**Balanced Accuracy**: Average recall per class
```python
from sklearn.metrics import balanced_accuracy_score
bal_acc = balanced_accuracy_score(y_true, y_pred)
```

**When to use**: Imbalanced datasets, ensures all classes matter equally

### Precision, Recall, F-Score

**Precision**: Of predicted positives, how many are actually positive
- Formula: TP / (TP + FP)
- Answers: "How reliable are positive predictions?"

**Recall** (Sensitivity): Of actual positives, how many are predicted positive
- Formula: TP / (TP + FN)
- Answers: "How complete is positive detection?"

**F1-Score**: Harmonic mean of precision and recall
- Formula: 2 * (precision * recall) / (precision + recall)
- Balanced measure when both precision and recall are important

```python
from sklearn.metrics import precision_recall_fscore_support, f1_score

precision, recall, f1, support = precision_recall_fscore_support(
    y_true, y_pred, average='weighted'
)

# Or individually
f1 = f1_score(y_true, y_pred, average='weighted')
```

**Averaging strategies for multiclass**:
- `binary`: Binary classification only
- `micro`: Calculate globally (total TP, FP, FN)
- `macro`: Calculate per class, unweighted mean (all classes equal)
- `weighted`: Calculate per class, weighted by support (class frequency)
- `samples`: For multilabel classification

**When to use**:
- `macro`: When all classes equally important (even rare ones)
- `weighted`: When class frequency matters
- `micro`: When overall performance across all samples matters

### Confusion Matrix

Shows true positives, false positives, true negatives, false negatives.

```python
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import matplotlib.pyplot as plt

cm = confusion_matrix(y_true, y_pred)
disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=['Class 0', 'Class 1'])
disp.plot()
plt.show()
```

### ROC Curve and AUC

**ROC (Receiver Operating Characteristic)**: Plot of true positive rate vs false positive rate at different thresholds

**AUC (Area Under Curve)**: Measures overall ability to discriminate between classes
- 1.0 = perfect classifier
- 0.5 = random classifier
- <0.5 = worse than random

```python
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt

# Requires probability predictions
y_proba = model.predict_proba(X_test)[:, 1]  # Probabilities for positive class

auc = roc_auc_score(y_true, y_proba)
fpr, tpr, thresholds = roc_curve(y_true, y_proba)

plt.plot(fpr, tpr, label=f'AUC = {auc:.3f}')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend()
plt.show()
```

**Multiclass ROC**: Use `multi_class='ovr'` (one-vs-rest) or `'ovo'` (one-vs-one)

```python
auc = roc_auc_score(y_true, y_proba, multi_class='ovr')
```

### Log Loss

Measures probability calibration quality.

```python
from sklearn.metrics import log_loss
loss = log_loss(y_true, y_proba)
```

**When to use**: When probability quality matters, not just class predictions
**Lower is better**: Perfect predictions have log loss of 0

### Classification Report

Comprehensive summary of precision, recall, f1-score per class.

```python
from sklearn.metrics import classification_report

print(classification_report(y_true, y_pred, target_names=['Class 0', 'Class 1']))
```

## Regression Metrics

### Mean Squared Error (MSE)
Average squared difference between predictions and true values.

```python
from sklearn.metrics import mean_squared_error
mse = mean_squared_error(y_true, y_pred)
rmse = mean_squared_error(y_true, y_pred, squared=False)  # Root MSE
```

**Characteristics**:
- Penalizes large errors heavily (squared term)
- Same units as target² (use RMSE for same units as target)
- Lower is better

### Mean Absolute Error (MAE)
Average absolute difference between predictions and true values.

```python
from sklearn.metrics import mean_absolute_error
mae = mean_absolute_error(y_true, y_pred)
```

**Characteristics**:
- More robust to outliers than MSE
- Same units as target
- More interpretable
- Lower is better

**MSE vs MAE**: Use MAE when outliers shouldn't dominate the metric

### R² Score (Coefficient of Determination)
Proportion of variance explained by the model.

```python
from sklearn.metrics import r2_score
r2 = r2_score(y_true, y_pred)
```

**Interpretation**:
- 1.0 = perfect predictions
- 0.0 = model as good as mean
- <0.0 = model worse than mean (possible!)
- Higher is better

**Note**: Can be negative for models that perform worse than predicting the mean.

### Mean Absolute Percentage Error (MAPE)
Percentage-based error metric.

```python
from sklearn.metrics import mean_absolute_percentage_error
mape = mean_absolute_percentage_error(y_true, y_pred)
```

**When to use**: When relative errors matter more than absolute errors
**Warning**: Undefined when true values are zero

### Median Absolute Error
Median of absolute errors (robust to outliers).

```python
from sklearn.metrics import median_absolute_error
med_ae = median_absolute_error(y_true, y_pred)
```

### Max Error
Maximum residual error.

```python
from sklearn.metrics import max_error
max_err = max_error(y_true, y_pred)
```

**When to use**: When worst-case performance matters

## Custom Scoring Functions

Create custom scorers for GridSearchCV and cross_val_score:

```python
from sklearn.metrics import make_scorer, fbeta_score

# F2 score (weights recall higher than precision)
f2_scorer = make_scorer(fbeta_score, beta=2)

# Custom function
def custom_metric(y_true, y_pred):
    # Your custom logic
    return score

custom_scorer = make_scorer(custom_metric, greater_is_better=True)

# Use in cross-validation or grid search
scores = cross_val_score(model, X, y, cv=5, scoring=custom_scorer)
```

## Scoring Parameter Options

Common scoring strings for `scoring` parameter:

**Classification**:
- `'accuracy'`, `'balanced_accuracy'`
- `'precision'`, `'recall'`, `'f1'` (add `_macro`, `_micro`, `_weighted` for multiclass)
- `'roc_auc'`, `'roc_auc_ovr'`, `'roc_auc_ovo'`
- `'log_loss'` (lower is better, negate for maximization)
- `'jaccard'` (Jaccard similarity)

**Regression**:
- `'r2'`
- `'neg_mean_squared_error'`, `'neg_root_mean_squared_error'`
- `'neg_mean_absolute_error'`
- `'neg_mean_absolute_percentage_error'`
- `'neg_median_absolute_error'`

**Note**: Many metrics are negated (neg_*) so GridSearchCV can maximize them.

## Validation Strategies

### Train-Test Split
Simple single split.

```python
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(
    X, y,
    test_size=0.2,
    random_state=42,
    stratify=y  # For classification with imbalanced classes
)
```

**When to use**: Large datasets, quick evaluation
**Parameters**:
- `test_size`: Proportion for test (typically 0.2-0.3)
- `stratify`: Preserves class proportions
- `random_state`: Reproducibility

### Train-Validation-Test Split
Three-way split for hyperparameter tuning.

```python
# First split: train+val and test
X_trainval, X_test, y_trainval, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

# Second split: train and validation
X_train, X_val, y_train, y_val = train_test_split(
    X_trainval, y_trainval, test_size=0.2, random_state=42
)

# Or use GridSearchCV with train+val, then evaluate on test
```

**When to use**: Model selection and final evaluation
**Strategy**:
1. Train: Model training
2. Validation: Hyperparameter tuning
3. Test: Final, unbiased evaluation (touch only once!)

### Learning Curves

Diagnose bias vs variance issues.

```python
from sklearn.model_selection import learning_curve
import matplotlib.pyplot as plt

train_sizes, train_scores, val_scores = learning_curve(
    model, X, y,
    cv=5,
    train_sizes=np.linspace(0.1, 1.0, 10),
    scoring='accuracy',
    n_jobs=-1
)

plt.plot(train_sizes, train_scores.mean(axis=1), label='Training score')
plt.plot(train_sizes, val_scores.mean(axis=1), label='Validation score')
plt.xlabel('Training set size')
plt.ylabel('Score')
plt.legend()
plt.show()
```

**Interpretation**:
- Large gap between train and validation: **Overfitting** (high variance)
- Both scores low: **Underfitting** (high bias)
- Scores converging but low: Need better features or more complex model
- Validation score still improving: More data would help

## Best Practices

### Metric Selection Guidelines

**Classification - Balanced classes**:
- Accuracy or F1-score

**Classification - Imbalanced classes**:
- Balanced accuracy
- F1-score (weighted or macro)
- ROC-AUC
- Precision-Recall curve

**Classification - Cost-sensitive**:
- Custom scorer with cost matrix
- Adjust threshold on probabilities

**Regression - Typical use**:
- RMSE (sensitive to outliers)
- R² (proportion of variance explained)

**Regression - Outliers present**:
- MAE (robust to outliers)
- Median absolute error

**Regression - Percentage errors matter**:
- MAPE

### Cross-Validation Guidelines

**Number of folds**:
- 5-10 folds typical
- More folds = more computation, less variance in estimate
- LeaveOneOut only for small datasets

**Stratification**:
- Always use for classification with imbalanced classes
- Use StratifiedKFold by default for classification

**Grouping**:
- Always use when samples are not independent
- Time series: Always use TimeSeriesSplit

**Nested cross-validation**:
- For unbiased performance estimate when doing hyperparameter tuning
- Outer loop: Performance estimation
- Inner loop: Hyperparameter selection

### Avoiding Common Pitfalls

1. **Data leakage**: Fit preprocessors only on training data within each CV fold (use Pipeline!)
2. **Test set leakage**: Never use test set for model selection
3. **Improper metric**: Use metrics appropriate for problem (balanced_accuracy for imbalanced data)
4. **Multiple testing**: More models evaluated = higher chance of random good results
5. **Temporal leakage**: For time series, use TimeSeriesSplit, not random splits
6. **Target leakage**: Features shouldn't contain information not available at prediction time
