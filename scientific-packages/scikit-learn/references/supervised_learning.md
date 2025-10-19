# Supervised Learning in scikit-learn

## Overview
Supervised learning algorithms learn patterns from labeled training data to make predictions on new data. Scikit-learn organizes supervised learning into 17 major categories.

## Linear Models

### Regression
- **LinearRegression**: Ordinary least squares regression
- **Ridge**: L2-regularized regression, good for multicollinearity
- **Lasso**: L1-regularized regression, performs feature selection
- **ElasticNet**: Combined L1/L2 regularization
- **LassoLars**: Lasso using Least Angle Regression algorithm
- **BayesianRidge**: Bayesian approach with automatic relevance determination

### Classification
- **LogisticRegression**: Binary and multiclass classification
- **RidgeClassifier**: Ridge regression for classification
- **SGDClassifier**: Linear classifiers with SGD training

**Use cases**: Baseline models, interpretable predictions, high-dimensional data, when linear relationships are expected

**Key parameters**:
- `alpha`: Regularization strength (higher = more regularization)
- `fit_intercept`: Whether to calculate intercept
- `solver`: Optimization algorithm ('lbfgs', 'saga', 'liblinear')

## Support Vector Machines (SVM)

- **SVC**: Support Vector Classification
- **SVR**: Support Vector Regression
- **LinearSVC**: Linear SVM using liblinear (faster for large datasets)
- **OneClassSVM**: Unsupervised outlier detection

**Use cases**: Complex non-linear decision boundaries, high-dimensional spaces, when clear margin of separation exists

**Key parameters**:
- `kernel`: 'linear', 'poly', 'rbf', 'sigmoid'
- `C`: Regularization parameter (lower = more regularization)
- `gamma`: Kernel coefficient ('scale', 'auto', or float)
- `degree`: Polynomial degree (for poly kernel)

**Performance tip**: SVMs don't scale well beyond tens of thousands of samples. Use LinearSVC for large datasets with linear kernel.

## Decision Trees

- **DecisionTreeClassifier**: Classification tree
- **DecisionTreeRegressor**: Regression tree
- **ExtraTreeClassifier/Regressor**: Extremely randomized tree

**Use cases**: Non-linear relationships, feature importance analysis, interpretable rules, handling mixed data types

**Key parameters**:
- `max_depth`: Maximum tree depth (controls overfitting)
- `min_samples_split`: Minimum samples to split a node
- `min_samples_leaf`: Minimum samples in leaf node
- `max_features`: Number of features to consider for splits
- `criterion`: 'gini', 'entropy' (classification); 'squared_error', 'absolute_error' (regression)

**Overfitting prevention**: Limit `max_depth`, increase `min_samples_split/leaf`, use pruning with `ccp_alpha`

## Ensemble Methods

### Random Forests
- **RandomForestClassifier**: Ensemble of decision trees
- **RandomForestRegressor**: Regression variant

**Use cases**: Robust general-purpose algorithm, reduces overfitting vs single trees, handles non-linear relationships

**Key parameters**:
- `n_estimators`: Number of trees (higher = better but slower)
- `max_depth`: Maximum tree depth
- `max_features`: Features per split ('sqrt', 'log2', int, float)
- `bootstrap`: Whether to use bootstrap samples
- `n_jobs`: Parallel processing (-1 uses all cores)

### Gradient Boosting
- **HistGradientBoostingClassifier/Regressor**: Histogram-based, fast for large datasets (>10k samples)
- **GradientBoostingClassifier/Regressor**: Traditional implementation, better for small datasets

**Use cases**: High-performance predictions, winning Kaggle competitions, structured/tabular data

**Key parameters**:
- `n_estimators`: Number of boosting stages
- `learning_rate`: Shrinks contribution of each tree
- `max_depth`: Maximum tree depth (typically 3-8)
- `subsample`: Fraction of samples per tree (enables stochastic gradient boosting)
- `early_stopping`: Stop when validation score stops improving

**Performance tip**: HistGradientBoosting is orders of magnitude faster for large datasets

### AdaBoost
- **AdaBoostClassifier/Regressor**: Adaptive boosting

**Use cases**: Boosting weak learners, less prone to overfitting than other methods

**Key parameters**:
- `estimator`: Base estimator (default: DecisionTreeClassifier with max_depth=1)
- `n_estimators`: Number of boosting iterations
- `learning_rate`: Weight applied to each classifier

### Bagging
- **BaggingClassifier/Regressor**: Bootstrap aggregating with any base estimator

**Use cases**: Reducing variance of unstable models, parallel ensemble creation

**Key parameters**:
- `estimator`: Base estimator to fit
- `n_estimators`: Number of estimators
- `max_samples`: Samples to draw per estimator
- `bootstrap`: Whether to use replacement

### Voting & Stacking
- **VotingClassifier/Regressor**: Combines different model types
- **StackingClassifier/Regressor**: Meta-learner trained on base predictions

**Use cases**: Combining diverse models, leveraging different model strengths

## Neural Networks

- **MLPClassifier**: Multi-layer perceptron classifier
- **MLPRegressor**: Multi-layer perceptron regressor

**Use cases**: Complex non-linear patterns, when gradient boosting is too slow, deep feature learning

**Key parameters**:
- `hidden_layer_sizes`: Tuple of hidden layer sizes (e.g., (100, 50))
- `activation`: 'relu', 'tanh', 'logistic'
- `solver`: 'adam', 'lbfgs', 'sgd'
- `alpha`: L2 regularization term
- `learning_rate`: Learning rate schedule
- `early_stopping`: Stop when validation score stops improving

**Important**: Feature scaling is critical for neural networks. Always use StandardScaler or similar.

## Nearest Neighbors

- **KNeighborsClassifier/Regressor**: K-nearest neighbors
- **RadiusNeighborsClassifier/Regressor**: Radius-based neighbors
- **NearestCentroid**: Classification using class centroids

**Use cases**: Simple baseline, irregular decision boundaries, when interpretability isn't critical

**Key parameters**:
- `n_neighbors`: Number of neighbors (typically 3-11)
- `weights`: 'uniform' or 'distance' (distance-weighted voting)
- `metric`: Distance metric ('euclidean', 'manhattan', 'minkowski')
- `algorithm`: 'auto', 'ball_tree', 'kd_tree', 'brute'

## Naive Bayes

- **GaussianNB**: Assumes Gaussian distribution of features
- **MultinomialNB**: For discrete counts (text classification)
- **BernoulliNB**: For binary/boolean features
- **CategoricalNB**: For categorical features
- **ComplementNB**: Adapted for imbalanced datasets

**Use cases**: Text classification, fast baseline, when features are independent, small training sets

**Key parameters**:
- `alpha`: Smoothing parameter (Laplace/Lidstone smoothing)
- `fit_prior`: Whether to learn class prior probabilities

## Linear/Quadratic Discriminant Analysis

- **LinearDiscriminantAnalysis**: Linear decision boundary with dimensionality reduction
- **QuadraticDiscriminantAnalysis**: Quadratic decision boundary

**Use cases**: When classes have Gaussian distributions, dimensionality reduction, when covariance assumptions hold

## Gaussian Processes

- **GaussianProcessClassifier**: Probabilistic classification
- **GaussianProcessRegressor**: Probabilistic regression with uncertainty estimates

**Use cases**: When uncertainty quantification is important, small datasets, smooth function approximation

**Key parameters**:
- `kernel`: Covariance function (RBF, Matern, RationalQuadratic, etc.)
- `alpha`: Noise level

**Limitation**: Doesn't scale well to large datasets (O(n³) complexity)

## Stochastic Gradient Descent

- **SGDClassifier**: Linear classifiers with SGD
- **SGDRegressor**: Linear regressors with SGD

**Use cases**: Very large datasets (>100k samples), online learning, when data doesn't fit in memory

**Key parameters**:
- `loss`: Loss function ('hinge', 'log_loss', 'squared_error', etc.)
- `penalty`: Regularization ('l2', 'l1', 'elasticnet')
- `alpha`: Regularization strength
- `learning_rate`: Learning rate schedule

## Semi-Supervised Learning

- **SelfTrainingClassifier**: Self-training with any base classifier
- **LabelPropagation**: Label propagation through graph
- **LabelSpreading**: Label spreading (modified label propagation)

**Use cases**: When labeled data is scarce but unlabeled data is abundant

## Feature Selection

- **VarianceThreshold**: Remove low-variance features
- **SelectKBest**: Select K highest scoring features
- **SelectPercentile**: Select top percentile of features
- **RFE**: Recursive feature elimination
- **RFECV**: RFE with cross-validation
- **SelectFromModel**: Select features based on importance
- **SequentialFeatureSelector**: Forward/backward feature selection

**Use cases**: Reducing dimensionality, removing irrelevant features, improving interpretability, reducing overfitting

## Probability Calibration

- **CalibratedClassifierCV**: Calibrate classifier probabilities

**Use cases**: When probability estimates are important (not just class predictions), especially with SVM and Naive Bayes

**Methods**:
- `sigmoid`: Platt scaling
- `isotonic`: Isotonic regression (more flexible, needs more data)

## Multi-Output Methods

- **MultiOutputClassifier**: Fit one classifier per target
- **MultiOutputRegressor**: Fit one regressor per target
- **ClassifierChain**: Models dependencies between targets
- **RegressorChain**: Regression variant

**Use cases**: Predicting multiple related targets simultaneously

## Specialized Regression

- **IsotonicRegression**: Monotonic regression
- **QuantileRegressor**: Quantile regression for prediction intervals

## Algorithm Selection Guidelines

**Start with**:
1. **Logistic Regression** (classification) or **LinearRegression/Ridge** (regression) as baseline
2. **RandomForestClassifier/Regressor** for general non-linear problems
3. **HistGradientBoostingClassifier/Regressor** when best performance is needed

**Consider dataset size**:
- Small (<1k samples): SVM, Gaussian Processes, any algorithm
- Medium (1k-100k): Random Forests, Gradient Boosting, Neural Networks
- Large (>100k): SGD, HistGradientBoosting, LinearSVC

**Consider interpretability needs**:
- High interpretability: Linear models, Decision Trees, Naive Bayes
- Medium: Random Forests (feature importance), Rule extraction
- Low (black box acceptable): Gradient Boosting, Neural Networks, SVM with RBF kernel

**Consider training time**:
- Fast: Linear models, Naive Bayes, Decision Trees
- Medium: Random Forests (parallelizable), SVM (small data)
- Slow: Gradient Boosting, Neural Networks, SVM (large data), Gaussian Processes
