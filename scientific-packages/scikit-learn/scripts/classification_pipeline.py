#!/usr/bin/env python3
"""
Complete classification pipeline with preprocessing, training, evaluation, and hyperparameter tuning.
Demonstrates best practices for scikit-learn workflows.
"""

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.impute import SimpleImputer
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score
import joblib


def create_preprocessing_pipeline(numeric_features, categorical_features):
    """
    Create preprocessing pipeline for mixed data types.

    Args:
        numeric_features: List of numeric column names
        categorical_features: List of categorical column names

    Returns:
        ColumnTransformer with appropriate preprocessing for each data type
    """
    numeric_transformer = Pipeline(steps=[
        ('imputer', SimpleImputer(strategy='median')),
        ('scaler', StandardScaler())
    ])

    categorical_transformer = Pipeline(steps=[
        ('imputer', SimpleImputer(strategy='constant', fill_value='missing')),
        ('onehot', OneHotEncoder(handle_unknown='ignore', sparse_output=True))
    ])

    preprocessor = ColumnTransformer(
        transformers=[
            ('num', numeric_transformer, numeric_features),
            ('cat', categorical_transformer, categorical_features)
        ])

    return preprocessor


def create_full_pipeline(preprocessor, classifier=None):
    """
    Create complete ML pipeline with preprocessing and classification.

    Args:
        preprocessor: Preprocessing ColumnTransformer
        classifier: Classifier instance (default: RandomForestClassifier)

    Returns:
        Complete Pipeline
    """
    if classifier is None:
        classifier = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)

    pipeline = Pipeline(steps=[
        ('preprocessor', preprocessor),
        ('classifier', classifier)
    ])

    return pipeline


def evaluate_model(pipeline, X_train, y_train, X_test, y_test, cv=5):
    """
    Evaluate model using cross-validation and test set.

    Args:
        pipeline: Trained pipeline
        X_train, y_train: Training data
        X_test, y_test: Test data
        cv: Number of cross-validation folds

    Returns:
        Dictionary with evaluation results
    """
    # Cross-validation on training set
    cv_scores = cross_val_score(pipeline, X_train, y_train, cv=cv, scoring='accuracy')

    # Test set evaluation
    y_pred = pipeline.predict(X_test)
    test_score = pipeline.score(X_test, y_test)

    # Get probabilities if available
    try:
        y_proba = pipeline.predict_proba(X_test)
        if len(np.unique(y_test)) == 2:
            # Binary classification
            auc = roc_auc_score(y_test, y_proba[:, 1])
        else:
            # Multiclass
            auc = roc_auc_score(y_test, y_proba, multi_class='ovr')
    except:
        auc = None

    results = {
        'cv_mean': cv_scores.mean(),
        'cv_std': cv_scores.std(),
        'test_score': test_score,
        'auc': auc,
        'classification_report': classification_report(y_test, y_pred),
        'confusion_matrix': confusion_matrix(y_test, y_pred)
    }

    return results


def tune_hyperparameters(pipeline, X_train, y_train, param_grid, cv=5):
    """
    Perform hyperparameter tuning using GridSearchCV.

    Args:
        pipeline: Pipeline to tune
        X_train, y_train: Training data
        param_grid: Dictionary of parameters to search
        cv: Number of cross-validation folds

    Returns:
        GridSearchCV object with best model
    """
    grid_search = GridSearchCV(
        pipeline,
        param_grid,
        cv=cv,
        scoring='f1_weighted',
        n_jobs=-1,
        verbose=1
    )

    grid_search.fit(X_train, y_train)

    print(f"Best parameters: {grid_search.best_params_}")
    print(f"Best CV score: {grid_search.best_score_:.3f}")

    return grid_search


def main():
    """
    Example usage of the classification pipeline.
    """
    # Load your data here
    # X, y = load_data()

    # Example with synthetic data
    from sklearn.datasets import make_classification
    X, y = make_classification(
        n_samples=1000,
        n_features=20,
        n_informative=15,
        n_redundant=5,
        random_state=42
    )

    # Convert to DataFrame for demonstration
    feature_names = [f'feature_{i}' for i in range(X.shape[1])]
    X = pd.DataFrame(X, columns=feature_names)

    # Split features into numeric and categorical (all numeric in this example)
    numeric_features = feature_names
    categorical_features = []

    # Split data (use stratify for imbalanced classes)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )

    # Create preprocessing pipeline
    preprocessor = create_preprocessing_pipeline(numeric_features, categorical_features)

    # Create full pipeline
    pipeline = create_full_pipeline(preprocessor)

    # Train model
    print("Training model...")
    pipeline.fit(X_train, y_train)

    # Evaluate model
    print("\nEvaluating model...")
    results = evaluate_model(pipeline, X_train, y_train, X_test, y_test)

    print(f"CV Accuracy: {results['cv_mean']:.3f} (+/- {results['cv_std']:.3f})")
    print(f"Test Accuracy: {results['test_score']:.3f}")
    if results['auc']:
        print(f"ROC-AUC: {results['auc']:.3f}")
    print("\nClassification Report:")
    print(results['classification_report'])

    # Hyperparameter tuning (optional)
    print("\nTuning hyperparameters...")
    param_grid = {
        'classifier__n_estimators': [100, 200],
        'classifier__max_depth': [10, 20, None],
        'classifier__min_samples_split': [2, 5]
    }

    grid_search = tune_hyperparameters(pipeline, X_train, y_train, param_grid)

    # Evaluate best model
    print("\nEvaluating tuned model...")
    best_pipeline = grid_search.best_estimator_
    y_pred = best_pipeline.predict(X_test)
    print(classification_report(y_test, y_pred))

    # Save model
    print("\nSaving model...")
    joblib.dump(best_pipeline, 'best_model.pkl')
    print("Model saved as 'best_model.pkl'")


if __name__ == "__main__":
    main()
