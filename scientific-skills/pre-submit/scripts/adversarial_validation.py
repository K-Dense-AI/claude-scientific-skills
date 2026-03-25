"""
Adversarial validation — detect distribution shift between train and test sets.

Usage:
    uv run python <skill-path>/scripts/adversarial_validation.py \
        --train data/train.csv --test data/test.csv --target target_column \
        [--n-folds 5] [--clf rf|lgbm|et]

Outputs:
    - Adversarial AUC (>0.55 = shift; >0.65 = investigate)
    - Top-20 leaking features (drop or transform these)
    - Sample weights file: artifacts/adversarial_weights.npy
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def run_adversarial_validation(
    train_path: str,
    test_path: str,
    target_col: str,
    n_folds: int = 5,
    save_weights: bool = True,
    clf: Any = None,
) -> dict:
    """
    Detect distribution shift between train and test sets.

    Parameters
    ----------
    train_path : str
        Path to the train CSV file.
    test_path : str
        Path to the test CSV file.
    target_col : str
        Target column name to drop from train before combining.
    n_folds : int
        Number of StratifiedKFold folds for AUC estimation.
    save_weights : bool
        If True, save sample weights to artifacts/adversarial_weights.npy.
    clf : sklearn-compatible classifier, optional
        Any object implementing fit(X, y) and predict_proba(X).
        Defaults to LGBMClassifier(n_estimators=300, random_state=42, verbose=-1).
    """
    try:
        from sklearn.metrics import roc_auc_score
        from sklearn.model_selection import StratifiedKFold
    except ImportError:
        raise SystemExit("uv add scikit-learn")

    if clf is None:
        try:
            from lightgbm import LGBMClassifier
        except ImportError:
            raise SystemExit("uv add lightgbm")
        clf = LGBMClassifier(n_estimators=300, random_state=42, verbose=-1)

    train = pd.read_csv(train_path)
    test = pd.read_csv(test_path)

    train_adv = train.drop(columns=[target_col], errors="ignore")

    combined = pd.concat(
        [train_adv.assign(_is_test=0), test.assign(_is_test=1)],
        ignore_index=True,
    )

    feats = (
        combined.select_dtypes(include=[np.number])
        .drop(columns=["_is_test"], errors="ignore")
        .columns.tolist()
    )
    X = combined[feats].fillna(-999)
    y = combined["_is_test"]

    cv = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=42)
    aucs: list[float] = []

    for tr, val in cv.split(X, y):
        clf.fit(X.iloc[tr], y.iloc[tr])
        aucs.append(roc_auc_score(y.iloc[val], clf.predict_proba(X.iloc[val])[:, 1]))

    # Final fit on full combined set for importance + weights
    clf.fit(X, y)
    auc = float(np.mean(aucs))

    importances = pd.Series(clf.feature_importances_, index=feats).sort_values(
        ascending=False
    )

    # AUC verdict
    if auc < 0.55:
        verdict = "✅ No shift — proceed normally"
    elif auc < 0.65:
        verdict = "⚠️  Mild shift — check top features; monitor LB-OOF gap"
    elif auc < 0.80:
        verdict = "❌ Moderate shift — drop or transform top leaking features"
    else:
        verdict = "🚨 Severe shift — likely ID/time leak — investigate immediately"

    print(f"\nAdversarial AUC: {auc:.4f}  →  {verdict}")
    print(f"\nTop-20 leaking features:\n{importances.head(20).to_string()}")

    result = {
        "auc": auc,
        "verdict": verdict,
        "top_features": importances.head(20).to_dict(),
    }

    if save_weights:
        # Weight train samples by P(is_test | features) — emphasise train rows
        # that "look like" test rows.
        train_X = X.iloc[: len(train_adv)]
        proba = clf.predict_proba(train_X)[:, 1]
        weights = np.clip(proba / (1 - proba + 1e-6), 0.1, 10.0)
        Path("artifacts").mkdir(exist_ok=True)
        np.save("artifacts/adversarial_weights.npy", weights)
        print("\nSample weights saved to artifacts/adversarial_weights.npy")
        result["weights_path"] = "artifacts/adversarial_weights.npy"

    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--train", required=True, help="Path to train CSV file")
    parser.add_argument("--test", required=True, help="Path to test CSV file")
    parser.add_argument("--target", required=True, help="Target column name")
    parser.add_argument("--n-folds", default=5, type=int)
    parser.add_argument(
        "--clf",
        default="lgbm",
        choices=["lgbm", "rf", "et"],
        help="Classifier: lgbm (default), rf (RandomForest), et (ExtraTrees)",
    )
    args = parser.parse_args()

    _clf = None
    if args.clf == "rf":
        from sklearn.ensemble import RandomForestClassifier
        _clf = RandomForestClassifier(n_estimators=300, random_state=42, n_jobs=-1)
    elif args.clf == "et":
        from sklearn.ensemble import ExtraTreesClassifier
        _clf = ExtraTreesClassifier(n_estimators=300, random_state=42, n_jobs=-1)
    # lgbm → None → function uses its default LGBMClassifier

    result = run_adversarial_validation(args.train, args.test, args.target, args.n_folds, clf=_clf)
    print("\nResult summary:")
    print(
        json.dumps({k: v for k, v in result.items() if k != "top_features"}, indent=2)
    )
