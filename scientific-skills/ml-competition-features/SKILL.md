---
name: ml-competition-features
description: "Build, review, and debug feature engineering and cross-validation pipelines for tabular ML competitions. Use when: implementing or caching engineered features; choosing and implementing the correct CV split strategy (GroupKFold, StratifiedKFold, TimeSeriesSplit, KFold); preventing target-leakage in fold-local encodings; structuring OOF accumulation arrays; diagnosing train/LB gap caused by leakage. NOT for model training, metrics, or tuning."
argument-hint: "Describe your task: e.g. 'implement target encoding without leakage', 'choose correct CV split for grouped data', 'debug OOF leakage', 'cache feature pipeline'"
license: MIT
metadata:
    skill-author: eak
---

# ML Competition — Features & Validation Strategy

## Overview

This skill covers the two data-layer concerns that must be correct before any model is worth running:

1. **Feature engineering** — encoding strategies, datetime features, aggregations, feature selection, and the versioned pkl cache discipline that prevents stale features from silently degrading results
2. **Validation strategy** — choosing the right CV split, accumulating OOF arrays correctly, preventing target leakage, and diagnosing leakage when train OOF ≫ LB

**Core invariants:**
- The feature cache version must be bumped on every change — stale cache is the #1 source of silent regressions
- Target encoding and fold-local statistics must always be computed **inside** cross-validation — computing them on the full dataset before splitting is leakage
- OOF arrays hold **training rows only**, never test rows

---

## Validation Strategy — Critical Rules

### CV split selection

| Data condition | Split to use |
|----------------|-------------|
| Natural integrity unit (user ID, entity ID, session ID) | `GroupKFold` |
| i.i.d. rows with class imbalance, no group column | `StratifiedKFold` |
| Rows are truly independent and target is balanced | `KFold` |
| Temporal ordering matters | `TimeSeriesSplit` or rolling-window split — **never GroupKFold** |

### OOF accumulation

```python
oof = np.zeros(len(train_df))           # allocate once, full training length
for fold, (tr_idx, va_idx) in enumerate(splits):
    model.fit(X[tr_idx], y[tr_idx])
    oof[va_idx] = model.predict_proba(X[va_idx])[:, 1]  # assign by index, not stack
```

- **Always collect full arrays, never fold-by-fold averages** — `oof[val] = preds`, not `oof_list.append(preds.mean())`
- **OOF arrays are train-rows only** — even if combined (train + auxiliary) data drives fold splits, only score `oof[va_idx[va_idx < n_train]]`

### Target encoding and fold-local statistics

```python
# ✅ CORRECT — computed inside fold loop, only on training portion
for fold, (tr_idx, va_idx) in enumerate(splits):
    te = train_df.iloc[tr_idx].groupby("cat_col")["target"].mean()
    train_df.loc[va_idx, "cat_te"] = train_df.loc[va_idx, "cat_col"].map(te)

# ❌ WRONG — computed before split, leaks validation targets into training features
train_df["cat_te"] = train_df.groupby("cat_col")["target"].transform("mean")
```

### Tuner vs final training folds

- Tuner folds: 3-fold for speed
- Final training: 5-fold × 5 seeds
- Both must use **identical CV strategy and same group column**

### Leakage checklist — before finalizing any feature

1. Is this feature computed using future information relative to the prediction time?
2. Does this feature use any information from the validation/test fold?
3. Is this a fold-local statistic computed after the split or before it?
4. Does OOF score dramatically exceed LB? (>+0.01 in AUC = suspect leakage)

---

## Feature Cache Discipline

```python
FEAT_CACHE = cfg.feat_cache   # e.g. "cache/features_v3.pkl"

def load_or_build_features(df):
    if os.path.exists(FEAT_CACHE):
        return pd.read_pickle(FEAT_CACHE)
    feats = engineer_features(df)
    feats.to_pickle(FEAT_CACHE)
    return feats
```

**Hard rules:**
- Bump the version number in `config.yaml` (`feat_cache: cache/features_v3.pkl`) on every change to `engineer_features()`
- Never patch a cached file in-place — delete old cache or bump version
- `build_model_matrices()` is called once per process and cached in a module-level variable; do not call it inside the fold loop

---

## Reference Files

| File | What it covers |
|------|----------------|
| [feature-engineering.md](./references/feature-engineering.md) | Encoding strategies, datetime features, aggregations, feature selection, cache discipline |
| [validation-strategy.md](./references/validation-strategy.md) | GroupKFold / TimeSeriesSplit, OOF accumulation, leakage prevention, leakage checklist |

---

## See Also

| Skill | When to use it instead |
|-------|------------------------|
| `ml-competition` | Full pipeline overview, task type decision guide, first-principles checklist |
| `ml-competition-setup` | Project structure, RunConfig, process management |
| `ml-competition-training` | Model training, competition metrics, correct output format |
| `ml-competition-tuning` | Optuna hyperparameter tuning |
| `ml-competition-advanced` | Pseudo-labeling, ensemble, post-processing, experiment tracking |
| `ml-competition-quality` | Coding rules, common pitfalls |
