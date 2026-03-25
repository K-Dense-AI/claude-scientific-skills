# Feature Engineering

## Overview

Feature engineering is the highest-leverage activity in tabular competitions — a good feature can outperform days of hyperparameter tuning. This file covers the architecture rules for `base/features.py`, all encoding strategies (ordinal, frequency, one-hot, target encoding), datetime and interaction features, aggregation patterns, and feature selection techniques.

**Three non-negotiable rules** govern every feature in this pipeline:
1. **Deterministic** — same inputs always produce the same outputs; no randomness in `engineer_features()`
2. **Versioned** — bump the cache file name (`features_v1.pkl` → `features_v2.pkl`) every time logic changes; stale cache is the most common source of unexplained OOF changes
3. **Leak-free** — target-based statistics (target encoding, aggregations over the target) must be computed fold-by-fold inside cross-validation, never on the full training set

**When to use:** When adding, modifying, or debugging features. Also consult before any cache-related confusion ("why did my OOF drop after I changed the feature?").

---

## Architecture Rule

Feature engineering belongs in `base/features.py`. It must be:
1. **Deterministic** — same inputs always produce same outputs
2. **Versioned** — bump the cache filename every time logic changes
3. **Leak-free** — aggregations over the target must be computed fold-by-fold

The output of `engineer_features()` is cached in `cache/features_vN.pkl`. The output of `build_model_matrices()` is an in-memory process-level cache (rebuilt each run).

---

## Encoding Strategies

### Ordinal / Label Encoding
```python
from sklearn.preprocessing import OrdinalEncoder

enc = OrdinalEncoder(handle_unknown="use_encoded_value", unknown_value=-1)
train[cat_cols] = enc.fit_transform(train[cat_cols])
test[cat_cols]  = enc.transform(test[cat_cols])
```
Use for: tree models (CatBoost handles natively; LGB/XGB need numeric indices).

### Frequency Encoding
```python
def freq_encode(train, test, col):
    freq = train[col].value_counts(normalize=True)
    train[f"{col}_freq"] = train[col].map(freq).fillna(0)
    test[f"{col}_freq"]  = test[col].map(freq).fillna(0)
```
Use for: high-cardinality categoricals where target correlation is weak.

### One-Hot Encoding
```python
train = pd.get_dummies(train, columns=low_card_cols, drop_first=True)
test  = test.reindex(columns=train.columns, fill_value=0)
```
Use for: low-cardinality (<15 categories) nominal features, NN inputs.

### Target Encoding (fold-by-fold — avoids leakage)
```python
def smooth_target_encode(train_fold, val_fold, col, target, smooth_k=20):
    global_mean = train_fold[target].mean()
    stats = train_fold.groupby(col)[target].agg(["mean", "count"])
    stats["encoded"] = (
        (stats["mean"] * stats["count"] + global_mean * smooth_k)
        / (stats["count"] + smooth_k)
    )
    return val_fold[col].map(stats["encoded"]).fillna(global_mean)

# Usage inside CV loop — NEVER computed globally
for fold, (tr_idx, va_idx) in enumerate(gkf.split(X, groups=groups)):
    train.loc[va_idx, f"{col}_te"] = smooth_target_encode(
        train.iloc[tr_idx], train.iloc[va_idx], col, target
    )
```

---

## DateTime Features

```python
def add_datetime_features(df, col):
    dt = pd.to_datetime(df[col])
    df[f"{col}_year"]       = dt.dt.year
    df[f"{col}_month"]      = dt.dt.month
    df[f"{col}_dayofweek"]  = dt.dt.dayofweek
    df[f"{col}_hour"]       = dt.dt.hour
    df[f"{col}_is_weekend"]  = (dt.dt.dayofweek >= 5).astype(int)
    df[f"{col}_quarter"]    = dt.dt.quarter
    # Cyclical encoding for periodic features (avoids discontinuity at boundaries)
    df[f"{col}_month_sin"]  = np.sin(2 * np.pi * dt.dt.month / 12)
    df[f"{col}_month_cos"]  = np.cos(2 * np.pi * dt.dt.month / 12)
    return df
```

**Age / elapsed time features:**
```python
df["days_since_event"] = (reference_date - pd.to_datetime(df["event_date"])).dt.days
df["account_age_days"] = (df["snapshot_date"] - df["signup_date"]).dt.days
```

---

## Aggregation Features (Group Statistics)

```python
def add_group_stats(train, test, group_col, num_cols):
    """Compute group-level statistics on train; apply to test."""
    stats = train.groupby(group_col)[num_cols].agg(["mean", "std", "min", "max", "count"])
    stats.columns = [f"{g}_{c}_{s}" for g, c in zip(
        [group_col] * len(stats.columns), stats.columns.get_level_values(0)
    ) for s in stats.columns.get_level_values(1)]
    stats = stats.reset_index()

    train = train.merge(stats, on=group_col, how="left")
    test  = test.merge(stats,  on=group_col, how="left")
    return train, test
```

**Deviation features:**
```python
df["val_minus_group_mean"] = df["val"] - df["group_mean_val"]
df["val_zscore_in_group"]  = (df["val"] - df["group_mean_val"]) / (df["group_std_val"] + 1e-6)
```

---

## Interaction Features

```python
# Ratio
df["ratio_a_b"] = df["a"] / (df["b"] + 1e-6)

# Product
df["product_a_b"] = df["a"] * df["b"]

# Binned numerical → crossed categorical
df["a_bin"] = pd.qcut(df["a"], q=10, labels=False, duplicates="drop")
df["b_bin"] = pd.qcut(df["b"], q=10, labels=False, duplicates="drop")
df["a_x_b"] = df["a_bin"].astype(str) + "_" + df["b_bin"].astype(str)
```

---

## Text Features

```python
from sklearn.feature_extraction.text import TfidfVectorizer
from scipy.sparse import hstack

# TF-IDF on a text column
tfidf = TfidfVectorizer(max_features=500, ngram_range=(1, 2), min_df=5)
X_text_train = tfidf.fit_transform(train["text_col"].fillna(""))
X_text_test  = tfidf.transform(test["text_col"].fillna(""))

# Concatenate with numeric features before training
X_train_full = hstack([X_num_train, X_text_train])
```

**Sentence embeddings (for NN):**
```python
# Use sentence-transformers for semantic embeddings
from sentence_transformers import SentenceTransformer
model_st = SentenceTransformer("all-MiniLM-L6-v2")
embeddings = model_st.encode(df["text_col"].tolist(), batch_size=256, show_progress_bar=True)
# Cache these — they're expensive to compute
```

---

## Imputation

```python
from sklearn.impute import SimpleImputer, KNNImputer

# Simple: fill with median (safe, fast)
imp = SimpleImputer(strategy="median")
X_train = imp.fit_transform(X_train)
X_test  = imp.transform(X_test)

# Flag missing values before imputing — the missingness pattern may be informative
for col in df.columns:
    if df[col].isnull().any():
        df[f"{col}_is_missing"] = df[col].isnull().astype(int)
```

---

## Feature Selection

### Importance-based (remove useless features)
```python
# After training a quick LGB model:
importances = pd.Series(model.feature_importance(), index=feature_names)
low_imp_cols = importances[importances < 1].index.tolist()
# Drop or keep — always validate on OOF before committing
```

### Correlation filter (remove near-duplicate features)
```python
corr_matrix = pd.DataFrame(X_train).corr().abs()
upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
to_drop = [col for col in upper.columns if any(upper[col] > 0.98)]
```

### Permutation importance (most honest)
```python
from sklearn.inspection import permutation_importance

result = permutation_importance(trained_model, X_val, y_val, n_repeats=5, random_state=42)
# Use result.importances_mean — negative means feature hurts
```

---

## Common Pitfalls

| Pitfall | Fix |
|---------|-----|
| Target encoding computed globally | Compute fold-by-fold inside CV loop only |
| Aggregations use test rows | Fit aggregations on train only; apply to test |
| DateTime features using future data | Only use data available at prediction time |
| Feature cache not bumped after change | Increment `vN` in `feat_cache` config immediately |
| Normalization fit on all data | `fit_transform` on train, `transform` on test/val |
| Interaction features before selection | Explosion of features — select first, interact after |

---

## See Also

| File | Why |
|------|-----|
| [validation-strategy.md](./validation-strategy.md) | Fold splits that govern where target encoding is computed |
| [model-training.md](./model-training.md) | `build_model_matrices()` consumes the output of `engineer_features()` |
| [project-structure.md](./project-structure.md) | Cache file naming conventions and `feat_cache` config key |
| [common-pitfalls.md](./common-pitfalls.md) | Pitfall #8 — feature cache version not bumped after logic change |
