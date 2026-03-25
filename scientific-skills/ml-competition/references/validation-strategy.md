# Validation Strategy

## Overview

Validation strategy is the foundation of every competition pipeline. The choice of cross-validation scheme controls whether your OOF score is an honest proxy for the leaderboard or a misleading in-sample estimate. A wrong split (e.g., `KFold` when a group column exists, or `StratifiedKFold` on time-ordered data) can produce OOF scores that are 0.01–0.05 higher than the actual LB — making every other improvement invisible.

This file covers: choosing the right CV method for your competition type; correct `GroupKFold`, `StratifiedKFold`, and `TimeSeriesSplit` setup; OOF array accumulation and saving; target encoding inside folds (leakage prevention); and the `leakage_check` utility.

**When to use:** Before writing a single line of training code. The CV scheme must be decided and locked before feature engineering, training, or tuning begins — changing it mid-competition invalidates all previous OOF comparisons.

---

## Choose the Right CV Strategy First

| Competition type | CV method | When to use |
|-----------------|-----------|-------------|
| Tabular — entity groups (user, farm, session) | `GroupKFold` | Rows from the same entity must not split across folds |
| Time-series — ordered data | `TimeSeriesSplit` / rolling window | Future leaks into past; temporal ordering matters |
| Tabular — no groups, rare positives | `StratifiedKFold` | Class imbalance but genuinely no groups |
| Tabular — no groups, balanced | `KFold` | Only when you are certain there are no groups |

**Default choice: GroupKFold.** When in doubt, look for an ID column. If one exists, use it as the group.

---

## GroupKFold Setup

```python
from sklearn.model_selection import GroupKFold

gkf = GroupKFold(n_splits=cfg.n_folds)

# groups must be the group column values aligned to the data being split
for fold, (tr_idx, va_idx) in enumerate(gkf.split(X, groups=groups)):
    # tr_idx, va_idx are indices into X
    x_tr, x_va = X[tr_idx], X[va_idx]
```

**Groups should be:**
- Stable identifiers: `user_id`, `farm_id`, `session_id`, `customer_id`
- Never derived from the target — that's leakage
- Ideally balanced: GroupKFold doesn't balance class distribution, only groups

---

## TimeSeriesSplit (temporal competitions)

```python
from sklearn.model_selection import TimeSeriesSplit

tscv = TimeSeriesSplit(n_splits=5, gap=0)   # gap: rows to skip between train/val

for fold, (tr_idx, va_idx) in enumerate(tscv.split(X)):
    x_tr, x_va = X[tr_idx], X[va_idx]
    # Training always uses past, val uses future — no temporal leakage
    # NOTE: early folds have much less training data; scores will be noisier
```

**Time-series pitfalls:**
- Never use GroupKFold on time-series data \u2014 it can mix past and future rows in the same fold
- Feature engineering must use only data available at prediction time (no future aggregates)
- Use expanding window (growing train set) for stable estimates; sliding window to simulate recency

---

## OOF Array Management

### Initialisation
```python
n_train = len(train)
oof = {t: np.zeros(n_train) for t in TARGETS}
tst = {t: np.zeros(len(test)) for t in TARGETS}
```

### Accumulation over seeds (not seeds × folds)
```python
# Each train row appears in exactly one fold's val set across folds.
# So each row gets predictions from exactly `len(seeds)` models.
oof[t][va_idx] += model.predict(x_va) / len(seeds)

# Test rows appear in ALL folds' train → average over folds × seeds
tst[t] += model.predict(x_test) / (n_folds * len(seeds))
```

### Saving OOF
```python
# base/common.py save_oof()
def save_oof(oof: dict, tst: dict, model_name: str, tag: str = "base"):
    path = os.path.join(cfg.oof_dir, f"{tag}_{model_name}.pkl")
    os.makedirs(cfg.oof_dir, exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump({"oof": oof, "test": tst}, f)
    logger.info(f"  [oof] saved → {path}")
```

**OOF file naming:**
- Base model: `oof/{tag}_{model}.pkl` → e.g. `oof/base_cat.pkl`
- Pseudo model: `oof/pseudo_{model}.pkl`
- Always save both `oof` (train preds) and `test` (test preds)

---

## Using Combined Data (train + prior / auxiliary)

When auxiliary data is available with known labels, it can augment training — but the validation set must stay clean.

```python
# Correct: split on combined, restrict OOF to train rows
n_train = len(train)
X_combined = np.vstack([X_train, X_aux])      # n_train + n_aux rows
groups_combined = np.concatenate([groups, groups_aux])

for fold, (tr_idx, va_idx) in enumerate(gkf.split(X_combined, groups=groups_combined)):
    # Train on combined
    x_tr = X_combined[tr_idx]
    y_tr = y_combined[target][tr_idx]
    
    # Val ONLY on original train rows
    va_train_idx = va_idx[va_idx < n_train]   # ← CRITICAL
    x_va = X_combined[va_train_idx]
    y_va = y_train[target][va_train_idx]
    
    # Accumulate OOF only on train rows
    oof[target][va_train_idx] += model.predict(x_va) / len(seeds)
```

**When to use combined data:**
- ✅ CatBoost with categorical features: robust to distribution shift
- ⚠️ LGB/XGB/NN: test carefully — prior data from a different distribution can **cause severe overfitting**. Always compare OOF score vs LB score before committing.

---

## Target Encoding — Leakage Prevention

Target encoding **must** be computed fold-by-fold, using only training fold labels.

```python
# WRONG — computed globally, leaks val targets
train["te_feature"] = train.groupby("cat_col")["target"].transform("mean")

# CORRECT — computed inside fold
for fold, (tr_idx, va_idx) in enumerate(gkf.split(X, groups=groups)):
    means = train.iloc[tr_idx].groupby("cat_col")["target"].mean()
    train.loc[train.index[va_idx], "te_feature"] = train.iloc[va_idx]["cat_col"].map(means)
```

With smoothing:
```python
def smooth_target_encode(train_fold, val_fold, col, target, smooth_k=20):
    global_mean = train_fold[target].mean()
    stats = train_fold.groupby(col)[target].agg(["mean", "count"])
    stats["encoded"] = (stats["mean"] * stats["count"] + global_mean * smooth_k) / (stats["count"] + smooth_k)
    return val_fold[col].map(stats["encoded"]).fillna(global_mean)
```

---

## Signs of Validation Leakage

| Symptom | Likely cause |
|---------|-------------|
| OOF >> LB by >0.005 | Target encoding computed globally; group boundary not respected |
| OOF improves but LB doesn't move | Val set shares rows with train (wrong KFold) |
| OOF >> LB after adding auxiliary data | Auxiliary data from different distribution; or auxiliary labels leak into val |
| Very low val loss but high LB loss | Loss function (BCE) ≠ competition metric — use custom eval metric |

---

## CV Configuration

| Setting | Tuning | Training |
|---------|--------|----------|
| `n_folds` | 3 | 5 |
| `seeds` | [42] | [42, 123, 777, 2024, 9999] |
| Group column | same | same |
| Metric | competition exact | competition exact |
| Early stopping | enabled | enabled |

---

## Diagnosis Workflow

When OOF > LB:
1. Check group column — are there within-group leaks?
2. Check target encoding — computed globally?
3. Check feature engineering — does it use future rows?
4. Check early stopping metric — is it the exact competition formula?
5. Check aux/prior data — does the distribution match the test set? Compare OOF without aux data.

---

## See Also

| File | Why |
|------|-----|
| [feature-engineering.md](./feature-engineering.md) | Target encoding must be computed fold-by-fold — the fold split is defined here |
| [model-training.md](./model-training.md) | Trainers consume the fold indices produced by the CV strategy |
| [experiment-tracking.md](./experiment-tracking.md) | OOF vs LB divergence diagnosis — the first place to look when CV and LB disagree |
| [common-pitfalls.md](./common-pitfalls.md) | Pitfalls #5 (combined split OOF indexing) and #12 (gating score leakage) |
