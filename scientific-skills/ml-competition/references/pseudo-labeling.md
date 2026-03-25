# Pseudo-Labeling

## Overview

Pseudo-labeling extends training data by treating high-confidence test predictions as labeled examples and retraining on the combined train + pseudo-test set. When applied correctly, it typically adds 0.001–0.003 to LB on competitions where the test distribution is similar to train. When applied incorrectly — too early, with low-confidence labels, or on distribution-shifted test data — it amplifies base model errors and degrades performance.

This file covers: preconditions that must hold before any pseudo-labeling attempt; pseudo label generation by task type (binary, regression, multiclass, multi-label); the fold loop pattern; confidence thresholds; sample weighting; and pitfalls.

**The single most important rule:** pseudo-labeling must only begin after all base models have converged and are individually tuned. Pseudo labels inherit all errors from the base models — running pseudo-labeling on an undertrained pipeline amplifies weaknesses rather than strengths.

**When to use:** After base models are fully trained (OOF stable, Optuna tuning complete, no pending architecture changes). After base ensemble is evaluated. Before adding the meta ensemble layer.

---

## When to Use

- ✅ All base models have converged (OOF stable, no pending architecture changes)
- ✅ Test set distribution is similar to train
- ✅ Pseudo label confidence is high (high agreement between base models)
- ❌ Never before base models are fully tuned — pseudo amplifies base model errors
- ❌ Never if base models disagree significantly on test predictions

---

## Core Pipeline

```python
# 1. Average base model test predictions
avg_test = {t: np.zeros(n_test) for t in TARGETS}
for model_name in label_models:
    d = load_oof(model_name, tag=base_tag)
    for t in TARGETS:
        avg_test[t] += d["test"][t] / len(label_models)

# 2. Hard pseudo labels at threshold 0.5
pseudo_y = {t: (avg_test[t] > PSEUDO_THRESHOLD).astype(float) for t in TARGETS}

# 3. Retrain each model on train_fold + pseudo_test
```

---

## Constants

```python
PSEUDO_THRESHOLD = 0.5   # binary classification only — hard labels at 0.5
PSEUDO_WEIGHT    = 0.5   # pseudo rows get half the sample weight of real rows (all tasks)
```

## Pseudo Label Generation by Task Type

The method for generating pseudo labels from averaged base model test predictions differs by task:

### Binary classification
```python
# Hard labels at threshold — cleaner than soft for tree models
pseudo_y = {t: (avg_test[t] > PSEUDO_THRESHOLD).astype(float) for t in TARGETS}
```

### Regression
```python
# Use soft labels directly — no thresholding
pseudo_y = {t: avg_test[t].copy() for t in TARGETS}   # continuous values
# Remove scale_pos_weight from XGB — not applicable for regression
```

### Multiclass
```python
# Hard labels via argmax; models receive integer class labels
# avg_test[t] shape: (n_test, n_classes) — avg softmax probabilities
pseudo_y = {t: avg_test[t].argmax(axis=1).astype(int) for t in TARGETS}
# Pass one-hot or integer labels depending on framework
```

### Multi-label / Multi-target binary
```python
# Each target is independent binary — same as single binary
pseudo_y = {t: (avg_test[t] > PSEUDO_THRESHOLD).astype(float) for t in TARGETS}
```

**Pseudo label quality check** (any task type):
```python
# High agreement between models = reliable pseudo labels
agreement = {t: (avg_test[t] > 0.8).mean() + (avg_test[t] < 0.2).mean() for t in TARGETS}
logger.info(f"Pseudo label confidence: {agreement}")  # aim for > 0.7 if possible
# If agreement < 0.5: base models disagree heavily — pseudo will amplify noise, skip it
```

---

## Fold Loop Pattern (all models)

```python
for seed in cfg.seeds:
    for fold, (tr_idx, va_idx) in enumerate(gkf.split(train, groups=groups)):
        # Augmented training: real fold rows + pseudo-labeled test
        X_fold = np.vstack([X_train[tr_idx], X_test])
        y_fold = np.concatenate([y_dict[t][tr_idx], pseudo_y[t]])
        w_fold = np.concatenate([
            np.ones(len(tr_idx)),           # real: full weight
            np.full(len(test), PSEUDO_WEIGHT),  # pseudo: half weight
        ])
        
        # Validation: real val rows only — no pseudo in val
        x_va = X_train[va_idx]
        y_va = y_dict[t][va_idx]
        
        # OOF: real rows only
        oof[t][va_idx] += model.predict(x_va) / len(cfg.seeds)
```

---

## CatBoost Pseudo

```python
# CB_FIT must use CatBoostLBMetric — NOT "Logloss"
CB_FIT = dict(
    eval_metric=CatBoostLBMetric(target),
    early_stopping_rounds=150,
    verbose=0,
)

pool_tr = cb.Pool(X_fold, y_fold, weight=w_fold, cat_features=cat_indices)
pool_va = cb.Pool(X_train_cb.iloc[va_idx], y_dict[t][va_idx], cat_features=cat_indices)
model = cb.CatBoostClassifier(**_cb_params(hp, seed + fold))
model.fit(pool_tr, eval_set=pool_va, verbose=0)
```

---

## LightGBM Pseudo

```python
params = dict(LGB_PARAMS)
params["seed"] = seed + fold

ds_tr = lgb_lib.Dataset(X_fold, y_fold, weight=w_fold)
ds_va = lgb_lib.Dataset(X_train[va_idx], y_dict[t][va_idx])

model = lgb_lib.train(
    params, ds_tr,
    num_boost_round=5000,
    valid_sets=[ds_va],
    feval=make_lgb_feval(t),
    callbacks=[
        lgb_lib.early_stopping(500, first_metric_only=True, verbose=False),
        lgb_lib.log_evaluation(period=-1),
    ],
)
```

---

## XGBoost Pseudo — scale_pos_weight

Compute from **real training fold labels only** — same as `xgb_trainer.py`. This is intentional: the pseudo rows are already down-weighted via `w_fold`; `scale_pos_weight` should reflect the real class imbalance, not the inflated pseudo negatives.

```python
n_neg = (y_dict[t][tr_idx] == 0).sum()
n_pos = max((y_dict[t][tr_idx] == 1).sum(), 1)
params = dict(XGB_PARAMS, scale_pos_weight=n_neg / n_pos, seed=seed + fold)
```

**NEVER remove `scale_pos_weight` entirely.** On a 1% positive-rate problem, omitting it collapses OOF by ~0.10 (from 0.928 to 0.825 observed). See [common-pitfalls.md](./common-pitfalls.md) #15.

---

## NN Pseudo

```python
x_num_aug = np.vstack([X_num_tr[tr_idx], X_num_te])
x_cat_aug = np.vstack([X_cat_tr[tr_idx], X_cat_te])
y_aug = {t: np.concatenate([y_dict[t][tr_idx], pseudo_y[t]]) for t in TARGETS}

# WeightedRandomSampler or sample_weight to downweight pseudo rows
w_aug = np.concatenate([np.ones(len(tr_idx)), np.full(len(test), PSEUDO_WEIGHT)])
```

---

## OOF Saving

```python
save_oof(oof_lgb, tst_lgb, "lgb", tag="pseudo")
# → oof/pseudo_lgb.pkl
```

Resume check in orchestrator:
```python
pseudo_exists = not (Path("oof") / f"pseudo_{m}.pkl").exists()
```

---

## Pitfalls

### Pseudo hurts if:
1. **Base models overfit** — pseudo amplifies bad test predictions back into training
2. **Prior/auxiliary data added simultaneously** — two changes at once; impossible to diagnose which caused regression
3. **Different label sets per model** — if you retrain only subset of models but use all models for pseudo labels, label source mismatch on resume
4. **Threshold too aggressive** — 0.3 or 0.7 creates very imbalanced pseudo label distribution
5. **Weight too high** — PSEUDO_WEIGHT=1.0 treats noisy pseudo rows same as real labels

### Diagnosis
- Compare base model OOF vs pseudo OOF — pseudo should be ≥ base by 0.001-0.003
- If pseudo < base: pseudo labels are too noisy (base models aren't good enough yet)
- If LB for pseudo << OOF for pseudo: pseudo-labeled test data shifted distribution of predictions

### Safe experiment protocol
1. Run base models → submit individual best
2. Run pseudo → check OOF gain ≥ 0.001 per model
3. Run meta on base only → submit
4. Run meta with pseudo → submit and compare
Never submit two changes at once.

---

## See Also

| File | Why |
|------|-----|
| [ensemble-meta.md](./ensemble-meta.md) | Ensemble is run on pseudo OOF outputs |
| [common-pitfalls.md](./common-pitfalls.md) | Pitfalls #6 (mismatched labels on partial resume), #15 (scale_pos_weight) |
| [validation-strategy.md](./validation-strategy.md) | Same fold splits must be used in pseudo retraining |
| [experiment-tracking.md](./experiment-tracking.md) | OOF gain threshold — only commit pseudo if OOF improves ≥ 0.001 |
