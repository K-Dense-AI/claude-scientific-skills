# Model Training Reference

## Overview

This file is the authoritative reference for implementing CatBoost, LightGBM, XGBoost, and Neural Network base models in a competition pipeline. It covers: the trainer architecture pattern (stateless fold engine); the critical distinction between training objective and eval metric; correct parameter sets for each task type (binary, regression, multiclass); training loops with GroupKFold and seed averaging; and auxiliary data handling.

**Before reading further**, confirm your task type in the SKILL.md Task Type Decision Guide. Several parameters documented here are **binary classification only** — `scale_pos_weight`, `is_unbalance`, `auto_class_weights: Balanced` — and must be removed for regression and multiclass. Applying them to the wrong task type produces silent score degradation.

**When to use:** When writing or reviewing any trainer file (`base/lgb_trainer.py`, `xgb_trainer.py`, `nn_trainer.py`) or model entrypoint (`train/cat.py`, `lgb.py`, etc.). Also consult when adapting a binary classifier for regression or multiclass.

---

## Trainer Architecture

Trainers are **stateless fold engines** — they receive data and params, return OOF + test preds. Entrypoints handle I/O, logging, and OOF saving.

> **Before writing any trainer:** check the Task Type Decision Guide in SKILL.md. Several parameters below are **binary classification only** and must be removed or changed for regression/multiclass.

```
entrypoint (train/cat.py)
  → _apply_tune(tune_dir)       # load + merge tuned params
  → load_or_build_matrices()    # get arrays from cache
  → trainer (base/lgb_trainer.py)  # fold loop, returns {oof, test, metrics}
  → print_scores() + save_oof() + make_submission()
```

---

## Training Objective vs Eval Metric

Every model has **two separate concerns**:

| | Tree models (CB/LGB/XGB) | Neural Network |
|-|--------------------------|----------------|
| **Training objective** | Built-in BCE / squared-error (internal to framework) | Custom loss: `FocalLoss`, `_SmoothBCE`, `MSELoss`, or `BCELoss` |
| **Eval metric** | Custom competition metric wrapper | Competition score computed per epoch |
| **Early stopping driven by** | Eval metric | Competition score (`competition_score`) |

For tree models you don't change the training objective — you only inject a custom eval metric for early stopping. For the NN you control both. See [competition-metrics.md](./competition-metrics.md) for full implementations.

---

## CatBoost

### Key params
```python
# ── Binary classification ─────────────────────────────────────
CB_BASE_PARAMS = {
    "iterations": 5000,
    "learning_rate": 0.05,
    "depth": 6,
    "l2_leaf_reg": 10.0,
    "min_data_in_leaf": 50,
    "auto_class_weights": "Balanced",  # BINARY IMBALANCED ONLY — remove for regression/multiclass
    "random_seed": seed,
    "eval_metric": CatBoostCompMetric(),  # competition metric wrapper — see competition-metrics.md
    "early_stopping_rounds": 150,         # noisy metric needs patience
    "verbose": 0,
    "cat_features": cat_indices,
}

# ── Regression: swap class → CatBoostRegressor, remove auto_class_weights ──
# CB_BASE_PARAMS = {"loss_function": "RMSE", "eval_metric": CatBoostCompMetric(), ...}

# ── Multiclass: loss_function="MultiClass", classes_count=N ──────────────
# CB_BASE_PARAMS = {"loss_function": "MultiClass", "classes_count": N, ...}
```

### Training loop
```python
gkf = GroupKFold(n_splits=cfg.n_folds)
for seed in cfg.seeds:
    for fold, (tr_idx, va_idx) in enumerate(gkf.split(X, groups=groups)):
        pool_tr = cb.Pool(X[tr_idx],  y[tr_idx],  cat_features=cat_idx)
        pool_va = cb.Pool(X[va_idx],  y[va_idx],  cat_features=cat_idx)

        model = cb.CatBoostClassifier(**params)
        model.fit(pool_tr, eval_set=pool_va, verbose=0)

        oof[t][va_idx] += model.predict_proba(X[va_idx])[:, 1] / len(cfg.seeds)
        tst[t] += model.predict_proba(X_test)[:, 1] / (cfg.n_folds * len(cfg.seeds))
```

> **Combined data pattern:** if training on `X_combined = vstack([X_train, X_aux])`, restrict val and OOF to train rows:
> `va_train_idx = va_idx[va_idx < n_train]` — index on combined but accumulate OOF only on original rows.

### Tuning key names → CB param names
```python
# Optuna search space uses shorthand; _build_params translates:
_BASE_KEYS = {"iters", "lr", "depth", "l2", "min_leaf"}
# iters → iterations, lr → learning_rate, l2 → l2_leaf_reg, min_leaf → min_data_in_leaf
```

---

## LightGBM

### Key params
```python
# ── Binary classification ─────────────────────────────────────
LGB_BASE_PARAMS = {
    "objective": "binary",      # regression: "regression" | multiclass: "multiclass" + num_class=N
    "metric": "None",           # REQUIRED — disables default metric for all task types
    "learning_rate": 0.005,
    "num_leaves": 31,
    "max_depth": 6,
    "min_child_samples": 30,
    "subsample": 0.8,
    "colsample_bytree": 0.6,
    "reg_alpha": 0.5,
    "reg_lambda": 3.0,
    "is_unbalance": True,        # BINARY IMBALANCED ONLY — remove for regression/multiclass
    "verbose": -1,
    "feature_pre_filter": False,
    "device_type": "gpu",
    "num_threads": -1,
}

# ── Regression: change objective, remove is_unbalance ────────────────────
# {"objective": "regression", "metric": "None", ...}   (no is_unbalance)

# ── Multiclass: add num_class ────────────────────────────────────────────
# {"objective": "multiclass", "num_class": N, "metric": "None", ...}
```

### Training call
```python
model = lgb.train(
    params,
    ds_tr,
    num_boost_round=5000,
    valid_sets=[ds_va],
    feval=make_lgb_feval(target),
    callbacks=[
        lgb.early_stopping(500, first_metric_only=True),
        lgb.log_evaluation(period=-1),
    ],
)
oof[t][va_idx] += model.predict(x_va) / len(seeds)
```

### Key translation map for Optuna keys
```python
_LGB_KEY_MAP = {
    "lr":                "learning_rate",
    "min_child_samples": "min_child_samples",  # passthrough
    # All other Optuna keys match LGB param names directly
}
```

---

## XGBoost

### Key params
```python
# ── Binary classification ─────────────────────────────────────
XGB_BASE_PARAMS = {
    "objective": "binary:logistic",  # regression: "reg:squarederror" | multiclass: "multi:softprob" + num_class=N
    "disable_default_eval_metric": 1,  # REQUIRED — suppress default eval for all task types
    "learning_rate": 0.01,
    "max_depth": 5,
    "min_child_weight": 30,
    "subsample": 0.7,
    "colsample_bytree": 0.5,
    "reg_alpha": 1.0,
    "reg_lambda": 5.0,
    "tree_method": "hist",
    "device": "cuda:0",
    "verbosity": 0,
    "nthread": -1,
}

# ── Regression: change objective, no scale_pos_weight ────────────────────
# {"objective": "reg:squarederror", "disable_default_eval_metric": 1, ...}

# ── Multiclass: add num_class ────────────────────────────────────────────
# {"objective": "multi:softprob", "num_class": N, "disable_default_eval_metric": 1, ...}
```

### scale_pos_weight — BINARY IMBALANCED ONLY
```python
# Binary classification with class imbalance: compute per fold, not globally
n_neg = (y_tr == 0).sum()
n_pos = max((y_tr == 1).sum(), 1)
params["scale_pos_weight"] = n_neg / n_pos

# DO NOT add scale_pos_weight for regression or multiclass — it has no effect or causes errors.
```

### Training call
```python
model = xgb.train(
    params, dtrain,
    num_boost_round=20000,
    evals=[(dval, "val")],
    custom_metric=make_xgb_eval(target),
    maximize=True,                  # REQUIRED with custom_metric
    early_stopping_rounds=200,
    verbose_eval=0,
)
oof[t][va_idx] += model.predict(dval) / len(seeds)
```

---

## Neural Network

### Architecture: Flexible tabular NN with embeddings + residual blocks
```python
class _FlexNN(nn.Module):
    def __init__(self, num_feat, cat_sizes, embed_dim=20, hidden_dim=256,
                 n_blocks=3, dropout=0.3, act_name="silu", head_hidden=64):
        # Embeddings for all categorical features
        # BN → in_proj → n_blocks of ResBlock → output head(s)
        pass

# Adapt the output head to task type:
#
# Binary classification:     nn.Linear(h, 1) + nn.Sigmoid()  → shape (n, 1), squeeze to (n,)
# Multi-target binary:       N × nn.Linear(h, 1) + nn.Sigmoid() → dict of (n,) arrays
# Regression:                nn.Linear(h, 1) + no activation  → shape (n, 1), squeeze
#                            (add nn.ReLU() if target is guaranteed non-negative)
# Multiclass (N classes):    nn.Linear(h, N) + nn.Softmax(dim=-1) → shape (n, N)
# Multi-label (N binary):    N × nn.Linear(h, 1) + nn.Sigmoid()  (same as multi-target)
```

### NN training loss by task type
```python
# Binary / multi-label / multi-target:  FocalLoss or BCEWithLogitsLoss
# Regression:                           nn.MSELoss() or nn.L1Loss() (MAE)
# Multiclass:                           nn.CrossEntropyLoss()
#
# For imbalanced binary: FocalLoss with gamma=2 suppresses easy negatives
# For noisy labels: SmoothBCE (label smoothing 0.05)
```

### Training loop pattern
```python
# checkpoint on competition score, not loss
if va_score > best_score:
    best_score = va_score
    best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
    no_improve = 0
else:
    no_improve += 1
if no_improve >= patience:
    break

model.load_state_dict(best_state)   # restore best before predicting
```

### Key training params
```python
NN_PARAMS = {
    "nn_epochs": 300,
    "nn_lr": 5e-4,
    "nn_hidden_dim": 256,
    "nn_n_blocks": 3,
    "nn_dropout": 0.30,
    "nn_loss_type": "focal",        # focal | smooth_bce | bce
    "nn_focal_gamma": 2.0,
    "nn_focal_alpha": 0.25,
    "nn_weight_decay": 1e-4,
    "nn_batch_size": 512,
    "nn_use_sampler": True,         # weighted sampler for class imbalance
    "nn_scheduler_type": "cosine_warm",
    "nn_patience": 20,
}
```

### Loss selection guidance

| Scenario | `loss_type` |
|----------|------------|
| Default — class-imbalanced, clean labels | `"focal"` |
| Pseudo-labeled rows mixed in (noisy labels) | `"smooth_bce"` |
| Balanced dataset, clean labels | `"bce"` |

- `focal`: downweights easy examples via `(1-pt)^gamma`; alpha=0.25 weights positives more
- `smooth_bce`: softens hard targets (1→0.975, 0→0.025 at s=0.05); reduces overconfidence from noisy labels
- Both work **in conjunction with** `use_sampler=True` which oversamples positives at the batch level

### GPU setup
```python
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# Always seed before model creation
torch.manual_seed(seed)
np.random.seed(seed)
torch.cuda.manual_seed_all(seed)
```

---

## OOF Accumulation Pattern (all models)

```python
# Initialise
oof = {t: np.zeros(n_train) for t in targets}
tst = {t: np.zeros(len(test)) for t in targets}

# Accumulate: divide by seeds (not seeds × folds) in fold loop,
# because each row appears in exactly one val fold
oof[t][va_idx] += model.predict(x_va) / len(seeds)

# Test: all folds × all seeds see every test row
tst[t] += model.predict(x_test) / (n_folds * len(seeds))
```

---

## GPU Consistency Rule

Whatever device you tune on **must match** the device you train on. If you tune on CUDA, both:
- `lgb_trainer.py` / `LGB_BASE_PARAMS` → `"device_type": "gpu"`
- `xgb_trainer.py` / `XGB_BASE_PARAMS` → `"device": "cuda:0"`

must match `tune_lgb.py` and `tune_xgb.py`.

---

## See Also

| File | Why |
|------|-----|
| [competition-metrics.md](./competition-metrics.md) | Full metric wrapper implementations for CB / LGB / XGB / NN |
| [hyperparameter-tuning.md](./hyperparameter-tuning.md) | Tuner objectives must mirror the exact trainer params shown here |
| [validation-strategy.md](./validation-strategy.md) | GroupKFold and OOF accumulation patterns used in trainer loops |
| [common-pitfalls.md](./common-pitfalls.md) | Pitfalls #2, #3 (wrong eval metric), #7 (scale_pos_weight), #15 (pseudo params) |
