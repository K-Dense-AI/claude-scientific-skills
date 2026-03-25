# Competition Metric Implementations

## Overview

This file is the single source of truth for all competition metric implementations. It covers: the fundamental distinction between training loss and eval metric; the generic `competition_score(y_true, y_pred)` pattern; per-framework metric injection for CatBoost (binary, regression, multiclass), LightGBM, XGBoost, and Neural Networks; NN training losses (FocalLoss, SmoothBCE, MSE, CrossEntropy); and shared scoring helpers in `base/common.py`.

**The most expensive bug in this pipeline:** using the wrong eval metric for early stopping. A model trained with `metric="binary_logloss"` when the competition uses a weighted AUC + LL blend will stop at the wrong iteration — typically 10× too early or too late — and every downstream experiment (tuning, ensemble, pseudo-labeling) will be built on a suboptimal base. This is silent: the model trains successfully and produces OOF scores; the OOF just does not track the leaderboard.

**How to use this file:**
1. Define `competition_score` in `base/metrics.py` — copy-adapt the template in Step 1 below
2. Copy the appropriate framework wrappers into `base/metrics.py` — one per framework you use
3. Import from `base/metrics.py` everywhere — never inline the formula in trainer or tuner files

---

## Two Distinct Concepts — Never Confuse Them

| Concept | What it is | Who uses it |
|---------|-----------|------------|
| **Training loss** | Differentiable objective the optimizer minimizes (BCE, Focal, SmoothBCE) | NN optimizer; tree models use objective internally |
| **Eval metric** | Competition score used for early stopping and model selection | All frameworks — must match the leaderboard formula exactly |

**These are completely independent.** Tree models (CatBoost/LGB/XGB) optimize their own built-in objective (binary cross-entropy internally) and separately evaluate on a custom metric for early stopping. NNs optimize a custom loss function AND monitor a custom eval score for checkpointing.

**The most common bug:** using a surrogate (plain logloss, accuracy) as the eval metric instead of the exact competition formula. This causes early stopping to fire at the wrong iteration.

---

## Single Source of Truth

Define metric wrappers once in `base/metrics.py`. Import into:
- All trainer files (early stopping)
- All tuner files (Optuna objective)
- The scoring helpers in `base/common.py`

Never inline the metric formula anywhere else.

---

## Step 1: Define Your Competition Metric

Define one function that takes `(y_true, y_pred)` and returns a float (higher = better). All competition-specific constants live here.

```python
# base/metrics.py — implement once for your competition
import numpy as np
from sklearn.metrics import roc_auc_score, log_loss, mean_squared_error, f1_score

TARGETS = [...]       # list of target column names — define per competition

# ── Choose the formula that matches the leaderboard exactly ──────────────
# Binary: AUC only
# def competition_score(y, p): return roc_auc_score(y, p)

# Binary: Log Loss only (negate → higher is better)
# def competition_score(y, p): return -log_loss(y, np.clip(p, 1e-15, 1-1e-15))

# Binary: Weighted AUC + LL blend
# NULL_LL = <null-model logloss — compute: log_loss(y_train, [y_train.mean()]*n)>
# W_AUC, W_LL = <auc_weight>, <ll_weight>  # competition-specific; fill in from problem statement
# def competition_score(y, p):
#     p = np.clip(p, 1e-15, 1-1e-15)
#     return roc_auc_score(y, p) * W_AUC + (1 - log_loss(y, p) / NULL_LL) * W_LL

# Regression: RMSE (negate for maximize direction)
# def competition_score(y, p): return -np.sqrt(mean_squared_error(y, p))

# Regression: RMSLE
# def competition_score(y, p): return -np.sqrt(((np.log1p(p) - np.log1p(y))**2).mean())

# Multi-class: Macro F1
# def competition_score(y, p): return f1_score(y, p.argmax(axis=1), average="macro")

def competition_score(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Single source of truth for the competition metric. HIGHER IS ALWAYS BETTER."""
    raise NotImplementedError("Implement this for your competition")
```

**Convention:** `competition_score` always maximizes. Negate RMSE/MAE/logloss if the competition is "lower is better".

---

## Common Competition Metrics Reference

| Competition type | Leaderboard metric | `competition_score` formula |
|-----------------|------------------|------------------------------|
| Binary clf | AUC-ROC | `roc_auc_score(y, p)` |
| Binary clf | Log Loss | `-log_loss(y, clip(p))` |
| Binary clf | Weighted AUC+LL blend | `w_auc * auc + w_ll * (1 - ll / null_ll)` |
| Regression | RMSE | `-sqrt(MSE(y, p))` |
| Regression | RMSLE | `-sqrt(mean((log1p(p)-log1p(y))²))` |
| Regression | MAE | `-mean(abs(y - p))` |
| Multi-class | Macro F1 | `f1_score(y, p.argmax(1), average="macro")` |
| Multi-class | Log Loss | `-log_loss(y, p)` |
| Multi-label | Mean AUC | `mean([roc_auc_score(y[:,i], p[:,i]) for i in range(n)])` |

**For multi-target competitions:** define one `competition_score` that averages or weights component scores across targets.

---

## CatBoost Custom Eval Metric

The `evaluate()` method receives **different input types depending on task**:
- **Binary classification**: `approxes[0]` = raw log-odds (logits) → **must apply sigmoid**
- **Regression**: `approxes[0]` = raw prediction values → **use directly, no activation**
- **Multiclass**: `approxes` = list of per-class log-odds arrays → apply softmax across classes

### Binary classification
```python
class CatBoostCompMetric:
    """
    For binary classification: CatBoost passes logits → must apply sigmoid.
    is_max_optimal=True → early stopping maximises this metric.
    """
    def is_max_optimal(self) -> bool:
        return True

    def evaluate(self, approxes, target, weight):
        p = 1.0 / (1.0 + np.exp(-np.array(approxes[0])))  # sigmoid: logit → prob
        p = np.clip(p, 1e-15, 1 - 1e-15)
        y = np.array(target)
        score = competition_score(y, p)
        return score, 1.0

    def get_final_error(self, error, weight):
        return error / (weight + 1e-38)
```

### Regression
```python
class CatBoostCompMetric:
    """For regression: approxes[0] is already the raw prediction — no activation."""
    def is_max_optimal(self) -> bool:
        return True   # set False if your competition metric is "lower is better" and you didn't negate

    def evaluate(self, approxes, target, weight):
        p = np.array(approxes[0])   # raw prediction, no sigmoid/softmax
        y = np.array(target)
        score = competition_score(y, p)
        return score, 1.0

    def get_final_error(self, error, weight):
        return error / (weight + 1e-38)
```

### Multiclass
```python
class CatBoostCompMetric:
    """For multiclass: approxes is a list of K arrays (one per class), logit-scaled."""
    def is_max_optimal(self) -> bool:
        return True

    def evaluate(self, approxes, target, weight):
        # Stack and softmax across classes
        logits = np.column_stack([np.array(a) for a in approxes])  # (n, K)
        logits -= logits.max(axis=1, keepdims=True)                 # numerical stability
        exp_l = np.exp(logits)
        p = exp_l / exp_l.sum(axis=1, keepdims=True)               # softmax
        y = np.array(target).astype(int)
        score = competition_score(y, p)   # competition_score receives (n,) labels + (n, K) probs
        return score, 1.0

    def get_final_error(self, error, weight):
        return error / (weight + 1e-38)
```

**Usage in trainer:**
```python
params = dict(
    eval_metric=CatBoostCompMetric(),
    early_stopping_rounds=150,   # noisy metrics need patience; tune up if metric oscillates
    # DO NOT set use_best_model=False — default True is correct
)
```

---

## LightGBM feval

LGB passes already-transformed outputs — probabilities for classification, raw values for regression.
For multiclass, `preds` is a **flattened** `(n * n_classes,)` array — must reshape.

```python
# Binary classification / regression
def make_lgb_feval():
    """
    Usage:
        params['metric'] = 'None'           # MUST disable default metric (all task types)
        lgb.early_stopping(500, first_metric_only=True)
        lgb.train(..., feval=make_lgb_feval())
    """
    def feval(preds: np.ndarray, eval_data) -> tuple:
        y = eval_data.get_label()
        p = np.clip(preds, 1e-15, 1 - 1e-15)  # safe for probabilities; harmless for regression
        score = competition_score(y, p)
        return "comp_score", score, True   # (name, score, higher_is_better)
    return feval

# Multiclass — preds is flattened (n * n_classes,), must reshape
def make_lgb_feval_multiclass(n_classes: int):
    def feval(preds: np.ndarray, eval_data) -> tuple:
        y = eval_data.get_label().astype(int)
        p = preds.reshape(len(y), n_classes)   # (n, K) probabilities
        score = competition_score(y, p)
        return "comp_score", score, True
    return feval
```

**Critical params:**
```python
params = {
    "metric": "None",   # REQUIRED — disables built-in logloss
    "verbose": -1,
    "feature_pre_filter": False,
}
lgb.train(
    params, ds_tr,
    valid_sets=[ds_va],
    feval=make_lgb_feval(),
    callbacks=[lgb.early_stopping(500, first_metric_only=True)],
)
```

---

## XGBoost custom_metric

XGBoost passes already-transformed outputs. For `multi:softprob`, `predt` is shape `(n * n_classes,)` — must reshape.

```python
# Binary classification / regression
def make_xgb_eval():
    """
    Usage:
        params['disable_default_eval_metric'] = 1   # MUST suppress default (all task types)
        xgb.train(..., custom_metric=make_xgb_eval(), maximize=True)
    """
    def eval_fn(predt: np.ndarray, dtrain) -> tuple:
        y = dtrain.get_label()
        p = np.clip(predt, 1e-15, 1 - 1e-15)  # safe for probabilities; harmless for regression
        score = competition_score(y, p)
        return "comp_score", score   # maximize=True handles direction
    return eval_fn

# Multiclass — predt is flattened (n * n_classes,), must reshape
def make_xgb_eval_multiclass(n_classes: int):
    def eval_fn(predt: np.ndarray, dtrain) -> tuple:
        y = dtrain.get_label().astype(int)
        n = len(y)
        p = predt.reshape(n, n_classes)   # (n, K) probabilities
        score = competition_score(y, p)
        return "comp_score", score
    return eval_fn
```

**Critical params:**
```python
params = {
    "disable_default_eval_metric": 1,   # REQUIRED
    "verbosity": 0,
}
xgb.train(
    params, dtrain,
    evals=[(dval, "val")],
    custom_metric=make_xgb_eval(),
    maximize=True,                       # REQUIRED — early-stop maximises
    early_stopping_rounds=200,
)
```

---

## Neural Network — Training Loss (optimizer objective)

The NN uses a **custom training loss**, separate from the competition metric used for checkpointing.

### FocalLoss (default)

Best for imbalanced binary classification. Downweights easy negatives so training focuses on hard examples.

```python
class FocalLoss(nn.Module):
    def __init__(self, alpha: float = 0.25, gamma: float = 2.0):
        super().__init__()
        self.alpha = alpha
        self.gamma = gamma

    def forward(self, pred, target):
        pred = pred.clamp(1e-7, 1 - 1e-7)
        bce = -target * torch.log(pred) - (1 - target) * torch.log(1 - pred)
        pt = target * pred + (1 - target) * (1 - pred)   # p_t = prob of correct class
        return (self.alpha * (1 - pt) ** self.gamma * bce).mean()
```

- `gamma=2.0`: standard focal weight — reduces loss for well-classified examples
- `alpha=0.25`: class-balance factor (positive class gets weight 0.25 vs 0.75 for negatives)
- Higher `gamma` → more focus on hard examples; too high → ignores easy positives

### Label-Smooth BCE

Alternative when focal is too aggressive. Softens targets: 1→(1-s/2), 0→(s/2).

```python
class _SmoothBCE(nn.Module):
    def __init__(self, s: float = 0.05):
        super().__init__()
        self.s = s

    def forward(self, p, y):
        p = p.clamp(1e-7, 1 - 1e-7)
        ys = y * (1 - self.s) + 0.5 * self.s   # smooth: 1→0.975, 0→0.025
        return (-ys * torch.log(p) - (1 - ys) * torch.log(1 - p)).mean()
```

- `s=0.05`: mild smoothing — prevents overconfident predictions
- Good when training data has label noise (e.g. pseudo-labeled rows)

### Which loss to use

| Situation | Recommended loss |
|-----------|-----------------|
| Strong class imbalance, clean labels (binary) | `focal` (default) |
| Pseudo-labeled data mixed in (binary) | `smooth_bce` — tolerant of noisy labels |
| Balanced classes, clean labels (binary) | `bce` (plain) |
| Regression | `mse` (MSELoss) or `mae` (L1Loss) or `huber` |
| Multi-class | `cross_entropy` (CrossEntropyLoss) |

### Loss selector in train_nn_fold
```python
if loss_type == "focal":
    criterion = FocalLoss(alpha=focal_alpha, gamma=focal_gamma)
elif loss_type == "smooth_bce":
    criterion = _SmoothBCE(s=label_smooth)
elif loss_type == "mse":
    criterion = nn.MSELoss()
elif loss_type == "mae":
    criterion = nn.L1Loss()
elif loss_type == "huber":
    criterion = nn.HuberLoss(delta=huber_delta)
elif loss_type == "cross_entropy":
    criterion = nn.CrossEntropyLoss()  # expects raw logits + integer targets
else:
    criterion = nn.BCELoss()
```

### Class imbalance: WeightedRandomSampler

For imbalanced binary/multi-label tasks, the NN uses a `WeightedRandomSampler` to oversample positives in each batch:

```python
if use_sampler:
    # Define positives based on your target(s). For single target:
    pos_mask = (y_tr_dict[targets[0]] == 1)
    # For multi-target / multi-label (positive if ANY target is 1):
    # pos_mask = np.any([y_tr_dict[t] == 1 for t in targets], axis=0)
    weights = np.ones(len(pos_mask))
    weights[pos_mask] = (1 - pos_mask.mean()) / max(pos_mask.mean(), 1e-6)
    sampler = WeightedRandomSampler(torch.FloatTensor(weights), len(weights), replacement=True)
```

This is complementary to `FocalLoss` — sampler fixes the batch distribution, focal loss refines gradient weighting. **Only for binary/multi-label tasks**; skip for regression and multiclass.

---

## Neural Network — Eval Metric (checkpointing)

The competition score is computed on the val set each epoch. The model state is checkpointed whenever it improves, and training uses early stopping based on it.

```python
# In epoch loop — checkpoint on competition score, NOT on loss
va_score = sum(competition_score(y_va[t], va_preds[i]) for i, t in enumerate(targets))
if va_score > best_score:
    best_score = va_score
    best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
    no_improve = 0
else:
    no_improve += 1
if no_improve >= patience:
    break

model.load_state_dict(best_state)   # ALWAYS restore best before predicting
```

**Key:** the loss goes down monotonically per epoch, but the competition score on val does not — it can oscillate. Early stopping with `patience=20` epochs without improvement is necessary.

---

## Full Framework Comparison

| Framework | Training objective | Early-stop metric | How to inject custom metric |
|-----------|-------------------|------------------|----------------------------|
| **CatBoost** | Built-in objective (Logloss / RMSE / MultiClass) | `CatBoostCompMetric` class | `eval_metric=CatBoostCompMetric()` |
| **LightGBM** | Built-in objective (binary / regression / multiclass) | `make_lgb_feval()` feval | `"metric": "None"` + `feval=` + `first_metric_only=True` |
| **XGBoost** | Built-in objective (binary:logistic / reg / softprob) | `make_xgb_eval()` custom_metric | `"disable_default_eval_metric": 1` + `custom_metric=` + `maximize=True` |
| **NN** | FocalLoss / SmoothBCE / MSELoss / CrossEntropyLoss | `competition_score()` per epoch | Compute in epoch loop; checkpoint manually |

**Critical API differences:**
- CatBoost binary: receives **raw log-odds** (logits) → apply sigmoid in `evaluate()`
- CatBoost regression: receives **raw prediction values** → use directly, no activation
- CatBoost multiclass: receives K logit arrays → apply softmax in `evaluate()`
- LGB/XGB: receive **already-transformed outputs** (probabilities or raw values) → no transformation needed
- NN: you control everything → apply appropriate activation in model forward pass (sigmoid / softmax / none)

---

## base/common.py — Shared Scoring Helpers

```python
from .metrics import competition_score, TARGETS

def compute_oof_score(y_dict: dict, oof_dict: dict) -> tuple[float, dict]:
    """
    Aggregate competition_score over all targets.
    Returns (mean_score, per_target_dict).
    Adapt aggregation (mean vs weighted sum) to match the competition formula.
    """
    per_target = {}
    for t in TARGETS:
        y_true = y_dict[t]
        y_pred = np.clip(oof_dict[t], 1e-15, 1 - 1e-15)
        per_target[t] = competition_score(y_true, y_pred)
    total = np.mean(list(per_target.values()))   # or weighted sum if competition specifies
    return total, per_target

def print_scores(y_dict, oof_dict, label=""):
    score, per = compute_oof_score(y_dict, oof_dict)
    for t, s in per.items():
        logger.info(f"    {t}: {s:.5f}")
    tag = f"  \u2190 {label}" if label else ""
    logger.info(f"    Score = {score:.6f}{tag}")
    return score
```

---

## See Also

| File | Why |
|------|-----|
| [model-training.md](./model-training.md) | Trainer parameter sets that reference the metric wrappers defined here |
| [hyperparameter-tuning.md](./hyperparameter-tuning.md) | Tuner objectives import `competition_score` directly from this file |
| [output-format.md](./output-format.md) | Metric name determines the prediction type required in the submission |
| [common-pitfalls.md](./common-pitfalls.md) | Pitfalls #2, #3, #10, #11 — all caused by wrong eval metric setup |
