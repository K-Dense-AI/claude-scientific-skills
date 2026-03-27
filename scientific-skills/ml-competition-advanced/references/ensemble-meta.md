# Ensemble & Meta-Learning

## Overview

Ensembling combines the predictions of multiple base models to improve generalization beyond any single model. This file documents the four ensemble levels used in this pipeline — weighted blend, LogReg stacking, dynamic gating, and weights+gating hybrid — and explains when each one is appropriate.

**A critical principle:** ensemble scores are only meaningful if computed on honest out-of-fold (OOF) predictions. Never evaluate an ensemble on in-sample predictions or on a subset of the val fold — the score will be inflated and will not reflect true generalization. All ensemble methods here use the same OOF arrays produced by the base trainer fold loops.

**Ensemble hierarchy by increasing complexity and overfit risk:** weighted blend → LogReg stacking → dynamic gating → weights+gating blend. On datasets with fewer than 20K rows, Level 1 (weighted blend) frequently beats more complex methods. Always start at Level 1 and move up only when OOF gain is confirmed.

**When to use:** After all base models are trained and their OOF pkl files saved. Do not attempt ensembling before base models are individually optimized — it will obscure rather than fix model-level bugs.

---

## Ensemble Hierarchy (ascending trust order)

```
Level 1 — Weighted blend (Nelder-Mead)          [A] weights
  → Fast, interpretable, no overfitting risk
  → OOF score = honest

Level 2 — LogReg stacking                       [B] stack
  → Learns non-linear weight combinations per target
  → OOF score = honest (cross-validated)
  → Often underperforms blend on small datasets

Level 3 — Dynamic gating (MoE)                  [C] gating
  → Per-sample routing: "which model is best for this input?"
  → OOF score = honest (GroupKFold on full OOF)
  → Most complex; most overfit risk; needs large dataset to help
  → On <20K rows, honest OOF often below Level 1

Level 4 — Weights + Gating blend                [D] weights_gating
  → Nelder-Mead over base OOFs + gating OOF as 9th input
  → Best of both: interpretable blend + gating as meta-feature
  → Usually the best or equal-best in practice
```

---

## Level 1: Nelder-Mead Weighted Blend

```python
from scipy.optimize import minimize
import numpy as np

def nelder_mead_blend(oofs: list[dict], y_dict: dict, score_fn, n_restarts=20):
    """
    Find softmax weights over models that maximise competition score.
    oofs: list of {target: np.array} — one per model
    """
    n_models = len(oofs)
    best_score = -np.inf
    best_w = None

    def objective(raw_w):
        w = np.exp(raw_w) / np.exp(raw_w).sum()   # softmax → sums to 1
        blended = {
            t: sum(w[i] * oofs[i][t] for i in range(n_models))
            for t in y_dict
        }
        sc, _ = score_fn(y_dict, blended)
        return -sc

    for _ in range(n_restarts):
        x0 = np.random.randn(n_models) * 0.1
        res = minimize(objective, x0, method="Nelder-Mead",
                       options={"maxiter": 5000, "xatol": 1e-6})
        if -res.fun > best_score:
            best_score = -res.fun
            best_w = np.exp(res.x) / np.exp(res.x).sum()

    return best_w, best_score
```

---

## Level 2: LogReg Stacking

```python
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold

def logistic_stacking(oofs: list[dict], y_dict: dict, groups: np.ndarray,
                      targets: list[str], n_folds: int = 5):
    """
    Train one LogReg per target with OOF features = stacked base model preds.
    Returns stacked OOF predictions.
    """
    n = len(list(y_dict.values())[0])
    meta_oof = {t: np.zeros(n) for t in targets}
    gkf = GroupKFold(n_splits=n_folds)

    for t in targets:
        X_meta = np.column_stack([oof[t] for oof in oofs])   # (n, n_models)
        y = y_dict[t]

        for tr_idx, va_idx in gkf.split(X_meta, groups=groups):
            lr = LogisticRegression(C=1.0, max_iter=1000)
            lr.fit(X_meta[tr_idx], y[tr_idx])
            meta_oof[t][va_idx] = lr.predict_proba(X_meta[va_idx])[:, 1]

    sc, _ = score_fn(y_dict, meta_oof)
    return meta_oof, sc
```

---

## Level 3: Dynamic Gating (Mixture of Experts)

```python
class GatingNetwork(nn.Module):
    """Learns which base model to trust per sample."""
    def __init__(self, n_features: int, n_models: int, hidden: int = 64):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_features, hidden),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(hidden, n_models),
            nn.Softmax(dim=-1),
        )
    def forward(self, x):
        return self.net(x)   # (batch, n_models) weights
```

### Honest OOF for gating
```python
def run_dynamic_gating_oof(oofs, y_dict, X_feat, groups, targets, n_folds=5):
    """
    Train gating network with proper GroupKFold OOF to avoid in-sample optimism.
    """
    n = len(list(y_dict.values())[0])
    gated_oof = {t: np.zeros(n) for t in targets}
    gkf = GroupKFold(n_splits=n_folds)

    # Inner ES split (~10%) from training fold to prevent gating overfit
    for fold, (tr_idx, va_idx) in enumerate(gkf.split(X_feat, groups=groups)):
        n_es = max(1, int(0.1 * len(tr_idx)))
        es_idx, pure_tr_idx = tr_idx[:n_es], tr_idx[n_es:]
        
        # Train gating on pure_tr, early-stop on es_idx, predict on va_idx
        # ... (train GatingNetwork, checkpoint on ES score) ...
        gated_oof[target][va_idx] = gated_predictions
    
    return gated_oof
```

**CRITICAL:** Never report in-sample gating score alongside honest OOF scores. They are not comparable.

---

## Submission Clipping Optimisation

Before submitting, find per-target clip ranges that maximise OOF score:

```python
def optimize_submission_clips(oof: dict, y_dict: dict) -> tuple[dict, dict]:
    """
    Find (lo, hi) clip per target that maximises competition score on OOF.
    Returns separate clips for AUC-driven and LL-driven targets.
    """
    from scipy.optimize import minimize_scalar

    clip_auc = {}
    clip_ll  = {}
    for t in TARGETS:
        def neg_score(hi):
            clipped = np.clip(oof[t], 1e-6, hi)
            return -_comp_score(y_dict[t], clipped, t)
        
        res = minimize_scalar(neg_score, bounds=(0.5, 1.0), method="bounded")
        clip_auc[t] = res.x
        clip_ll[t]  = 1.0 - res.x   # symmetric low clip for LL benefit

    return clip_auc, clip_ll
```

---

## Meta run() Interface

```python
def run(
    use_pseudo: bool = False,
    base_tags: dict | None = None,
    models: list | None = None,
    pseudo_exclude: list | None = None,
    steps: list | None = None,   # None = all; or e.g. ["gating", "weights_gating"]
    force: bool = False,         # True = recompute even if OOF cache exists
) -> dict[str, float]:
    """
    Returns dict of honest OOF scores keyed by sub-step name.
    Keys: "weights", "stack", "gating", "weights_gating"
    Only sub-steps in `steps` that ran will appear in the dict.
    """
    # [A] Nelder-Mead blend of all base OOFs      → scores["weights"]
    # [B] LogReg stacking (GroupKFold OOF)        → scores["stack"]
    # [C] Dynamic gating (honest GroupKFold OOF)  → scores["gating"]
    # [D] Nelder-Mead blend of base OOFs + gating OOF as 9th input
    #     → scores["weights_gating"]
    return scores
```

### Continuation / Sub-step Caching

Each sub-step checks `oof/meta_meta_{name}.pkl` on disk before running:

```python
# Skip if cache exists and force=False
if not force and Path(f"oof/meta_meta_{name}.pkl").exists():
    oof, tst, score = _load_from_cache(name)
    continue

# Useful CLI patterns:
# Only recompute gating + blend (loads weights from cache for [D]):
#   python train.py --model meta --meta_steps gating weights_gating
# Force-recompute just the blend:
#   python train.py --model meta --meta_steps weights_gating --meta_force
```

### The weights_gating Blend (Level 4)

[D] blends **all base OOFs + the honest gating OOF** as an additional 9th input via Nelder-Mead:

```
Inputs to [D]: cat, lgb, xgb, nn, pseudo_cat, pseudo_lgb, pseudo_xgb, pseudo_nn, gating
```

Key insight: the gating OOF is treated as a *meta-feature* in the final blend, not a standalone prediction. The optimiser assigns it whatever weight truly improves OOF. In practice, when the gating network doesn't generalise well (small dataset, ~13K rows), it gets low weight (~0.17) and the blend falls back to something close to the plain `weights` blend.

**Note on reporting:** Never put in-sample gating scores next to honest OOF scores — see [common-pitfalls.md](./common-pitfalls.md) #12.

---

## Loading OOF Files

```python
def load_oof(model_name: str, tag: str = "base") -> dict:
    """Returns dict with keys 'oof' and 'test'."""
    path = os.path.join(cfg.oof_dir, f"{tag}_{model_name}.pkl")
    with open(path, "rb") as f:
        return pickle.load(f)

# Usage
d = load_oof("cat", tag="base")
oof_cat = d["oof"]   # {target: np.array(n_train,)}
tst_cat = d["test"]  # {target: np.array(n_test,)}
```

---

## See Also

| File | Why |
|------|-----|
| [pseudo-labeling.md](./pseudo-labeling.md) | Pseudo OOF files are consumed as inputs to the ensemble |
| [experiment-tracking.md](./experiment-tracking.md) | Logging ensemble OOF scores and deciding which variant to submit |
| [submission-postprocessing.md](./submission-postprocessing.md) | Post-process the ensemble output before generating the submission file |
| [common-pitfalls.md](./common-pitfalls.md) | Pitfall #12 — in-sample gating score mixed with honest OOF |
