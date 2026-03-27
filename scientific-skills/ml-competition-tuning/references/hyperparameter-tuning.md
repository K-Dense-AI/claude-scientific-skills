# Hyperparameter Tuning with Optuna

## Overview

Hyperparameter tuning with Optuna systematically searches the parameter space to find configurations that maximize the competition metric. This file covers: the tuning architecture (`tune.py` → `run_study()` → per-model objectives); the `tune_config.yaml` format; the `load_tuned_params` contract (the single most common source of tuning bugs); and search spaces for CatBoost, LightGBM, XGBoost, and NN.

**The most common tuning bug:** `load_tuned_params` returns the `params` sub-dict instead of the full JSON. The caller then calls `params[t].update(tuned)`, which silently merges `"value"`, `"trial"`, and `"timestamp"` keys into the model params — causing either a runtime error or, worse, silently wrong training. See the [load_tuned_params contract](#load_tuned_params-contract--critical) section and [common-pitfalls.md](./common-pitfalls.md) #1 for the fix.

**Two non-negotiable rules for tuning correctness:**
1. Tuner folds must use **identical CV splits** as training (same `GroupKFold`, same group column)
2. Tuner must use the **exact same competition metric** as training — never a surrogate like log-loss when the competition uses a weighted AUC+LL blend

**When to use:** After base models are implemented and producing reasonable OOF scores. Run tuning before pseudo-labeling and ensemble — tuned base models produce better pseudo labels and ensemble inputs.

---

## Architecture

```
tune.py  (orchestrator)
  → loads tune_config.yaml
  → builds features once (shared across all runs)
  → calls run_study() per model
       → make_{model}_objective() returns Optuna objective fn
       → saves {model}_best.json + {model}_trials.csv
```

---

## tune_config.yaml

```yaml
feat_cache: cache/features_v1.pkl

runs:
  - name: cat_full_lb
    model: cat
    trials: 60
    objective: full_lb
    out_dir: tuning/full_lb

  - name: lgb_full_lb
    model: lgb
    trials: 60
    objective: full_lb
    out_dir: tuning/full_lb
```

---

## run_study() Pattern

```python
def run_study(label, objective_fn, n_trials, pruner=None, out_dir="tuning"):
    sampler = TPESampler(seed=SEED_TUNE)
    study = optuna.create_study(
        direction="maximize",
        sampler=sampler,
        pruner=pruner or optuna.pruners.NopPruner(),
    )
    study.optimize(objective_fn, n_trials=n_trials, show_progress_bar=True)
    
    # Save best config
    best = study.best_trial
    best_data = {
        "model": label,
        "trial": best.number,
        "value": best.value,
        "params": best.params,
        "n_trials": len(study.trials),
        "timestamp": datetime.now().strftime("%Y%m%d_%H%M%S"),
    }
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    with open(f"{out_dir}/{label}_best.json", "w") as f:
        json.dump(best_data, f, indent=2)
    
    # Save all trials to CSV
    df = study.trials_dataframe()
    df.to_csv(f"{out_dir}/{label}_trials.csv", index=False)
    
    return [best_data]   # top-5 list for results.json
```

---

## load_tuned_params Contract — CRITICAL

```python
# base/common.py — return the FULL json dict, not just params sub-dict
def load_tuned_params(model_name: str, tune_dir: str) -> dict | None:
    path = os.path.join(tune_dir, f"{model_name}_best.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)   # returns {"value": 0.959, "params": {...}, ...}

# Caller (train/lgb.py, train/cat.py, etc.) — always extract params separately
tuned = load_tuned_params("lgb", tune_dir)
if not tuned:
    return LGB_BASE_PARAMS
hp = tuned.get("params", {})             # ← extract params sub-dict
for k, v in hp.items():
    params[_LGB_KEY_MAP.get(k, k)] = v   # apply with key translation
logger.info(f"loaded (value={tuned.get('value', '?')})")  # ← logs the score
```

**Bug to never repeat:**  
If `load_tuned_params` returns only `data.get("params", {})`, then `tuned.get("value")` returns `None` and you log `value=?`. The fix is to return the full dict and let callers do `tuned["params"]`.

---

## Optuna Objective Structure

```python
def make_lgb_objective(matrices, train, score_fn):
    """Return an Optuna objective function for LightGBM."""
    X = matrices["X_num_train"]
    y_dict = matrices["y_dict"]
    groups = matrices["groups"]
    gkf = GroupKFold(n_splits=N_FOLDS_TUNE)   # 3-fold for speed

    def objective(trial: optuna.Trial) -> float:
        params = {
            "learning_rate": trial.suggest_float("lr", 1e-3, 0.1, log=True),
            "num_leaves":    trial.suggest_int("num_leaves", 20, 200),
            "max_depth":     trial.suggest_int("max_depth", 3, 12),
            # ... more params
            "metric": "None",          # SAME as training
            "device_type": "gpu",      # SAME as training
        }
        
        oof = {t: np.zeros(len(train)) for t in TARGETS}
        for seed in [SEED_TUNE]:
            for fold, (tr_idx, va_idx) in enumerate(gkf.split(X, groups=groups)):
                for t in TARGETS:
                    # ... train fold, predict val
                    oof[t][va_idx] += predictions
        
        score, _ = score_fn(y_dict, oof)
        return score

    return objective
```

### Tuning folds vs training folds
- Tuning: 3 folds × 1 seed (speed)
- Training: 5 folds × 5 seeds (quality)
- This discrepancy is acceptable — tuning finds good hyperparameter region; seeds add variance reduction

---

## CatBoost Tuning — Optuna Key Names

Use shorthand that matches the `_build_params()` translator exactly:

```python
params_trial = {
    "depth":           trial.suggest_int("depth", 4, 10),
    "lr":              trial.suggest_float("lr", 1e-3, 0.1, log=True),
    "l2":              trial.suggest_float("l2", 1.0, 30.0, log=True),
    "min_leaf":        trial.suggest_int("min_leaf", 10, 300),
    "rsm":             trial.suggest_float("rsm", 0.5, 1.0),
    "random_strength": trial.suggest_float("random_strength", 0.1, 5.0),
    "bootstrap_type":  trial.suggest_categorical("bootstrap_type", ["MVS", "Bayesian"]),
    "subsample":       trial.suggest_float("subsample", 0.5, 1.0),
    "iters": 3000,     # fixed at tuning time
}
```

The `_build_params(hp, seed, cat_idx, target)` in `train/cat.py` translates:
- `iters` → `iterations`
- `lr` → `learning_rate`
- `l2` → `l2_leaf_reg`
- `min_leaf` → `min_data_in_leaf`
- All other keys pass through as-is

---

## NN Tuning — Pruning

```python
from optuna.pruners import MedianPruner

pruner = MedianPruner(n_startup_trials=8, n_warmup_steps=15)
# n_startup_trials: don't prune in first 8 trials
# n_warmup_steps: don't prune before epoch 15

# In objective, report intermediate values for pruning
for epoch in range(epochs):
    # ... train epoch
    trial.report(val_score, epoch)
    if trial.should_prune():
        raise optuna.TrialPruned()
```

---

## Results File Format

```json
// {model}_best.json
{
  "model": "lgb",
  "trial": 42,
  "value": 0.9527,
  "params": {
    "lr": 0.0209,
    "num_leaves": 57,
    "max_depth": 10,
    "min_child_samples": 114
  },
  "n_trials": 60,
  "timestamp": "20260226_155633"
}
```

```json
// results.json — dict of model → top-5 list (for resuming multi-run studies)
{
  "cat": [{"model": "cat", "trial": 48, "value": 0.9593, "params": {...}}],
  "lgb": [{"model": "lgb", "trial": 42, "value": 0.9527, "params": {...}}]
}
```

---

## See Also

| File | Why |
|------|-----|
| [model-training.md](./model-training.md) | Trainer parameter sets that the tuner must mirror |
| [competition-metrics.md](./competition-metrics.md) | Metric wrappers reused in tuner objectives — must be identical to training |
| [common-pitfalls.md](./common-pitfalls.md) | Pitfall #1 (`load_tuned_params` full dict contract), #9 (value=? symptom) |
| [project-structure.md](./project-structure.md) | `tuning/` directory conventions and `tune_config.yaml` location |
