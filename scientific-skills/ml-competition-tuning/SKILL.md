---
name: ml-competition-tuning
description: "Set up and run Optuna hyperparameter tuning for tabular ML competition models (CatBoost, LightGBM, XGBoost, Neural Network). Use when: writing or debugging Optuna objective functions; configuring tune_config.yaml; implementing run_study() and per-model search spaces; fixing the load_tuned_params contract bug (most common tuning error); ensuring tuner folds match training folds. NOT for base model training, feature engineering, or ensemble."
argument-hint: "Describe your task: e.g. 'set up Optuna for LGB multiclass', 'fix load_tuned_params returning wrong value', 'implement CatBoost search space', 'tune NN learning rate and depth'"
license: MIT
metadata:
    skill-author: eak
---

# ML Competition — Hyperparameter Tuning

## Overview

This skill covers Optuna-based hyperparameter tuning for all model families used in tabular competitions.

**The single most common tuning bug:** `load_tuned_params` returns the `params` sub-dict instead of the full JSON. The caller then calls `params[t].update(tuned)`, which silently merges `"value"`, `"trial"`, and `"timestamp"` keys into model params — causing either a runtime error or, worse, silently wrong training.

**Two non-negotiable rules:**
1. Tuner folds must use **identical CV splits** as training (same `GroupKFold`, same group column)
2. Tuner must use the **exact same competition metric** as training — never a surrogate like log-loss when the competition uses a weighted AUC+LL blend

**When to run tuning:** After base models produce reasonable OOF scores. Run before pseudo-labeling and ensemble — tuned base models produce better pseudo labels and ensemble inputs.

---

## load_tuned_params Contract — Critical

```python
# common.py — always return the FULL json dict
def load_tuned_params(model_name, tune_dir):
    path = os.path.join(tune_dir, f"{model_name}_best.json")
    data = json.load(open(path))
    return data          # {"value": 0.959, "params": {...}, "trial": 48, ...}

# caller — extract params separately; log value
tuned = load_tuned_params("cat", tune_dir)
hp = tuned.get("params", {})
params[t].update(hp)
logger.info(f"loaded tuned params (value={tuned.get('value', '?')})")

# ❌ BUG — merges "value", "trial", "timestamp" into model params
# params[t].update(tuned)
```

The `_best.json` format must always be `{"value": float, "params": {...}, "trial": int}`. Never save just the params dict.

---

## Tuner Fold Rules

- Tuner: **3-fold** for speed (acceptable approximation)
- Final training: **5-fold × 5 seeds** (full quality)
- Both must use the **same CV strategy**: if training uses `GroupKFold(group_col="entity_id")`, tuner must too
- Tuner metric must be **identical** to the training eval metric — same `competition_score` function, same `GroupKFold` split

---

## Reference Files

| File | What it covers |
|------|----------------|
| [hyperparameter-tuning.md](./references/hyperparameter-tuning.md) | Optuna setup, `run_study()`, `load_tuned_params` contract, per-model search spaces (CB/LGB/XGB/NN) |

---

## See Also

| Skill | When to use it instead |
|-------|------------------------|
| `ml-competition` | Full pipeline overview, task type decision guide, first-principles checklist |
| `ml-competition-setup` | Project structure, RunConfig, process management |
| `ml-competition-training` | Base model training, metrics, output format |
| `ml-competition-features` | Feature engineering, validation strategy |
| `ml-competition-advanced` | Pseudo-labeling, ensemble, post-processing, experiment tracking |
| `ml-competition-quality` | Coding rules, common pitfalls |
