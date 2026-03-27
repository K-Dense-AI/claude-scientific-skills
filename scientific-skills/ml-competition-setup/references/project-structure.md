# Project Structure & Config

## Overview

A well-structured project is the prerequisite for every other improvement. Without consistent file locations, a `RunConfig` singleton, and a YAML-driven orchestrator, every experiment requires manual coordination — making results non-reproducible and pipeline bugs hard to isolate.

This file documents: the canonical package layout; the `RunConfig` dataclass singleton pattern; the YAML orchestrator with resume-by-existence logic; and naming conventions for OOF files, submission files, and tuning outputs.

**Two principles govern the structure:**
1. **Single source of truth for config** — all user-tunable knobs live in `config.yaml` and are loaded into `RunConfig` at startup; no magic constants scattered across trainer files
2. **Resume by existence** — the orchestrator skips any step whose output file already exists; this means you can interrupt and rerun `train.py` at any point without duplicating work or corrupting results

**When to use:** When starting a new competition project, when onboarding code from a notebook to a pipeline, or when debugging "why did my train.py re-run step X when it already produced output?"

---

## Package Layout

```
project/
├── train.py                  # YAML-driven orchestrator
├── tune.py                   # Optuna orchestrator
├── config.yaml               # All user-facing knobs
├── cache/
│   └── features_vN.pkl       # Versioned feature cache
├── oof/
│   └── {tag}_{model}.pkl     # OOF + test preds per model
├── submissions/
│   └── {tag}_{model}_{ts}.csv
├── tuning/
│   ├── {model}_best.json     # {"value": 0.959, "params": {...}}
│   └── {model}_trials.csv
└── package/
    ├── base/
    │   ├── config.py         # RunConfig dataclass + singleton
    │   ├── common.py         # TARGETS, LB constants, load/save helpers
    │   ├── features.py       # engineer_features() + build_model_matrices()
    │   ├── metrics.py        # Competition metric wrappers (CB/LGB/XGB/NN)
    │   ├── lgb_trainer.py    # Stateless LGB fold engine
    │   ├── xgb_trainer.py    # Stateless XGB fold engine
    │   └── nn_trainer.py     # Stateless NN fold engine
    ├── train/
    │   ├── cat.py            # CatBoost entrypoint
    │   ├── lgb.py            # LGB entrypoint
    │   ├── xgb.py            # XGB entrypoint
    │   ├── nn.py             # NN entrypoint
    │   ├── pseudo.py         # Pseudo-label retraining
    │   ├── meta.py           # Ensemble meta (weights / stack / gating / weights_gating)
    │   ├── meta_gating.py    # Dynamic gating network (honest GroupKFold OOF)
    │   └── calibrate.py      # Post-hoc isotonic calibration on any meta OOF
    └── tune/
        ├── common.py         # run_study(), make_scorer(), load_results()
        ├── tune_cat.py       # make_cb_objective()
        ├── tune_lgb.py       # make_lgb_objective()
        ├── tune_xgb.py       # make_xgb_objective()
        └── tune_nn.py        # make_nn_objective()
```

---

## RunConfig Singleton

**Golden rule**: import `cfg` everywhere; call `init_config()` exactly once at startup in `train.py` / `tune.py`.

```python
# base/config.py
from dataclasses import dataclass, field

@dataclass
class RunConfig:
    data_dir:   str  = "data/raw"
    oof_dir:    str  = "oof"
    sub_dir:    str  = "submissions"
    feat_cache: str  = "cache/features_v1.pkl"
    n_folds:    int  = 5
    seeds: list = field(default_factory=lambda: [42, 123, 777, 2024, 9999])

cfg: RunConfig = RunConfig()  # singleton

def init_config(yaml_dict: dict) -> RunConfig:
    """Mutate singleton in-place — all importers see updated values."""
    cfg.data_dir   = yaml_dict.get("data_dir",   cfg.data_dir)
    cfg.oof_dir    = yaml_dict.get("oof_dir",     cfg.oof_dir)
    cfg.feat_cache = yaml_dict.get("feat_cache",  cfg.feat_cache)
    cfg.n_folds    = yaml_dict.get("n_folds",     cfg.n_folds)
    cfg.seeds      = yaml_dict.get("seeds",       cfg.seeds)
    return cfg
```

Usage in any module:
```python
from package.base.config import cfg
gkf = GroupKFold(n_splits=cfg.n_folds)
```

---

## config.yaml

```yaml
tag: base
tune_dir: tuning/run1
# feat_cache: cache/features_v1.pkl   ← bump vN on every feature change

# Per-model tune dirs (override tune_dir for specific models)
tune_dirs:
  cat: tuning/full_lb
  lgb: tuning/full_lb

# CV settings
n_folds: 5
seeds: [42, 123, 777, 2024, 9999]
feat_cache: cache/features_v1.pkl

train: [cat, lgb, xgb, nn]

pseudo:
  enabled: true
  models: [cat, lgb, xgb, nn]

meta:
  enabled: true
  use_pseudo: false
  models: [cat, lgb, xgb, nn]
  pseudo_exclude: []
  steps: null       # null = all sub-steps; or list e.g. [gating, weights_gating]
  force: false      # true = recompute even if OOF cache exists

calibrate:
  enabled: false              # off by default — verify on OOF before enabling; can hurt when ensemble already smooths scores
  source_model: meta_weights
  source_tag: meta

auto:
  enabled: false
  min_score: 0.94   # drop models below this tuning score
  top_n: null       # keep only top N models
```

---

## YAML Orchestrator (train.py) Pattern

Key design decisions:
1. `_deep_merge(base, override)` — YAML over defaults, CLI over YAML
2. `_init_config(cfg)` called **after** all overrides are applied
3. Resume by existence: `oof/{tag}_{model}.pkl` → skip and score from disk
4. `_model_tune_dir(model)` → per-model dir → fallback to global `tune_dir`
5. Steps resolved canonically: `["cat","lgb","xgb","nn","pseudo","meta"]` order always

```python
# Resume pattern — always load + score so final summary is complete
oof_path = Path(cfg["oof_dir"]) / f"{cfg['tag']}_{step}.pkl"
if oof_path.exists():
    _d = pickle.load(open(oof_path, "rb"))
    _sc, _ = weighted_lb_score(_M["y_dict"], _d["oof"])
    scores[step] = _sc
    continue
```

**Known limitation**: `_deep_merge` skips explicit YAML `null` values because of `elif v is not None` guard — explicit nulls cannot override defaults to None via YAML.

---

## base/common.py — Constants and Data Loading

Define all competition-specific constants here, import them everywhere else.

```python
# base/common.py

# ── Competition constants ─────────────────────────────────────
TARGETS = [...]          # list of target column names — define per competition
GROUP_COL = "..."        # group column for GroupKFold

# For metric scaling (AUC+LL blend type):
# NULL_LL_SCALING = <null-model log-loss — compute: log_loss(y_train, [y_mean]*n)>

# ── Data loading ──────────────────────────────────────────────
def load_data():
    """Load train / test DataFrames. Add auxiliary sources as needed."""
    train = pd.read_csv(f"{cfg.data_dir}/train.csv")
    test  = pd.read_csv(f"{cfg.data_dir}/test.csv")
    # Optional: parse date columns, load auxiliary data
    # aux = pd.read_csv(f"{cfg.data_dir}/aux.csv")
    return train, test   # or return train, test, aux
```

**Conventions:**
- TARGETS is a list even for single-target competitions — keeps all training code uniform
- For multiclass: OOF shape is `{target: np.zeros((n_train, n_classes))}` instead of `{target: np.zeros(n_train)}`
- GROUP_COL is `None` for competitions with no natural group; use `StratifiedKFold` in that case
- Auxiliary data is loaded in `load_data()` and passed explicitly — never inferred inside trainers
- NULL_LL_SCALING is only needed for AUC+LL blend metrics; remove it for pure RMSE/AUC/F1 competitions

---

## Feature Cache Discipline

- Bump `features_vN` version **every time** `engineer_features()` changes
- Process-level in-memory matrix cache: `_matrices_cache: dict | None = None` global
- `build_model_matrices()` is called once per process, returns all arrays for all model families
- Feature cache is shared between `train.py` and `tune.py` — same pkl path

```python
# common.py — process-level matrix cache
_matrices_cache = None

def load_or_build_matrices(verbose=True):
    global _matrices_cache
    if _matrices_cache is not None:
        return _matrices_cache    # zero recomputation on subsequent calls
    _matrices_cache = build_model_matrices(...)
    return _matrices_cache
```

---

## See Also

| File | Why |
|------|-----|
| [model-training.md](./model-training.md) | Trainer files that live in `package/base/` |
| [feature-engineering.md](./feature-engineering.md) | Feature cache versioning and `feat_cache` config key |
| [hyperparameter-tuning.md](./hyperparameter-tuning.md) | `tuning/` directory layout and YAML orchestrator integration |
| [experiment-tracking.md](./experiment-tracking.md) | `oof/` and `submissions/` directory naming conventions |
