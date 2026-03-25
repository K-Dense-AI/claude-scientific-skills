---
name: ml-competition
description: "Build, debug, review and improve ML competition pipelines for tabular competitions (Kaggle, Zindi, etc.) covering binary classification, multiclass, regression, multi-label, and ranking tasks. Use when: identifying competition type and setting up the metric; structuring a competition codebase from scratch; reviewing training code for correctness bugs; implementing CatBoost/LightGBM/XGBoost/NN base models; debugging metric leakage or train-LB gaps; adding pseudo-labeling; setting up Optuna hyperparameter tuning with YAML config; building ensemble (weighted blend / LogReg stacking / dynamic gating); managing OOF files; aligning early-stopping metrics with competition objective; adding auxiliary/prior data safely; engineering and selecting features; post-processing and calibrating predictions. NOT for NLP/CV competitions."
argument-hint: "Describe your task: e.g. 'review training code for bugs', 'implement pseudo-labeling for regression', 'set up Optuna for LGB multiclass', 'debug LB gap'"
license: MIT
metadata:
    skill-author: eak
---

# ML Competition Pipeline Skill

## Overview

This is the core skill for tabular ML competitions (Kaggle, Zindi, and equivalents). It covers the full pipeline from data ingestion to final submission: feature engineering, model training (CatBoost / LightGBM / XGBoost / NN), Optuna hyperparameter tuning, pseudo-labeling, ensemble meta-learning, and submission post-processing.

**Critical principle: every layer of the pipeline owns exactly one thing.** Config drives knobs. Trainers are stateless fold engines. Metrics are defined once and reused everywhere. Violating these boundaries is the single most common source of silent bugs.

This skill is **not** for NLP or CV competitions — those require separate skill sets for tokenization, transformer architectures, and image augmentation pipelines.

## When to Use This Skill

Use this skill when:
- Starting or structuring a new competition codebase from scratch
- Reviewing any training or tuning code before wasting GPU hours
- Implementing or modifying a model family (CatBoost, LGB, XGB, NN)
- Debugging a train OOF vs LB submission gap
- Setting up or fixing Optuna tuning — metric alignment, `best.json` format
- Implementing pseudo-labeling or ensemble meta-learners
- Safely incorporating auxiliary or prior-period data
- Engineering, selecting, or caching features
- Post-processing, calibrating, or clipping predictions before submission

---

## Task Type Decision Guide

**Identify your task type first** — it determines which parameters, objectives, and techniques apply.

| Task Type | Framework objectives | Binary-only params to remove |
|-----------|---------------------|------------------------------|
| **Binary classification** | XGB: `binary:logistic` · LGB: `binary` · CB: `CatBoostClassifier` | — (all apply) |
| **Regression** | XGB: `reg:squarederror` · LGB: `regression` · CB: `CatBoostRegressor` | `scale_pos_weight`, `is_unbalance`, `auto_class_weights: Balanced`, threshold=0.5 |
| **Multiclass** | XGB: `multi:softprob` · LGB: `multiclass` · CB: `MultiClass` loss | `scale_pos_weight`, `is_unbalance` |
| **Multi-label** | N independent binary models OR single NN with N sigmoid heads | — (each target is binary) |
| **Ranking** | XGB: `rank:pairwise` · LGB: `lambdarank` · CB: `YetiRank` | all imbalance params |

**TARGETS pattern:**
- Single target: `TARGETS = ["target"]` — dict-of-arrays still works, simplest to keep uniform
- Multi-target / multi-label: `TARGETS = ["target_a", "target_b", "target_c"]` — dict-of-arrays everywhere, one model per target per fold
- Multiclass: `TARGETS = ["label"]` — predictions are 2D `(n, n_classes)`; OOF shape changes

**CatBoost API difference by task:**
- Binary/multiclass: `approxes[0]` = raw **logits** → must apply sigmoid (binary) or softmax (multiclass)
- Regression: `approxes[0]` = raw **prediction values** → no activation needed

---

## Separation of Concerns — Non-Negotiable

Every layer owns exactly one thing. Never blur boundaries.

| Layer | Owns | Key files |
|-------|------|-----------|
| **Config** | All user-tunable knobs | `base/config.py` (RunConfig dataclass + singleton) |
| **Features** | Raw → engineered DataFrames, versioned pkl cache | `base/features.py`, `cache/features_vN.pkl` |
| **Matrices** | NumPy/DataFrame arrays per model family | `base/features.py::build_model_matrices()` |
| **Metrics** | Competition metric, framework-specific wrappers | `base/metrics.py` |
| **Trainers** | Stateless fold-loop engines, injectable params | `base/lgb_trainer.py`, `xgb_trainer.py`, `nn_trainer.py` |
| **Entrypoints** | Load → call trainer → score → save OOF | `train/cat.py`, `lgb.py`, `xgb.py`, `nn.py` |
| **Tuning** | Optuna objectives, run_study, save JSON | `tune/tune_*.py`, `tune.py` |
| **Pseudo** | Retrain on train + pseudo-labeled test | `train/pseudo.py` |
| **Meta** | Weighted blend / stacking / gating on OOF | `train/meta.py`, `train/meta_gating.py` |
| **Orchestrator** | YAML-driven sequential runner + resume logic | `train.py` |

---

## Coding Standards — Non-Negotiable

Apply to every `src/*.py` and `scripts/*.py` file. These are hard quality gates, not suggestions.

| Rule | What it means |
|------|---------------|
| **No dead code** | No unused imports, variables, parameters, or uncalled private helpers (`_name`) |
| **Clear contracts** | Every public function has explicit input/output docs; optional params state their default behavior |
| **Single responsibility** | Functions do one job — if a function loads data AND computes metrics AND logs, split it |
| **Explicit types and names** | Type hints on all signatures; use `train_df`, `class_weights` not `tmp`, `d`, `x1` |
| **Predictable data handling** | Validate required columns up front and `raise` with a specific message; never silently mutate |
| **Structured logging** | Use `logger.*`; no `print`, no commented-out debug blocks, no stale TODOs |

**Quick review gate** — before committing any Python change, verify:
1. No unused imports, variables, or parameters
2. No uncalled private helper functions
3. Function signatures match actual behavior (not aspirational behavior)
4. Error messages are specific and actionable
5. Logs are concise and useful for iteration debugging

See [coding-rules.md](./references/coding-rules.md) for anti-patterns and extended examples.

---

## Process Management — Non-Negotiable

Training scripts can run for hours. The most expensive mistake is launching a duplicate process.

### Before every training launch — pre-flight check

```bash
# 1. Is training already running?
RUNNING_PIDS=$(pgrep -f "python scripts/train.py" 2>/dev/null)
[ -n "$RUNNING_PIDS" ] && echo "⚠️  Already running — PIDs: $RUNNING_PIDS" && exit 0

# 2. Are artifacts already fresh? (< 5 min old)
[ -f "artifacts/oof.npy" ] && \
  ARTIFACT_AGE=$(( $(date +%s) - $(stat -c %Y artifacts/oof.npy) )) && \
  [ $ARTIFACT_AGE -lt 300 ] && echo "✅ Artifacts are fresh — skip retraining"
```

### Correct launch pattern

```bash
# TRAIN_PID=$! MUST be on its own line immediately after & — not on the same line
nohup uv run python scripts/train.py > train.log 2>&1 &
TRAIN_PID=$!
echo "Training started — PID: $TRAIN_PID"

# Wait loop
while kill -0 $TRAIN_PID 2>/dev/null; do sleep 60; tail -5 train.log; done
echo "✅ Done"; tail -50 train.log
```

### Hard rules

- **Never launch without the pre-flight check** — two processes writing the same `oof.npy` corrupt results silently
- **`TRAIN_PID=$!` must be on its own line** — assigning on the same compound line as `&` captures an empty string
- **Never pipe `train.py` to `head`** — `head` closes the pipe and kills the process; always redirect to a log file
- **Fast exit ≠ failure** — check artifact timestamps and CPU load before concluding anything went wrong
- **Kill before relaunch** — if code changed and an old process is running, kill it first, then relaunch

See [process-management.md](./references/process-management.md) for all four workflow patterns (pre-flight, launch/wait, kill/relaunch, diagnose fast exit).

---

## Competition Metrics — Key Rules

Two completely independent concepts — never confuse them:

| Concept | Who uses it | What it is |
|---------|------------|------------|
| **Training loss** | Optimizer | Differentiable objective: BCE, Focal, SmoothBCE, MSELoss, CrossEntropyLoss |
| **Eval metric** | Early stopping | **Exact competition formula** — must match the leaderboard precisely |

**Getting the eval metric wrong wastes all training** — early stopping fires at the wrong iteration.

### Framework injection summary

| Framework | Metric injection | Input type in callback |
|-----------|-----------------|------------------------|
| **CatBoost binary** | `eval_metric=CatBoostCompMetric()` | Raw logits → apply `sigmoid` |
| **CatBoost regression** | `eval_metric=CatBoostCompMetric()` | Raw predictions → use directly |
| **CatBoost multiclass** | `eval_metric=CatBoostCompMetric()` | K logit arrays → apply `softmax` |
| **LightGBM** | `"metric": "None"` + `feval=make_lgb_feval()` | Already-transformed probs/values |
| **XGBoost** | `"disable_default_eval_metric": 1` + `custom_metric=make_xgb_eval()` + `maximize=True` | Already-transformed probs/values |
| **NN** | Compute in epoch loop; `model.load_state_dict(best_state)` | You control activation in `forward()` |

**Single source of truth**: define `competition_score(y_true, y_pred) -> float` once in `base/metrics.py`; never inline the formula elsewhere. `competition_score` always maximizes — negate RMSE/MAE/logloss when the leaderboard is "lower is better".

See [competition-metrics.md](./references/competition-metrics.md) for full implementations of every wrapper.

---

## Output Format — Submission Prediction Type

> **The metric determines the prediction type. The target column values in `train.csv` do NOT.**

The #1 source of silent score destruction: submitting class labels when the metric expects probabilities. A 0.91 AUC model submitting "Yes/No" strings scores ~0.5 — indistinguishable from random.

| Metric | Prediction type | `predict` call |
|--------|----------------|----------------|
| AUC-ROC, PR-AUC, Log Loss, Brier | **float [0, 1]** | `model.predict_proba(X)[:, 1]` |
| Accuracy, F1, Cohen's Kappa, QWK | **class label** | `model.predict(X)` |
| RMSE, MAE, RMSLE, MAPE | **continuous value** | `model.predict(X)` |
| Multi-class log loss, MCLL | **prob matrix (n, K)** | `model.predict_proba(X)` |
| MAP@K, NDCG@K | **ranked list** | task-specific |

### Identifying the metric from competition text

| Phrase in README | Metric | Type |
|-----------------|--------|------|
| "AUC", "ROC AUC", "AUROC" | AUC-ROC | probability |
| "Log Loss", "Logloss", "Cross-Entropy" | Log loss | probability |
| "Average Precision", "PR-AUC" | AP | probability |
| "Brier Score" | Brier | probability |
| "Accuracy", "F1", "F-beta" | Accuracy / F1 | class label |
| "Quadratic Weighted Kappa", "QWK" | QWK | ordinal label |
| "RMSE", "MAE", "RMSLE", "MAPE" | Regression | value |
| "NDCG", "MAP@K" | Ranking | score / ranked list |

### Hard rules

- **Always derive submission format from `sample_submission.csv`**, never from `train.csv` target dtype — train targets may be strings ("Yes"/"No") while submission requires floats
- **RMSLE**: clip predictions to ≥ 0 before saving; do NOT log-transform before saving
- **Multiclass log-loss**: column order in submission must match `sample_submission.csv` exactly
- **OOF arrays**: always collect full arrays, never fold-by-fold averages — `oof[val] = model.predict_proba(X[val])[:, 1]`

### Scout checklist — before writing any submission code

1. What is the **exact metric name**? (Quote from README.)
2. Is it a **probability metric**? → float predictions required.
3. What are the **exact column names** in `sample_submission.csv`?
4. What **value types** do those columns contain? (Check `sample_submission.csv`, not `train.csv`.)
5. Copy 2–3 **example rows** verbatim from `sample_submission.csv`.

See [output-format.md](./references/output-format.md) for metric-by-metric sklearn calls and OOF collection patterns.

---

## First-Principles Checklist (new competition)

Complete in this exact order — each depends on the previous:

1. [ ] Identify **task type** (binary/regression/multiclass/multi-label/ranking) and read the leaderboard formula carefully
2. [ ] Implement **competition metric** exactly — single function, unit-tested against sample submission scores
3. [ ] Wrap metric for every framework (CatBoost class / LGB feval / XGB custom_metric / NN monitor)
4. [ ] Identify **group column**: use `GroupKFold` when a natural integrity unit (entity ID, session ID) exists; use `StratifiedKFold` for i.i.d. data with class imbalance; use `KFold` only when rows are truly independent and balanced; use `TimeSeriesSplit` for temporal competitions.
5. [ ] Build **feature pipeline** with versioned cache; bump version on any change
6. [ ] Build **model matrices** — single `build_model_matrices()`, process-level cache
7. [ ] Create **RunConfig dataclass** singleton; `init_config(yaml_dict)` at startup
8. [ ] Write **base model entrypoints** with `tune_dir=None` param
9. [ ] Set up **load_tuned_params** returning full JSON dict (see [hyperparameter-tuning.md](./references/hyperparameter-tuning.md))
10. [ ] Set up **Optuna** with identical folds/metrics to training
11. [ ] Write **OOF save/load** with versioned tags
12. [ ] Add **pseudo-labeling** only after base models converge
13. [ ] Add **meta ensemble** after pseudo converges
14. [ ] Wire everything in **YAML orchestrator** with resume-by-existence logic
15. [ ] Add **calibration step** as an optional post-meta pass — but verify on OOF first. For ensembles of 4+ models the ensemble already softens extremes and calibration often hurts (e.g. −0.006 OOF observed on a 4-model binary blend). Keep `enabled: false` by default; enable only when OOF gain is confirmed. See [submission-postprocessing.md](./references/submission-postprocessing.md).

---

## Critical Rules (hard constraints)

### Training loss vs eval metric — two separate things
- **Training loss** (what the optimizer minimizes): tree models use their built-in objective; NNs use `FocalLoss`, `SmoothBCE`, `MSELoss`, or `BCELoss` depending on task type
- **Eval metric** (what drives early stopping): must be the **exact competition formula** for every model
- These are completely independent — getting the eval metric wrong wastes all your training

See [competition-metrics.md](./references/competition-metrics.md) for the generic `competition_score` pattern and correct per-framework implementation.

**Critical API difference (task-type dependent):**
- CatBoost **classification**: `approxes[0]` = raw logits → must apply sigmoid (binary) or softmax (multiclass)
- CatBoost **regression**: `approxes[0]` = raw prediction → no activation, use directly
- LGB and XGB: always pass already-transformed outputs (probabilities for classification, raw values for regression) → no activation needed

### Validation hygiene
- `GroupKFold` when a natural integrity unit exists (user ID, entity ID, session ID)
- `StratifiedKFold` for i.i.d. rows with class imbalance (no group column)
- `KFold` only when rows are truly independent and the target is balanced
- For time-series: use `TimeSeriesSplit` or a rolling-window split — never GroupKFold when temporal ordering matters
- OOF arrays are **train-rows only** — even if combined data drives fold splits
- Tuner folds can be 3-fold for speed; final training uses 5-fold × 5 seeds
- Target encoding and any statistics must be computed fold-by-fold inside cross-validation

### load_tuned_params contract
```python
# common.py — always return the FULL json dict
def load_tuned_params(model_name, tune_dir):
    data = json.load(open(path))
    return data          # {"value": 0.959, "params": {...}, "trial": 48, ...}

# caller — extract params separately; log value
tuned = load_tuned_params("cat", tune_dir)
hp = tuned.get("params", {})
params[t].update(hp)
logger.info(f"loaded (value={tuned.get('value', '?')})")
# BUG if you do: params[t].update(tuned)  ← merges "value", "trial" etc. into params
```

### Prior / auxiliary data
- CatBoost handles heterogeneous cardinality and distribution shift better than other frameworks → auxiliary data often helps
- LGB / XGB / NN **may overfit severely** when auxiliary data comes from a different distribution than test; always verify OOF vs LB gap before committing
- When using combined splits: fold on combined, but score OOF on `va_idx[va_idx < n_train]` only
- **Never add two changes at once** — add auxiliary data alone, submit, then decide

---

## Reference Files

Each file begins with an **Overview** explaining its scope and when to use it, and ends with a **See Also** table of directly related files.

| File | What it covers |
|------|----------------|
| [project-structure.md](./references/project-structure.md) | Package layout, `RunConfig` singleton, YAML orchestrator, resume-by-existence logic |
| [model-training.md](./references/model-training.md) | CB/LGB/XGB/NN params by task type, trainer architecture, training objective vs eval metric |
| [competition-metrics.md](./references/competition-metrics.md) | `competition_score` pattern, per-framework metric wrappers (CB/LGB/XGB/NN), training losses |
| [hyperparameter-tuning.md](./references/hyperparameter-tuning.md) | Optuna setup, `run_study()`, `load_tuned_params` contract, per-model search spaces |
| [validation-strategy.md](./references/validation-strategy.md) | GroupKFold / TimeSeriesSplit, OOF accumulation, leakage prevention, leakage checklist |
| [feature-engineering.md](./references/feature-engineering.md) | Encoding strategies, datetime features, aggregations, feature selection, cache discipline |
| [ensemble-meta.md](./references/ensemble-meta.md) | Weighted blend (Nelder-Mead), LogReg stacking, dynamic gating, weights+gating hybrid |
| [pseudo-labeling.md](./references/pseudo-labeling.md) | When/how/weight, per-task label generation, confidence check, pitfalls |
| [submission-postprocessing.md](./references/submission-postprocessing.md) | Calibration (Platt/isotonic), OOF-optimised clipping, domain constraints, YAML toggle |
| [output-format.md](./references/output-format.md) | Metric → prediction type table, submission format by task, OOF collection patterns, scout checklist |
| [experiment-tracking.md](./references/experiment-tracking.md) | Score ledger, OOF as LB proxy, OOF vs LB divergence diagnosis, submission decision logic |
| [common-pitfalls.md](./references/common-pitfalls.md) | 16 production bugs with ❌ / ✅ patterns — read before finalizing any component |
| [coding-rules.md](./references/coding-rules.md) | No dead code, clear contracts, single responsibility, explicit types, structured logging |
| [process-management.md](./references/process-management.md) | Pre-flight checks, PID tracking, launch/wait/kill patterns, fast-exit diagnosis |
