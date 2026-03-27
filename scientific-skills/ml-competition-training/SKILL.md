---
name: ml-competition-training
description: "Implement and review base model training for tabular ML competitions: CatBoost, LightGBM, XGBoost, and Neural Networks. Use when: writing or reviewing trainer files and model entrypoints; implementing per-framework competition metric wrappers (CB/LGB/XGB/NN); aligning early-stopping metric with the competition objective; identifying correct prediction type for submission (probability vs class label vs value); debugging train OOF vs LB gaps caused by wrong metric or output format; handling auxiliary/prior data. NOT for hyperparameter tuning, ensemble, or pseudo-labeling."
argument-hint: "Describe your task: e.g. 'implement CatBoost binary trainer', 'fix early stopping metric for multiclass', 'debug why submission scores ~0.5 despite good OOF', 'add auxiliary data to LGB'"
license: MIT
metadata:
    skill-author: eak
---

# ML Competition — Model Training, Metrics & Output Format

## Overview

This skill covers three tightly coupled concerns for base model training:

1. **Model training** — stateless trainer architecture for CatBoost/LGB/XGB/NN; correct params by task type; training objective vs eval metric separation; auxiliary data safety rules
2. **Competition metrics** — `competition_score` single source of truth; per-framework wrappers; the CatBoost logit vs prediction API difference that silently destroys scores
3. **Output format** — metric determines prediction type; deriving submission format from `sample_submission.csv`; OOF collection patterns

**The two most expensive silent bugs:**
- Wrong eval metric in early stopping → model stops at wrong iteration → wasted training
- Submitting class labels when metric expects probabilities → 0.91 AUC model scores ~0.5

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

---

## Training Objective vs Eval Metric — Two Separate Things

| | Tree models (CB/LGB/XGB) | Neural Network |
|---|---|---|
| **Training objective** | Built-in BCE / squared-error (internal to framework) | Custom loss: `FocalLoss`, `SmoothBCE`, `MSELoss`, `BCELoss` |
| **Eval metric** | Custom competition metric wrapper | Competition score computed per epoch |
| **Early stopping driven by** | Eval metric | Competition score |

Getting the eval metric wrong wastes all training — early stopping fires at the wrong iteration.

### Framework injection summary

| Framework | Metric injection | Input type in callback |
|-----------|-----------------|------------------------|
| **CatBoost binary** | `eval_metric=CatBoostCompMetric()` | Raw logits → apply `sigmoid` |
| **CatBoost regression** | `eval_metric=CatBoostCompMetric()` | Raw predictions → use directly |
| **CatBoost multiclass** | `eval_metric=CatBoostCompMetric()` | K logit arrays → apply `softmax` |
| **LightGBM** | `"metric": "None"` + `feval=make_lgb_feval()` | Already-transformed probs/values |
| **XGBoost** | `"disable_default_eval_metric": 1` + `custom_metric=make_xgb_eval()` + `maximize=True` | Already-transformed probs/values |
| **NN** | Compute in epoch loop; `model.load_state_dict(best_state)` | You control activation in `forward()` |

**Single source of truth**: define `competition_score(y_true, y_pred) -> float` once in `base/metrics.py`. `competition_score` always maximizes — negate RMSE/MAE/logloss when the leaderboard is "lower is better".

**Critical CatBoost API difference:**
- `CatBoostClassifier`: `approxes[0]` = raw **logits** → must apply `sigmoid` (binary) or `softmax` (multiclass)
- `CatBoostRegressor`: `approxes[0]` = raw **prediction values** → no activation needed

---

## Output Format — Submission Prediction Type

> **The metric determines the prediction type. The target column values in `train.csv` do NOT.**

| Metric | Prediction type | `predict` call |
|--------|----------------|----------------|
| AUC-ROC, PR-AUC, Log Loss, Brier | **float [0, 1]** | `model.predict_proba(X)[:, 1]` |
| Accuracy, F1, Cohen's Kappa, QWK | **class label** | `model.predict(X)` |
| RMSE, MAE, RMSLE, MAPE | **continuous value** | `model.predict(X)` |
| Multi-class log loss | **prob matrix (n, K)** | `model.predict_proba(X)` |
| MAP@K, NDCG@K | **ranked list** | task-specific |

### Hard rules

- **Always derive submission format from `sample_submission.csv`**, not `train.csv` target dtype
- **RMSLE**: clip predictions to ≥ 0 before saving; do NOT log-transform before saving
- **Multiclass log-loss**: column order in submission must match `sample_submission.csv` exactly
- **OOF arrays**: always collect full arrays — `oof[val] = model.predict_proba(X[val])[:, 1]`, never fold-by-fold averages

### Scout checklist — before writing any submission code

1. What is the **exact metric name**? (Quote from README.)
2. Is it a **probability metric**? → float predictions required.
3. What are the **exact column names** in `sample_submission.csv`?
4. What **value types** do those columns contain? (Check `sample_submission.csv`, not `train.csv`.)
5. Copy 2–3 **example rows** verbatim from `sample_submission.csv`.

---

## Prior / Auxiliary Data Safety Rules

- CatBoost handles heterogeneous cardinality and distribution shift better → auxiliary data often helps
- LGB / XGB / NN **may overfit severely** when auxiliary data comes from a different distribution; verify OOF vs LB gap before committing
- When using combined splits: fold on combined data, but score OOF on `va_idx[va_idx < n_train]` only
- **Never add two changes at once** — add auxiliary data alone, submit, then decide

---

## Reference Files

| File | What it covers |
|------|----------------|
| [model-training.md](./references/model-training.md) | CB/LGB/XGB/NN params by task type, trainer architecture, training objective vs eval metric |
| [competition-metrics.md](./references/competition-metrics.md) | `competition_score` pattern, per-framework metric wrappers (CB/LGB/XGB/NN), training losses |
| [output-format.md](./references/output-format.md) | Metric → prediction type table, submission format by task, OOF collection patterns, scout checklist |

---

## See Also

| Skill | When to use it instead |
|-------|------------------------|
| `ml-competition` | Full pipeline overview, task type decision guide, first-principles checklist |
| `ml-competition-setup` | Project structure, RunConfig, process management |
| `ml-competition-features` | Feature engineering, validation strategy |
| `ml-competition-tuning` | Optuna hyperparameter tuning |
| `ml-competition-advanced` | Pseudo-labeling, ensemble, post-processing, experiment tracking |
| `ml-competition-quality` | Coding rules, common pitfalls |
