# Common Pitfalls (Hard-Won Lessons)

## Overview

This file catalogs bugs and mistakes discovered in production competition pipelines. Each entry follows the same format: **Bug** (what went wrong and why it is silent), **Fix** (the correct pattern with code), and a reference to the canonical source for the topic. Every entry here represents a real score loss — read them before starting any new pipeline component.

**How to use:** Before finalizing any trainer, tuner, or ensemble component, scan this file. If your code matches a ❌ pattern, fix it before running.

---

These are bugs and mistakes discovered in production. Never repeat them.

---

## 1. load_tuned_params returns params sub-dict instead of full JSON

**Bug:**
```python
# WRONG — returns only {"lr": 0.02, "depth": 5, ...}
def load_tuned_params(model_name, tune_dir):
    return json.load(open(path)).get("params", {})

# In caller:
tuned = load_tuned_params("cat", tune_dir)
logger.info(f"value={tuned.get('value', '?')}")  # logs "value=?" always
params[t].update(tuned)  # accidentally merges "value", "trial" keys into params
```

**Fix:**
```python
# CORRECT — return full dict; caller extracts params separately
def load_tuned_params(model_name, tune_dir):
    return json.load(open(path))   # {"value": 0.959, "params": {...}}

# Caller:
tuned = load_tuned_params("cat", tune_dir)
hp = tuned.get("params", {})
params[t].update(hp)
logger.info(f"value={tuned.get('value', '?')}")  # now shows 0.959
```

---

## 2. Early stopping on wrong metric

**Bug:** Model trains with `metric="binary_logloss"` (default) while competition metric is weighted AUC + LL. Early stopping at wrong iteration, often 10x too early or too late.

**Fix:** See [competition-metrics.md](./competition-metrics.md). Every framework requires explicit metric suppression + custom metric injection.

- CatBoost: `eval_metric=CatBoostCompMetric()` — **not** `"Logloss"`
- LGB: `"metric": "None"` + `feval=make_lgb_feval()`
- XGB: `"disable_default_eval_metric": 1` + `custom_metric=make_xgb_eval()` + `maximize=True`

---

## 3. CatBoost eval_metric as string vs custom class

**Bug:** `early_stopping_rounds` with string `eval_metric="Logloss"` in pseudo training but `CatBoostCompMetric` in base training. Different metrics → different best iterations → inconsistent ensemble.

**Fix:** Use `CatBoostCompMetric()` everywhere — base training AND pseudo training. Never use `"Logloss"` as early stop metric when the competition uses a custom metric.

---

## 4. Adding auxiliary / prior-period data to LGB/XGB/NN

**Bug:** Adding large auxiliary rows (e.g. prior-period data) to the training set for LGB/XGB/NN causes severe overfitting when the auxiliary data comes from a different distribution than test. OOF may improve while LB drops significantly.

**Root cause:** Auxiliary data from a different time period, source, or group composition. Tree models overfit to the combined distribution which does not match the test set distribution. CatBoost is more robust to distribution shift due to its oblivious tree structure and categorical handling.

**Fix:** Verify OOF vs LB gap before committing any data addition. For LGB/XGB/NN, use train-only rows as the default. If auxiliary data helps OOF but hurts LB, remove it.

**Rule:** Never add two changes at once. Add auxiliary data alone, submit, then decide.

---

## 5. OOF accumulated wrong when using combined splits

**Bug:**
```python
# WRONG — oof indexed on combined (larger) array
oof[t][va_idx] += predictions / len(seeds)
# va_idx can be > n_train when splitting on combined data → IndexError or wrong OOF
```

**Fix:**
```python
# CORRECT
va_train_idx = va_idx[va_idx < n_train]  # restrict to actual train rows
oof[t][va_train_idx] += predictions / len(seeds)
```

---

## 6. Pseudo partial resume with mismatched labels

**Scenario:** pseudo_cat.pkl exists, pseudo_lgb.pkl doesn't. If base_cat was retrained since pseudo_cat was saved, the label ensemble (avg of 4 base models' test preds) now uses new cat predictions but pseudo_cat.pkl used old ones.

**Risk:** pseudo_cat OOF is from a different label set than pseudo_lgb OOF → meta ensemble degrades.

**Accepted tech debt — workaround:** Delete all `pseudo_*.pkl` files whenever any base model OOF changes, then re-run pseudo fully:
```bash
rm oof/pseudo_*.pkl
python train.py --model pseudo meta
```

---

## 7. Scale_pos_weight computed from full data instead of fold

**Bug:**
```python
# WRONG — uses full dataset class balance, not fold-specific balance
n_pos = (y_all == 1).sum()
n_neg = (y_all == 0).sum()
params["scale_pos_weight"] = n_neg / n_pos
```

**Fix:**
```python
# CORRECT — per-fold balance
n_pos = max((y_tr == 1).sum(), 1)
n_neg = (y_tr == 0).sum()
params["scale_pos_weight"] = n_neg / n_pos
```

---

## 8. Feature cache version not bumped

**Bug:** `features_v3.pkl` used from cache after `engineer_features()` was changed. All models trained on stale features.

**Fix:** Increment `vN` in `feat_cache` config every time feature engineering logic changes. The cache path is set in `config.yaml` → `feat_cache: cache/features_v4.pkl`.

---

## 9. Tuner score value=? in logs

**Symptom:** `loaded params from 'tuning/cat_best.json' (value=?)`

**Cause:** See bug #1 — `load_tuned_params` returning only `params` sub-dict.

---

## 10. XGB maximize=True missing

**Bug:** `xgb.train(..., custom_metric=make_xgb_eval(t))` without `maximize=True`. XGBoost treats the custom metric as minimise-by-default → early stops at wrong direction → model stops too early.

**Fix:** Always pair `custom_metric` with `maximize=True`.

---

## 11. LGB metric="None" missing

**Bug:** `params` has `feval` set but `"metric"` not set to `"None"`. LGB evaluates BOTH the default `binary_logloss` AND the feval. `first_metric_only=True` in early stopping then watches the default logloss, not the competition metric.

**Fix:** `params["metric"] = "None"` is mandatory when using custom feval.

---

## 12. in-sample gating score mixed with honest OOF

**Bug:** `sc_g` from fitting the gating model on the full training set (not cross-validated) placed in `scores` dict alongside honest OOF scores like `sc_w=0.963` and `sc_s=0.964`.

**Risk:** Misleading — in-sample gating will appear better than honest scores. You may submit a miscalibrated gating model.

**Fix:** Always use `run_dynamic_gating_oof()` which computes honest GroupKFold OOF. Mark in-sample scores clearly if you compute them at all.

---

## 13. --tune_dir CLI clobbers per-model tune_dirs

**Bug:**
```python
# In train.py CLI override:
if args.tune_dir is not None:
    cfg["tune_dir"] = args.tune_dir
    for _m in BASE_MODELS:
        cfg["tune_dirs"][_m] = args.tune_dir   # overwrites per-model YAML settings
```

If `config.yaml` has `tune_dirs: {cat: tuning/v2, lgb: tuning/v1}` and you pass `--tune_dir tuning/v3`, all per-model dirs get silently overwritten.

**Accepted tech debt — workaround:** Never use `--tune_dir` when `tune_dirs` is set per-model in config.yaml. Use per-model YAML keys instead, or pass per-model tags explicitly.

---

## 15. Removing scale_pos_weight from XGBoost pseudo retraining

**Bug:** Attempting to "clean up" XGBoost pseudo training by removing `scale_pos_weight` computation because pseudo labels add negatives at a different ratio:

```python
# WRONG — removing scale_pos_weight entirely
params = dict(XGB_PARAMS, seed=seed + fold)  # no scale_pos_weight
```

**Result:** On a 1% positive-rate competition, this collapses pseudo_xgb OOF from ~0.928 to **0.825** — a catastrophic drop. The model loses its ability to detect the minority class.

**Root cause:** XGBoost's internal loss function treats all samples equally without `scale_pos_weight`. With 99% negatives, it learns to predict near-zero for everything.

**Fix:** Always compute `scale_pos_weight` from the **real training fold labels only** (not the pseudo rows) — same as `xgb_trainer.py`:
```python
n_neg = (y_dict[t][tr_idx] == 0).sum()
n_pos = max((y_dict[t][tr_idx] == 1).sum(), 1)
params = dict(XGB_PARAMS, scale_pos_weight=n_neg / n_pos, seed=seed + fold)
```

**Rule:** Never touch parameters that exist in the base trainer without first verifying they're intentionally different. Compare pseudo training against `xgb_trainer.py` before any param changes.

---

## 16. _deep_merge cannot process explicit YAML null

**Bug:**
```python
def _deep_merge(base, override):
    for k, v in override.items():
        elif v is not None:   # ← YAML null (None) silently ignored
            result[k] = v
```

If you want to override `tune_dir: "tuning/v1"` back to `null` in YAML, it won't work.

**Workaround:** Use a sentinel string like `"none"` and handle it in the resolver.

---

## See Also

| File | Why |
|------|-----|
| [competition-metrics.md](./competition-metrics.md) | Correct metric injection patterns that prevent pitfalls #2, #3, #10, #11 |
| [validation-strategy.md](./validation-strategy.md) | OOF accumulation patterns that prevent pitfalls #5, #12 |
| [hyperparameter-tuning.md](./hyperparameter-tuning.md) | Full `load_tuned_params` contract (pitfall #1 and #9) |
| [pseudo-labeling.md](./pseudo-labeling.md) | Correct pseudo loop that prevents pitfalls #6 and #15 |
| [feature-engineering.md](./feature-engineering.md) | Cache versioning rules that prevent pitfall #8 |
