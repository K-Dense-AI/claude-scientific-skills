# Experiment Tracking

## Overview

Experiment tracking is the practice of recording every OOF evaluation and every LB submission in a structured log. Without it, competition work degrades into guesswork: you lose track of which feature version or ensemble weight produced a given LB score, you submit stale predictions, and you cannot tell whether a new change is genuinely better than a previous one.

This file covers: the score ledger format; using OOF as a real-time leaderboard proxy; diagnosing OOF vs LB divergence; the decision framework for what to submit; and the reproducibility checklist.

**The core discipline:** OOF score is computed from the saved `.pkl` file — never from training logs or in-memory arrays. Training logs may represent a mid-epoch or mid-fold state. The `.pkl` always represents the final, post-averaged fold result.

**When to use:** Throughout the competition, for every experiment. Refer to the divergence table whenever your OOF and LB scores do not move together — this is the earliest diagnostic signal for leakage, metric bugs, or distribution shift.

---

## Why Tracking Matters

ML competitions involve many parallel experiments: different models, feature versions, ensemble weights, post-processing settings. Without structured tracking you will lose results, submit stale predictions, or forget which version produced your best LB score.

---

## Score Ledger (the only source of truth)

Keep a running `experiments.md` or `experiments.csv` in the project root. Log **every submission** and **every OOF evaluation** you intend to submit.

```
| Date       | Tag            | OOF Score | LB Score | Notes                              |
|------------|----------------|-----------|----------|------------------------------------|
| YYYY-MM-DD | base_cat       | 0.XXXX    | 0.XXXX   | baseline CatBoost, no tuning       |
| YYYY-MM-DD | base_lgb       | 0.XXXX    | 0.XXXX   | LGB v1, tuned 60 trials            |
| YYYY-MM-DD | ensemble_v1    | 0.XXXX    | 0.XXXX   | Nelder-Mead blend cat+lgb+xgb+nn   |
| YYYY-MM-DD | pseudo_cat     | 0.XXXX    | 0.XXXX   | pseudo threshold=0.5, weight=0.5   |
| YYYY-MM-DD | ensemble_v2    | 0.XXXX    | 0.XXXX   | ensemble on pseudo models          |
```

**Fields to always log:**
- `Tag`: reproducible name matching the OOF pkl file (e.g. `base_cat`, `pseudo_lgb`)
- `OOF Score`: computed from saved pkl — not from training logs (logs may be mid-epoch)
- `LB Score`: filled in after submission
- `Notes`: what changed vs previous — one concise sentence

---

## OOF Score Is Your Real-Time LB Proxy

You have limited LB submissions per day. Use OOF as the primary signal:

```python
# Always compute OOF score from the saved pkl, not from memory
from package.base.common import load_oof, compute_oof_score, print_scores

d = load_oof("cat", tag="base")
score, per = compute_oof_score(y_dict, d["oof"])
print_scores(y_dict, d["oof"], label="base_cat")
```

**Rule:** Only submit to LB when OOF shows a meaningful gain (≥ 0.001 for most competitions; ≥ 0.0003 for high-variance metrics like AUC).

---

## When OOF and LB Diverge

| Pattern | Diagnosis |
|---------|-----------|
| OOF >> LB consistently (≥ 0.005 gap) | Validation leakage — check group column, target encoding |
| OOF improves but LB flat | Feature helps only in val distribution — likely leakage |
| OOF flat but LB improves | OOF metric implementation is wrong — recheck formula |
| OOF and LB improve together | Changes are genuine — trust this signal |
| Adding auxiliary data: OOF ↑ but LB ↓ | Auxiliary data from the wrong distribution — revert || Meta gating OOF >> meta weights OOF by ≥ 0.003 | In-sample gating score mixed with honest OOF — see [common-pitfalls.md](./common-pitfalls.md) #12. Always use `run_dynamic_gating_oof()`. |
| Calibration: OOF drops after calibration step | Expected for 4+ model ensembles — ensemble already softens extremes. Keep `calibrate.enabled: false` by default; enable only when OOF gain is confirmed. |
**Diagnosis workflow:**
1. Compute OOF on each model individually — compare to known-good baseline
2. Check `OOF - LB` gap is stable across experiments (a sudden jump indicates a bug)
3. If LB moves and OOF doesn't: your OOF formula diverges from the actual leaderboard formula — recheck `competition_score`

---

## Deciding What to Submit

Follow this priority order:

1. **Single best model** — establish a non-ensemble LB baseline first
2. **Best base ensemble** — Nelder-Mead blend of all base models
3. **Best pseudo ensemble** — only after pseudo OOF gain ≥ 0.001
4. **Calibrated / clipped variant** — submit the calibrated version only if OOF gain justifies a submission slot

**At end of competition (final 2 submissions):**
- Submit 1: best single model (hedge against ensemble overfit)
- Submit 2: best ensemble (maximum expected performance)

---

## Reproducibility Checklist

Before closing a run that produced a promising submission:

```
[ ] OOF pkl saved: oof/{tag}_{model}.pkl
[ ] Submission CSV saved: submissions/{tag}_{model}_{timestamp}.csv
[ ] Config YAML snapshot saved (or committed to git)
[ ] Random seeds logged (seeds list in config)
[ ] Feature cache version noted (feat_cache: cache/features_vN.pkl)
[ ] Score logged in experiments.md with LB score
```

---

## Git-Based Experiment Branches (optional but recommended)

```bash
# Create experiment branch from main
git checkout -b exp/pseudo-cat-v2

# Train, evaluate, log score
# If it's worse, discard:
git checkout main
git branch -D exp/pseudo-cat-v2

# If it's better, merge or cherry-pick the pkl and config:
git checkout main
git checkout exp/pseudo-cat-v2 -- oof/pseudo_cat.pkl config.yaml
```

---

## Resume Logic (YAML Orchestrator)

The orchestrator skips completed steps using the OOF pkl as the existence check:

```python
oof_path = Path(cfg["oof_dir"]) / f"{cfg['tag']}_{step}.pkl"
if oof_path.exists():
    d = pickle.load(open(oof_path, "rb"))
    sc, _ = compute_oof_score(y_dict, d["oof"])
    scores[step] = sc
    logger.info(f"  [skip] {step} (score={sc:.6f})")
    continue
```

**To re-run a step:** delete its pkl file and re-run. NEVER edit pkl files manually.

**Danger: stale pseudo pkls.** If any base model is retrained after pseudo pkls exist, the pseudo labels no longer correspond to the current base model OOF. Delete all `pseudo_*.pkl` files before re-running pseudo.

---

## Logging Pattern

Use structured logging in every training script:

```python
from loguru import logger

logger.add("logs/{time}.log", level="INFO", rotation="50 MB")

# At start of run
logger.info(f"[{step}] tag={cfg.tag}  tune_dir={tune_dir}  seeds={cfg.seeds}")

# After each fold
logger.info(f"  fold={fold}  val_score={va_score:.6f}  best_iter={best_iter}")

# After full model
logger.info(f"  [{step}] OOF = {oof_score:.6f}")
```

Keep logs for every experiment — they contain the per-fold scores needed to diagnose variance issues.

---

## See Also

| File | Why |
|------|-----|
| [common-pitfalls.md](./common-pitfalls.md) | Pitfall #12 (in-sample gating inflates OOF score) |
| [validation-strategy.md](./validation-strategy.md) | The CV strategy determines whether OOF is an honest LB proxy |
| [submission-postprocessing.md](./submission-postprocessing.md) | Final step before filing; post-processed OOF score should be logged here |
