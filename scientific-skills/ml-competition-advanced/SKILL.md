---
name: ml-competition-advanced
description: "Implement and debug late-stage improvements for tabular ML competitions: pseudo-labeling, ensemble meta-learning, submission post-processing, and experiment tracking. Use when: adding pseudo-labeling after base models converge; building ensemble (weighted blend / LogReg stacking / dynamic gating); calibrating predictions (Platt / isotonic); clipping or constraining output values; tracking OOF scores and diagnosing OOF vs LB divergence; deciding which submission to make final. NOT for base model training, feature engineering, or hyperparameter tuning."
argument-hint: "Describe your task: e.g. 'implement pseudo-labeling for regression', 'build weighted blend ensemble', 'calibrate multiclass predictions', 'diagnose OOF vs LB gap', 'decide submission'"
license: MIT
metadata:
    skill-author: eak
---

# ML Competition — Advanced: Pseudo-Labeling, Ensemble, Post-Processing & Tracking

## Overview

This skill covers the four late-stage pipeline components applied after base models converge:

1. **Pseudo-labeling** — when/how to use test predictions as additional training labels; per-task label generation; confidence gating; weight decay over rounds; pitfalls
2. **Ensemble meta-learning** — weighted blend (Nelder-Mead), LogReg stacking, dynamic gating, and the weights+gating hybrid; honest OOF discipline
3. **Submission post-processing** — calibration (Platt/isotonic), OOF-optimized clipping, domain constraints; YAML toggle; when calibration hurts vs helps
4. **Experiment tracking** — score ledger, using OOF as LB proxy, diagnosing OOF vs LB divergence, submission decision logic

**Order of operations**: base models → hyperparameter tuning → pseudo-labeling → ensemble → post-processing. Violating this order contaminates the OOF used for ensemble training.

---

## Prior / Auxiliary Data Rules

- CatBoost handles heterogeneous cardinality and distribution shift better → auxiliary data often helps
- LGB / XGB / NN **may overfit severely** when auxiliary data comes from a different distribution; always verify OOF vs LB gap before committing
- When using combined splits: fold on combined data, but score OOF on `va_idx[va_idx < n_train]` only
- **Never add two changes at once** — add auxiliary data alone, submit, then decide

---

## Calibration — When It Helps vs Hurts

Add calibration only as an optional post-meta pass, **disabled by default**. Enable only after OOF gain is confirmed.

**Key danger**: For ensembles of 4+ models, the ensemble already softens extremes and calibration often hurts (−0.006 OOF observed on a 4-model binary blend). Always verify on OOF before enabling.

```yaml
# config.yaml
calibration:
  enabled: false   # only enable after confirming OOF gain
  method: isotonic  # or platt
```

---

## First-Principles Checklist — Late-Stage Steps

These depend on all earlier pipeline steps being complete and stable:

11. [ ] Write **OOF save/load** with versioned tags
12. [ ] Add **pseudo-labeling** only after base models converge
13. [ ] Add **meta ensemble** after pseudo converges
14. [ ] Wire everything in **YAML orchestrator** with resume-by-existence logic
15. [ ] Add **calibration step** as optional post-meta pass — verify OOF gain, keep `enabled: false` by default

---

## Reference Files

| File | What it covers |
|------|----------------|
| [pseudo-labeling.md](./references/pseudo-labeling.md) | When/how/weight, per-task label generation, confidence check, pitfalls |
| [ensemble-meta.md](./references/ensemble-meta.md) | Weighted blend (Nelder-Mead), LogReg stacking, dynamic gating, weights+gating hybrid |
| [submission-postprocessing.md](./references/submission-postprocessing.md) | Calibration (Platt/isotonic), OOF-optimised clipping, domain constraints, YAML toggle |
| [experiment-tracking.md](./references/experiment-tracking.md) | Score ledger, OOF as LB proxy, OOF vs LB divergence diagnosis, submission decision logic |

---

## See Also

| Skill | When to use it instead |
|-------|------------------------|
| `ml-competition` | Full pipeline overview, task type decision guide, first-principles checklist |
| `ml-competition-setup` | Project structure, RunConfig, process management |
| `ml-competition-features` | Feature engineering, validation strategy |
| `ml-competition-training` | Base model training, competition metrics, output format |
| `ml-competition-tuning` | Optuna hyperparameter tuning |
| `ml-competition-quality` | Coding rules, common pitfalls |
