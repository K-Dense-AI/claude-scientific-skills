# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

"""genomic_cv_python.py

Cross-validation methods for genomic prediction (Python).

Includes:
- K-fold CV
- Stratified K-fold CV (by subpopulation)
- Leave-one-out CV
- Forward (chronological) validation
- Relatedness-aware CV (kinship-threshold grouping)
- GroupKFold by known family IDs (recommended when pedigree/family labels exist)
- RR-BLUP via marker ridge regression
- GBLUP via kernel ridge regression on the genomic relationship matrix (GRM)

Run:
  python genomic_cv_python.py

Dependencies:
  pip install numpy pandas scikit-learn

Optional (for larger/alternative GRM computations): scipy
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from collections import deque

from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import (
    BaseCrossValidator,
    GroupKFold,
    KFold,
    LeaveOneOut,
    StratifiedKFold,
)


def _regression_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict:
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    corr = np.corrcoef(y_true, y_pred)[0, 1]
    rmse = float(np.sqrt(mean_squared_error(y_true, y_pred)))
    r2 = float(r2_score(y_true, y_pred))
    return {"correlation": float(corr), "rmse": rmse, "r2": r2}


# =============================================================================
# Simulated genotype/phenotype data
# =============================================================================


def simulate_genomic_data(
    n_individuals: int = 400,
    n_markers: int = 2000,
    n_qtl: int = 50,
    n_subpopulations: int = 3,
    n_families: int = 60,
    family_freq_sd: float = 0.03,
    n_years: int = 6,
    year_effect_sd: float = 0.5,
    heritability: float = 0.5,
    seed: int = 42,
) -> dict:
    """Simulate genomic data with population structure, family relatedness, and time.

    Population structure:
      Subpopulations have different allele frequencies.

    Relatedness:
      Individuals are grouped into families; family members share family-specific
      allele frequencies around their subpopulation frequencies.

    Time:
      Each individual has an integer time label (e.g., year/generation) with an
      additive year effect on phenotype.

    Returns a dict containing standardized genotype matrix X, phenotype y,
    subpopulation labels, family IDs, time labels, and a VanRaden-style GRM.
    """

    rng = np.random.default_rng(seed)

    if not (1 <= n_families <= n_individuals):
        raise ValueError("n_families must be in [1, n_individuals]")
    if not (0 < heritability < 1):
        raise ValueError("heritability must be in (0, 1)")
    if n_subpopulations < 1:
        raise ValueError("n_subpopulations must be >= 1")
    if n_years < 2:
        raise ValueError("n_years must be >= 2")

    # Assign each family to a subpopulation (roughly balanced)
    family_subpops = np.repeat(np.arange(n_subpopulations), n_families // n_subpopulations)
    rem = n_families % n_subpopulations
    if rem:
        family_subpops = np.concatenate([family_subpops, np.arange(rem)])
    family_subpops = family_subpops[:n_families]
    rng.shuffle(family_subpops)

    # Family sizes that sum to n_individuals
    base = n_individuals // n_families
    family_sizes = np.full(n_families, base, dtype=int)
    family_sizes[: (n_individuals - base * n_families)] += 1
    rng.shuffle(family_sizes)

    family_ids = np.concatenate([np.full(family_sizes[f], f, dtype=int) for f in range(n_families)])
    rng.shuffle(family_ids)
    subpopulations = family_subpops[family_ids]

    # Base allele frequencies
    base_freqs = rng.uniform(0.1, 0.9, n_markers)
    subpop_freqs = []
    for _ in range(n_subpopulations):
        deviation = rng.normal(0, 0.10, n_markers)
        subpop_freqs.append(np.clip(base_freqs + deviation, 0.05, 0.95))

    # Genotypes coded 0/1/2
    M = np.zeros((n_individuals, n_markers), dtype=float)
    for f in range(n_families):
        sp = int(family_subpops[f])
        freqs = subpop_freqs[sp]
        family_freqs = np.clip(freqs + rng.normal(0, family_freq_sd, n_markers), 0.05, 0.95)
        idx = np.where(family_ids == f)[0]
        # sample each member independently from the family allele frequencies
        M[idx, :] = rng.binomial(2, family_freqs, size=(idx.size, n_markers))

    # Standardize markers for ridge (RR-BLUP on markers)
    M_centered = M - M.mean(axis=0)
    sd = M_centered.std(axis=0, ddof=1)
    sd[sd == 0] = 1.0
    X = M_centered / sd

    # Choose QTL and simulate additive genetic value
    qtl_idx = rng.choice(n_markers, size=n_qtl, replace=False)
    qtl_effects = rng.normal(0, 1, n_qtl)
    g = X[:, qtl_idx] @ qtl_effects

    var_g = float(np.var(g))
    var_e = var_g * (1 - heritability) / heritability

    time = np.repeat(np.arange(n_years), n_individuals // n_years)
    rem2 = n_individuals % n_years
    if rem2:
        time = np.concatenate([time, np.arange(rem2)])
    time = time[:n_individuals]
    rng.shuffle(time)
    year_effect = rng.normal(0, year_effect_sd, n_years)

    e = rng.normal(0, np.sqrt(var_e), n_individuals)
    y = g + year_effect[time] + e

    # VanRaden GRM (Method 1): G = (M - 2p)(M - 2p)' / (2 * sum p(1-p))
    p = M.mean(axis=0) / 2.0
    denom = 2.0 * float(np.sum(p * (1 - p)))
    W = (M - 2.0 * p) / np.sqrt(denom)
    G = (W @ W.T) / n_markers

    return {
        "genotypes": X,
        "phenotypes": y,
        "subpopulations": subpopulations,
        "family_ids": family_ids,
        "time": time,
        "kinship": G,
        "qtl_indices": qtl_idx,
    }


# =============================================================================
# CV splitters
# =============================================================================


class RelatednessAwareCV(BaseCrossValidator):
    """Kinship-threshold grouping CV.

    1) Build a graph where i~j if kinship[i,j] > threshold
    2) Connected components define "related groups" (e.g., families/clusters)
    3) Split groups into folds; keep each group entirely in train or test

    This is a simple, practical leakage-avoidance splitter when you do NOT have
    explicit family/pedigree IDs.
    """

    def __init__(self, kinship_matrix: np.ndarray, threshold: float = 0.125, n_splits: int = 5, random_state: int = 42):
        self.kinship_matrix = np.asarray(kinship_matrix)
        self.threshold = float(threshold)
        self.n_splits = int(n_splits)
        self.random_state = int(random_state)

    def get_n_splits(self, X=None, y=None, groups=None):
        return self.n_splits

    def _connected_components(self) -> list[list[int]]:
        K = self.kinship_matrix
        n = K.shape[0]
        adj = K > self.threshold
        np.fill_diagonal(adj, False)

        visited = np.zeros(n, dtype=bool)
        comps: list[list[int]] = []

        for i in range(n):
            if visited[i]:
                continue
            comp: list[int] = []
            q = deque([i])
            visited[i] = True
            while q:
                cur = q.popleft()
                comp.append(int(cur))
                neigh = np.where(adj[cur])[0]
                for j in neigh:
                    if not visited[j]:
                        visited[j] = True
                        q.append(int(j))
            comps.append(comp)

        return comps

    def split(self, X, y=None, groups=None):
        rng = np.random.default_rng(self.random_state)
        comps = self._connected_components()
        n_groups = len(comps)

        group_order = np.arange(n_groups)
        rng.shuffle(group_order)

        # Assign groups to folds as evenly as possible
        folds = np.array_split(group_order, self.n_splits)

        n = len(X)
        for fold_groups in folds:
            test_idx = np.concatenate([np.array(comps[g], dtype=int) for g in fold_groups])
            mask = np.ones(n, dtype=bool)
            mask[test_idx] = False
            train_idx = np.where(mask)[0]
            yield train_idx, test_idx


# =============================================================================
# Models
# =============================================================================


def rrblup_ridge_fit_predict(X_train, y_train, X_test, alpha=1.0):
    """RR-BLUP via ridge regression on markers."""
    model = Ridge(alpha=float(alpha))
    model.fit(X_train, y_train)
    return model.predict(X_test)


class GBLUP:
    """GBLUP as kernel ridge regression on the GRM.

    y_hat_test = mu + G_test,train @ (G_train + lambda I)^-1 @ (y_train - mu)

    Here lambda_param acts like sigma_e^2 / sigma_g^2.
    """

    def __init__(self, lambda_param: float = 1.0):
        self.lambda_param = float(lambda_param)
        self.mu_ = None
        self.alpha_ = None

    def fit(self, y_train: np.ndarray, G_train: np.ndarray):
        y_train = np.asarray(y_train)
        G_train = np.asarray(G_train)
        if G_train.shape[0] != G_train.shape[1]:
            raise ValueError("G_train must be square")
        if G_train.shape[0] != y_train.shape[0]:
            raise ValueError("y_train length must match G_train")

        n = y_train.shape[0]
        mu = float(np.mean(y_train))
        A = G_train + self.lambda_param * np.eye(n)
        alpha = np.linalg.solve(A, y_train - mu)

        self.mu_ = mu
        self.alpha_ = alpha
        return self

    def predict(self, G_test_train: np.ndarray) -> np.ndarray:
        if self.mu_ is None or self.alpha_ is None:
            raise RuntimeError("GBLUP model not fit")
        return self.mu_ + np.asarray(G_test_train) @ self.alpha_


# =============================================================================
# CV runners
# =============================================================================


def kfold_cv_rrblup(X, y, n_splits=5, random_state=42):
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    rows = []
    for fold, (tr, te) in enumerate(kf.split(X), 1):
        y_pred = rrblup_ridge_fit_predict(X[tr], y[tr], X[te], alpha=1.0)
        m = _regression_metrics(y[te], y_pred)
        rows.append({"fold": fold, "train_size": tr.size, "test_size": te.size, **m})
    return pd.DataFrame(rows)


def stratified_kfold_cv_rrblup(X, y, strata, n_splits=5, random_state=42):
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    rows = []
    for fold, (tr, te) in enumerate(skf.split(X, strata), 1):
        y_pred = rrblup_ridge_fit_predict(X[tr], y[tr], X[te], alpha=1.0)
        m = _regression_metrics(y[te], y_pred)
        rows.append({"fold": fold, "train_size": tr.size, "test_size": te.size, **m})
    return pd.DataFrame(rows)


def group_kfold_cv_rrblup(X, y, groups, n_splits=5):
    """GroupKFold is the cleanest leakage avoidance when groups are known."""
    gkf = GroupKFold(n_splits=n_splits)
    rows = []
    for fold, (tr, te) in enumerate(gkf.split(X, y, groups=groups), 1):
        y_pred = rrblup_ridge_fit_predict(X[tr], y[tr], X[te], alpha=1.0)
        m = _regression_metrics(y[te], y_pred)
        rows.append({"fold": fold, "train_size": tr.size, "test_size": te.size, **m})
    return pd.DataFrame(rows)


def relatedness_aware_cv_rrblup(X, y, kinship, threshold=0.125, n_splits=5, random_state=42):
    cv = RelatednessAwareCV(kinship, threshold=threshold, n_splits=n_splits, random_state=random_state)
    rows = []
    for fold, (tr, te) in enumerate(cv.split(X), 1):
        y_pred = rrblup_ridge_fit_predict(X[tr], y[tr], X[te], alpha=1.0)
        m = _regression_metrics(y[te], y_pred)
        rows.append({"fold": fold, "train_size": tr.size, "test_size": te.size, **m})
    return pd.DataFrame(rows)


def gblup_kfold_cv(y, kinship, n_splits=5, random_state=42, lambda_param=1.0):
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    rows = []
    n = y.shape[0]
    for fold, (tr, te) in enumerate(kf.split(np.arange(n)), 1):
        G_train = kinship[np.ix_(tr, tr)]
        G_test_train = kinship[np.ix_(te, tr)]
        model = GBLUP(lambda_param=lambda_param).fit(y_train=y[tr], G_train=G_train)
        y_pred = model.predict(G_test_train)
        m = _regression_metrics(y[te], y_pred)
        rows.append({"fold": fold, "train_size": tr.size, "test_size": te.size, **m})
    return pd.DataFrame(rows)


def forward_validation_rrblup(X, y, time, min_train_timepoints=2):
    """Forward (chronological) validation.

    Train on all time <= t and test on the next time point.
    """
    time = np.asarray(time)
    uniq = np.sort(np.unique(time))
    if uniq.size < (min_train_timepoints + 1):
        raise ValueError("Not enough time points for forward validation")

    rows = []
    for i in range(min_train_timepoints - 1, uniq.size - 1):
        train_t = uniq[i]
        test_t = uniq[i + 1]
        tr = np.where(time <= train_t)[0]
        te = np.where(time == test_t)[0]
        if tr.size == 0 or te.size == 0:
            continue
        y_pred = rrblup_ridge_fit_predict(X[tr], y[tr], X[te], alpha=1.0)
        m = _regression_metrics(y[te], y_pred)
        rows.append({
            "split": len(rows) + 1,
            "train_time_max": int(train_t),
            "test_time": int(test_t),
            "train_size": tr.size,
            "test_size": te.size,
            **m,
        })

    return pd.DataFrame(rows)


def loo_cv_rrblup(X, y):
    loo = LeaveOneOut()
    y_pred = np.zeros_like(y, dtype=float)
    for tr, te in loo.split(X):
        y_pred[te[0]] = rrblup_ridge_fit_predict(X[tr], y[tr], X[te], alpha=1.0)[0]
    m = _regression_metrics(y, y_pred)
    return pd.DataFrame([{"n_samples": int(y.shape[0]), **m}])


# =============================================================================
# Demo
# =============================================================================


def _print_block(title: str, df: pd.DataFrame):
    print(f"\n{title}")
    print("-" * len(title))
    print(df.to_string(index=False))
    if "correlation" in df.columns and len(df) > 1:
        print(f"\nMean correlation: {df['correlation'].mean():.4f} (±{df['correlation'].std():.4f})")


def run_demo():
    data = simulate_genomic_data(
        n_individuals=300,
        n_markers=800,
        n_qtl=40,
        n_subpopulations=3,
        n_families=60,
        n_years=6,
        heritability=0.6,
        seed=42,
    )

    X = data["genotypes"]
    y = data["phenotypes"]
    subp = data["subpopulations"]
    fam = data["family_ids"]
    time = data["time"]
    G = data["kinship"]

    print("=" * 80)
    print("Genomic prediction CV demo (simulated data)")
    print("=" * 80)
    print(f"Individuals: {X.shape[0]} | Markers: {X.shape[1]}")
    print(f"Subpops: {len(np.unique(subp))} | Families: {len(np.unique(fam))} | Time points: {len(np.unique(time))}")

    _print_block("1) Standard 5-fold CV (RR-BLUP via Ridge)", kfold_cv_rrblup(X, y, n_splits=5, random_state=42))
    _print_block("2) Stratified 5-fold CV by subpopulation", stratified_kfold_cv_rrblup(X, y, subp, n_splits=5, random_state=42))
    _print_block("3) GroupKFold by family_id (leakage avoidance)", group_kfold_cv_rrblup(X, y, fam, n_splits=5))
    _print_block("4) Kinship-threshold CV (connected components)", relatedness_aware_cv_rrblup(X, y, G, threshold=0.125, n_splits=5, random_state=42))
    _print_block("5) GBLUP with GRM (5-fold CV)", gblup_kfold_cv(y, G, n_splits=5, random_state=42, lambda_param=1.0))
    _print_block("6) Forward validation (chronological)", forward_validation_rrblup(X, y, time=time, min_train_timepoints=2))

    # LOOCV is included but not run by default (can be slow)
    print("\nLOOCV is available via: loo_cv_rrblup(X, y)")


if __name__ == "__main__":
    run_demo()
