# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# Licensed under the Apache License, Version 2.0.

# genomic_cv_r.R
#
# Cross-validation methods for genomic prediction (R).
#
# Includes:
# - Standard K-fold CV
# - Stratified K-fold CV (by subpopulation)
# - GroupKFold by known family IDs (recommended)
# - Relatedness-aware CV from kinship threshold groups
# - GBLUP CV using the genomic relationship matrix
# - Forward (chronological) validation
#
# Install deps:
#   install.packages(c("rrBLUP", "caret"))
#
# Run:
#   source("genomic_cv_r.R")

suppressPackageStartupMessages({
  if (!requireNamespace("rrBLUP", quietly = TRUE)) {
    stop("Package 'rrBLUP' is required. Install via install.packages('rrBLUP')")
  }
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package 'caret' is required. Install via install.packages('caret')")
  }
})

library(rrBLUP)
library(caret)


# =============================================================================
# Utilities
# =============================================================================

regression_metrics <- function(y_true, y_pred) {
  y_true <- as.numeric(y_true)
  y_pred <- as.numeric(y_pred)
  corr <- suppressWarnings(cor(y_true, y_pred))
  rmse <- sqrt(mean((y_true - y_pred)^2))
  r2 <- 1 - sum((y_true - y_pred)^2) / sum((y_true - mean(y_true))^2)
  list(correlation = corr, rmse = rmse, r2 = r2)
}


# =============================================================================
# Simulated genotype/phenotype data
# =============================================================================

simulate_genomic_data <- function(
  n_individuals = 400,
  n_markers = 2000,
  n_qtl = 50,
  n_subpopulations = 3,
  n_families = 60,
  family_freq_sd = 0.03,
  n_years = 6,
  year_effect_sd = 0.5,
  heritability = 0.5,
  seed = 42
) {
  set.seed(seed)

  if (n_families < 1 || n_families > n_individuals) {
    stop("n_families must be in [1, n_individuals]")
  }
  if (heritability <= 0 || heritability >= 1) {
    stop("heritability must be in (0, 1)")
  }
  if (n_years < 2) {
    stop("n_years must be >= 2")
  }

  # Assign families to subpops (roughly balanced)
  family_subpops <- rep(0:(n_subpopulations - 1), length.out = n_families)
  family_subpops <- sample(family_subpops)

  # Family sizes that sum to n_individuals
  base <- floor(n_individuals / n_families)
  family_sizes <- rep(base, n_families)
  family_sizes[1:(n_individuals - base * n_families)] <- family_sizes[1:(n_individuals - base * n_families)] + 1
  family_sizes <- sample(family_sizes)

  family_ids <- unlist(lapply(1:n_families, function(f) rep(f, family_sizes[f])))
  family_ids <- sample(family_ids)
  subpopulations <- family_subpops[family_ids]

  # Base allele freqs
  base_freqs <- runif(n_markers, 0.1, 0.9)
  subpop_freqs <- vector("list", n_subpopulations)
  for (sp in 1:n_subpopulations) {
    deviation <- rnorm(n_markers, 0, 0.10)
    subpop_freqs[[sp]] <- pmin(0.95, pmax(0.05, base_freqs + deviation))
  }

  # Genotypes M (0/1/2)
  M <- matrix(0, nrow = n_individuals, ncol = n_markers)
  for (f in 1:n_families) {
    sp <- family_subpops[f] + 1
    freqs <- subpop_freqs[[sp]]
    family_freqs <- pmin(0.95, pmax(0.05, freqs + rnorm(n_markers, 0, family_freq_sd)))
    idx <- which(family_ids == f)
    M[idx, ] <- matrix(rbinom(length(idx) * n_markers, 2, rep(family_freqs, each = length(idx))), nrow = length(idx), byrow = TRUE)
  }

  # Standardize X for ridge on markers
  M_centered <- scale(M, center = TRUE, scale = FALSE)
  sdv <- apply(M_centered, 2, sd)
  sdv[sdv == 0] <- 1
  X <- sweep(M_centered, 2, sdv, "/")

  # QTL effects
  qtl_idx <- sample(1:n_markers, n_qtl)
  qtl_effects <- rnorm(n_qtl, 0, 1)
  g <- as.numeric(X[, qtl_idx] %*% qtl_effects)

  var_g <- var(g)
  var_e <- var_g * (1 - heritability) / heritability

  time <- rep(0:(n_years - 1), length.out = n_individuals)
  time <- sample(time)
  year_effect <- rnorm(n_years, 0, year_effect_sd)

  e <- rnorm(n_individuals, 0, sqrt(var_e))
  y <- g + year_effect[time + 1] + e

  # VanRaden GRM
  p <- colMeans(M) / 2
  denom <- 2 * sum(p * (1 - p))
  W <- sweep(M, 2, 2 * p, "-") / sqrt(denom)
  G <- (W %*% t(W)) / n_markers

  list(
    genotypes = X,
    phenotypes = y,
    subpopulations = subpopulations,
    family_ids = family_ids,
    time = time,
    kinship = G,
    qtl_indices = qtl_idx
  )
}


# =============================================================================
# RR-BLUP / GBLUP helpers
# =============================================================================

rrblup_fit_predict <- function(X_train, y_train, X_test) {
  fit <- mixed.solve(y = y_train, Z = X_train, method = "REML")
  as.numeric(X_test %*% fit$u + fit$beta)
}

# GBLUP via kernel ridge on GRM with lambda = Ve/Vu estimated from training

gblup_fit_predict <- function(y_train, K_train, K_test_train) {
  fit <- mixed.solve(y = y_train, K = K_train, method = "REML")
  mu <- as.numeric(fit$beta)
  lambda <- as.numeric(fit$Ve / fit$Vu)
  n <- length(y_train)
  alpha <- solve(K_train + diag(lambda, n), y_train - mu)
  as.numeric(mu + K_test_train %*% alpha)
}


# =============================================================================
# CV methods
# =============================================================================

standard_kfold_cv_rrblup <- function(X, y, n_folds = 5, seed = 42) {
  set.seed(seed)
  folds <- sample(rep(1:n_folds, length.out = length(y)))
  out <- data.frame()
  for (f in 1:n_folds) {
    te <- which(folds == f)
    tr <- which(folds != f)
    pred <- rrblup_fit_predict(X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE])
    m <- regression_metrics(y[te], pred)
    out <- rbind(out, data.frame(fold = f, train_size = length(tr), test_size = length(te),
                                 correlation = m$correlation, rmse = m$rmse, r2 = m$r2))
  }
  out
}

stratified_kfold_cv_rrblup <- function(X, y, strata, n_folds = 5, seed = 42) {
  set.seed(seed)
  strata <- as.factor(strata)
  folds <- rep(NA_integer_, length(y))
  for (s in levels(strata)) {
    idx <- which(strata == s)
    folds[idx] <- sample(rep(1:n_folds, length.out = length(idx)))
  }
  out <- data.frame()
  for (f in 1:n_folds) {
    te <- which(folds == f)
    tr <- which(folds != f)
    pred <- rrblup_fit_predict(X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE])
    m <- regression_metrics(y[te], pred)
    out <- rbind(out, data.frame(fold = f, train_size = length(tr), test_size = length(te),
                                 correlation = m$correlation, rmse = m$rmse, r2 = m$r2))
  }
  out
}

family_group_kfold_cv_rrblup <- function(X, y, family_ids, n_folds = 5, seed = 42) {
  set.seed(seed)
  idx_list <- groupKFold(family_ids, k = n_folds)
  out <- data.frame()
  for (f in 1:n_folds) {
    tr <- idx_list[[f]]
    te <- setdiff(seq_along(y), tr)
    pred <- rrblup_fit_predict(X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE])
    m <- regression_metrics(y[te], pred)
    out <- rbind(out, data.frame(fold = f, train_size = length(tr), test_size = length(te),
                                 correlation = m$correlation, rmse = m$rmse, r2 = m$r2))
  }
  out
}

find_related_groups <- function(K, threshold = 0.125) {
  n <- nrow(K)
  adj <- K > threshold
  diag(adj) <- FALSE

  visited <- rep(FALSE, n)
  groups <- list()

  for (i in 1:n) {
    if (visited[i]) next
    q <- c(i)
    visited[i] <- TRUE
    comp <- c()

    while (length(q) > 0) {
      cur <- q[1]
      q <- q[-1]
      comp <- c(comp, cur)
      neigh <- which(adj[cur, ])
      for (j in neigh) {
        if (!visited[j]) {
          visited[j] <- TRUE
          q <- c(q, j)
        }
      }
    }

    groups[[length(groups) + 1]] <- comp
  }

  groups
}

relatedness_aware_cv_rrblup <- function(X, y, K, threshold = 0.125, n_folds = 5, seed = 42) {
  set.seed(seed)
  groups <- find_related_groups(K, threshold)
  group_order <- sample(seq_along(groups))
  folds <- split(group_order, cut(seq_along(group_order), breaks = n_folds, labels = FALSE))

  out <- data.frame()
  for (f in 1:n_folds) {
    te <- unlist(groups[folds[[f]]])
    tr <- setdiff(seq_along(y), te)
    pred <- rrblup_fit_predict(X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE])
    m <- regression_metrics(y[te], pred)
    out <- rbind(out, data.frame(fold = f, train_size = length(tr), test_size = length(te),
                                 correlation = m$correlation, rmse = m$rmse, r2 = m$r2))
  }
  out
}

gblup_kfold_cv <- function(y, K, n_folds = 5, seed = 42) {
  set.seed(seed)
  folds <- sample(rep(1:n_folds, length.out = length(y)))
  out <- data.frame()
  for (f in 1:n_folds) {
    te <- which(folds == f)
    tr <- which(folds != f)

    K_train <- K[tr, tr, drop = FALSE]
    K_test_train <- K[te, tr, drop = FALSE]

    pred <- gblup_fit_predict(y[tr], K_train, K_test_train)
    m <- regression_metrics(y[te], pred)

    out <- rbind(out, data.frame(fold = f, train_size = length(tr), test_size = length(te),
                                 correlation = m$correlation, rmse = m$rmse, r2 = m$r2))
  }
  out
}

forward_validation_rrblup <- function(X, y, time, min_train_timepoints = 2) {
  uniq <- sort(unique(time))
  if (length(uniq) < (min_train_timepoints + 1)) stop("Not enough time points")

  out <- data.frame()
  split_id <- 0
  for (i in seq(min_train_timepoints - 1, length(uniq) - 2)) {
    train_t <- uniq[i + 1]
    test_t <- uniq[i + 2]
    tr <- which(time <= train_t)
    te <- which(time == test_t)
    if (length(tr) == 0 || length(te) == 0) next

    pred <- rrblup_fit_predict(X[tr, , drop = FALSE], y[tr], X[te, , drop = FALSE])
    m <- regression_metrics(y[te], pred)

    split_id <- split_id + 1
    out <- rbind(out, data.frame(split = split_id, train_time_max = train_t, test_time = test_t,
                                 train_size = length(tr), test_size = length(te),
                                 correlation = m$correlation, rmse = m$rmse, r2 = m$r2))
  }

  out
}


# =============================================================================
# Demo
# =============================================================================

run_demo <- function() {
  data <- simulate_genomic_data(
    n_individuals = 300,
    n_markers = 800,
    n_qtl = 40,
    n_subpopulations = 3,
    n_families = 60,
    n_years = 6,
    heritability = 0.6,
    seed = 42
  )

  X <- data$genotypes
  y <- data$phenotypes
  subp <- data$subpopulations
  fam <- data$family_ids
  time <- data$time
  K <- data$kinship

  cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat("Genomic prediction CV demo (simulated data)\n")
  cat(paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat(sprintf("Individuals: %d | Markers: %d\n", nrow(X), ncol(X)))
  cat(sprintf("Subpops: %d | Families: %d | Time points: %d\n", length(unique(subp)), length(unique(fam)), length(unique(time))))

  print(standard_kfold_cv_rrblup(X, y, n_folds = 5, seed = 42))
  print(stratified_kfold_cv_rrblup(X, y, subp, n_folds = 5, seed = 42))
  print(family_group_kfold_cv_rrblup(X, y, fam, n_folds = 5, seed = 42))
  print(relatedness_aware_cv_rrblup(X, y, K, threshold = 0.125, n_folds = 5, seed = 42))
  print(gblup_kfold_cv(y, K, n_folds = 5, seed = 42))
  print(forward_validation_rrblup(X, y, time, min_train_timepoints = 2))

  invisible(data)
}

if (interactive() || sys.nframe() == 0) {
  run_demo()
}
