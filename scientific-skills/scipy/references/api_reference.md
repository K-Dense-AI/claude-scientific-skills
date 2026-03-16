# SciPy API Quick Reference

## scipy.optimize — Optimization & Root Finding

| Function | Purpose | Key Parameters |
|----------|---------|----------------|
| `minimize(fun, x0)` | General minimization | `method`, `bounds`, `constraints`, `tol` |
| `minimize_scalar(fun)` | 1D minimization | `method`, `bounds` |
| `curve_fit(f, xdata, ydata)` | Non-linear least squares | `p0`, `sigma`, `bounds`, `absolute_sigma` |
| `least_squares(fun, x0)` | Robust least squares | `bounds`, `loss`, `method` |
| `root(fun, x0)` | Find roots of vector function | `method`, `tol` |
| `brentq(f, a, b)` | Bracketed scalar root | `xtol`, `rtol` |
| `linprog(c)` | Linear programming | `A_ub`, `b_ub`, `A_eq`, `b_eq`, `bounds` |
| `differential_evolution(func, bounds)` | Global optimization | `strategy`, `maxiter`, `tol`, `seed` |
| `dual_annealing(func, bounds)` | Simulated annealing | `maxiter`, `seed` |
| `shgo(func, bounds)` | Simplicial homology global | `n`, `iters`, `sampling_method` |

## scipy.stats — Statistics

### Hypothesis Tests

| Function | Purpose | Returns |
|----------|---------|---------|
| `ttest_ind(a, b)` | Independent two-sample t-test | `statistic`, `pvalue` |
| `ttest_rel(a, b)` | Paired t-test | `statistic`, `pvalue` |
| `ttest_1samp(a, popmean)` | One-sample t-test | `statistic`, `pvalue` |
| `f_oneway(*groups)` | One-way ANOVA | `statistic`, `pvalue` |
| `kruskal(*groups)` | Kruskal-Wallis H-test | `statistic`, `pvalue` |
| `mannwhitneyu(x, y)` | Mann-Whitney U test | `statistic`, `pvalue` |
| `wilcoxon(x, y)` | Wilcoxon signed-rank | `statistic`, `pvalue` |
| `chi2_contingency(table)` | Chi-squared independence | `chi2`, `p`, `dof`, `expected` |
| `shapiro(x)` | Shapiro-Wilk normality | `statistic`, `pvalue` |
| `normaltest(x)` | D'Agostino-Pearson normality | `statistic`, `pvalue` |
| `kstest(x, cdf)` | Kolmogorov-Smirnov | `statistic`, `pvalue` |
| `pearsonr(x, y)` | Pearson correlation | `r`, `pvalue` |
| `spearmanr(x, y)` | Spearman correlation | `rho`, `pvalue` |
| `kendalltau(x, y)` | Kendall's tau | `tau`, `pvalue` |

### Distribution Methods (all 80+ distributions)

| Method | Purpose |
|--------|---------|
| `.rvs(size=n)` | Generate random samples |
| `.pdf(x)` / `.pmf(x)` | Probability density / mass |
| `.cdf(x)` | Cumulative distribution function |
| `.ppf(q)` | Percent point / quantile function |
| `.sf(x)` | Survival function (1 - CDF) |
| `.fit(data)` | MLE parameter estimation |
| `.interval(confidence)` | Confidence interval |

### Common Distributions

`norm`, `t`, `chi2`, `f`, `gamma`, `beta`, `expon`, `lognorm`, `weibull_min`, `uniform`, `poisson`, `binom`, `nbinom`, `bernoulli`, `multinomial`

## scipy.signal — Signal Processing

| Function | Purpose | Key Parameters |
|----------|---------|----------------|
| `butter(N, Wn)` | Butterworth filter design | `btype`, `fs`, `output` |
| `cheby1(N, rp, Wn)` | Chebyshev Type I filter | `btype`, `fs` |
| `filtfilt(b, a, x)` | Zero-phase filtering | `padlen`, `axis` |
| `sosfiltfilt(sos, x)` | Zero-phase SOS filtering | `padlen` |
| `welch(x, fs)` | Power spectral density | `nperseg`, `noverlap`, `window` |
| `spectrogram(x, fs)` | Time-frequency spectrogram | `nperseg`, `noverlap`, `mode` |
| `find_peaks(x)` | Detect peaks | `height`, `distance`, `prominence`, `width` |
| `savgol_filter(x, N, p)` | Savitzky-Golay smoothing | `window_length`, `polyorder` |
| `resample(x, num)` | Resample signal | `t`, `axis` |
| `hilbert(x)` | Analytic signal | `N` |
| `cwt(data, wavelet, widths)` | Continuous wavelet transform | `wavelet`, `widths` |

## scipy.linalg — Linear Algebra

| Function | Purpose | Key Parameters |
|----------|---------|----------------|
| `solve(A, b)` | Solve Ax = b | `assume_a`, `overwrite_a` |
| `inv(A)` | Matrix inverse | `overwrite_a` |
| `det(A)` | Determinant | — |
| `eig(A)` | Eigenvalues + vectors | `left`, `right` |
| `eigh(A)` | Symmetric eigendecomposition | `eigvals_only`, `subset_by_value` |
| `svd(A)` | Singular value decomposition | `full_matrices`, `compute_uv` |
| `lu(A)` | LU decomposition | `permute_l` |
| `cholesky(A)` | Cholesky decomposition | `lower` |
| `qr(A)` | QR decomposition | `mode` |
| `expm(A)` | Matrix exponential | — |
| `norm(A)` | Matrix/vector norm | `ord` |

## scipy.interpolate — Interpolation

| Function | Purpose | Key Parameters |
|----------|---------|----------------|
| `CubicSpline(x, y)` | Cubic spline (recommended) | `bc_type`, `extrapolate` |
| `interp1d(x, y)` | 1D interpolation | `kind`, `fill_value`, `bounds_error` |
| `UnivariateSpline(x, y)` | Smoothing spline | `s`, `k` |
| `RBFInterpolator(y, d)` | RBF N-D interpolation | `kernel`, `epsilon` |
| `griddata(points, values, xi)` | Scattered data interpolation | `method` ('linear', 'cubic', 'nearest') |
| `RegularGridInterpolator(pts, vals)` | Regular grid interpolation | `method`, `bounds_error` |

## scipy.integrate — Integration & ODEs

| Function | Purpose | Key Parameters |
|----------|---------|----------------|
| `quad(func, a, b)` | Definite integral | `args`, `limit`, `epsabs` |
| `dblquad(func, a, b, gfun, hfun)` | Double integral | `epsabs`, `epsrel` |
| `nquad(func, ranges)` | N-dimensional integral | `opts` |
| `solve_ivp(fun, t_span, y0)` | Initial value ODE | `method`, `t_eval`, `rtol`, `atol`, `events` |
| `trapezoid(y, x)` | Trapezoidal rule | `dx`, `axis` |
| `simpson(y, x)` | Simpson's rule | `dx`, `axis` |

### ODE Solver Methods for `solve_ivp`

| Method | Type | Best For |
|--------|------|----------|
| `RK45` | Explicit | General non-stiff problems (default) |
| `RK23` | Explicit | Low-accuracy, non-stiff |
| `DOP853` | Explicit | High-accuracy, non-stiff |
| `Radau` | Implicit | Stiff problems |
| `BDF` | Implicit | Stiff problems |
| `LSODA` | Auto | Automatically switches stiff/non-stiff |

## scipy.sparse — Sparse Matrices

| Function | Purpose |
|----------|---------|
| `csr_matrix(data)` | Compressed Sparse Row |
| `csc_matrix(data)` | Compressed Sparse Column |
| `coo_matrix(data)` | Coordinate format (fast construction) |
| `eye(n)` | Sparse identity matrix |
| `diags(diagonals, offsets)` | Sparse diagonal matrix |
| `block_diag(*mats)` | Sparse block diagonal |
| `sparse.linalg.spsolve(A, b)` | Sparse linear solve |
| `sparse.linalg.eigsh(A, k)` | Sparse symmetric eigenvalues |
| `sparse.linalg.svds(A, k)` | Sparse SVD |
| `sparse.linalg.gmres(A, b)` | Iterative solver (GMRES) |
| `sparse.linalg.cg(A, b)` | Conjugate gradient solver |

## scipy.spatial — Spatial Data Structures

| Function | Purpose | Key Parameters |
|----------|---------|----------------|
| `KDTree(data)` | K-D tree for fast queries | `leafsize` |
| `cKDTree(data)` | Cython K-D tree (faster) | `leafsize` |
| `ConvexHull(points)` | Convex hull | — |
| `Delaunay(points)` | Delaunay triangulation | — |
| `Voronoi(points)` | Voronoi diagram | — |
| `distance.pdist(X)` | Pairwise distances | `metric` |
| `distance.cdist(XA, XB)` | Cross distances | `metric` |
| `distance.squareform(X)` | Condensed ↔ square matrix | — |
