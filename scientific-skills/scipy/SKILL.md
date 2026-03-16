---
name: scipy
description: Scientific computing with SciPy. Use for optimization (minimize, curve fitting, linear programming), statistics (hypothesis testing, distributions, descriptive stats), signal processing (filtering, FFT, spectrograms), linear algebra (decompositions, eigenvalues), interpolation, numerical integration, ODE solvers, sparse matrices, and spatial data structures.
license: BSD-3-Clause license
metadata:
    skill-author: Adem Usta
---

# SciPy: Scientific Computing

## Overview

SciPy is the foundational Python library for scientific and technical computing. Built on NumPy, it provides efficient, well-tested implementations of algorithms for optimization, integration, interpolation, eigenvalue problems, algebraic equations, differential equations, statistics, signal processing, and more. It underpins nearly every scientific Python workflow.

## Installation

```bash
# Install SciPy using uv
uv pip install scipy

# Commonly used with
uv pip install numpy matplotlib pandas
```

## When to Use This Skill

Use the SciPy skill when:

- Performing numerical optimization (minimization, curve fitting, root finding)
- Running statistical tests (t-tests, chi-squared, ANOVA, KS tests)
- Processing signals (filtering, FFT, spectrograms, peak detection)
- Solving linear algebra problems (eigenvalues, decompositions, linear systems)
- Interpolating data (1D/2D splines, griddata)
- Numerically integrating functions or solving ODEs
- Working with sparse matrices for large-scale problems
- Computing spatial distances, KD-trees, or convex hulls

Consider alternatives when:
- You need GPU-accelerated computation → Use **PyTorch** or **CuPy**
- You need symbolic mathematics → Use **sympy**
- You need deep learning → Use **pytorch-lightning** or **transformers**
- You need statistical modeling with formulas → Use **statsmodels**

## Quick Start

### Optimization: Curve Fitting

```python
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Define model function
def exponential_decay(x, a, b, c):
    return a * np.exp(-b * x) + c

# Generate noisy data
x_data = np.linspace(0, 10, 100)
y_data = exponential_decay(x_data, 2.5, 0.5, 0.3) + 0.1 * np.random.randn(100)

# Fit the model
popt, pcov = curve_fit(exponential_decay, x_data, y_data, p0=[1, 1, 0])
print(f"Fitted parameters: a={popt[0]:.3f}, b={popt[1]:.3f}, c={popt[2]:.3f}")

# Parameter uncertainties (1-sigma)
perr = np.sqrt(np.diag(pcov))
print(f"Uncertainties: ±{perr[0]:.3f}, ±{perr[1]:.3f}, ±{perr[2]:.3f}")

# Plot
plt.scatter(x_data, y_data, s=5, label='Data')
plt.plot(x_data, exponential_decay(x_data, *popt), 'r-', label='Fit')
plt.legend()
plt.savefig('curve_fit.png', dpi=150)
```

### Statistical Hypothesis Testing

```python
from scipy import stats

# Two-sample t-test
group_a = [23.1, 25.4, 22.8, 24.9, 26.1, 23.7, 25.0]
group_b = [28.3, 27.1, 29.5, 26.8, 30.2, 27.9, 28.7]

t_stat, p_value = stats.ttest_ind(group_a, group_b)
print(f"t-statistic: {t_stat:.4f}, p-value: {p_value:.6f}")

# Normality test (Shapiro-Wilk)
stat, p = stats.shapiro(group_a)
print(f"Shapiro-Wilk: stat={stat:.4f}, p={p:.4f}")

# Non-parametric alternative (Mann-Whitney U)
u_stat, p_value = stats.mannwhitneyu(group_a, group_b, alternative='two-sided')
print(f"Mann-Whitney U: stat={u_stat:.4f}, p={p_value:.6f}")
```

## Core Capabilities

### 1. Optimization (`scipy.optimize`)

Functions for finding minima, roots, and fitting models to data.

```python
from scipy.optimize import minimize, minimize_scalar, root, linprog

# Unconstrained minimization
def rosenbrock(x):
    return (1 - x[0])**2 + 100 * (x[1] - x[0]**2)**2

result = minimize(rosenbrock, x0=[0, 0], method='L-BFGS-B')
print(f"Minimum at: {result.x}, value: {result.fun:.6f}")

# Constrained optimization
constraints = [
    {'type': 'ineq', 'fun': lambda x: x[0] + x[1] - 1},  # x+y >= 1
]
bounds = [(0, None), (0, None)]  # x, y >= 0
result = minimize(rosenbrock, x0=[0.5, 0.5], method='SLSQP',
                  bounds=bounds, constraints=constraints)

# Root finding
from scipy.optimize import brentq
root_val = brentq(lambda x: x**3 - 2*x - 5, 1, 3)
print(f"Root: {root_val:.6f}")

# Linear programming
# Minimize c^T x subject to A_ub @ x <= b_ub
c = [-1, -2]  # Coefficients to minimize
A_ub = [[1, 1], [2, 1], [1, 0]]
b_ub = [10, 16, 8]
result = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=[(0, None), (0, None)])
print(f"Optimal: x={result.x}, value={result.fun:.2f}")
```

**Key functions:**
- `minimize()`: General-purpose minimization (Nelder-Mead, BFGS, L-BFGS-B, SLSQP, trust-constr)
- `curve_fit()`: Non-linear least squares curve fitting
- `least_squares()`: Robust non-linear least squares with bounds
- `root()` / `brentq()`: Root finding for equations
- `linprog()`: Linear programming
- `differential_evolution()`: Global optimization

**See:** `references/api_reference.md` for full parameter details.

### 2. Statistics (`scipy.stats`)

Comprehensive statistical functions: distributions, tests, and descriptive statistics.

```python
from scipy import stats
import numpy as np

# Fit a distribution to data
data = np.random.gamma(shape=2.0, scale=1.5, size=1000)
shape, loc, scale = stats.gamma.fit(data)
print(f"Fitted gamma: shape={shape:.3f}, loc={loc:.3f}, scale={scale:.3f}")

# Kolmogorov-Smirnov test (goodness of fit)
ks_stat, ks_p = stats.kstest(data, 'gamma', args=(shape, loc, scale))
print(f"KS test: stat={ks_stat:.4f}, p={ks_p:.4f}")

# One-way ANOVA
f_stat, p_val = stats.f_oneway(group1, group2, group3)

# Chi-squared test of independence
contingency_table = np.array([[10, 20, 30], [6, 14, 25]])
chi2, p, dof, expected = stats.chi2_contingency(contingency_table)

# Pearson and Spearman correlations
r, p = stats.pearsonr(x, y)
rho, p = stats.spearmanr(x, y)

# Descriptive statistics
result = stats.describe(data)
print(f"Mean={result.mean:.3f}, Variance={result.variance:.3f}")
print(f"Skewness={result.skewness:.3f}, Kurtosis={result.kurtosis:.3f}")
```

**80+ continuous distributions** (norm, gamma, beta, lognorm, weibull, etc.) all with:
- `pdf()`, `cdf()`, `ppf()`, `sf()` — density, CDF, quantile, survival
- `rvs()` — random sampling
- `fit()` — MLE parameter estimation

### 3. Signal Processing (`scipy.signal`)

Filter design, spectral analysis, and signal manipulation.

```python
from scipy import signal
import numpy as np

# Design a Butterworth bandpass filter
fs = 1000  # Sampling frequency (Hz)
lowcut, highcut = 20, 100  # Band: 20-100 Hz
b, a = signal.butter(4, [lowcut, highcut], btype='band', fs=fs)

# Apply the filter
t = np.linspace(0, 1, fs, endpoint=False)
noisy_signal = np.sin(2*np.pi*50*t) + 0.5*np.random.randn(fs)
filtered = signal.filtfilt(b, a, noisy_signal)

# Spectrogram
f, t_spec, Sxx = signal.spectrogram(noisy_signal, fs=fs,
                                      nperseg=128, noverlap=64)

# Peak detection
peaks, properties = signal.find_peaks(filtered, height=0.5,
                                       distance=10, prominence=0.3)
print(f"Found {len(peaks)} peaks at indices: {peaks}")

# FFT-based power spectral density
f_psd, Pxx = signal.welch(noisy_signal, fs=fs, nperseg=256)
```

### 4. Linear Algebra (`scipy.linalg`)

Matrix decompositions, solvers, and specialized operations.

```python
from scipy import linalg
import numpy as np

A = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 10]])
b = np.array([1, 2, 3])

# Solve linear system Ax = b
x = linalg.solve(A, b)

# LU decomposition
P, L, U = linalg.lu(A)

# SVD (Singular Value Decomposition)
U, s, Vt = linalg.svd(A)
print(f"Singular values: {s}")

# Eigenvalues and eigenvectors
eigenvalues, eigenvectors = linalg.eig(A)

# Matrix exponential (useful for ODEs)
expm_A = linalg.expm(A)

# Cholesky decomposition (positive definite matrices)
M = np.array([[4, 2], [2, 3]])
L = linalg.cholesky(M, lower=True)

# Condition number
cond = np.linalg.cond(A)
print(f"Condition number: {cond:.2f}")
```

### 5. Interpolation (`scipy.interpolate`)

Fit smooth functions through data points.

```python
from scipy.interpolate import (
    interp1d, CubicSpline, RBFInterpolator, griddata
)
import numpy as np

# 1D interpolation
x = np.array([0, 1, 2, 3, 4, 5])
y = np.array([0, 0.8, 0.9, 0.1, -0.8, -1.0])

# Cubic spline (recommended)
cs = CubicSpline(x, y)
x_fine = np.linspace(0, 5, 100)
y_fine = cs(x_fine)
y_deriv = cs(x_fine, 1)  # First derivative

# 2D scattered data interpolation
points = np.random.rand(100, 2)
values = np.sin(points[:, 0] * np.pi) * np.cos(points[:, 1] * np.pi)
grid_x, grid_y = np.mgrid[0:1:50j, 0:1:50j]
grid_z = griddata(points, values, (grid_x, grid_y), method='cubic')

# RBF interpolation (N-dimensional)
rbf = RBFInterpolator(points, values, kernel='thin_plate_spline')
query_points = np.random.rand(10, 2)
result = rbf(query_points)
```

### 6. Integration (`scipy.integrate`)

Numerical integration and ODE solvers.

```python
from scipy.integrate import quad, dblquad, solve_ivp
import numpy as np

# Definite integral
result, error = quad(lambda x: np.exp(-x**2), 0, np.inf)
print(f"∫exp(-x²)dx from 0 to ∞ = {result:.6f} (error: {error:.2e})")

# Double integral
result, error = dblquad(lambda y, x: np.exp(-(x**2 + y**2)),
                         0, 1, 0, 1)

# Solve ODE: dy/dt = -ky (exponential decay)
def decay(t, y, k=0.5):
    return -k * y

sol = solve_ivp(decay, t_span=[0, 10], y0=[1.0],
                t_eval=np.linspace(0, 10, 100), method='RK45')
# sol.t = time points, sol.y = solution values

# Lotka-Volterra predator-prey system
def lotka_volterra(t, z, a=1.5, b=1.0, c=3.0, d=1.0):
    x, y = z
    return [a*x - b*x*y, -c*y + d*x*y]

sol = solve_ivp(lotka_volterra, [0, 20], [10, 5],
                t_eval=np.linspace(0, 20, 500), method='RK45')
```

### 7. Sparse Matrices (`scipy.sparse`)

Efficient storage and computation for large sparse matrices.

```python
from scipy import sparse
from scipy.sparse.linalg import spsolve, eigsh
import numpy as np

# Create sparse matrices
row = [0, 0, 1, 2, 2]
col = [0, 2, 1, 0, 2]
data = [1, 3, 2, 4, 5]
A = sparse.csr_matrix((data, (row, col)), shape=(3, 3))

# Convert formats
A_csc = A.tocsc()  # Column-compressed
A_dense = A.toarray()

# Sparse linear solve
b = np.array([1, 2, 3])
x = spsolve(A.tocsc(), b)

# Sparse eigenvalues (largest k)
eigenvalues, eigenvectors = eigsh(A.toarray(), k=2)

# Build sparse identity + diagonals
I = sparse.eye(100, format='csr')
D = sparse.diags([1, -2, 1], [-1, 0, 1], shape=(100, 100), format='csr')
```

### 8. Spatial Data Structures (`scipy.spatial`)

Distance computation, KD-trees, convex hulls, and Voronoi diagrams.

```python
from scipy.spatial import (
    KDTree, ConvexHull, Voronoi, distance, Delaunay
)
import numpy as np

# KD-Tree for fast nearest neighbor queries
points = np.random.rand(10000, 3)
tree = KDTree(points)

# Query nearest neighbors
query = np.array([0.5, 0.5, 0.5])
dist, idx = tree.query(query, k=5)  # 5 nearest neighbors
print(f"Nearest indices: {idx}, distances: {dist}")

# Pairwise distances
from scipy.spatial.distance import pdist, squareform, cdist
dists = pdist(points[:100], metric='euclidean')
dist_matrix = squareform(dists)

# Cross-distance between two sets
dists_ab = cdist(points[:50], points[50:100], metric='cosine')

# Convex hull
hull = ConvexHull(points[:, :2])  # 2D projection
print(f"Hull area: {hull.area:.4f}, volume: {hull.volume:.4f}")
```

## Common Scientific Workflows

### Fitting Experimental Data

```python
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2

def model(x, a, b, c):
    return a * np.exp(-b * x) + c

# Fit with uncertainties
popt, pcov = curve_fit(model, x_data, y_data, sigma=y_errors,
                        absolute_sigma=True)

# Reduced chi-squared
residuals = y_data - model(x_data, *popt)
chi2_red = np.sum((residuals / y_errors)**2) / (len(y_data) - len(popt))
print(f"Reduced χ² = {chi2_red:.3f}")

# Confidence intervals
from scipy.stats import t as t_dist
dof = len(y_data) - len(popt)
t_val = t_dist.ppf(0.975, dof)
ci = t_val * np.sqrt(np.diag(pcov))
for i, (p, c) in enumerate(zip(popt, ci)):
    print(f"param[{i}] = {p:.4f} ± {c:.4f} (95% CI)")
```

### Power Spectral Analysis

```python
from scipy import signal
import numpy as np

# Generate multi-frequency signal
fs = 2000  # Hz
t = np.arange(0, 5, 1/fs)
sig = (np.sin(2*np.pi*50*t) + 0.5*np.sin(2*np.pi*120*t)
       + 0.2*np.random.randn(len(t)))

# Welch's PSD estimate
f, Pxx = signal.welch(sig, fs=fs, nperseg=1024)

# Detect dominant frequencies
peaks, props = signal.find_peaks(10*np.log10(Pxx), height=-20, prominence=5)
dominant_freqs = f[peaks]
print(f"Dominant frequencies: {dominant_freqs} Hz")
```

### Bootstrap Confidence Intervals

```python
from scipy.stats import bootstrap
import numpy as np

data = np.random.exponential(scale=2.0, size=100)
result = bootstrap(
    (data,),
    statistic=np.mean,
    n_resamples=10000,
    confidence_level=0.95,
    method='percentile'
)
print(f"95% CI for mean: [{result.confidence_interval.low:.3f}, "
      f"{result.confidence_interval.high:.3f}]")
```

## Key Parameters to Adjust

### Optimization
- `method`: Algorithm choice ('Nelder-Mead', 'BFGS', 'L-BFGS-B', 'SLSQP', 'trust-constr')
- `tol`: Convergence tolerance
- `maxiter`: Maximum iterations
- `bounds`: Parameter bounds for bounded optimization

### Signal Processing
- `fs`: Sampling frequency (critical — always set explicitly)
- `nperseg`: Segment length for spectral analysis (trade-off: frequency vs time resolution)
- `order`: Filter order (higher = sharper cutoff but more ringing)

### Integration
- `method`: ODE solver ('RK45', 'RK23', 'BDF', 'Radau', 'LSODA')
- `rtol`, `atol`: Relative and absolute error tolerances

## Common Pitfalls and Best Practices

1. **Always check convergence**: Verify `result.success` after optimization
2. **Provide good initial guesses**: `x0` for `minimize()`, `p0` for `curve_fit()`
3. **Use `filtfilt` not `lfilter`**: `filtfilt` applies zero-phase filtering (no time delay)
4. **Set `absolute_sigma=True`** in `curve_fit()` when you have known error bars
5. **Choose the right sparse format**: CSR for row operations, CSC for column operations
6. **Use `solve_ivp` over `odeint`**: Modern API with better error control
7. **Check matrix condition number** before solving linear systems
8. **Use appropriate statistical tests**: Check normality before using parametric tests

## Additional Resources

- **Official Documentation**: https://docs.scipy.org/doc/scipy/
- **Tutorial**: https://docs.scipy.org/doc/scipy/tutorial/
- **API Reference**: https://docs.scipy.org/doc/scipy/reference/
- **Cookbook**: https://scipy-cookbook.readthedocs.io/

## Reference Documentation

- **API Reference**: See [references/api_reference.md](references/api_reference.md) for module-by-module function listing
