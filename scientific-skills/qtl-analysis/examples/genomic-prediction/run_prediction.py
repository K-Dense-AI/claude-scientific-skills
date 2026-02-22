#!/usr/bin/env python3
"""
Example: Genomic Prediction (GBLUP)

This example demonstrates genomic prediction using GBLUP (Genomic Best Linear Unbiased Prediction).
Predicts breeding values for individuals using genome-wide markers.

Equivalent to QTLmax: "Genomic prediction"
https://open.qtlmax.com/guide/index.php/2025/07/12/genomic-prediction/

Auto-installs: pandas, numpy, scikit-learn, matplotlib
"""

import subprocess
import sys
import os


def install_packages():
    """Install required packages without root"""
    packages = ["pandas", "numpy", "scikit-learn", "matplotlib"]
    for pkg in packages:
        try:
            __import__(pkg.replace("-", "_"))
        except ImportError:
            print(f"Installing {pkg}...")
            subprocess.check_call(
                [sys.executable, "-m", "pip", "install", "--user", "-q", pkg]
            )


def generate_prediction_data(n_train=200, n_test=50, n_snps=500):
    """Generate training and testing data for genomic prediction"""
    import numpy as np

    np.random.seed(42)

    # Generate genotypes
    train_geno = np.random.binomial(2, 0.3, (n_train, n_snps))
    test_geno = np.random.binomial(2, 0.3, (n_test, n_snps))

    # True breeding values (sum of random SNP effects)
    snp_effects = np.random.normal(0, 0.3, n_snps)

    # Training phenotypes = TBV + noise
    train_tbv = train_geno @ snp_effects
    train_pheno = train_tbv + np.random.normal(0, 1, n_train)

    # Test phenotypes (for validation)
    test_tbv = test_geno @ snp_effects
    test_pheno = test_tbv + np.random.normal(0, 1, n_test)

    return (
        train_geno,
        test_geno,
        train_pheno,
        test_pheno,
        train_tbv,
        test_tbv,
        snp_effects,
    )


def run_genomic_prediction(train_geno, test_geno, train_pheno, output_dir):
    """Run GBLUP genomic prediction"""
    import numpy as np
    import pandas as pd
    from sklearn.linear_model import Ridge

    print("Running genomic prediction (GBLUP)...")

    # Build genomic relationship matrix (G matrix)
    # G = ZZ' / sum(2*p*(1-p))
    p = train_geno.mean(axis=0) / 2
    Z = train_geno - 2 * p
    G = Z @ Z.T / (2 * np.sum(p * (1 - p)))

    # Fit prediction model (Ridge regression as GBLUP proxy)
    model = Ridge(alpha=1.0)
    model.fit(train_geno, train_pheno)

    # Predict for test individuals
    test_pred = model.predict(test_geno)

    # Calculate accuracy
    from scipy.stats import pearsonr

    # We'll store predictions for later validation
    results = pd.DataFrame({"predicted": test_pred, "index": range(len(test_pred))})
    results.to_csv(f"{output_dir}/gebv_predictions.csv", index=False)

    return model, test_pred


def create_prediction_plot(test_pred, test_tbv, output_dir):
    """Create genomic prediction plots"""
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import pearsonr

    print("Creating prediction plots...")

    # Calculate accuracy
    r, _ = pearsonr(test_pred, test_tbv)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # 1. Predicted vs True GEBV
    ax1 = axes[0]
    ax1.scatter(test_tbv, test_pred, alpha=0.6, s=50)
    ax1.plot(
        [test_tbv.min(), test_tbv.max()],
        [test_tbv.min(), test_tbv.max()],
        "r--",
        linewidth=2,
    )
    ax1.set_xlabel("True Breeding Value")
    ax1.set_ylabel("Predicted GEBV")
    ax1.set_title(f"GEBV Prediction\nr = {r:.3f}")
    ax1.grid(True, alpha=0.3)

    # 2. Distribution of GEBVs
    ax2 = axes[1]
    ax2.hist(test_tbv, bins=20, alpha=0.5, label="True", color="blue")
    ax2.hist(test_pred, bins=20, alpha=0.5, label="Predicted", color="red")
    ax2.set_xlabel("Breeding Value")
    ax2.set_ylabel("Frequency")
    ax2.set_title("Distribution of GEBVs")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # 3. Prediction accuracy by rank
    ax3 = axes[2]
    sorted_idx = np.argsort(test_tbv)[::-1]
    ax3.bar(range(len(test_pred)), test_pred[sorted_idx], alpha=0.7, label="Predicted")
    ax3.plot(
        range(len(test_pred)), test_tbv[sorted_idx], "r-", linewidth=2, label="True"
    )
    ax3.set_xlabel("Rank")
    ax3.set_ylabel("Breeding Value")
    ax3.set_title("Top Individuals Ranking")
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/genomic_prediction.png", dpi=150, bbox_inches="tight")
    print(f"Prediction plot saved: {output_dir}/genomic_prediction.png")

    return r


def main():
    print("=" * 60)
    print("Example: Genomic Prediction (GBLUP)")
    print("=" * 60)

    # Install packages
    print("\n[1/5] Installing dependencies...")
    install_packages()

    # Setup
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    # Generate data
    print("\n[2/5] Generating prediction data...")
    (
        train_geno,
        test_geno,
        train_pheno,
        test_pheno,
        train_tbv,
        test_tbv,
        snp_effects,
    ) = generate_prediction_data()
    print(f"  Training: {train_geno.shape[0]} individuals")
    print(f"  Test: {test_geno.shape[0]} individuals")
    print(f"  Markers: {train_geno.shape[1]}")

    # Run prediction
    print("\n[3/5] Running genomic prediction...")
    model, test_pred = run_genomic_prediction(
        train_geno, test_geno, train_pheno, output_dir
    )

    # Calculate accuracy
    from scipy.stats import pearsonr

    r, pval = pearsonr(test_pred, test_tbv)
    print(f"  Prediction accuracy (r): {r:.3f}")
    print(f"  p-value: {pval:.2e}")

    # Create plots
    print("\n[4/5] Creating visualizations...")
    r_value = create_prediction_plot(test_pred, test_tbv, output_dir)

    # Summary
    print("\n[5/5] Summary")
    print("=" * 40)
    print(f"Training samples: {train_geno.shape[0]}")
    print(f"Test samples: {test_geno.shape[0]}")
    print(f"SNP markers: {train_geno.shape[1]}")
    print(f"Prediction accuracy: r = {r:.3f}")
    print(f"\nOutputs:")
    print(f"  - {output_dir}/gebv_predictions.csv")
    print(f"  - {output_dir}/genomic_prediction.png")
    print("\n✅ Genomic prediction example complete!")
    print("\nIn QTLmax: Genomic Prediction → GBLUP")


if __name__ == "__main__":
    main()
