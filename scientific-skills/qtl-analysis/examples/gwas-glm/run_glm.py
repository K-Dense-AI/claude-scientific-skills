#!/usr/bin/env python3
# Copyright 2026 Clayton Young (borealBytes / Superior Byte Works, LLC)
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Author: Clayton Young <Clayton@SuperiorByteWorks.com>
# LinkedIn: https://linkedin.com/in/claytoneyoung/
# GitHub: https://github.com/borealBytes

#!/usr/bin/env python3
"""
Example: GWAS with Generalized Linear Model (GLM)

This example demonstrates GWAS using GLM (logistic/linear regression) without
the kinship correction. This is the faster alternative to LMM for initial screens.

Equivalent to QTLmax: "GWAS procedure (Generalized linear model)"
https://open.qtlmax.com/guide/index.php/2025/07/11/gwas-procedure-generalized-linear-model/

Auto-installs: pandas, numpy, scipy, statsmodels, matplotlib
"""

import subprocess
import sys
import os


def install_packages():
    """Install required packages without root"""
    packages = ["pandas", "numpy", "scipy", "statsmodels", "matplotlib", "scikit-learn"]
    for pkg in packages:
        try:
            __import__(pkg.replace("-", "_"))
        except ImportError:
            print(f"Installing {pkg}...")
            subprocess.check_call(
                [sys.executable, "-m", "pip", "install", "--user", "-q", pkg]
            )


def generate_synthetic_data(n_individuals=200, n_snps=5000):
    """Generate synthetic GWAS data"""
    import numpy as np
    import pandas as pd

    np.random.seed(42)

    # Generate genotypes
    genotypes = np.random.binomial(2, 0.3, (n_individuals, n_snps))

    # Generate binary phenotype (disease status)
    # Some SNPs are causal
    causal_effects = np.random.normal(0, 0.8, 10)
    causal_snps = np.random.choice(n_snps, 10, replace=False)

    logit_p = -0.5 + np.dot(genotypes[:, causal_snps], causal_effects)
    probability = 1 / (1 + np.exp(-logit_p))
    phenotype = (np.random.random(n_individuals) < probability).astype(int)

    # Create dataframes
    snp_ids = [f"rs{i}" for i in range(n_snps)]
    sample_ids = [f"Sample{i}" for i in range(n_individuals)]

    # PLINK-style PED file
    ped_data = []
    for i, (sample_id, pheno) in enumerate(zip(sample_ids, phenotype)):
        row = [
            sample_id,
            sample_id,
            "0",
            "0",
            "0",
            "2",
        ]  # Family ID, ID, dad, mom, sex, pheno
        row.extend(genotypes[i, :].tolist())
        ped_data.append(row)

    # MAP file
    map_data = []
    for i, snp_id in enumerate(snp_ids):
        chrom = (i // 500) + 1
        map_data.append([chrom, snp_id, "0", i * 10000])

    return ped_data, map_data, phenotype, sample_ids, snp_ids, causal_snps


def run_glm_gwas(ped_data, map_data, phenotype, sample_ids, snp_ids, output_dir):
    """Run GLM-based GWAS"""
    import numpy as np
    import pandas as pd
    from scipy import stats
    import statsmodels.api as sm

    print("Running GLM-based GWAS...")

    # Convert PED to array
    genotypes = np.array([row[6:] for row in ped_data])

    results = []
    for i in range(genotypes.shape[1]):
        geno = genotypes[:, i]

        # Logistic regression for binary phenotype
        try:
            X = sm.add_constant(geno)
            model = sm.Logit(phenotype, X).fit(disp=0, method="newton")
            pval = model.pvalues[1]
            beta = model.params[1]
            se = model.bse[1]
        except:
            pval = 1.0
            beta = 0.0
            se = 0.0

        chrom = (i // 500) + 1
        pos = i * 10000

        results.append(
            {
                "CHR": chrom,
                "SNP": snp_ids[i],
                "POS": pos,
                "A1": "A",
                "A2": "G",
                "BETA": beta,
                "SE": se,
                "P": pval,
                "OR": np.exp(beta),
            }
        )

    results_df = pd.DataFrame(results)
    results_df.to_csv(f"{output_dir}/glm_gwas_results.csv", index=False)

    return results_df


def create_manhattan_plot(results_df, output_file):
    """Create Manhattan plot"""
    import numpy as np
    import matplotlib.pyplot as plt

    print("Creating Manhattan plot...")

    fig, ax = plt.subplots(figsize=(14, 6))

    # -log10(p) transform
    results_df["NEG_LOG10_P"] = -np.log10(results_df["P"].clip(lower=1e-300))

    # Plot by chromosome
    chroms = sorted(results_df["CHR"].unique())
    x_pos = 0
    tick_positions = []
    tick_labels = []

    for i, chrom in enumerate(chroms):
        chrom_data = results_df[results_df["CHR"] == chrom].sort_values("POS")
        color = "#1f77b4" if i % 2 == 0 else "#ff7f0e"

        ax.scatter(
            range(x_pos, x_pos + len(chrom_data)),
            chrom_data["NEG_LOG10_P"],
            c=color,
            s=8,
            alpha=0.6,
        )

        tick_positions.append(x_pos + len(chrom_data) / 2)
        tick_labels.append(str(chrom))
        x_pos += len(chrom_data)

    # Significance threshold
    ax.axhline(
        y=-np.log10(5e-8),
        color="red",
        linestyle="--",
        linewidth=1.5,
        label="Genome-wide significance (5e-8)",
    )
    ax.axhline(
        y=-np.log10(1e-5),
        color="blue",
        linestyle="--",
        linewidth=1,
        label="Suggestive (1e-5)",
    )

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    ax.set_xlabel("Chromosome", fontsize=12)
    ax.set_ylabel("-log10(p-value)", fontsize=12)
    ax.set_title(
        "GLM-GWAS Manhattan Plot\n(Logistic Regression)", fontsize=14, fontweight="bold"
    )
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    print(f"Manhattan plot saved: {output_file}")


def create_qq_plot(results_df, output_file):
    """Create QQ plot"""
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy import stats

    print("Creating QQ plot...")

    pvals = results_df["P"].dropna()
    n = len(pvals)

    # Calculate lambda GC
    chi2 = stats.chi2.ppf(1 - pvals, df=1)
    lambda_gc = np.median(chi2) / stats.chi2.ppf(0.5, df=1)

    # Expected vs observed
    expected = -np.log10(np.linspace(1 / n, 1, n))
    observed = -np.log10(sorted(pvals))

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(expected, observed, s=5, alpha=0.5, c="#1f77b4")
    ax.plot([0, max(expected)], [0, max(expected)], "r--", linewidth=2, label="y=x")

    ax.set_xlabel("Expected -log10(p)", fontsize=12)
    ax.set_ylabel("Observed -log10(p)", fontsize=12)
    ax.set_title(f"QQ Plot (λ GC = {lambda_gc:.3f})", fontsize=14, fontweight="bold")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    print(f"QQ plot saved: {output_file}")


def main():
    print("=" * 60)
    print("Example: GWAS with Generalized Linear Model (GLM)")
    print("=" * 60)

    # Install packages
    print("\n[1/5] Installing dependencies...")
    install_packages()

    # Setup
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    # Generate data
    print("\n[2/5] Generating synthetic GWAS data...")
    ped_data, map_data, phenotype, sample_ids, snp_ids, causal_snps = (
        generate_synthetic_data()
    )
    print(f"  Generated: {len(sample_ids)} samples, {len(snp_ids)} SNPs")
    print(f"  Causal SNPs: {causal_snps[:5]}... (10 total)")

    # Run GWAS
    print("\n[3/5] Running GLM-GWAS...")
    results_df = run_glm_gwas(
        ped_data, map_data, phenotype, sample_ids, snp_ids, output_dir
    )

    # Find significant SNPs
    sig_snps = results_df[results_df["P"] < 5e-8]
    print(f"  Significant SNPs (p < 5e-8): {len(sig_snps)}")

    # Create plots
    print("\n[4/5] Creating visualizations...")
    create_manhattan_plot(results_df, f"{output_dir}/glm_manhattan.png")
    create_qq_plot(results_df, f"{output_dir}/glm_qqplot.png")

    # Summary
    print("\n[5/5] Summary")
    print("=" * 40)
    print(f"Samples: {len(sample_ids)}")
    print(f"SNPs tested: {len(snp_ids)}")
    print(f"Significant (p < 5e-8): {len(sig_snps)}")
    import numpy as np
    pvals = results_df['P'].values
    pvals_valid = pvals[pvals > 0]
    if len(pvals_valid) > 0:
        genomic_inflation = np.median(-np.log10(pvals_valid)) / 0.455
        print(f"Genomic inflation λ: {genomic_inflation:.3f}")
    else:
        print(f"Genomic inflation λ: N/A (no valid p-values)")
    print(f"\nOutputs:")
    print(f"  - {output_dir}/glm_gwas_results.csv")
    print(f"  - {output_dir}/glm_manhattan.png")
    print(f"  - {output_dir}/glm_qqplot.png")
    print("\n✅ GLM-GWAS example complete!")
    print("\nIn QTLmax: GWAS → Generalized linear model procedure")


if __name__ == "__main__":
    main()
