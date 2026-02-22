#!/usr/bin/env python3
"""
Unified CLI for QTL Analysis

Provides a single entry point for all QTL analysis workflows:
- gwas: GWAS with GEMMA or PLINK
- eqtl: eQTL mapping with tensorQTL
- lodscan: Classical QTL with R/qtl2
- manhattan: Create Manhattan plots
- qqplot: Create QQ plots
- kinship: Calculate kinship matrices
"""

import argparse
import sys
import subprocess
import json
from pathlib import Path


def run_gwas(args):
    """Run GWAS with GEMMA or PLINK"""
    if args.method == "lmm":
        print(f"Running LMM-GWAS with GEMMA...")
        # Step 1: Calculate kinship
        subprocess.run(
            [
                "gemma",
                "-bfile",
                args.geno.replace(".bed", ""),
                "-gk",
                "1",
                "-o",
                "kinship",
            ],
            check=True,
        )

        # Step 2: Run GWAS
        subprocess.run(
            [
                "gemma",
                "-bfile",
                args.geno.replace(".bed", ""),
                "-k",
                "output/kinship.sXX.txt",
                "-lmm",
                "4",
                "-p",
                args.pheno,
                "-o",
                "gwas",
            ],
            check=True,
        )

        print(f"✅ GWAS complete. Results: output/gwas.assoc.txt")

    elif args.method == "glm":
        print(f"Running GLM-GWAS with PLINK...")
        subprocess.run(
            [
                "plink2",
                "--bfile",
                args.geno.replace(".bed", ""),
                "--pheno",
                args.pheno,
                "--glm",
                "allow-no-covars",
                "--out",
                "gwas",
            ],
            check=True,
        )

        print(f"✅ GWAS complete. Results: gwas.PHENO1.glm.linear")


def run_eqtl(args):
    """Run eQTL mapping with tensorQTL"""
    print(f"Running {args.mode}-eQTL with tensorQTL...")

    script = f"""
import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis

# Load data
genotype_df, variant_df = genotypeio.load_genotypes('{args.geno}')
expression_df = pd.read_parquet('{args.expr}')
covariates_df = pd.read_csv('{args.covariates}', index_col=0) if '{args.covariates}' else None

# Run eQTL mapping
results = cis.map_cis(
    genotype_df, expression_df,
    covariates_df=covariates_df,
    window={args.window},
    nperm=10000
)

results.to_csv('{args.output}/cis_eqtl_results.csv')
print(f"Found {{len(results[results['pval'] < 5e-8])}} significant eQTLs")
"""

    subprocess.run(["python", "-c", script], check=True)
    print(f"✅ eQTL mapping complete. Results: {args.output}/cis_eqtl_results.csv")


def run_lodscan(args):
    """Run classical QTL LOD scan with R/qtl2"""
    print(f"Running LOD scan with R/qtl2...")

    r_script = f"""
library(qtl2)

# Load cross data
cross <- read_cross2("{args.cross}")

# Insert pseudomarkers
map <- insert_pseudomarkers(cross$gmap, step=1)

# Calculate genotype probabilities
pr <- calc_genoprob(cross, map, error_prob=0.002)

# Perform genome scan
out <- scan1(pr, cross$pheno)

# Permutation test for significance
if ({args.perms} > 0) {{
    operm <- scan1perm(pr, cross$pheno, n_perm={args.perms})
    threshold <- summary(operm, alpha=c(0.05, 0.1))
    write.csv(threshold, "{args.output}/thresholds.csv")
}}

# Write results
write.csv(out, "{args.output}/lod_scores.csv")

# Find peaks
peaks <- find_peaks(out, map, threshold=3, drop=1.5)
write.csv(peaks, "{args.output}/peaks.csv")

cat("✅ LOD scan complete\\n")
cat("Found", nrow(peaks), "QTL peaks\\n")
"""

    subprocess.run(["Rscript", "-e", r_script], check=True)
    print(f"✅ LOD scan complete. Results in {args.output}/")


def plot_manhattan(args):
    """Create Manhattan plot from GWAS results"""
    import pandas as pd
    import matplotlib.pyplot as plt

    print(f"Creating Manhattan plot...")

    # Load results (GEMMA format)
    df = pd.read_csv(args.input, sep="\t")

    # Calculate -log10(p)
    df["neg_log10_p"] = -df["p_lrt"].apply(lambda x: -(10**x))

    # Create plot
    fig, ax = plt.subplots(figsize=(14, 6))

    # Color by chromosome
    colors = {0: "#1f77b4", 1: "#ff7f0e"}
    df["color"] = df["chr"] % 2

    # Plot points
    for color_val, color_code in colors.items():
        subset = df[df["color"] == color_val]
        ax.scatter(subset["pos"], subset["neg_log10_p"], c=color_code, s=5, alpha=0.6)

    # Add significance lines
    if args.significance:
        ax.axhline(
            y=-np.log10(args.significance), color="red", linestyle="--", linewidth=1.5
        )

    ax.set_xlabel("Chromosomal Position")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title("Manhattan Plot")

    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print(f"✅ Manhattan plot saved: {args.output}")


def plot_qq(args):
    """Create QQ plot from p-values"""
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy import stats

    print(f"Creating QQ plot...")

    # Load p-values
    df = pd.read_csv(args.input, sep="\t")
    pvals = df["p_lrt"].dropna()

    # Calculate expected p-values
    n = len(pvals)
    expected = -np.log10(np.linspace(1 / n, 1, n))
    observed = -np.log10(sorted(pvals))

    # Calculate lambda GC
    chi2 = stats.chi2.ppf(1 - pvals, df=1)
    lambda_gc = np.median(chi2) / stats.chi2.ppf(0.5, df=1)

    # Create plot
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(expected, observed, s=5, alpha=0.6)
    ax.plot([0, max(expected)], [0, max(expected)], "r--", linewidth=1.5, label="y=x")

    ax.set_xlabel("Expected -log10(p)")
    ax.set_ylabel("Observed -log10(p)")
    ax.set_title(f"QQ Plot (λ GC = {lambda_gc:.3f})")
    ax.legend()

    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print(f"✅ QQ plot saved: {args.output}")
    print(f"Genomic inflation factor (λ GC): {lambda_gc:.3f}")


def calc_kinship(args):
    """Calculate kinship matrix"""
    print(f"Calculating {args.method} kinship matrix...")

    if args.method == "vanraden":
        subprocess.run(
            [
                "gemma",
                "-bfile",
                args.geno.replace(".bed", ""),
                "-gk",
                "1",
                "-o",
                args.output,
            ],
            check=True,
        )
    elif args.method == "sample":
        subprocess.run(
            [
                "gemma",
                "-bfile",
                args.geno.replace(".bed", ""),
                "-gk",
                "2",
                "-o",
                args.output,
            ],
            check=True,
        )
    elif args.method == "dominant":
        # Calculate using Python/numpy for dominant model
        import numpy as np
        import pandas as pd

        # Load genotypes (simplified - assumes PLINK format)
        print("Loading genotypes...")
        subprocess.run(
            [
                "plink",
                "--bfile",
                args.geno.replace(".bed", ""),
                "--recode",
                "A",
                "--out",
                "temp_geno",
            ],
            check=True,
        )

        # Calculate dominant kinship
        geno = pd.read_csv("temp_geno.raw", sep="\s+", header=0)
        geno_mat = geno.iloc[:, 6:].values  # Skip first 6 columns
        geno_mat = np.nan_to_num(geno_mat)

        # VanRaden method
        p = np.mean(geno_mat, axis=0) / 2
        W = geno_mat - 2 * p
        K = np.dot(W, W.T) / (2 * np.sum(p * (1 - p)))

        # Save
        np.savetxt(f"{args.output}.txt", K, delimiter="\t")
        print(f"✅ Kinship matrix saved: {args.output}.txt")

    print(f"✅ Kinship calculation complete")


def main():
    parser = argparse.ArgumentParser(description="Unified CLI for QTL Analysis")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # GWAS command
    gwas_parser = subparsers.add_parser("gwas", help="Run GWAS")
    gwas_parser.add_argument("--geno", required=True, help="Genotype file prefix")
    gwas_parser.add_argument("--pheno", required=True, help="Phenotype file")
    gwas_parser.add_argument(
        "--method", choices=["lmm", "glm"], default="lmm", help="GWAS method"
    )
    gwas_parser.add_argument("--output", default="gwas_results", help="Output prefix")

    # eQTL command
    eqtl_parser = subparsers.add_parser("eqtl", help="Run eQTL mapping")
    eqtl_parser.add_argument("--geno", required=True, help="Genotype VCF")
    eqtl_parser.add_argument("--expr", required=True, help="Expression BED")
    eqtl_parser.add_argument("--mode", choices=["cis", "trans"], default="cis")
    eqtl_parser.add_argument("--covariates", help="Covariates file")
    eqtl_parser.add_argument(
        "--window", type=int, default=1000000, help="cis window (bp)"
    )
    eqtl_parser.add_argument(
        "--output", default="eqtl_results", help="Output directory"
    )

    # LOD scan command
    lod_parser = subparsers.add_parser("lodscan", help="Run classical QTL LOD scan")
    lod_parser.add_argument("--cross", required=True, help="Cross data JSON")
    lod_parser.add_argument(
        "--method",
        default="haley-knott",
        choices=["haley-knott", "imp", "em"],
        help="Analysis method",
    )
    lod_parser.add_argument(
        "--perms", type=int, default=1000, help="Number of permutations"
    )
    lod_parser.add_argument("--output", default="lod_results", help="Output directory")

    # Manhattan plot command
    manhattan_parser = subparsers.add_parser("manhattan", help="Create Manhattan plot")
    manhattan_parser.add_argument("--input", required=True, help="GWAS results file")
    manhattan_parser.add_argument("--output", required=True, help="Output PNG")
    manhattan_parser.add_argument(
        "--significance", type=float, default=5e-8, help="Significance threshold"
    )

    # QQ plot command
    qq_parser = subparsers.add_parser("qqplot", help="Create QQ plot")
    qq_parser.add_argument("--input", required=True, help="GWAS results file")
    qq_parser.add_argument("--output", required=True, help="Output PNG")

    # Kinship command
    kinship_parser = subparsers.add_parser("kinship", help="Calculate kinship matrix")
    kinship_parser.add_argument("--geno", required=True, help="Genotype file prefix")
    kinship_parser.add_argument(
        "--method",
        choices=["vanraden", "sample", "dominant"],
        default="vanraden",
        help="Kinship method",
    )
    kinship_parser.add_argument("--output", default="kinship", help="Output prefix")

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    # Dispatch to appropriate handler
    command_map = {
        "gwas": run_gwas,
        "eqtl": run_eqtl,
        "lodscan": run_lodscan,
        "manhattan": plot_manhattan,
        "qqplot": plot_qq,
        "kinship": calc_kinship,
    }

    try:
        command_map[args.command](args)
    except subprocess.CalledProcessError as e:
        print(f"❌ Error: Command failed with return code {e.returncode}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
