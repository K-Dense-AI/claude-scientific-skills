#!/usr/bin/env python3
"""
Generate all example visualizations for the QTL Analysis skill
Creates output PNGs and GIFs for every example
"""

import os
import sys
import subprocess


def install_packages(packages):
    """Install packages without root"""
    for pkg in packages:
        try:
            __import__(pkg.replace("-", "_"))
        except ImportError:
            print(f"Installing {pkg}...")
            subprocess.check_call(
                [sys.executable, "-m", "pip", "install", "--user", "-q", pkg]
            )


def generate_all_outputs():
    """Generate all example outputs"""
    print("=" * 70)
    print("Generating QTL Analysis Example Visualizations")
    print("=" * 70)

    examples_dir = "scientific-skills/qtl-analysis/examples"
    examples = [
        "gwas-lmm",
        "gwas-glm",
        "eqtl-cis",
        "classical-qtl",
        "ld-decay",
        "population-structure",
        "admixture",
        "kmeans-clustering",
        "genomic-prediction",
        "marker-selection",
        "vcf-validation",
        "snp-filtering",
        "phenotype-plots",
        "imputation",
        "blup",
    ]

    success = []
    failed = []

    for example in examples:
        print(f"\n{'=' * 70}")
        print(f"Processing: {example}")
        print("=" * 70)

        example_path = os.path.join(examples_dir, example)
        if not os.path.exists(example_path):
            print(f"  ⚠️  Directory not found: {example_path}")
            failed.append(example)
            continue

        # Create output directory
        output_dir = os.path.join(example_path, "output")
        os.makedirs(output_dir, exist_ok=True)

        # Find Python files
        py_files = [f for f in os.listdir(example_path) if f.endswith(".py")]

        if not py_files:
            print(f"  ⚠️  No Python files found")
            failed.append(example)
            continue

        # Run the main script
        main_script = py_files[0]  # Use first .py file
        script_path = os.path.join(example_path, main_script)

        print(f"  Running: {main_script}")
        try:
            # Change to example directory and run
            original_dir = os.getcwd()
            os.chdir(example_path)
            result = subprocess.run(
                [sys.executable, main_script],
                capture_output=True,
                text=True,
                timeout=120,
            )
            os.chdir(original_dir)

            if result.returncode == 0:
                print(f"  ✅ Success")
                success.append(example)
            else:
                print(f"  ❌ Failed: {result.stderr[:200]}")
                failed.append(example)
        except Exception as e:
            os.chdir(original_dir)
            print(f"  ❌ Error: {str(e)[:200]}")
            failed.append(example)

    # Summary
    print(f"\n{'=' * 70}")
    print("SUMMARY")
    print("=" * 70)
    print(f"✅ Successful: {len(success)}/{len(examples)}")
    print(f"❌ Failed: {len(failed)}")

    if failed:
        print(f"\nFailed examples: {', '.join(failed)}")

    return len(success), len(failed)


if __name__ == "__main__":
    success, failed = generate_all_outputs()
    sys.exit(0 if failed == 0 else 1)
