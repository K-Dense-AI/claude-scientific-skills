#!/usr/bin/env python3
"""
Basic GRN inference script using arboreto GRNBoost2.

This script demonstrates the standard workflow for gene regulatory network inference:
1. Load expression data
2. Optionally load transcription factor names
3. Run GRNBoost2 inference
4. Save results

Usage:
    python basic_grn_inference.py <expression_file> [options]

Example:
    python basic_grn_inference.py expression_data.tsv -t tf_names.txt -o network.tsv
"""

import argparse
import pandas as pd
from arboreto.algo import grnboost2
from arboreto.utils import load_tf_names


def main():
    parser = argparse.ArgumentParser(
        description='Infer gene regulatory network using GRNBoost2'
    )
    parser.add_argument(
        'expression_file',
        help='Path to expression data file (TSV/CSV format)'
    )
    parser.add_argument(
        '-t', '--tf-file',
        help='Path to file containing transcription factor names (one per line)',
        default=None
    )
    parser.add_argument(
        '-o', '--output',
        help='Output file path for network results',
        default='network_output.tsv'
    )
    parser.add_argument(
        '-s', '--seed',
        type=int,
        help='Random seed for reproducibility',
        default=42
    )
    parser.add_argument(
        '--sep',
        help='Separator for input file (default: tab)',
        default='\t'
    )
    parser.add_argument(
        '--transpose',
        action='store_true',
        help='Transpose the expression matrix (use if genes are rows)'
    )

    args = parser.parse_args()

    # Load expression data
    print(f"Loading expression data from {args.expression_file}...")
    expression_data = pd.read_csv(args.expression_file, sep=args.sep, index_col=0)

    # Transpose if needed
    if args.transpose:
        print("Transposing expression matrix...")
        expression_data = expression_data.T

    print(f"Expression data shape: {expression_data.shape}")
    print(f"  Observations (rows): {expression_data.shape[0]}")
    print(f"  Genes (columns): {expression_data.shape[1]}")

    # Load TF names if provided
    tf_names = None
    if args.tf_file:
        print(f"Loading transcription factor names from {args.tf_file}...")
        tf_names = load_tf_names(args.tf_file)
        print(f"  Found {len(tf_names)} transcription factors")
    else:
        print("No TF file provided. Using all genes as potential regulators.")

    # Run GRNBoost2
    print("\nRunning GRNBoost2 inference...")
    print("  (This may take a while depending on dataset size)")

    network = grnboost2(
        expression_data=expression_data,
        tf_names=tf_names,
        seed=args.seed
    )

    print(f"\nInference complete!")
    print(f"  Total regulatory links inferred: {len(network)}")
    print(f"  Unique TFs: {network['TF'].nunique()}")
    print(f"  Unique targets: {network['target'].nunique()}")

    # Save results
    print(f"\nSaving results to {args.output}...")
    network.to_csv(args.output, sep='\t', index=False)

    # Display top 10 predictions
    print("\nTop 10 predicted regulatory relationships:")
    print(network.head(10).to_string(index=False))

    print("\nDone!")


if __name__ == '__main__':
    main()
