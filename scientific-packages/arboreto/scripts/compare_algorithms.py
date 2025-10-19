#!/usr/bin/env python3
"""
Compare GRNBoost2 and GENIE3 algorithms on the same dataset.

This script runs both algorithms on the same expression data and compares:
- Runtime
- Number of predicted links
- Top predicted relationships
- Overlap between predictions

Usage:
    python compare_algorithms.py <expression_file> [options]

Example:
    python compare_algorithms.py expression_data.tsv -t tf_names.txt
"""

import argparse
import time
import pandas as pd
from arboreto.algo import grnboost2, genie3
from arboreto.utils import load_tf_names


def compare_networks(network1, network2, name1, name2, top_n=100):
    """Compare two inferred networks."""
    print(f"\n{'='*60}")
    print("Network Comparison")
    print(f"{'='*60}")

    # Basic statistics
    print(f"\n{name1} Statistics:")
    print(f"  Total links: {len(network1)}")
    print(f"  Unique TFs: {network1['TF'].nunique()}")
    print(f"  Unique targets: {network1['target'].nunique()}")
    print(f"  Importance range: [{network1['importance'].min():.3f}, {network1['importance'].max():.3f}]")

    print(f"\n{name2} Statistics:")
    print(f"  Total links: {len(network2)}")
    print(f"  Unique TFs: {network2['TF'].nunique()}")
    print(f"  Unique targets: {network2['target'].nunique()}")
    print(f"  Importance range: [{network2['importance'].min():.3f}, {network2['importance'].max():.3f}]")

    # Compare top predictions
    print(f"\nTop {top_n} Predictions Overlap:")

    # Create edge sets for top N predictions
    top_edges1 = set(
        zip(network1.head(top_n)['TF'], network1.head(top_n)['target'])
    )
    top_edges2 = set(
        zip(network2.head(top_n)['TF'], network2.head(top_n)['target'])
    )

    # Calculate overlap
    overlap = top_edges1 & top_edges2
    only_net1 = top_edges1 - top_edges2
    only_net2 = top_edges2 - top_edges1

    overlap_pct = (len(overlap) / top_n) * 100

    print(f"  Shared edges: {len(overlap)} ({overlap_pct:.1f}%)")
    print(f"  Only in {name1}: {len(only_net1)}")
    print(f"  Only in {name2}: {len(only_net2)}")

    # Show some example overlapping edges
    if overlap:
        print(f"\nExample overlapping predictions:")
        for i, (tf, target) in enumerate(list(overlap)[:5], 1):
            print(f"  {i}. {tf} -> {target}")


def main():
    parser = argparse.ArgumentParser(
        description='Compare GRNBoost2 and GENIE3 algorithms'
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
        '--grnboost2-output',
        help='Output file path for GRNBoost2 results',
        default='grnboost2_network.tsv'
    )
    parser.add_argument(
        '--genie3-output',
        help='Output file path for GENIE3 results',
        default='genie3_network.tsv'
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
    parser.add_argument(
        '--top-n',
        type=int,
        help='Number of top predictions to compare (default: 100)',
        default=100
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
    print("\n" + "="*60)
    print("Running GRNBoost2...")
    print("="*60)
    start_time = time.time()

    grnboost2_network = grnboost2(
        expression_data=expression_data,
        tf_names=tf_names,
        seed=args.seed
    )

    grnboost2_time = time.time() - start_time
    print(f"GRNBoost2 completed in {grnboost2_time:.2f} seconds")

    # Save GRNBoost2 results
    grnboost2_network.to_csv(args.grnboost2_output, sep='\t', index=False)
    print(f"Results saved to {args.grnboost2_output}")

    # Run GENIE3
    print("\n" + "="*60)
    print("Running GENIE3...")
    print("="*60)
    start_time = time.time()

    genie3_network = genie3(
        expression_data=expression_data,
        tf_names=tf_names,
        seed=args.seed
    )

    genie3_time = time.time() - start_time
    print(f"GENIE3 completed in {genie3_time:.2f} seconds")

    # Save GENIE3 results
    genie3_network.to_csv(args.genie3_output, sep='\t', index=False)
    print(f"Results saved to {args.genie3_output}")

    # Compare runtimes
    print("\n" + "="*60)
    print("Runtime Comparison")
    print("="*60)
    print(f"GRNBoost2: {grnboost2_time:.2f} seconds")
    print(f"GENIE3: {genie3_time:.2f} seconds")
    speedup = genie3_time / grnboost2_time
    print(f"Speedup: {speedup:.2f}x (GRNBoost2 is {speedup:.2f}x faster)")

    # Compare networks
    compare_networks(
        grnboost2_network,
        genie3_network,
        "GRNBoost2",
        "GENIE3",
        top_n=args.top_n
    )

    print("\n" + "="*60)
    print("Comparison complete!")
    print("="*60)


if __name__ == '__main__':
    main()
