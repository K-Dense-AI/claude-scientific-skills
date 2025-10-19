#!/usr/bin/env python3
"""
Distributed GRN inference script using arboreto with custom Dask configuration.

This script demonstrates how to use arboreto with a custom Dask LocalCluster
for better control over computational resources.

Usage:
    python distributed_inference.py <expression_file> [options]

Example:
    python distributed_inference.py expression_data.tsv -t tf_names.txt -w 8 -m 4GB
"""

import argparse
import pandas as pd
from dask.distributed import Client, LocalCluster
from arboreto.algo import grnboost2
from arboreto.utils import load_tf_names


def main():
    parser = argparse.ArgumentParser(
        description='Distributed GRN inference using GRNBoost2 with custom Dask cluster'
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
        '-w', '--workers',
        type=int,
        help='Number of Dask workers',
        default=4
    )
    parser.add_argument(
        '-m', '--memory-limit',
        help='Memory limit per worker (e.g., "4GB", "2000MB")',
        default='4GB'
    )
    parser.add_argument(
        '--threads',
        type=int,
        help='Threads per worker',
        default=2
    )
    parser.add_argument(
        '--dashboard-port',
        type=int,
        help='Port for Dask dashboard (default: 8787)',
        default=8787
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

    # Set up Dask cluster
    print(f"\nSetting up Dask LocalCluster...")
    print(f"  Workers: {args.workers}")
    print(f"  Threads per worker: {args.threads}")
    print(f"  Memory limit per worker: {args.memory_limit}")
    print(f"  Dashboard: http://localhost:{args.dashboard_port}")

    cluster = LocalCluster(
        n_workers=args.workers,
        threads_per_worker=args.threads,
        memory_limit=args.memory_limit,
        diagnostics_port=args.dashboard_port
    )
    client = Client(cluster)

    print(f"\nDask cluster ready!")
    print(f"  Dashboard available at: {client.dashboard_link}")

    # Run GRNBoost2
    print("\nRunning GRNBoost2 inference with distributed computation...")
    print("  (Monitor progress via the Dask dashboard)")

    try:
        network = grnboost2(
            expression_data=expression_data,
            tf_names=tf_names,
            seed=args.seed,
            client_or_address=client
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

    finally:
        # Clean up Dask resources
        print("\nClosing Dask cluster...")
        client.close()
        cluster.close()


if __name__ == '__main__':
    main()
