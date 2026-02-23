#!/usr/bin/env python3
"""
Create animated Manhattan plot for GWAS-LMM example
Shows SNPs appearing one by one to visualize the scan process
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import os


def create_manhattan_animation(output_file="output/manhattan_animation.gif"):
    """Create animated Manhattan plot"""

    # Generate GWAS data
    np.random.seed(42)
    n_snps = 200
    chromosomes = np.repeat([1, 2, 3, 4], n_snps // 4)
    positions = np.concatenate([np.arange(50) * 100000 for _ in range(4)])

    # P-values (most non-significant, a few significant)
    pvals = np.random.uniform(0.01, 0.99, n_snps)
    # Add significant hits
    pvals[25] = 5e-10  # Chr 2
    pvals[120] = 2e-8  # Chr 4

    neg_log_p = -np.log10(pvals)

    # Setup figure
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.set_xlim(0, n_snps)
    ax.set_ylim(0, max(neg_log_p) * 1.1)
    ax.axhline(
        y=-np.log10(5e-8),
        color="red",
        linestyle="--",
        linewidth=2,
        label="Genome-wide significance (5e-8)",
    )
    ax.set_xlabel("Chromosome", fontsize=12)
    ax.set_ylabel("-log10(p-value)", fontsize=12)
    ax.set_title(
        "GWAS-LMM Manhattan Plot Animation\nSNPs being tested sequentially",
        fontsize=14,
        fontweight="bold",
    )
    ax.legend(loc="upper right")

    # Create scatter points (initially empty)
    colors = ["#1f77b4" if c % 2 == 1 else "#ff7f0e" for c in chromosomes]
    scatter = ax.scatter([], [], c=[], s=20, alpha=0.7)

    # Add chromosome labels
    tick_positions = [25, 75, 125, 175]
    tick_labels = ["Chr 1", "Chr 2", "Chr 3", "Chr 4"]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)

    # Animation update function
    def update(frame):
        # Gradually reveal SNPs
        n_show = min((frame + 1) * 5, n_snps)

        x = list(range(n_show))
        y = neg_log_p[:n_show]
        c = colors[:n_show]

        scatter.set_offsets(np.column_stack([x, y]))
        scatter.set_color(c)

        # Add annotation for significant hit
        if n_show > 25 and frame < 50:
            ax.annotate(
                f"Significant hit!\np=5e-10",
                xy=(25, neg_log_p[25]),
                xytext=(50, max(neg_log_p) * 0.8),
                bbox=dict(boxstyle="round", facecolor="yellow", alpha=0.8),
                arrowprops=dict(arrowstyle="->", color="red", lw=2),
                fontsize=10,
                fontweight="bold",
            )

        return (scatter,)

    # Create animation
    n_frames = n_snps // 5
    anim = FuncAnimation(fig, update, frames=n_frames, interval=100, blit=False)

    # Save as GIF
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    writer = PillowWriter(fps=10)
    anim.save(output_file, writer=writer)

    plt.close()
    print(f"Animation saved: {output_file}")
    return output_file


if __name__ == "__main__":
    create_manhattan_animation()
