"""Plotting functions for meta-analysis."""

import logging
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from ..models import MetaAnalysisResult

logger = logging.getLogger(__name__)


def create_forest_plot(
    result: MetaAnalysisResult,
    output_path: Path,
    title: str | None = None,
    effect_label: str = "Odds Ratio",
) -> None:
    """
    Create forest plot for meta-analysis results.

    Args:
        result: Meta-analysis result
        output_path: Output file path
        title: Plot title
        effect_label: Label for effect measure
    """
    study_ids = result.model_details.get("study_ids", [])
    log_effects = result.model_details.get("log_effects", [])
    variances = result.model_details.get("variances", [])

    if not study_ids:
        logger.warning("No study data available for forest plot")
        return

    n_studies = len(study_ids)

    # Convert log effects to original scale
    effects = [math.exp(le) for le in log_effects]
    ci_lowers = [math.exp(le - 1.96 * math.sqrt(v)) for le, v in zip(log_effects, variances)]
    ci_uppers = [math.exp(le + 1.96 * math.sqrt(v)) for le, v in zip(log_effects, variances)]

    # Pooled effect
    pooled = math.exp(result.pooled_effect)
    pooled_lower = math.exp(result.pooled_ci_lower)
    pooled_upper = math.exp(result.pooled_ci_upper)

    # Create figure
    fig, ax = plt.subplots(figsize=(10, max(6, n_studies * 0.4 + 2)))

    # Plot individual studies
    y_positions = np.arange(n_studies, 0, -1)

    for i, (study_id, effect, ci_low, ci_high) in enumerate(
        zip(study_ids, effects, ci_lowers, ci_uppers)
    ):
        y = y_positions[i]

        # Error bar (CI)
        ax.plot([ci_low, ci_high], [y, y], "k-", linewidth=1.5)

        # Point estimate
        weight = 1.0 / variances[i] if variances[i] > 0 else 1.0
        size = min(200, 50 + weight * 10)
        ax.scatter(effect, y, s=size, c="black", marker="s", zorder=3)

        # Study label
        ax.text(-0.1, y, study_id, ha="right", va="center", fontsize=9, transform=ax.get_yaxis_transform())

        # Effect size text
        ax.text(
            1.02,
            y,
            f"{effect:.2f} [{ci_low:.2f}, {ci_high:.2f}]",
            ha="left",
            va="center",
            fontsize=8,
            transform=ax.get_yaxis_transform(),
        )

    # Pooled estimate (diamond)
    y_pooled = 0
    diamond_y = [y_pooled, y_pooled - 0.15, y_pooled, y_pooled + 0.15, y_pooled]
    diamond_x = [pooled_lower, pooled, pooled_upper, pooled, pooled_lower]
    ax.fill(diamond_x, diamond_y, color="navy", alpha=0.7, zorder=4)
    ax.plot(diamond_x, diamond_y, color="navy", linewidth=2, zorder=4)

    # Pooled label
    ax.text(
        -0.1,
        y_pooled,
        "Pooled",
        ha="right",
        va="center",
        fontsize=10,
        fontweight="bold",
        transform=ax.get_yaxis_transform(),
    )
    ax.text(
        1.02,
        y_pooled,
        f"{pooled:.2f} [{pooled_lower:.2f}, {pooled_upper:.2f}]",
        ha="left",
        va="center",
        fontsize=9,
        fontweight="bold",
        transform=ax.get_yaxis_transform(),
    )

    # Null line
    ax.axvline(x=1, color="gray", linestyle="--", linewidth=1, zorder=1)

    # Formatting
    ax.set_ylim(-0.5, n_studies + 0.5)
    ax.set_xlabel(effect_label, fontsize=11)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    # Log scale for x-axis
    ax.set_xscale("log")
    ax.set_xlim(
        min(min(ci_lowers) * 0.8, 0.5), max(max(ci_uppers) * 1.2, 2.0)
    )

    # Title
    if title is None:
        title = f"Forest Plot: {result.risk_factor}"
    ax.set_title(title, fontsize=13, fontweight="bold", pad=20)

    # Heterogeneity stats
    stats_text = (
        f"Heterogeneity: I² = {result.i_squared:.1f}%, "
        f"τ² = {result.tau_squared:.3f}, "
        f"Q = {result.q_statistic:.2f} (p = {result.q_p_value:.3f})"
    )
    ax.text(
        0.5,
        -0.05,
        stats_text,
        ha="center",
        va="top",
        fontsize=9,
        transform=ax.transAxes,
    )

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    logger.info(f"Forest plot saved to {output_path}")


def create_funnel_plot(
    result: MetaAnalysisResult, output_path: Path, title: str | None = None
) -> None:
    """
    Create funnel plot for publication bias assessment.

    Args:
        result: Meta-analysis result
        output_path: Output file path
        title: Plot title
    """
    log_effects = result.model_details.get("log_effects", [])
    variances = result.model_details.get("variances", [])

    if not log_effects:
        logger.warning("No study data available for funnel plot")
        return

    # Standard errors
    se = [math.sqrt(v) for v in variances]

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))

    # Scatter plot
    ax.scatter(log_effects, se, s=80, alpha=0.6, edgecolors="black", linewidth=1)

    # Pooled estimate line
    ax.axvline(x=result.pooled_effect, color="navy", linestyle="--", linewidth=2, label="Pooled Effect")

    # Pseudo 95% confidence interval funnel
    se_min = 0
    se_max = max(se) * 1.1

    se_range = np.linspace(se_min, se_max, 100)
    ci_lower = result.pooled_effect - 1.96 * se_range
    ci_upper = result.pooled_effect + 1.96 * se_range

    ax.plot(ci_lower, se_range, "gray", linestyle=":", linewidth=1, label="95% CI")
    ax.plot(ci_upper, se_range, "gray", linestyle=":", linewidth=1)

    # Invert y-axis (precision increases upward)
    ax.invert_yaxis()

    # Labels
    ax.set_xlabel("Log Odds Ratio", fontsize=11)
    ax.set_ylabel("Standard Error", fontsize=11)

    if title is None:
        title = f"Funnel Plot: {result.risk_factor}"
    ax.set_title(title, fontsize=13, fontweight="bold")

    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)

    # Add Egger test result if available
    if result.egger_p_value is not None:
        ax.text(
            0.05,
            0.95,
            f"Egger's test: p = {result.egger_p_value:.3f}",
            transform=ax.transAxes,
            verticalalignment="top",
            fontsize=9,
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
        )

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    logger.info(f"Funnel plot saved to {output_path}")


def create_sensitivity_plot(
    sensitivity_results: list[dict[str, Any]], output_path: Path, title: str = "Leave-One-Out Analysis"
) -> None:
    """
    Create plot for leave-one-out sensitivity analysis.

    Args:
        sensitivity_results: List of sensitivity analysis results
        output_path: Output file path
        title: Plot title
    """
    if not sensitivity_results:
        logger.warning("No sensitivity data available")
        return

    excluded = [r["excluded_study"] for r in sensitivity_results]
    pooled = [math.exp(r["pooled_effect"]) for r in sensitivity_results]
    ci_lower = [math.exp(r["pooled_ci_lower"]) for r in sensitivity_results]
    ci_upper = [math.exp(r["pooled_ci_upper"]) for r in sensitivity_results]

    fig, ax = plt.subplots(figsize=(10, max(6, len(excluded) * 0.3)))

    y_positions = np.arange(len(excluded), 0, -1)

    for i, (study, effect, ci_low, ci_high) in enumerate(zip(excluded, pooled, ci_lower, ci_upper)):
        y = y_positions[i]

        # Error bar
        ax.plot([ci_low, ci_high], [y, y], "k-", linewidth=1.5)

        # Point estimate
        ax.scatter(effect, y, s=100, c="black", marker="o", zorder=3)

        # Study label
        ax.text(-0.05, y, f"Excluding: {study}", ha="right", va="center", fontsize=8, transform=ax.get_yaxis_transform())

    # Null line
    ax.axvline(x=1, color="gray", linestyle="--", linewidth=1)

    ax.set_ylim(0, len(excluded) + 1)
    ax.set_xlabel("Odds Ratio", fontsize=11)
    ax.set_yticks([])
    ax.set_xscale("log")
    ax.set_title(title, fontsize=13, fontweight="bold")

    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    logger.info(f"Sensitivity plot saved to {output_path}")
