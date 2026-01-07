"""PRISMA flow diagram generation."""

import logging
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

logger = logging.getLogger(__name__)


class PRISMAFlowData:
    """Data for PRISMA flow diagram."""

    def __init__(self) -> None:
        """Initialize PRISMA flow data."""
        self.identification_databases: dict[str, int] = {}
        self.identification_other: int = 0
        self.records_after_dedup: int = 0
        self.duplicates_removed: int = 0
        self.records_screened: int = 0
        self.records_excluded: int = 0
        self.fulltext_assessed: int = 0
        self.fulltext_excluded: int = 0
        self.exclusion_reasons: dict[str, int] = {}
        self.studies_included: int = 0

    def from_screening_stats(self, stats: dict[str, Any]) -> None:
        """Populate from screening database statistics."""
        # This would be populated from actual screening data
        pass


def generate_prisma_diagram(
    flow_data: PRISMAFlowData, output_path: Path, title: str = "PRISMA Flow Diagram"
) -> None:
    """
    Generate PRISMA 2020 flow diagram.

    Args:
        flow_data: PRISMA flow data
        output_path: Output file path
        title: Diagram title
    """
    fig, ax = plt.subplots(figsize=(12, 16))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 20)
    ax.axis("off")

    # Title
    ax.text(5, 19.5, title, ha="center", va="top", fontsize=16, fontweight="bold")

    # Box style
    box_style = "round,pad=0.1"
    box_color = "#E8F4F8"
    excluded_color = "#FFE8E8"

    # Identification
    y = 18
    _add_box(
        ax,
        2.5,
        y,
        2,
        1.5,
        f"Records identified from:\nDatabases (n = {sum(flow_data.identification_databases.values())})",
        box_color,
    )
    _add_box(
        ax,
        5.5,
        y,
        2,
        1.5,
        f"Records identified from:\nOther sources (n = {flow_data.identification_other})",
        box_color,
    )

    # Database breakdown (small text)
    if flow_data.identification_databases:
        db_text = "\n".join(
            [f"{db}: {count}" for db, count in flow_data.identification_databases.items()]
        )
        ax.text(3.5, y - 0.8, db_text, ha="center", va="center", fontsize=7)

    # Duplicates removed
    y = 15.5
    _add_box(
        ax,
        2.5,
        y,
        5,
        1.2,
        f"Records after duplicates removed\n(n = {flow_data.records_after_dedup})",
        box_color,
    )
    _add_box(
        ax,
        8,
        y,
        1.5,
        1.2,
        f"Duplicates removed\n(n = {flow_data.duplicates_removed})",
        excluded_color,
    )

    _add_arrow(ax, 5, 16.5, 5, 16.2)

    # Screening
    y = 13.5
    _add_box(
        ax, 2.5, y, 5, 1.2, f"Records screened\n(n = {flow_data.records_screened})", box_color
    )
    _add_box(
        ax,
        8,
        y,
        1.5,
        1.2,
        f"Records excluded\n(n = {flow_data.records_excluded})",
        excluded_color,
    )

    _add_arrow(ax, 5, 14.7, 5, 14.4)

    # Full-text assessment
    y = 11.5
    _add_box(
        ax,
        2.5,
        y,
        5,
        1.2,
        f"Full-text articles assessed\n(n = {flow_data.fulltext_assessed})",
        box_color,
    )

    # Exclusion reasons box
    exclusion_text = "Full-text articles excluded:\n"
    for reason, count in flow_data.exclusion_reasons.items():
        exclusion_text += f"  {reason}: {count}\n"
    exclusion_text += f"Total: {flow_data.fulltext_excluded}"

    _add_box(ax, 8, y, 1.5, 2, exclusion_text, excluded_color)

    _add_arrow(ax, 5, 12.7, 5, 12.4)

    # Included studies
    y = 9.5
    _add_box(
        ax,
        2.5,
        y,
        5,
        1.5,
        f"Studies included in review\n(n = {flow_data.studies_included})",
        "#E8F8E8",
    )

    # Section labels
    ax.text(0.5, 18, "Identification", ha="left", va="center", fontsize=12, fontweight="bold")
    ax.text(0.5, 13.5, "Screening", ha="left", va="center", fontsize=12, fontweight="bold")
    ax.text(0.5, 11.5, "Included", ha="left", va="center", fontsize=12, fontweight="bold")

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    logger.info(f"PRISMA diagram saved to {output_path}")


def _add_box(
    ax: plt.Axes, x: float, y: float, width: float, height: float, text: str, color: str
) -> None:
    """Add a box to the diagram."""
    box = FancyBboxPatch(
        (x, y - height / 2),
        width,
        height,
        boxstyle="round,pad=0.1",
        facecolor=color,
        edgecolor="black",
        linewidth=1.5,
    )
    ax.add_patch(box)
    ax.text(x + width / 2, y, text, ha="center", va="center", fontsize=9, wrap=True)


def _add_arrow(ax: plt.Axes, x1: float, y1: float, x2: float, y2: float) -> None:
    """Add an arrow between boxes."""
    arrow = FancyArrowPatch(
        (x1, y1),
        (x2, y2),
        arrowstyle="->",
        linewidth=2,
        color="black",
        mutation_scale=20,
    )
    ax.add_patch(arrow)


def generate_prisma_table(flow_data: PRISMAFlowData, output_path: Path) -> None:
    """
    Generate PRISMA flow counts table as CSV.

    Args:
        flow_data: PRISMA flow data
        output_path: Output CSV path
    """
    import csv

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Stage", "Count"])

        # Identification
        writer.writerow(["Records identified from databases", sum(flow_data.identification_databases.values())])
        for db, count in flow_data.identification_databases.items():
            writer.writerow([f"  - {db}", count])
        writer.writerow(["Records identified from other sources", flow_data.identification_other])

        # Deduplication
        writer.writerow(["Duplicates removed", flow_data.duplicates_removed])
        writer.writerow(["Records after deduplication", flow_data.records_after_dedup])

        # Screening
        writer.writerow(["Records screened (title/abstract)", flow_data.records_screened])
        writer.writerow(["Records excluded at screening", flow_data.records_excluded])

        # Full-text
        writer.writerow(["Full-text articles assessed", flow_data.fulltext_assessed])
        writer.writerow(["Full-text articles excluded", flow_data.fulltext_excluded])
        for reason, count in flow_data.exclusion_reasons.items():
            writer.writerow([f"  - {reason}", count])

        # Included
        writer.writerow(["Studies included in review", flow_data.studies_included])

    logger.info(f"PRISMA table saved to {output_path}")
