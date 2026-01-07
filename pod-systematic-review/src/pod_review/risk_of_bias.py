"""Risk of bias assessment tools (NOS and RoB2)."""

import csv
import logging
from pathlib import Path
from typing import Any

from .models import NOSAssessment, RoB2Assessment

logger = logging.getLogger(__name__)


def create_nos_template(output_path: Path) -> None:
    """
    Create CSV template for NOS assessment.

    Args:
        output_path: Output path for template
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "citation_id",
        "assessor_id",
        # Selection (max 4 stars)
        "representativeness",  # 0 or 1
        "selection_nonexposed",  # 0 or 1
        "ascertainment_exposure",  # 0 or 1
        "outcome_not_present",  # 0 or 1
        # Comparability (max 2 stars)
        "comparability",  # 0, 1, or 2
        # Outcome (max 3 stars)
        "assessment_outcome",  # 0 or 1
        "followup_duration",  # 0 or 1
        "followup_completeness",  # 0 or 1
        "notes",
    ]

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        # Add example row
        writer.writerow(
            {
                "citation_id": "EXAMPLE_2024",
                "assessor_id": "R1",
                "representativeness": "1",
                "selection_nonexposed": "1",
                "ascertainment_exposure": "1",
                "outcome_not_present": "1",
                "comparability": "2",
                "assessment_outcome": "1",
                "followup_duration": "1",
                "followup_completeness": "1",
                "notes": "Good quality study, 9/9 stars",
            }
        )

    logger.info(f"NOS template created: {output_path}")


def create_rob2_template(output_path: Path) -> None:
    """
    Create CSV template for RoB2 assessment.

    Args:
        output_path: Output path for template
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "citation_id",
        "assessor_id",
        "randomization",  # Low/Some concerns/High
        "deviations",
        "missing_data",
        "measurement",
        "selection_reported",
        "overall",  # Overall judgment
        "notes",
    ]

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        # Add example row
        writer.writerow(
            {
                "citation_id": "EXAMPLE_RCT_2024",
                "assessor_id": "R1",
                "randomization": "Low",
                "deviations": "Low",
                "missing_data": "Some concerns",
                "measurement": "Low",
                "selection_reported": "Low",
                "overall": "Some concerns",
                "notes": "Missing data for 8% of participants",
            }
        )

    logger.info(f"RoB2 template created: {output_path}")


def load_nos_assessments(csv_path: Path) -> list[NOSAssessment]:
    """
    Load NOS assessments from CSV.

    Args:
        csv_path: Path to CSV file

    Returns:
        List of NOS assessments
    """
    assessments = []

    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)

        for row in reader:
            # Skip example
            if row["citation_id"] == "EXAMPLE_2024":
                continue

            assessment = NOSAssessment(
                citation_id=row["citation_id"],
                assessor_id=row["assessor_id"],
                representativeness=int(row["representativeness"]),
                selection_nonexposed=int(row["selection_nonexposed"]),
                ascertainment_exposure=int(row["ascertainment_exposure"]),
                outcome_not_present=int(row["outcome_not_present"]),
                comparability=int(row["comparability"]),
                assessment_outcome=int(row["assessment_outcome"]),
                followup_duration=int(row["followup_duration"]),
                followup_completeness=int(row["followup_completeness"]),
                notes=row.get("notes", ""),
            )
            assessments.append(assessment)

    logger.info(f"Loaded {len(assessments)} NOS assessments from {csv_path}")
    return assessments


def load_rob2_assessments(csv_path: Path) -> list[RoB2Assessment]:
    """
    Load RoB2 assessments from CSV.

    Args:
        csv_path: Path to CSV file

    Returns:
        List of RoB2 assessments
    """
    assessments = []

    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)

        for row in reader:
            # Skip example
            if row["citation_id"].startswith("EXAMPLE"):
                continue

            assessment = RoB2Assessment(
                citation_id=row["citation_id"],
                assessor_id=row["assessor_id"],
                randomization=row["randomization"],
                deviations=row["deviations"],
                missing_data=row["missing_data"],
                measurement=row["measurement"],
                selection_reported=row["selection_reported"],
                overall=row["overall"],
                notes=row.get("notes", ""),
            )
            assessments.append(assessment)

    logger.info(f"Loaded {len(assessments)} RoB2 assessments from {csv_path}")
    return assessments


def generate_nos_summary_table(assessments: list[NOSAssessment], output_path: Path) -> None:
    """
    Generate NOS summary table.

    Args:
        assessments: List of NOS assessments
        output_path: Output CSV path
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        fieldnames = ["Citation ID", "Total Stars", "Quality Rating", "Selection", "Comparability", "Outcome"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for assessment in assessments:
            selection_stars = (
                assessment.representativeness
                + assessment.selection_nonexposed
                + assessment.ascertainment_exposure
                + assessment.outcome_not_present
            )
            outcome_stars = (
                assessment.assessment_outcome
                + assessment.followup_duration
                + assessment.followup_completeness
            )

            writer.writerow(
                {
                    "Citation ID": assessment.citation_id,
                    "Total Stars": assessment.total_stars,
                    "Quality Rating": assessment.quality_rating,
                    "Selection": f"{selection_stars}/4",
                    "Comparability": f"{assessment.comparability}/2",
                    "Outcome": f"{outcome_stars}/3",
                }
            )

    logger.info(f"NOS summary table saved to {output_path}")


def generate_rob2_summary_table(assessments: list[RoB2Assessment], output_path: Path) -> None:
    """
    Generate RoB2 summary table.

    Args:
        assessments: List of RoB2 assessments
        output_path: Output CSV path
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "Citation ID",
            "Randomization",
            "Deviations",
            "Missing Data",
            "Measurement",
            "Selection Reported",
            "Overall",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for assessment in assessments:
            writer.writerow(
                {
                    "Citation ID": assessment.citation_id,
                    "Randomization": assessment.randomization,
                    "Deviations": assessment.deviations,
                    "Missing Data": assessment.missing_data,
                    "Measurement": assessment.measurement,
                    "Selection Reported": assessment.selection_reported,
                    "Overall": assessment.overall,
                }
            )

    logger.info(f"RoB2 summary table saved to {output_path}")
