"""Data extraction module for systematic review."""

import csv
import json
import logging
from pathlib import Path
from typing import Any

from .models import ExtractionRecord, RiskFactor, StudyDesign, EffectType

logger = logging.getLogger(__name__)


def create_extraction_template(output_path: Path) -> None:
    """
    Create CSV template for data extraction.

    Args:
        output_path: Output path for template
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "citation_id",
        "extractor_id",
        "study_design",
        "country",
        "setting",
        "recruitment_period",
        "sample_size",
        "mean_age",
        "sd_age",
        "age_range",
        "female_n",
        "female_percent",
        "surgery_type",
        "primary_revision",
        "pod_definition",
        "pod_assessment_tool",
        "pod_assessment_timing",
        "pod_incidence_n",
        "pod_incidence_percent",
        # Risk factors (repeatable section - instruction only)
        "risk_factor_name",
        "risk_factor_category",
        "comparison",
        "effect_type",
        "effect_size",
        "ci_lower",
        "ci_upper",
        "p_value",
        "adjustment",
        "covariates",
        "subgroup",
        "notes",
    ]

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        # Add example row
        writer.writerow(
            {
                "citation_id": "EXAMPLE_2024",
                "extractor_id": "R1",
                "study_design": "prospective_cohort",
                "country": "USA",
                "setting": "multicenter",
                "recruitment_period": "2020-2022",
                "sample_size": "450",
                "mean_age": "72.5",
                "sd_age": "5.2",
                "age_range": "65-89",
                "female_n": "270",
                "female_percent": "60",
                "surgery_type": "THA",
                "primary_revision": "primary",
                "pod_definition": "DSM-5",
                "pod_assessment_tool": "CAM",
                "pod_assessment_timing": "POD 1-7",
                "pod_incidence_n": "68",
                "pod_incidence_percent": "15.1",
                "risk_factor_name": "Age (per year)",
                "risk_factor_category": "demographic",
                "comparison": "per year increase",
                "effect_type": "OR",
                "effect_size": "1.08",
                "ci_lower": "1.03",
                "ci_upper": "1.14",
                "p_value": "0.002",
                "adjustment": "adjusted",
                "covariates": "sex; ASA score; surgery duration",
                "subgroup": "",
                "notes": "Adjusted for multiple covariates",
            }
        )

    logger.info(f"Extraction template created: {output_path}")


def create_extraction_schema() -> dict[str, Any]:
    """
    Create JSON schema for data extraction validation.

    Returns:
        JSON schema dictionary
    """
    schema = {
        "$schema": "http://json-schema.org/draft-07/schema#",
        "title": "ExtractionRecord",
        "type": "object",
        "required": [
            "citation_id",
            "extractor_id",
            "study_design",
            "pod_definition",
            "pod_assessment_tool",
        ],
        "properties": {
            "citation_id": {"type": "string"},
            "extractor_id": {"type": "string"},
            "study_design": {
                "type": "string",
                "enum": [
                    "prospective_cohort",
                    "retrospective_cohort",
                    "case_control",
                    "rct_secondary",
                ],
            },
            "country": {"type": "string"},
            "setting": {"type": "string"},
            "recruitment_period": {"type": "string"},
            "sample_size": {"type": ["integer", "null"]},
            "mean_age": {"type": ["number", "null"]},
            "sd_age": {"type": ["number", "null"]},
            "age_range": {"type": "string"},
            "female_n": {"type": ["integer", "null"]},
            "female_percent": {"type": ["number", "null"]},
            "surgery_type": {"type": "string"},
            "primary_revision": {"type": "string"},
            "pod_definition": {"type": "string"},
            "pod_assessment_tool": {"type": "string"},
            "pod_assessment_timing": {"type": "string"},
            "pod_incidence_n": {"type": ["integer", "null"]},
            "pod_incidence_percent": {"type": ["number", "null"]},
            "risk_factors": {
                "type": "array",
                "items": {
                    "type": "object",
                    "required": ["name", "effect_type", "effect_size"],
                    "properties": {
                        "name": {"type": "string"},
                        "category": {"type": "string"},
                        "comparison": {"type": "string"},
                        "effect_type": {"type": "string", "enum": ["OR", "RR", "HR"]},
                        "effect_size": {"type": "number"},
                        "ci_lower": {"type": ["number", "null"]},
                        "ci_upper": {"type": ["number", "null"]},
                        "p_value": {"type": ["number", "null"]},
                        "adjustment": {"type": "string"},
                        "covariates": {"type": "array", "items": {"type": "string"}},
                        "subgroup": {"type": "string"},
                    },
                },
            },
            "notes": {"type": "string"},
        },
    }
    return schema


def load_extraction_csv(csv_path: Path) -> list[ExtractionRecord]:
    """
    Load extraction records from CSV file.

    Args:
        csv_path: Path to CSV file

    Returns:
        List of ExtractionRecords
    """
    records: dict[str, ExtractionRecord] = {}

    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)

        for row in reader:
            citation_id = row["citation_id"]

            # Skip example row
            if citation_id == "EXAMPLE_2024":
                continue

            # Create or get existing record
            if citation_id not in records:
                records[citation_id] = ExtractionRecord(
                    citation_id=citation_id,
                    extractor_id=row["extractor_id"],
                    study_design=StudyDesign(row["study_design"]),
                    country=row.get("country", ""),
                    setting=row.get("setting", ""),
                    recruitment_period=row.get("recruitment_period", ""),
                    sample_size=int(row["sample_size"]) if row.get("sample_size") else None,
                    mean_age=float(row["mean_age"]) if row.get("mean_age") else None,
                    sd_age=float(row["sd_age"]) if row.get("sd_age") else None,
                    age_range=row.get("age_range", ""),
                    female_n=int(row["female_n"]) if row.get("female_n") else None,
                    female_percent=float(row["female_percent"])
                    if row.get("female_percent")
                    else None,
                    surgery_type=row.get("surgery_type", ""),
                    primary_revision=row.get("primary_revision", ""),
                    pod_definition=row.get("pod_definition", ""),
                    pod_assessment_tool=row.get("pod_assessment_tool", ""),
                    pod_assessment_timing=row.get("pod_assessment_timing", ""),
                    pod_incidence_n=int(row["pod_incidence_n"])
                    if row.get("pod_incidence_n")
                    else None,
                    pod_incidence_percent=float(row["pod_incidence_percent"])
                    if row.get("pod_incidence_percent")
                    else None,
                    notes=row.get("notes", ""),
                )

            # Add risk factor if present
            if row.get("risk_factor_name"):
                risk_factor = RiskFactor(
                    name=row["risk_factor_name"],
                    category=row.get("risk_factor_category", ""),
                    comparison=row.get("comparison", ""),
                    effect_type=EffectType(row["effect_type"]),
                    effect_size=float(row["effect_size"]),
                    ci_lower=float(row["ci_lower"]) if row.get("ci_lower") else None,
                    ci_upper=float(row["ci_upper"]) if row.get("ci_upper") else None,
                    p_value=float(row["p_value"]) if row.get("p_value") else None,
                    adjustment=row.get("adjustment", "crude"),
                    covariates=[c.strip() for c in row.get("covariates", "").split(";")]
                    if row.get("covariates")
                    else [],
                    subgroup=row.get("subgroup", ""),
                )
                records[citation_id].risk_factors.append(risk_factor)

    logger.info(f"Loaded {len(records)} extraction records from {csv_path}")
    return list(records.values())


def save_extraction_json(records: list[ExtractionRecord], output_path: Path) -> None:
    """
    Save extraction records to JSON file.

    Args:
        records: List of extraction records
        output_path: Output path
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    data = [record.model_dump() for record in records]

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, default=str)

    logger.info(f"Saved {len(records)} extraction records to {output_path}")


def validate_extraction_record(record: ExtractionRecord) -> list[str]:
    """
    Validate extraction record for completeness.

    Args:
        record: Extraction record

    Returns:
        List of validation errors (empty if valid)
    """
    errors = []

    # Required fields
    if not record.citation_id:
        errors.append("Missing citation_id")
    if not record.extractor_id:
        errors.append("Missing extractor_id")
    if not record.study_design:
        errors.append("Missing study_design")
    if not record.pod_definition:
        errors.append("Missing POD definition")
    if not record.pod_assessment_tool:
        errors.append("Missing POD assessment tool")

    # Sample size validation
    if record.sample_size is not None and record.sample_size <= 0:
        errors.append("Sample size must be positive")

    # Female percentage validation
    if record.female_percent is not None and not (0 <= record.female_percent <= 100):
        errors.append("Female percentage must be between 0 and 100")

    # Risk factors validation
    if not record.risk_factors:
        errors.append("No risk factors extracted")

    for idx, rf in enumerate(record.risk_factors):
        if not rf.name:
            errors.append(f"Risk factor {idx + 1}: Missing name")
        if rf.effect_size <= 0:
            errors.append(f"Risk factor {idx + 1}: Effect size must be positive")
        if rf.ci_lower is not None and rf.ci_upper is not None:
            if rf.ci_lower >= rf.ci_upper:
                errors.append(f"Risk factor {idx + 1}: CI lower must be < CI upper")

    return errors
