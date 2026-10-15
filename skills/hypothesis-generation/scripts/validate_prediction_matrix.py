#!/usr/bin/env python3
"""Validate a prediction/rival-hypothesis CSV without ranking candidates."""

from __future__ import annotations

import argparse
import sys
import re
from collections import Counter
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _common import (  # noqa: E402
    ValidationError,
    error_exit,
    issue,
    read_csv_records,
    read_json,
    require_identifier,
    require_text,
    require_unique,
    split_identifiers,
    write_json_report,
)
from validate_hypothesis_schema import load_hypothesis_record, validate_record  # noqa: E402

FIELDS = (
    "prediction_id",
    "hypothesis_id",
    "rival_hypothesis_ids",
    "conditions",
    "observable",
    "expected_if_focal",
    "expected_if_rivals",
    "falsifier",
    "indeterminate_result",
    "boundary_conditions",
    "measurement_ids",
    "negative_control_ids",
    "analysis_ids",
    "uncertainty",
)


def load_matrix(raw_path: str) -> list[dict[str, Any]]:
    rows = read_csv_records(raw_path, fields=FIELDS)
    parsed: list[dict[str, Any]] = []
    prediction_ids: list[str] = []
    for line_number, row in enumerate(rows, start=2):
        context = f"matrix row {line_number}"
        prediction_id = require_identifier(
            row["prediction_id"], f"{context}.prediction_id"
        )
        prediction_ids.append(prediction_id)
        parsed.append(
            {
                "prediction_id": prediction_id,
                "hypothesis_id": require_identifier(
                    row["hypothesis_id"], f"{context}.hypothesis_id"
                ),
                "rival_hypothesis_ids": split_identifiers(
                    row["rival_hypothesis_ids"],
                    f"{context}.rival_hypothesis_ids",
                ),
                "conditions": require_text(
                    row["conditions"], f"{context}.conditions", minimum=5
                ),
                "observable": require_text(
                    row["observable"], f"{context}.observable", minimum=5
                ),
                "expected_if_focal": require_text(
                    row["expected_if_focal"],
                    f"{context}.expected_if_focal",
                    minimum=5,
                ),
                "expected_if_rivals": require_text(
                    row["expected_if_rivals"],
                    f"{context}.expected_if_rivals",
                    minimum=5,
                ),
                "falsifier": require_text(
                    row["falsifier"], f"{context}.falsifier", minimum=10
                ),
                "indeterminate_result": require_text(
                    row["indeterminate_result"],
                    f"{context}.indeterminate_result",
                    minimum=10,
                ),
                "boundary_conditions": require_text(
                    row["boundary_conditions"],
                    f"{context}.boundary_conditions",
                    minimum=5,
                ),
                "measurement_ids": split_identifiers(
                    row["measurement_ids"], f"{context}.measurement_ids"
                ),
                "negative_control_ids": split_identifiers(
                    row["negative_control_ids"],
                    f"{context}.negative_control_ids",
                ),
                "analysis_ids": split_identifiers(
                    row["analysis_ids"], f"{context}.analysis_ids"
                ),
                "uncertainty": require_text(
                    row["uncertainty"], f"{context}.uncertainty", minimum=10
                ),
            }
        )
    require_unique(prediction_ids, "prediction matrix")
    return parsed


def _normalized_expectation(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", value.casefold()).strip()


def _unknown(
    values: list[str], allowed: set[str], code: str, field: str
) -> list[dict[str, str]]:
    return [
        issue(code, f"{field}:{value}")
        for value in sorted(set(values) - allowed)
    ]


def validate_matrix(
    rows: list[dict[str, Any]], record: dict[str, Any] | None = None
) -> dict[str, Any]:
    errors: list[dict[str, str]] = []
    warnings: list[dict[str, str]] = []
    hypothesis_counts: Counter[str] = Counter()

    record_ids: dict[str, set[str]] | None = None
    if record is not None:
        record_report = validate_record(record)
        if not record_report["valid"]:
            raise ValidationError(
                "optional hypothesis record must pass schema validation first"
            )
        record_ids = {
            "hypotheses": {
                item["hypothesis_id"] for item in record["hypotheses"]
            },
            "predictions": {
                item["prediction_id"] for item in record["predictions"]
            },
            "measurements": {
                item["measurement_id"] for item in record["operationalizations"]
            },
            "controls": {
                item["control_id"] for item in record["negative_controls"]
            },
            "analyses": {
                item["analysis_id"] for item in record["analysis_plan"]["analyses"]
            },
        }

    for row in rows:
        prediction_id = row["prediction_id"]
        hypothesis_id = row["hypothesis_id"]
        hypothesis_counts[hypothesis_id] += 1
        if hypothesis_id in row["rival_hypothesis_ids"]:
            errors.append(issue("FOCAL_HYPOTHESIS_LISTED_AS_RIVAL", prediction_id))
        if (
            _normalized_expectation(row["expected_if_focal"])
            == _normalized_expectation(row["expected_if_rivals"])
        ):
            errors.append(issue("FOCAL_AND_RIVAL_EXPECTATIONS_IDENTICAL", prediction_id))
        if not row["negative_control_ids"]:
            warnings.append(issue("NEGATIVE_CONTROL_LINK_MISSING", prediction_id))

        if record_ids is not None:
            errors.extend(
                _unknown(
                    [hypothesis_id],
                    record_ids["hypotheses"],
                    "UNKNOWN_HYPOTHESIS_ID",
                    prediction_id,
                )
            )
            errors.extend(
                _unknown(
                    row["rival_hypothesis_ids"],
                    record_ids["hypotheses"],
                    "UNKNOWN_RIVAL_HYPOTHESIS_ID",
                    prediction_id,
                )
            )
            errors.extend(
                _unknown(
                    [prediction_id],
                    record_ids["predictions"],
                    "UNKNOWN_PREDICTION_ID",
                    prediction_id,
                )
            )
            errors.extend(
                _unknown(
                    row["measurement_ids"],
                    record_ids["measurements"],
                    "UNKNOWN_MEASUREMENT_ID",
                    prediction_id,
                )
            )
            errors.extend(
                _unknown(
                    row["negative_control_ids"],
                    record_ids["controls"],
                    "UNKNOWN_CONTROL_ID",
                    prediction_id,
                )
            )
            errors.extend(
                _unknown(
                    row["analysis_ids"],
                    record_ids["analyses"],
                    "UNKNOWN_ANALYSIS_ID",
                    prediction_id,
                )
            )
            record_prediction = next(
                (
                    item
                    for item in record["predictions"]
                    if item["prediction_id"] == prediction_id
                ),
                None,
            )
            if (
                record_prediction is not None
                and record_prediction["hypothesis_id"] != hypothesis_id
            ):
                errors.append(
                    issue("PREDICTION_HYPOTHESIS_MISMATCH", prediction_id)
                )

    return {
        "schema_version": "2.0",
        "valid": not errors,
        "status": "INVALID_MATRIX" if errors else "VALID_FOR_HUMAN_REVIEW",
        "errors": errors,
        "warnings": warnings,
        "prediction_count": len(rows),
        "hypothesis_prediction_counts": dict(sorted(hypothesis_counts.items())),
        "prediction_ids": sorted(row["prediction_id"] for row in rows),
        "record_cross_check_performed": record is not None,
        "notice": (
            "This report checks CSV structure, identifiers, declared contrasts, "
            "and optional cross-links. It cannot determine whether a prediction "
            "is scientifically discriminating, sufficiently precise, or likely, "
            "and it never ranks or selects a hypothesis."
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate a bounded local prediction/rival CSV and emit identifiers "
            "and rule codes without scientific scoring."
        )
    )
    parser.add_argument("matrix", help="Local prediction/rival matrix CSV")
    parser.add_argument(
        "--record", help="Optional local hypothesis record JSON for cross-checks"
    )
    parser.add_argument("-o", "--output", help="Optional local JSON report")
    parser.add_argument(
        "--force", action="store_true", help="Replace an existing output file"
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        record = (
            load_hypothesis_record(read_json(args.record)) if args.record else None
        )
        report = validate_matrix(load_matrix(args.matrix), record)
        write_json_report(report, args.output, force=args.force)
        return 0 if report["valid"] else 1
    except ValidationError as exc:
        return error_exit(exc)


if __name__ == "__main__":
    raise SystemExit(main())
