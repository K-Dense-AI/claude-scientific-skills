#!/usr/bin/env python3
"""Audit Markdown citation keys against a local reference CSV without network use."""

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
    read_markdown,
    require_doi,
    require_enum,
    require_identifier,
    require_text,
    require_unique,
    require_url,
    write_json_report,
)

REFERENCE_FIELDS = (
    "reference_id",
    "title",
    "authors",
    "year",
    "doi",
    "url",
    "verification_status",
)
VERIFICATION_STATUSES = {
    "verified_primary",
    "verified_secondary",
    "not_verified",
}
CITATION_GROUP_RE = re.compile(r"\[([^\]\n]{1,1000})\]")
CITATION_KEY_RE = re.compile(r"@([A-Za-z][A-Za-z0-9._:-]{0,95})")
YEAR_RE = re.compile(r"^(?:1[5-9]\d{2}|20\d{2}|2100)$")


def load_references(raw_path: str) -> list[dict[str, str]]:
    rows = read_csv_records(
        raw_path,
        required_fields=REFERENCE_FIELDS,
        allowed_fields=REFERENCE_FIELDS,
    )
    parsed: list[dict[str, str]] = []
    ids: list[str] = []
    for line_number, row in enumerate(rows, start=2):
        context = f"references row {line_number}"
        reference_id = require_identifier(
            row["reference_id"], f"{context}.reference_id"
        )
        ids.append(reference_id)
        year = require_text(row["year"], f"{context}.year", maximum=4)
        if not YEAR_RE.fullmatch(year):
            raise ValidationError(f"{context}.year must be a four-digit year")
        parsed.append(
            {
                "reference_id": reference_id,
                "title": require_text(
                    row["title"], f"{context}.title", minimum=3, maximum=2_000
                ),
                "authors": require_text(
                    row["authors"], f"{context}.authors", minimum=2, maximum=2_000
                ),
                "year": year,
                "doi": require_doi(
                    row["doi"], f"{context}.doi", allow_empty=True
                ),
                "url": require_url(
                    row["url"], f"{context}.url", allow_empty=True
                ),
                "verification_status": require_enum(
                    row["verification_status"],
                    VERIFICATION_STATUSES,
                    f"{context}.verification_status",
                ),
            }
        )
    require_unique(ids, "reference IDs")
    return parsed


def extract_citations(markdown: str) -> tuple[dict[str, list[int]], list[int]]:
    citations: dict[str, list[int]] = {}
    malformed_lines: list[int] = []
    for line_number, line in enumerate(markdown.splitlines(), start=1):
        matched_starts = 0
        for group in CITATION_GROUP_RE.finditer(line):
            content = group.group(1)
            if "@" not in content:
                continue
            matched_starts += content.count("@")
            keys = CITATION_KEY_RE.findall(content)
            if not keys:
                malformed_lines.append(line_number)
                continue
            for key in keys:
                citations.setdefault(key, []).append(line_number)
        if line.count("[@") > matched_starts:
            malformed_lines.append(line_number)
    return citations, sorted(set(malformed_lines))


def audit(markdown: str, references: list[dict[str, str]]) -> dict[str, Any]:
    citations, malformed_lines = extract_citations(markdown)
    by_id = {row["reference_id"]: row for row in references}
    cited_ids = set(citations)
    reference_ids = set(by_id)
    missing_ids = sorted(cited_ids - reference_ids)
    uncited_ids = sorted(reference_ids - cited_ids)
    unverified_ids = sorted(
        reference_id
        for reference_id in cited_ids.intersection(reference_ids)
        if by_id[reference_id]["verification_status"] == "not_verified"
    )
    no_locator_ids = sorted(
        row["reference_id"]
        for row in references
        if not row["doi"] and not row["url"]
    )
    errors: list[dict[str, str]] = [
        issue("CITATION_WITHOUT_REFERENCE", reference_id)
        for reference_id in missing_ids
    ]
    errors.extend(
        issue("MALFORMED_CITATION_SYNTAX", f"line:{line_number}")
        for line_number in malformed_lines
    )
    warnings: list[dict[str, str]] = [
        issue("REFERENCE_NOT_CITED", reference_id) for reference_id in uncited_ids
    ]
    warnings.extend(
        issue("CITED_REFERENCE_NOT_VERIFIED", reference_id)
        for reference_id in unverified_ids
    )
    warnings.extend(
        issue("REFERENCE_HAS_NO_PERSISTENT_LOCATOR", reference_id)
        for reference_id in no_locator_ids
    )
    verification_counts = Counter(
        row["verification_status"] for row in references
    )
    return {
        "schema_version": "2.0",
        "valid": not errors,
        "status": "VALID" if not errors else "CITATION_INCONSISTENCIES",
        "errors": errors,
        "warnings": warnings,
        "citation_occurrence_count": sum(len(lines) for lines in citations.values()),
        "cited_reference_count": len(cited_ids),
        "reference_count": len(references),
        "missing_reference_ids": missing_ids,
        "missing_reference_line_numbers": {
            reference_id: sorted(set(citations[reference_id]))
            for reference_id in missing_ids
        },
        "uncited_reference_ids": uncited_ids,
        "unverified_cited_reference_ids": unverified_ids,
        "references_without_persistent_locator": no_locator_ids,
        "verification_counts": dict(sorted(verification_counts.items())),
        "malformed_citation_line_numbers": malformed_lines,
        "notice": (
            "This local audit checks structured citation-key consistency and "
            "identifier format only. It does not query registries, verify that a "
            "source exists, confirm that a citation supports a claim, or echo "
            "manuscript/reference prose."
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Audit Pandoc-style Markdown citations such as [@ref-id] against a "
            "strict local reference CSV. No network calls are made."
        )
    )
    parser.add_argument("manuscript", help="Local Markdown file")
    parser.add_argument("references", help="Local reference CSV")
    parser.add_argument("-o", "--output", help="Optional local JSON report")
    parser.add_argument(
        "--force", action="store_true", help="Replace an existing output file"
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        report = audit(
            read_markdown(args.manuscript),
            load_references(args.references),
        )
        write_json_report(report, args.output, force=args.force)
        return 0 if report["valid"] else 1
    except ValidationError as exc:
        return error_exit(exc)


if __name__ == "__main__":
    raise SystemExit(main())
