#!/usr/bin/env python3
"""Lexically lint causal versus associational claims in bounded Markdown."""

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
    read_markdown,
    write_json_report,
)

TAG_RE = re.compile(
    r"\[(claim|estimand|identification|confounding|selection|collider|"
    r"reverse-causation):([A-Za-z0-9_.:-]+)\]",
    re.IGNORECASE,
)
CAUSAL_RE = re.compile(
    r"\b(?:causes?|caused|causal\s+(?:effect|impact)|effect\s+of|"
    r"leads?\s+to|results?\s+in|increases?|reduces?|prevents?|improves?|"
    r"worsens?|drives?|mediates?|produces?)\b",
    re.IGNORECASE,
)
ASSOCIATION_RE = re.compile(
    r"\b(?:associated\s+with|association|correlates?\s+with|correlation|"
    r"co-var(?:y|ies|ied)\s+with|predicts?|linked\s+to)\b",
    re.IGNORECASE,
)
CLAIM_TYPES = {"causal", "associational", "descriptive", "predictive", "mechanistic"}
IDENTIFICATION_TYPES = {
    "randomized",
    "quasi_experimental",
    "observational_assumption_dependent",
    "mechanistic_experiment",
    "other_assumption_dependent",
}
RISK_STATES = {"assessed", "unresolved", "not_applicable"}
REQUIRED_CAUSAL_TAGS = (
    "estimand",
    "identification",
    "confounding",
    "selection",
    "collider",
    "reverse-causation",
)


def _finding(code: str, line_number: int) -> dict[str, Any]:
    return {"code": code, "line": line_number}


def lint(markdown: str) -> dict[str, Any]:
    errors: list[dict[str, Any]] = []
    warnings: list[dict[str, Any]] = []
    claim_counts: Counter[str] = Counter()
    causal_trigger_lines = 0
    association_trigger_lines = 0
    in_fence = False

    for line_number, raw_line in enumerate(markdown.splitlines(), start=1):
        stripped = raw_line.strip()
        if stripped.startswith("```") or stripped.startswith("~~~"):
            in_fence = not in_fence
            continue
        if in_fence or not stripped or stripped.startswith("#"):
            continue

        matches = [(key.casefold(), value) for key, value in TAG_RE.findall(stripped)]
        key_counts = Counter(key for key, _ in matches)
        for key, count in key_counts.items():
            if count > 1:
                errors.append(_finding(f"DUPLICATE_{key.upper()}_TAG", line_number))
        tags = {key: value for key, value in matches}
        claim_type = tags.get("claim", "").casefold()
        has_causal_language = CAUSAL_RE.search(stripped) is not None
        has_association_language = ASSOCIATION_RE.search(stripped) is not None

        if claim_type:
            if claim_type not in CLAIM_TYPES:
                errors.append(_finding("INVALID_CLAIM_TYPE", line_number))
            else:
                claim_counts[claim_type] += 1

        if has_causal_language:
            causal_trigger_lines += 1
            if not claim_type:
                errors.append(_finding("UNMARKED_CAUSAL_LANGUAGE", line_number))
            elif claim_type not in {"causal", "mechanistic"}:
                errors.append(
                    _finding("CAUSAL_LANGUAGE_IN_NONCAUSAL_CLAIM", line_number)
                )

        if has_association_language:
            association_trigger_lines += 1
            if not claim_type:
                warnings.append(
                    _finding("UNMARKED_ASSOCIATIONAL_LANGUAGE", line_number)
                )
            elif claim_type == "causal" and not has_causal_language:
                warnings.append(
                    _finding("CAUSAL_TAG_WITH_ASSOCIATION_ONLY_LANGUAGE", line_number)
                )

        if claim_type == "causal":
            for required_tag in REQUIRED_CAUSAL_TAGS:
                if required_tag not in tags:
                    code_tag = required_tag.replace("-", "_").upper()
                    errors.append(
                        _finding(
                            f"CAUSAL_CLAIM_MISSING_{code_tag}_TAG",
                            line_number,
                        )
                    )
            identification = tags.get("identification", "").casefold()
            if identification and identification not in IDENTIFICATION_TYPES:
                errors.append(_finding("INVALID_IDENTIFICATION_TAG", line_number))
            for risk_tag in (
                "confounding",
                "selection",
                "collider",
                "reverse-causation",
            ):
                state = tags.get(risk_tag, "").casefold()
                code_tag = risk_tag.replace("-", "_").upper()
                if state and state not in RISK_STATES:
                    errors.append(
                        _finding(f"INVALID_{code_tag}_STATE", line_number)
                    )
                if state == "unresolved":
                    warnings.append(
                        _finding(f"UNRESOLVED_{code_tag}_RISK", line_number)
                    )
        elif claim_type and "estimand" in tags:
            warnings.append(_finding("ESTIMAND_TAG_ON_NONCAUSAL_CLAIM", line_number))

    if in_fence:
        errors.append(issue("UNCLOSED_MARKDOWN_FENCE", "document"))

    return {
        "schema_version": "2.0",
        "valid": not errors,
        "status": "INVALID_CLAIM_MARKUP" if errors else "VALID_LEXICAL_LINT",
        "errors": errors,
        "warnings": warnings,
        "claim_type_counts": dict(sorted(claim_counts.items())),
        "causal_trigger_line_count": causal_trigger_lines,
        "associational_trigger_line_count": association_trigger_lines,
        "notice": (
            "This deterministic lexical lint does not determine whether language "
            "is scientifically causal, whether an estimand is well defined, or "
            "whether identification assumptions hold. Review every flagged and "
            "unflagged claim manually."
        ),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Lint bounded local Markdown for causal/associational claim annotations "
            "and emit only rule codes and line numbers."
        )
    )
    parser.add_argument("document", help="Local Markdown document")
    parser.add_argument("-o", "--output", help="Optional local JSON report")
    parser.add_argument(
        "--force", action="store_true", help="Replace an existing output file"
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        report = lint(read_markdown(args.document))
        write_json_report(report, args.output, force=args.force)
        return 0 if report["valid"] else 1
    except ValidationError as exc:
        return error_exit(exc)


if __name__ == "__main__":
    raise SystemExit(main())
