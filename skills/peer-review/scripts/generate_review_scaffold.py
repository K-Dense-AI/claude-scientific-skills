#!/usr/bin/env python3
"""Generate a local structured peer-review draft scaffold from validated intake."""

from __future__ import annotations

import argparse
import sys
import re
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _common import (  # noqa: E402
    ValidationError,
    error_exit,
    read_markdown,
    write_markdown,
)
from validate_review_intake import validate_intake  # noqa: E402

TEMPLATE_PATH = (
    Path(__file__).resolve().parents[1] / "assets" / "review_scaffold_template.md"
)
PLACEHOLDER_RE = re.compile(r"\{\{[A-Z0-9_]+\}\}")


def generate(payload: Any, template_path: Path = TEMPLATE_PATH) -> str:
    """Render only safe intake identifiers; never interpolate manuscript prose."""
    report = validate_intake(payload)
    if not report["valid"]:
        codes = sorted({item["code"] for item in report["errors"]})
        raise ValidationError(
            "intake is blocked; resolve these controls first: " + ", ".join(codes)
        )
    template = read_markdown(template_path)
    replacements = {
        "{{REVIEW_ID}}": report["review_id"],
        "{{REVIEWER_CAPACITY}}": report["normalized_scope"]["capacity"],
        "{{PEER_REVIEW_MODEL}}": report["normalized_scope"]["peer_review_model"],
        "{{AI_PLAN}}": report["normalized_scope"]["ai_plan"],
    }
    rendered = template
    for placeholder, value in replacements.items():
        rendered = rendered.replace(placeholder, value)
    unresolved = sorted(set(PLACEHOLDER_RE.findall(rendered)))
    if unresolved:
        raise ValidationError(
            "scaffold template contains unresolved placeholders: "
            + ", ".join(unresolved)
        )
    return rendered


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a private local Markdown review scaffold after intake "
            "authorization and conflict controls pass."
        )
    )
    parser.add_argument("intake", help="Local review intake JSON")
    parser.add_argument("-o", "--output", required=True, help="Output Markdown path")
    parser.add_argument(
        "--force", action="store_true", help="Replace an existing output file"
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        from _common import read_json

        rendered = generate(read_json(args.intake))
        destination = write_markdown(rendered, args.output, force=args.force)
        print(f"Created local review scaffold: {destination}")
        return 0
    except ValidationError as exc:
        return error_exit(exc)


if __name__ == "__main__":
    raise SystemExit(main())
