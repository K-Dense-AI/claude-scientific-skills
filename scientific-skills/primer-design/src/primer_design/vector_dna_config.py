"""Vector .dna file path configuration.

Maps vector registry names to local SnapGene .dna file paths.
Used by the MCP server for in-silico cloning (recombinant plasmid generation).

Duet vectors with multiple MCS sites share the same physical .dna file.

Paths are resolved relative to the `vectors/` directory bundled with this
skill repo (primer-design/vectors/), so the skill works on any machine after
cloning without additional path configuration.
"""

from __future__ import annotations

from pathlib import Path

# vectors/ directory bundled with the skill repo
# __file__ -> .../primer-design/src/primer_design/vector_dna_config.py
_VECTORS = Path(__file__).parent.parent.parent / "vectors"

VECTOR_DNA_PATHS: dict[str, Path] = {
    "pET-21a(+)":       _VECTORS / "pET-21a(+).dna",
    "pET-28a(+)":       _VECTORS / "pET-28a(+).dna",
    "pMAL-c6T":         _VECTORS / "pMAL-c6T.dna",
    "pETDuet-1:MCS1":   _VECTORS / "pETDuet-1.dna",
    "pETDuet-1:MCS2":   _VECTORS / "pETDuet-1.dna",
    "pACYCDuet-1:MCS1": _VECTORS / "pACYCDuet-1.dna",
    "pACYCDuet-1:MCS2": _VECTORS / "pACYCDuet-1.dna",
}


def get_vector_dna_path(vector_name: str) -> Path | None:
    """Look up the .dna file path for a vector name.

    Supports exact match and fuzzy alias resolution.
    Returns None if no .dna file is configured or the file doesn't exist.
    """
    # Direct match
    path = VECTOR_DNA_PATHS.get(vector_name)
    if path and path.exists():
        return path

    # Fuzzy match (same normalization as vector_registry.get_vector)
    def _norm(s: str) -> str:
        return s.lower().replace(" ", "").replace("-", "").replace("(", "").replace(")", "").replace("+", "")

    target = _norm(vector_name)
    for canonical, p in VECTOR_DNA_PATHS.items():
        if _norm(canonical) == target and p.exists():
            return p

    return None
