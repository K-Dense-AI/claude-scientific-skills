"""Vector .dna file path configuration.

Maps vector registry names to local SnapGene .dna file paths.
Used by the MCP server for in-silico cloning (recombinant plasmid generation).

Duet vectors with multiple MCS sites share the same physical .dna file.
"""

from __future__ import annotations

from pathlib import Path

# Base directory for common vectors
_KU = Path(r"C:\Users\Jahyun\OneDrive - 고려대학교\저장소")
_HOSEO = Path(
    r"C:\Users\Jahyun\OneDrive - 호서대학교"
    r"\바탕 화면 [Labtop]\USB backup\google drive\Data\Sequences\Vectors"
)

VECTOR_DNA_PATHS: dict[str, Path] = {
    "pET-21a(+)": _KU / "4. API Synthesis" / "pET-21a(+).dna",
    "pET-28a(+)": _HOSEO / "pET-28a(+).dna",
    "pMAL-c6T": _KU / "Sequences" / "pMAL-c6T.dna",
    "pETDuet-1:MCS1": _KU / "8. D-Gal to MA and D-tagatose" / "Genes"
                       / "Fusion protein" / "Module 2" / "pETDuet-1.dna",
    "pETDuet-1:MCS2": _KU / "8. D-Gal to MA and D-tagatose" / "Genes"
                       / "Fusion protein" / "Module 2" / "pETDuet-1.dna",
    "pACYCDuet-1:MCS1": _KU / "8. D-Gal to MA and D-tagatose" / "Genes"
                         / "pACYCDuet-1.dna",
    "pACYCDuet-1:MCS2": _KU / "8. D-Gal to MA and D-tagatose" / "Genes"
                         / "pACYCDuet-1.dna",
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
