#!/usr/bin/env python3
"""
Primer Design MCP Server
=========================
FastMCP 기반 로컬 MCP 서버.
RE cloning, colony PCR, expression analysis, order sheet 등
primer design 기능을 Claude Code에서 tool로 직접 호출.

IMPORTANT: stdio transport 사용 -stdout으로의 print() 절대 금지.
           모든 로깅은 stderr로만 출력.
"""

from __future__ import annotations

import logging
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

from mcp.server.fastmcp import FastMCP

# stderr-only logging (stdout은 JSON-RPC 메시지 전용)
logging.basicConfig(
    stream=sys.stderr,
    level=logging.INFO,
    format="%(asctime)s [primer-design] %(levelname)s: %(message)s",
)
logger = logging.getLogger("primer-design")

# ── Imports from primer_design package ────────────────────────────────────────

from .restriction_cloning_mode import RestrictionCloningDesigner
from .colony_pcr_mode import ColonyPCRDesigner
from .expression_analyzer import ExpressionAnalyzer
from .order_sheet import PrimerOrderSheet
from .vector_registry import (
    EXPRESSION_VECTORS,
    RESTRICTION_ENZYMES,
    check_reading_frame,
    format_frame_report,
)

# ── Singleton instances ───────────────────────────────────────────────────────

_re_designer = RestrictionCloningDesigner()
_colony_designer = ColonyPCRDesigner()
_expression_analyzer = ExpressionAnalyzer()

# ── FastMCP Server ────────────────────────────────────────────────────────────

mcp = FastMCP("primer-design")


# ── Tool 1: design_re_cloning_primers ─────────────────────────────────────────

@mcp.tool()
def design_re_cloning_primers(
    insert_seq: str,
    re_5prime: str,
    re_3prime: str,
    vector_name: str | None = None,
    include_start_codon: bool = True,
    include_stop_codon: bool = False,
    target_tm: float = 62.0,
    gene_name: str = "Insert",
    gene_description: str = "",
    output_dir: str | None = None,
    generate_report_png: bool = False,
) -> dict:
    """Design RE cloning primers for inserting a CDS into an expression vector.

    Automatically generates a SnapGene .dna construct file (output_dir or CWD).
    Optionally generates a PNG report when generate_report_png=True.

    Args:
        insert_seq: Insert CDS sequence (DNA, ATG to stop). Spaces are stripped.
        re_5prime: 5' restriction enzyme name (e.g. "BamHI", "NdeI", "NheI").
        re_3prime: 3' restriction enzyme name (e.g. "XhoI", "NotI", "HindIII").
        vector_name: Optional vector name for reading frame check (e.g. "pET-28a(+)").
        include_start_codon: Include ATG start codon in forward primer (default True).
        include_stop_codon: Include stop codon in reverse primer (default False).
        target_tm: Target annealing Tm in Celsius (default 62.0).
        gene_name: Gene name for labeling report and SnapGene features (default "Insert").
        gene_description: Gene product description (e.g. "beta-galactosidase").
            Added as /note qualifier on CDS in SnapGene .dna file.
        output_dir: Output directory for generated files.
            Defaults to current working directory.
        generate_report_png: Generate a visual PNG report (default False).

    Returns:
        dict with F/R primer sequences, Tm, GC%, QC results, frame check, warnings,
        snapgene_path, and optionally report_image_path.
    """
    clean_seq = insert_seq.replace(" ", "")
    logger.info(
        "design_re_cloning_primers: %s/%s, vector=%s, insert=%d bp",
        re_5prime, re_3prime, vector_name, len(clean_seq),
    )
    result = _re_designer.design(
        insert_seq=insert_seq,
        re_5prime=re_5prime,
        re_3prime=re_3prime,
        vector_name=vector_name,
        target_tm=target_tm,
        include_start_codon=include_start_codon,
        include_stop_codon=include_stop_codon,
    )

    # ── Output directory (always generate files) ──────────────────────────
    out = Path(output_dir) if output_dir else Path.cwd()
    out.mkdir(parents=True, exist_ok=True)

    safe_name = gene_name.replace(" ", "_")
    if vector_name:
        safe_vec = (vector_name
                    .replace("(", "").replace(")", "")
                    .replace(":", "-").replace(" ", "_"))
        base = f"{safe_name}_{safe_vec}_{re_5prime}_{re_3prime}"
    else:
        base = f"{safe_name}_{re_5prime}_{re_3prime}"

    # ── Colony PCR primer suggestion (for .dna annotation) ──────────────
    colony_pcr_result = None
    if vector_name:
        try:
            colony_pcr_result = _colony_designer.suggest(
                vector_name=vector_name,
                insert_length_bp=len(clean_seq),
            )
            result["colony_pcr"] = colony_pcr_result
            logger.info(
                "Colony PCR: %s / %s (band: %d bp)",
                colony_pcr_result["f_name"],
                colony_pcr_result["r_name"],
                colony_pcr_result["expected_band_with_insert"],
            )
        except Exception as exc:
            logger.warning("Colony PCR suggest skipped: %s", exc)

    # ── Always generate SnapGene .dna construct file ──────────────────────
    try:
        from .snapgene_writer import write_cloning_construct
        from .vector_dna_config import get_vector_dna_path

        vec_dna = get_vector_dna_path(vector_name) if vector_name else None
        if vec_dna:
            logger.info("Using base vector: %s", vec_dna)
        else:
            logger.info("No vector .dna file found; generating PCR product only")

        dna_path = write_cloning_construct(
            design_result=result,
            insert_seq=clean_seq,
            output_path=out / f"{base}_construct.dna",
            gene_name=gene_name,
            vector_dna_path=vec_dna,
            gene_description=gene_description,
            colony_pcr=colony_pcr_result,
        )
        result["snapgene_path"] = str(dna_path)
        logger.info("SnapGene file: %s", dna_path)
    except Exception as exc:
        logger.error("SnapGene generation failed: %s", exc)
        result["snapgene_error"] = str(exc)

    # ── Optional PNG report ───────────────────────────────────────────────
    if generate_report_png:
        try:
            from .cloning_report import generate_cloning_report

            report_path = generate_cloning_report(
                design_result=result,
                output_path=out / f"{base}_report.png",
                gene_name=gene_name,
            )
            result["report_image_path"] = str(report_path)
            logger.info("Report image: %s", report_path)
        except Exception as exc:
            logger.error("Report generation failed: %s", exc)
            result["report_image_error"] = str(exc)

    # ── Expression viability check ──────────────────────────────────────
    try:
        expr_check = _check_expression_viability(
            insert_seq=clean_seq,
            gene_name=gene_name,
            re_5prime=re_5prime,
            re_3prime=re_3prime,
            design_result=result,
        )
        result["expression_check"] = expr_check
        logger.info("Expression check: %s", expr_check.get("verdict", "?"))
    except Exception as exc:
        logger.error("Expression check failed: %s", exc)
        result["expression_check_error"] = str(exc)

    return result


def _check_expression_viability(
    insert_seq: str,
    gene_name: str,
    re_5prime: str,
    re_3prime: str,
    design_result: dict,
) -> dict:
    """Comprehensive expression viability check for the cloned construct.

    Analyzes the full vector context including promoter, RBS, start codon,
    fusion protein topology, and insert quality.

    Checks:
    1. Vector expression context (promoter, RBS, start codon, direction)
    2. Fusion protein topology and reading frame
    3. Full expressed protein sequence and MW
    4. Internal RE sites in the insert
    5. Premature stop codons within the CDS
    6. Insert CDS quality (CAI, rare codons, strain recommendation)
    """
    from Bio.Seq import Seq

    warnings: list[str] = []
    seq = insert_seq.upper()
    frame_check = design_result.get("frame_check", {})
    vector_name = frame_check.get("vector_name", "")

    # ── 1. Vector expression context ──────────────────────────────────
    # Determine promoter, RBS, expression system from vector name
    vector_context = _get_vector_expression_context(vector_name)

    # ── 2. Fusion protein topology from frame_check ───────────────────
    topology = frame_check.get("topology", "")
    in_frame_5 = frame_check.get("in_frame_5prime", None)
    in_frame_3 = frame_check.get("in_frame_3prime", None)
    linker_5_aa = frame_check.get("linker_5prime_aa", "")
    linker_3_aa = frame_check.get("linker_3prime_aa", "")

    # Determine if insert's native stop codon terminates translation
    insert_has_native_stop = seq[-3:] in ("TAA", "TAG", "TGA")
    include_stop = design_result.get("include_stop_codon", False)
    has_stop = insert_has_native_stop or include_stop

    # ── 3. Build full expressed protein ───────────────────────────────
    # Translate insert (exclude terminal stop from protein)
    insert_protein = ""
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa = {"TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L",
              "CTA":"L","CTG":"L","ATT":"I","ATC":"I","ATA":"I","ATG":"M",
              "GTT":"V","GTC":"V","GTA":"V","GTG":"V","TCT":"S","TCC":"S",
              "TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
              "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A",
              "GCA":"A","GCG":"A","TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
              "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","AAT":"N","AAC":"N",
              "AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
              "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R",
              "CGA":"R","CGG":"R","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
              "GGT":"G","GGC":"G","GGA":"G","GGG":"G"}.get(codon, "X")
        if aa == "*":
            break
        insert_protein += aa

    # Full fusion protein
    if linker_5_aa and in_frame_5:
        fusion_protein = linker_5_aa + insert_protein
    else:
        fusion_protein = insert_protein

    if not has_stop and linker_3_aa and in_frame_3:
        fusion_protein += linker_3_aa

    # Fusion MW calculation
    avg_aa_mw = 110.0
    water = 18.015
    fusion_mw_kda = (
        len(fusion_protein) * avg_aa_mw
        - (len(fusion_protein) - 1) * water + water
    ) / 1000.0 if fusion_protein else 0.0

    # ── 4. Expression direction and start codon verification ──────────
    # Check if vector ATG or insert ATG is used
    start_codon_source = "unknown"
    if linker_5_aa:
        # Vector provides N-terminal region including ATG
        start_codon_source = "vector"
    elif seq[:3] == "ATG":
        start_codon_source = "insert"
    else:
        warnings.append(
            "No ATG start codon found at insert start, "
            "and no vector N-terminal fusion detected"
        )

    # Reading frame warnings
    if in_frame_5 is False:
        warnings.append(
            f"5' reading frame MISMATCH: insert is NOT in-frame with "
            f"vector N-terminal tags through {re_5prime}"
        )
    if not has_stop and in_frame_3 is False:
        warnings.append(
            f"3' reading frame MISMATCH: insert is NOT in-frame with "
            f"vector C-terminal tags through {re_3prime}"
        )

    # Active tags
    active_tags = []
    if linker_5_aa and in_frame_5:
        # Parse N-terminal tags from topology
        for tag in ["N-His6", "N-His8", "Thrombin", "T7-tag", "MBP",
                     "Factor Xa", "S-tag", "TEV"]:
            if tag in topology:
                active_tags.append({"name": tag, "position": "N-terminal"})
    if not has_stop and linker_3_aa and in_frame_3:
        for tag in ["C-His6", "C-His8", "S-tag"]:
            if tag in topology:
                active_tags.append({"name": tag, "position": "C-terminal"})

    # C-terminal tag blocked by stop codon
    if has_stop and "C-His" in topology:
        c_tag_blocked = True
    else:
        c_tag_blocked = False

    # ── 5. Internal RE site check ─────────────────────────────────────
    re5_site = design_result.get("re_5prime_site", "")
    re3_site = design_result.get("re_3prime_site", "")

    internal_cuts = []
    for re_name, re_site in [(re_5prime, re5_site), (re_3prime, re3_site)]:
        if not re_site:
            continue
        for check_site, label in [
            (re_site, re_site),
            (str(Seq(re_site).reverse_complement()), f"RC({re_site})"),
        ]:
            if check_site == re_site or check_site != re_site:
                count = seq.count(check_site)
                if count > 0:
                    positions = []
                    start = 0
                    while True:
                        idx = seq.find(check_site, start)
                        if idx == -1:
                            break
                        positions.append(idx + 1)
                        start = idx + 1
                    internal_cuts.append({
                        "enzyme": re_name,
                        "site": label,
                        "count": count,
                        "positions_bp": positions,
                    })
                    warnings.append(
                        f"CRITICAL: {re_name} ({label}) cuts within the insert "
                        f"at {count} position(s): bp {positions}"
                    )

    # ── 6. Premature stop codon check ─────────────────────────────────
    stop_codons = {"TAA", "TAG", "TGA"}
    premature_stops = []
    for i in range(0, len(seq) - 5, 3):  # exclude last codon
        codon = seq[i:i + 3]
        if codon in stop_codons:
            premature_stops.append({
                "codon": codon,
                "position_codon": i // 3 + 1,
                "position_bp": i + 1,
            })
    if premature_stops:
        warnings.append(
            f"Premature stop codon(s) at codon position(s): "
            f"{[s['position_codon'] for s in premature_stops]}. "
            f"Protein will be truncated!"
        )

    # ── 7. CDS length check ───────────────────────────────────────────
    if len(seq) % 3 != 0:
        warnings.append(
            f"Insert length ({len(seq)} bp) is not a multiple of 3. "
            f"Reading frame will be disrupted!"
        )

    # ── 8. Insert expression quality ──────────────────────────────────
    expr_result = _expression_analyzer.analyze(cds_seq=seq)

    # ── 9. Construct verdict ──────────────────────────────────────────
    has_critical = (
        bool(internal_cuts)
        or bool(premature_stops)
        or (len(seq) % 3 != 0)
        or in_frame_5 is False
    )

    if has_critical:
        verdict = "FAIL"
        verdict_detail = "Critical issues found - construct will NOT express correctly"
    elif expr_result.get("cai", 1.0) < 0.2:
        verdict = "WARNING"
        verdict_detail = "Very low CAI - expression may fail"
    elif len(warnings) > 0:
        verdict = "WARNING"
        verdict_detail = "Minor issues found - review before proceeding"
    else:
        verdict = "PASS"
        verdict_detail = "No issues detected - construct should express correctly"

    return {
        "verdict": verdict,
        "verdict_detail": verdict_detail,
        # Vector context
        "vector_context": vector_context,
        "start_codon_source": start_codon_source,
        # Fusion protein
        "topology": topology,
        "active_tags": active_tags,
        "c_terminal_tag_blocked": c_tag_blocked,
        "fusion_protein_length_aa": len(fusion_protein),
        "fusion_protein_mw_kda": round(fusion_mw_kda, 2),
        "insert_protein_length_aa": len(insert_protein),
        "n_terminal_linker_aa": linker_5_aa,
        "c_terminal_linker_aa": linker_3_aa if not has_stop else "(blocked by stop codon)",
        # Reading frame
        "in_frame_5prime": in_frame_5,
        "in_frame_3prime": in_frame_3,
        "insert_has_native_stop": insert_has_native_stop,
        # Insert quality
        "insert_cds_length_bp": len(seq),
        "cds_length_valid": len(seq) % 3 == 0,
        "gc_content": expr_result["basic_info"]["gc_content"],
        "cai": expr_result["cai"],
        "rare_codon_pct": expr_result["rare_codons"]["rare_codon_pct"],
        "rare_codon_clusters": len(expr_result.get("rare_codon_clusters", [])),
        "strain_recommendation": expr_result["strain_recommendation"]["primary_strain"],
        "map_removal": expr_result["map_removal"]["will_be_removed"],
        # Issues
        "internal_re_sites": internal_cuts,
        "premature_stop_codons": premature_stops,
        "warnings": warnings,
    }


# ── Vector expression context lookup ─────────────────────────────────────────

_VECTOR_EXPRESSION_INFO: dict[str, dict] = {
    "pET-21a(+)": {
        "promoter": "T7/lac",
        "inducer": "IPTG",
        "rbs": "vector-provided (T7 gene 10 leader)",
        "host_requirement": "BL21(DE3) or DE3 lysogen",
        "expression_type": "cytoplasmic, high-level",
        "notes": "C-terminal His-tag only (no N-terminal tag)",
    },
    "pET-28a(+)": {
        "promoter": "T7/lac",
        "inducer": "IPTG",
        "rbs": "vector-provided (T7 gene 10 leader)",
        "host_requirement": "BL21(DE3) or DE3 lysogen",
        "expression_type": "cytoplasmic, high-level",
        "notes": "N-terminal His6 + thrombin + T7-tag; optional C-terminal His6",
    },
    "pMAL-c6T": {
        "promoter": "Ptac",
        "inducer": "IPTG",
        "rbs": "malE translation initiation signals",
        "host_requirement": "any E. coli strain",
        "expression_type": "cytoplasmic, MBP fusion",
        "notes": "MBP-TEV-insert fusion; amylose resin purification",
    },
    "pETDuet-1:MCS1": {
        "promoter": "T7/lac (MCS1)",
        "inducer": "IPTG",
        "rbs": "vector-provided",
        "host_requirement": "BL21(DE3) or DE3 lysogen",
        "expression_type": "cytoplasmic, co-expression",
        "notes": "Duet vector MCS1; N-terminal His6; co-expression with MCS2",
    },
    "pETDuet-1:MCS2": {
        "promoter": "T7/lac (MCS2)",
        "inducer": "IPTG",
        "rbs": "vector-provided",
        "host_requirement": "BL21(DE3) or DE3 lysogen",
        "expression_type": "cytoplasmic, co-expression",
        "notes": "Duet vector MCS2; optional S-tag; co-expression with MCS1",
    },
    "pACYCDuet-1:MCS1": {
        "promoter": "T7/lac (MCS1)",
        "inducer": "IPTG",
        "rbs": "vector-provided",
        "host_requirement": "BL21(DE3) or DE3 lysogen",
        "expression_type": "cytoplasmic, co-expression",
        "notes": "CmR; compatible with pET/pETDuet (ColA ori); N-terminal His6",
    },
    "pACYCDuet-1:MCS2": {
        "promoter": "T7/lac (MCS2)",
        "inducer": "IPTG",
        "rbs": "vector-provided",
        "host_requirement": "BL21(DE3) or DE3 lysogen",
        "expression_type": "cytoplasmic, co-expression",
        "notes": "CmR; compatible with pET/pETDuet (ColA ori); optional S-tag",
    },
}


def _get_vector_expression_context(vector_name: str) -> dict:
    """Return expression context info for a vector."""
    if not vector_name:
        return {"promoter": "unknown", "notes": "No vector specified"}

    # Direct match
    info = _VECTOR_EXPRESSION_INFO.get(vector_name)
    if info:
        return info

    # Fuzzy match
    vn = vector_name.lower().replace(" ", "").replace("-", "").replace("(", "").replace(")", "").replace("+", "")
    for name, data in _VECTOR_EXPRESSION_INFO.items():
        norm = name.lower().replace(" ", "").replace("-", "").replace("(", "").replace(")", "").replace("+", "")
        if vn == norm:
            return data

    return {"promoter": "unknown", "notes": f"Vector '{vector_name}' not in expression database"}


# ── Tool 2: recommend_re_pair ─────────────────────────────────────────────────

@mcp.tool()
def recommend_re_pair(
    insert_seq: str,
    vector_name: str,
    prefer_hf: bool = True,
) -> list[dict]:
    """Recommend optimal RE pairs for cloning an insert into a vector.

    Scores each RE combination based on buffer compatibility, reading frame,
    and HF variant availability.

    Args:
        insert_seq: Insert CDS sequence (DNA).
        vector_name: Vector name (e.g. "pET-28a(+)", "pET-21a(+)").
        prefer_hf: Prefer HF (High Fidelity) enzyme variants (default True).

    Returns:
        List of RE pair recommendations sorted by score (highest first).
        Each entry has re_5prime, re_3prime, score, buffer, frame_ok, reason.
    """
    logger.info(
        "recommend_re_pair: vector=%s, insert=%d bp",
        vector_name, len(insert_seq.replace(" ", "")),
    )
    return _re_designer.recommend_re_pair(
        insert_seq=insert_seq,
        vector_name=vector_name,
        prefer_hf=prefer_hf,
    )


# ── Tool 3: suggest_colony_pcr ────────────────────────────────────────────────

@mcp.tool()
def suggest_colony_pcr(
    vector_name: str,
    insert_length_bp: int,
) -> dict:
    """Suggest colony PCR primers and expected band sizes for a vector+insert.

    Uses universal primer database (Macrogen standard primers) matched to
    the vector type.

    Args:
        vector_name: Vector name (e.g. "pET-28a(+)", "pETDuet-1:MCS1").
            Fuzzy matching supported (e.g. "pET28a").
        insert_length_bp: Insert length in base pairs.

    Returns:
        dict with F/R primer name, sequence, Tm, expected band sizes
        (with insert and empty vector), and annealing temperature.
    """
    logger.info(
        "suggest_colony_pcr: vector=%s, insert=%d bp",
        vector_name, insert_length_bp,
    )
    return _colony_designer.suggest(
        vector_name=vector_name,
        insert_length_bp=insert_length_bp,
    )


# ── Tool 4: analyze_expression ────────────────────────────────────────────────

@mcp.tool()
def analyze_expression(cds_seq: str) -> dict:
    """Analyze a CDS sequence for E. coli expression optimization.

    Performs comprehensive analysis including:
    - Protein info (length, MW, GC content)
    - Rare codon frequency and cluster detection (E. coli K12)
    - CAI (Codon Adaptation Index)
    - Signal peptide prediction (rule-based)
    - MAP (Met aminopeptidase) removal prediction
    - Expression strain recommendation (BL21/Rosetta/Rosetta 2)

    Args:
        cds_seq: CDS DNA sequence (ATG to stop codon, inclusive).

    Returns:
        dict with basic_info, rare_codons, cai, signal_peptide,
        strain_recommendation, and warnings.
    """
    logger.info(
        "analyze_expression: cds=%d bp",
        len(cds_seq.replace(" ", "").replace("\n", "")),
    )
    return _expression_analyzer.analyze(cds_seq=cds_seq)


# ── Tool 5: check_reading_frame ───────────────────────────────────────────────

@mcp.tool()
def check_reading_frame_tool(
    vector_name: str,
    re_5prime: str,
    re_3prime: str,
    insert_has_atg: bool = True,
    insert_has_stop: bool = False,
    insert_cds_bp: int | None = None,
) -> dict:
    """Verify reading frame compatibility for an RE cloning strategy.

    Checks whether the insert will be in-frame with the vector's N-terminal
    and C-terminal tags after ligation.

    Args:
        vector_name: Vector name (e.g. "pET-28a(+)").
        re_5prime: 5' restriction enzyme name.
        re_3prime: 3' restriction enzyme name.
        insert_has_atg: Whether the insert has its own start codon (default True).
        insert_has_stop: Whether the insert has its own stop codon (default False).
        insert_cds_bp: Insert CDS length in bp (optional, for length validation).

    Returns:
        dict with in_frame_5prime, in_frame_3prime, topology diagram,
        linker sequences, and warnings. Also includes a human-readable
        frame_report string.
    """
    logger.info(
        "check_reading_frame: %s %s/%s",
        vector_name, re_5prime, re_3prime,
    )
    result = check_reading_frame(
        vector_name=vector_name,
        re_5prime=re_5prime,
        re_3prime=re_3prime,
        insert_has_atg=insert_has_atg,
        insert_has_stop=insert_has_stop,
        insert_cds_bp=insert_cds_bp,
    )
    result["frame_report"] = format_frame_report(result)
    return result


# ── Tool 6: list_vectors ──────────────────────────────────────────────────────

@mcp.tool()
def list_vectors() -> list[dict]:
    """List all available expression vectors and their RE sites.

    Returns:
        List of vectors, each with name, RE sites, tag info, and aliases.
    """
    vectors = []
    for name, data in EXPRESSION_VECTORS.items():
        tags = {}
        for tag_name, tag_info in data.get("tags", {}).items():
            tags[tag_name] = tag_info["type"]

        vectors.append({
            "name": name,
            "re_sites": list(data["re_sites"].keys()),
            "tags": tags,
            "aliases": data.get("aliases", []),
            "notes": data.get("notes", ""),
        })
    return vectors


# ── Tool 7: list_restriction_enzymes ──────────────────────────────────────────

@mcp.tool()
def list_restriction_enzymes() -> list[dict]:
    """List all available restriction enzymes with recognition sites and cut info.

    Returns:
        List of enzymes, each with name, recognition site, overhang type,
        overhang length, and cut positions.
    """
    enzymes = []
    for name, info in sorted(RESTRICTION_ENZYMES.items()):
        enzymes.append({
            "name": name,
            "recognition": info["recognition"],
            "overhang": info["overhang"],
            "overhang_len": info["overhang_len"],
            "cut_top": info["cut_top"],
            "cut_bottom": info["cut_bottom"],
        })
    return enzymes


# ── Tool 8: generate_macrogen_order ───────────────────────────────────────────

@mcp.tool()
def generate_macrogen_order(
    primers: list[dict],
    project_name: str = "primer_order",
    output_dir: str | None = None,
) -> dict:
    """Generate a Macrogen-compatible primer order sheet (XLSX).

    Creates an XLSX file with Macrogen Oligo Order format
    (No., Oligo Name, Sequence, Amount, Purification).

    Args:
        primers: List of primer dicts, each with "name" (str) and "sequence" (str).
            Example: [{"name": "gudD_BamHI_F", "sequence": "GCGCGGATCC..."}]
        project_name: Project name for the output file (default "primer_order").
        output_dir: Output directory path. Defaults to current working directory.

    Returns:
        dict with the generated file path and order summary.
    """
    logger.info(
        "generate_macrogen_order: %d primers, project=%s",
        len(primers), project_name,
    )
    sheet = PrimerOrderSheet(project_name=project_name)

    for p in primers:
        sheet.add_custom_primer(
            name=p["name"],
            sequence=p["sequence"],
        )

    out_dir = Path(output_dir) if output_dir else Path.cwd()
    xlsx_path = sheet.to_macrogen_oligo(output_path=None)

    # If output_dir specified, move to that directory
    if output_dir:
        out_dir.mkdir(parents=True, exist_ok=True)
        target = out_dir / xlsx_path.name
        xlsx_path.rename(target)
        xlsx_path = target

    summary = sheet.summary()

    return {
        "file_path": str(xlsx_path),
        "total_primers": summary["total_primers"],
        "total_length_nt": summary["total_length_nt"],
        "estimated_cost_krw": summary["estimated_cost_krw"],
    }


# ── E. coli K12 codon usage (best codon per amino acid) ──────────────────────

_ECOLI_BEST_CODON: dict[str, str] = {
    "F": "TTC", "L": "CTG", "I": "ATT", "M": "ATG", "V": "GTG",
    "S": "AGC", "P": "CCG", "T": "ACC", "A": "GCG", "Y": "TAC",
    "H": "CAC", "Q": "CAG", "N": "AAC", "K": "AAA", "D": "GAT",
    "E": "GAA", "C": "TGC", "W": "TGG", "R": "CGC", "G": "GGC",
    "*": "TAA",
}


def _codon_optimize_for_ecoli(protein_seq: str) -> str:
    """Back-translate protein to DNA using E. coli K12 optimal codons."""
    codons = []
    for aa in protein_seq.upper():
        best = _ECOLI_BEST_CODON.get(aa)
        if best is None:
            raise ValueError(f"Unknown amino acid: '{aa}'")
        codons.append(best)
    return "".join(codons)


def _search_ncbi_gene(gene_name: str, organism: str) -> int | None:
    """Search NCBI Gene database and return Gene ID."""
    from Bio import Entrez

    queries = [
        f'{gene_name}[Gene Name] AND "{organism}"[Organism]',
        f'{gene_name}[All Fields] AND "{organism}"[Organism]',
    ]
    for query in queries:
        logger.info("NCBI Gene search: %s", query)
        try:
            handle = Entrez.esearch(db="gene", term=query, retmax=5)
            record = Entrez.read(handle)
            handle.close()
        except Exception as exc:
            logger.warning("Entrez esearch failed: %s", exc)
            time.sleep(0.5)
            continue

        id_list = record.get("IdList", [])
        if id_list:
            gene_id = int(id_list[0])
            logger.info("Found Gene ID: %d (%d results)", gene_id, len(id_list))
            return gene_id
        time.sleep(0.4)
    return None


def _fetch_cds_from_gene_id(gene_id: int, gene_name_hint: str = "") -> dict | None:
    """Fetch CDS nucleotide sequence from NCBI Gene ID."""
    from Bio import Entrez, SeqIO

    # 1. Gene summary for metadata
    logger.info("Fetching Gene summary for ID %d", gene_id)
    official_name = ""
    description = ""
    try:
        handle = Entrez.efetch(db="gene", id=str(gene_id), rettype="xml")
        xml_data = handle.read()
        handle.close()
        root = ET.fromstring(xml_data)
        for el in root.iter("Gene-ref_locus"):
            official_name = el.text or ""
            break
        for el in root.iter("Gene-ref_desc"):
            description = el.text or ""
            break
    except Exception as exc:
        logger.warning("Gene XML parse: %s", exc)
    time.sleep(0.4)

    # 2. Link Gene → Nucleotide (RefSeq)
    nuc_ids: list[str] = []
    for link_name in ("gene_nuccore_refseqrna", "gene_nuccore_refseqgene", "gene_nuccore"):
        try:
            handle = Entrez.elink(dbfrom="gene", db="nuccore", id=str(gene_id), linkname=link_name)
            link_results = Entrez.read(handle)
            handle.close()
        except Exception as exc:
            logger.warning("Elink %s failed: %s", link_name, exc)
            time.sleep(0.4)
            continue

        for linkset in link_results:
            for link_db in linkset.get("LinkSetDb", []):
                for link in link_db.get("Link", []):
                    nuc_ids.append(link["Id"])
        if nuc_ids:
            logger.info("Found %d nucleotide IDs via %s", len(nuc_ids), link_name)
            break
        time.sleep(0.4)

    if not nuc_ids:
        logger.warning("No linked nucleotide records for Gene %d", gene_id)
        return None

    # 3. Fetch GenBank and extract CDS
    # Build target name set for case-insensitive matching
    target_names: set[str] = set()
    if official_name:
        target_names.add(official_name.lower())
    if gene_name_hint:
        target_names.add(gene_name_hint.lower())

    def _matches_gene(feat_gene: str, feat_locus: str) -> bool:
        """Check if CDS gene/locus_tag matches any target name."""
        if not target_names:
            return False
        return (feat_gene.lower() in target_names
                or feat_locus.lower() in target_names)

    fallback_cds = None  # gene name 매칭 실패 시 첫 CDS를 fallback으로 보관

    for nuc_id in nuc_ids[:10]:
        time.sleep(0.4)
        try:
            handle = Entrez.efetch(db="nuccore", id=nuc_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
        except Exception as exc:
            logger.warning("Nucleotide fetch %s failed: %s", nuc_id, exc)
            continue

        for feature in record.features:
            if feature.type != "CDS":
                continue
            feat_gene = feature.qualifiers.get("gene", [""])[0]
            feat_locus = feature.qualifiers.get("locus_tag", [""])[0]
            cds_seq = str(feature.location.extract(record.seq)).upper()

            # Accept alternative start codons (ATG, GTG, TTG) common in prokaryotes
            if cds_seq[:3] not in ("ATG", "GTG", "TTG") or len(cds_seq) % 3 != 0:
                continue

            protein = feature.qualifiers.get("translation", [None])[0]
            if protein is None:
                try:
                    protein = str(feature.location.extract(record.seq).translate(to_stop=True))
                except Exception:
                    continue

            product = feature.qualifiers.get("product", [""])[0]
            accession = record.id or record.name

            cds_info = {
                "cds_seq": cds_seq,
                "protein_seq": protein,
                "source": accession,
                "description": product or description,
                "gene_name_official": feat_gene or official_name,
                "locus_tag": feat_locus,
            }

            # Gene name 매칭: gene qualifier 또는 locus_tag 비교
            if _matches_gene(feat_gene, feat_locus):
                logger.info(
                    "Matched CDS by gene/locus '%s'/'%s': %d bp, %d aa",
                    feat_gene, feat_locus, len(cds_seq), len(protein),
                )
                return cds_info

            # Single-CDS record (mRNA/RefSeq) → gene name 무관하게 사용
            cds_count = sum(1 for f in record.features if f.type == "CDS")
            if cds_count == 1:
                logger.info(
                    "Single-CDS record: gene=%s, %d bp, %d aa",
                    feat_gene or feat_locus, len(cds_seq), len(protein),
                )
                return cds_info

            # Fallback: 첫 유효 CDS 저장 (gene name 매칭이 안 될 경우 대비)
            if fallback_cds is None:
                fallback_cds = cds_info

    # Gene name 매칭 실패 시 fallback 사용하지 않음 (잘못된 gene 반환 방지)
    if fallback_cds:
        logger.warning(
            "No CDS matched gene name(s) %s; discarding fallback '%s'",
            target_names, fallback_cds.get("gene_name_official"),
        )
    return None


# ── Tool 9: fetch_gene_sequence ───────────────────────────────────────────────

@mcp.tool()
def fetch_gene_sequence(
    gene_name: str,
    organism: str = "Escherichia coli",
    gene_id: int | None = None,
    codon_optimize: bool = False,
) -> dict:
    """Fetch a gene's CDS nucleotide sequence from NCBI for primer design.

    Searches the NCBI Gene database, retrieves the CDS (ATG to stop),
    and returns DNA/protein sequences ready for design_re_cloning_primers.

    Default: returns native genomic sequence (gDNA PCR용).
    codon_optimize=True일 때만 E. coli K12 최적 코돈으로 back-translate.

    Args:
        gene_name: Gene name or symbol (e.g. "gudD", "lacZ", "malE").
        organism: Organism name for NCBI search (default "Escherichia coli").
        gene_id: Optional NCBI Gene ID. If provided, skips search step.
        codon_optimize: Back-translate using E. coli K12 optimal codons
            (default False). Only enable when explicitly requested.

    Returns:
        dict with gene_name, organism, gene_id, cds_seq, protein_seq,
        cds_length_bp, protein_length_aa, source, description.
    """
    from Bio import Entrez
    Entrez.email = "jahyunlee00299@gmail.com"
    Entrez.tool = "primer-design-mcp"

    logger.info(
        "fetch_gene_sequence: gene=%s, organism=%s, gene_id=%s, optimize=%s",
        gene_name, organism, gene_id, codon_optimize,
    )

    # Step 1: Resolve Gene ID
    if gene_id is None:
        gene_id = _search_ncbi_gene(gene_name, organism)
        if gene_id is None:
            return {
                "error": f"Gene '{gene_name}' not found for '{organism}'. "
                         f"Try a different name or provide gene_id directly.",
                "gene_name": gene_name,
                "organism": organism,
            }

    # Step 2: Fetch CDS
    try:
        cds_result = _fetch_cds_from_gene_id(gene_id, gene_name_hint=gene_name)
    except Exception as exc:
        return {
            "error": f"Failed to fetch CDS for Gene ID {gene_id}: {exc}",
            "gene_name": gene_name, "gene_id": gene_id,
        }

    if cds_result is None:
        return {
            "error": f"No valid CDS for Gene ID {gene_id} ('{gene_name}'). "
                     f"Provide the sequence directly.",
            "gene_name": gene_name, "gene_id": gene_id,
        }

    cds_seq = cds_result["cds_seq"]
    protein_seq = cds_result["protein_seq"]

    # Step 3: Codon optimization (optional, only when explicitly requested)
    if codon_optimize:
        logger.info("Codon-optimizing %d aa for E. coli K12", len(protein_seq))
        try:
            optimized = _codon_optimize_for_ecoli(protein_seq)
            if cds_seq[-3:] in ("TAA", "TAG", "TGA"):
                optimized += _ECOLI_BEST_CODON["*"]
            cds_seq = optimized
        except Exception as exc:
            return {
                "error": f"Codon optimization failed: {exc}",
                "gene_name": gene_name, "gene_id": gene_id,
                "native_cds_seq": cds_result["cds_seq"],
                "protein_seq": protein_seq,
            }

    # Detect native start codon (ATG/GTG/TTG)
    native_start = cds_result["cds_seq"][:3]
    start_codon_note = ""
    if native_start != "ATG":
        start_codon_note = (
            f"Native start codon is {native_start} (not ATG). "
            f"For recombinant expression, the RE cloning primer will "
            f"replace it with ATG automatically."
        )

    result = {
        "gene_name": cds_result.get("gene_name_official") or gene_name,
        "organism": organism,
        "gene_id": gene_id,
        "cds_seq": cds_seq,
        "protein_seq": protein_seq,
        "cds_length_bp": len(cds_seq),
        "protein_length_aa": len(protein_seq),
        "source": cds_result["source"],
        "description": cds_result["description"],
        "codon_optimized": codon_optimize,
        "native_start_codon": native_start,
    }
    if start_codon_note:
        result["start_codon_note"] = start_codon_note
    if cds_result.get("locus_tag"):
        result["locus_tag"] = cds_result["locus_tag"]

    logger.info(
        "fetch_gene_sequence OK: %s, %d bp, %d aa, source=%s",
        result["gene_name"], len(cds_seq), len(protein_seq), result["source"],
    )
    return result


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    logger.info("Starting primer-design MCP server (stdio transport)")
    mcp.run()
