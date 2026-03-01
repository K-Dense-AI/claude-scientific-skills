#!/usr/bin/env python3
"""
Stress test: fetch random E. coli genes and run through
fetch_gene_sequence → design_re_cloning_primers pipeline.

Tests random genes to find bugs and edge cases.
Results are written to stress_test_results.json.

Usage:
  python stress_test_genes.py              # full test (1500 genes)
  python stress_test_genes.py --quick 50   # quick test (50 genes)
"""
from __future__ import annotations

import argparse
import json
import os
import random
import sys
import tempfile
import time
import traceback
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))


def log(msg: str):
    print(msg, flush=True)


# ── Step 1: Get E. coli K12 gene list from NCBI ─────────────────────────────

def fetch_ecoli_gene_list(n: int = 2000) -> list[dict]:
    """Fetch a list of E. coli K12 genes from NCBI Gene database."""
    from Bio import Entrez
    Entrez.email = "jahyunlee00299@gmail.com"
    Entrez.tool = "primer-design-stress-test"

    log(f"[STEP 1] Fetching E. coli K12 gene list (up to {n} genes)...")
    handle = Entrez.esearch(
        db="gene",
        term='"Escherichia coli str. K-12 substr. MG1655"[Organism] AND alive[property]',
        retmax=n,
        sort="relevance",
    )
    record = Entrez.read(handle)
    handle.close()
    id_list = record.get("IdList", [])
    log(f"  Found {len(id_list)} gene IDs")

    # Fetch gene summaries in batches
    genes = []
    batch_size = 200
    for i in range(0, len(id_list), batch_size):
        batch = id_list[i:i + batch_size]
        time.sleep(0.5)
        try:
            handle = Entrez.esummary(db="gene", id=",".join(batch))
            summaries = Entrez.read(handle)
            handle.close()
        except Exception as exc:
            log(f"  WARN: batch {i} esummary failed: {exc}")
            continue

        doc_sums = summaries.get("DocumentSummarySet", {}).get("DocumentSummary", [])
        for doc in doc_sums:
            name = doc.get("Name", "")
            desc = doc.get("Description", "")
            gene_id = doc.attributes.get("uid", "") if hasattr(doc, "attributes") else ""
            if name:
                genes.append({
                    "gene_id": int(gene_id) if gene_id else None,
                    "gene_name": name,
                    "description": desc,
                })
        log(f"  Fetched {len(genes)} gene summaries...")

    return genes


# ── Step 2: Test pipeline ────────────────────────────────────────────────────

VECTOR_COMBOS = [
    ("NdeI", "XhoI", "pET-28a(+)"),
    ("NheI", "HindIII", "pET-21a(+)"),
    ("BamHI", "HindIII", "pMAL-c6T"),
    ("NdeI", "BamHI", "pET-28a(+)"),
    ("NcoI", "XhoI", "pET-28a(+)"),
    ("BamHI", "XhoI", "pET-28a(+)"),
]


def test_single_gene(gene_info: dict, test_idx: int, total: int) -> dict:
    """Test fetch + design pipeline for a single gene."""
    from primer_design.mcp_server import fetch_gene_sequence, design_re_cloning_primers

    gene_name = gene_info["gene_name"]
    gene_id = gene_info.get("gene_id")
    result = {
        "gene_name": gene_name,
        "gene_id": gene_id,
        "fetch_ok": False,
        "design_ok": False,
        "fetch_error": None,
        "design_error": None,
        "cds_length": None,
        "protein_length": None,
        "vector_combo": None,
    }

    tag = f"[{test_idx + 1}/{total}]"

    # Fetch gene sequence
    try:
        fetch_result = fetch_gene_sequence(
            gene_name=gene_name,
            organism="Escherichia coli",
            gene_id=gene_id,
        )
    except Exception as exc:
        result["fetch_error"] = f"EXCEPTION: {type(exc).__name__}: {exc}"
        log(f"  {tag} {gene_name}: FETCH EXCEPTION - {type(exc).__name__}: {exc}")
        return result

    if "error" in fetch_result:
        result["fetch_error"] = fetch_result["error"]
        # Short log for common errors
        err = fetch_result["error"]
        if "not found" in err.lower() or "no valid" in err.lower():
            log(f"  {tag} {gene_name}: no CDS found")
        else:
            log(f"  {tag} {gene_name}: fetch error - {err[:100]}")
        return result

    result["fetch_ok"] = True
    cds_seq = fetch_result.get("cds_seq", "")
    result["cds_length"] = len(cds_seq)
    result["protein_length"] = fetch_result.get("protein_length_aa")

    # Skip very short genes for primer design
    if len(cds_seq) < 30:
        result["design_error"] = f"CDS too short ({len(cds_seq)} bp)"
        log(f"  {tag} {gene_name}: {len(cds_seq)} bp (too short, skip design)")
        return result

    # Design RE cloning primers with random vector combo
    re5, re3, vec = random.choice(VECTOR_COMBOS)
    result["vector_combo"] = f"{re5}/{re3}/{vec}"

    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            design_result = design_re_cloning_primers(
                insert_seq=cds_seq,
                re_5prime=re5,
                re_3prime=re3,
                vector_name=vec,
                gene_name=gene_name,
                gene_description=fetch_result.get("description", ""),
                output_dir=tmpdir,
            )
        except Exception as exc:
            tb_lines = traceback.format_exc().strip().split("\n")
            short_tb = "\n".join(tb_lines[-5:])
            result["design_error"] = f"EXCEPTION: {type(exc).__name__}: {exc}\n{short_tb}"
            log(f"  {tag} {gene_name}: DESIGN EXCEPTION - {type(exc).__name__}: {exc}")
            return result

    if "error" in design_result:
        result["design_error"] = str(design_result.get("error", design_result))
        log(f"  {tag} {gene_name}: design error - {result['design_error'][:100]}")
        return result

    result["design_ok"] = True
    log(f"  {tag} {gene_name}: OK ({len(cds_seq)} bp, {result['vector_combo']})")
    return result


# ── Results ──────────────────────────────────────────────────────────────────

def save_results(out_path, results, errors_fetch, errors_design,
                 fetch_ok, design_ok, total):
    fetch_error_types: dict[str, int] = {}
    for e in errors_fetch:
        err = e["fetch_error"] or ""
        if "EXCEPTION" in err:
            key = "Exception: " + err.split(":")[1].strip() if ":" in err else "Unknown Exception"
        elif "No valid CDS" in err:
            key = "No valid CDS"
        elif "not found" in err.lower():
            key = "Gene not found"
        elif "No linked nucleotide" in err:
            key = "No linked nucleotide"
        else:
            key = err[:80]
        fetch_error_types[key] = fetch_error_types.get(key, 0) + 1

    design_error_types: dict[str, int] = {}
    for e in errors_design:
        err = e["design_error"] or ""
        if "EXCEPTION" in err:
            parts = err.split(":")
            key = "Exception: " + parts[1].strip() if len(parts) > 1 else "Unknown Exception"
        elif "too short" in err:
            key = "CDS too short"
        else:
            key = err[:100]
        design_error_types[key] = design_error_types.get(key, 0) + 1

    summary = {
        "total_tested": total,
        "fetch_ok": fetch_ok,
        "fetch_ok_pct": round(100 * fetch_ok / max(total, 1), 1),
        "design_ok": design_ok,
        "design_ok_pct": round(100 * design_ok / max(fetch_ok, 1), 1),
        "fetch_errors": len(errors_fetch),
        "design_errors": len(errors_design),
        "fetch_error_types": dict(sorted(fetch_error_types.items(), key=lambda x: -x[1])),
        "design_error_types": dict(sorted(design_error_types.items(), key=lambda x: -x[1])),
        "design_exception_samples": [
            {
                "gene_name": e["gene_name"],
                "gene_id": e["gene_id"],
                "cds_length": e["cds_length"],
                "vector_combo": e["vector_combo"],
                "error": (e["design_error"] or "")[:500],
            }
            for e in errors_design
            if (e.get("design_error") or "").startswith("EXCEPTION")
        ][:30],
        "fetch_exception_samples": [
            {
                "gene_name": e["gene_name"],
                "gene_id": e["gene_id"],
                "error": (e["fetch_error"] or "")[:300],
            }
            for e in errors_fetch
            if (e.get("fetch_error") or "").startswith("EXCEPTION")
        ][:10],
    }

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False, default=str)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", type=int, default=None,
                        help="Quick test with N genes (default: full 1500)")
    args = parser.parse_args()

    out_path = Path(__file__).parent / "stress_test_results.json"
    n_test = args.quick or 1500

    # 1. Get gene list
    genes = fetch_ecoli_gene_list(n=2000)
    if not genes:
        log("ERROR: No genes fetched")
        sys.exit(1)

    # 2. Shuffle and take sample
    random.seed(42)
    random.shuffle(genes)
    sample = genes[:n_test]

    log(f"\n{'='*60}")
    log(f"  Starting stress test: {len(sample)} genes")
    log(f"{'='*60}\n")

    # 3. Test each gene
    results = []
    errors_fetch = []
    errors_design = []
    fetch_ok = 0
    design_ok = 0

    for idx, gene_info in enumerate(sample):
        r = test_single_gene(gene_info, idx, len(sample))
        results.append(r)

        if r["fetch_ok"]:
            fetch_ok += 1
        elif r["fetch_error"]:
            errors_fetch.append(r)

        if r["design_ok"]:
            design_ok += 1
        elif r["design_error"] and r["fetch_ok"]:
            errors_design.append(r)

        # Progress every 50
        if (idx + 1) % 50 == 0:
            log(f"\n--- Progress: {idx + 1}/{len(sample)} | "
                f"fetch={fetch_ok} design={design_ok} "
                f"fetch_err={len(errors_fetch)} design_err={len(errors_design)} ---\n")
            save_results(out_path, results, errors_fetch, errors_design,
                         fetch_ok, design_ok, len(results))

    # 4. Save final results
    save_results(out_path, results, errors_fetch, errors_design,
                 fetch_ok, design_ok, len(results))

    log(f"\n{'='*60}")
    log(f"  STRESS TEST COMPLETE")
    log(f"{'='*60}")
    log(f"Total tested: {len(results)}")
    log(f"Fetch OK:  {fetch_ok}/{len(results)} ({100*fetch_ok/len(results):.1f}%)")
    log(f"Design OK: {design_ok}/{fetch_ok} fetched ({100*design_ok/max(fetch_ok,1):.1f}%)")
    log(f"Fetch errors:  {len(errors_fetch)}")
    log(f"Design errors: {len(errors_design)}")
    log(f"Results: {out_path}")


if __name__ == "__main__":
    main()
