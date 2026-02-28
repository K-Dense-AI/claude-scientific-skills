#!/usr/bin/env python3
"""
PsXR (PsXyl1) Cofactor Specificity Mutagenesis - iPCR Primer Designer
======================================================================
Khattab et al. (2011) BBRC 404:634-637 기반
E223A/S271A (AA) 및 E223D/S271A (DA) 이중 변이체 제작용 inverse PCR primer 설계

Target constructs:
  1. pETduet-1_PsXyl1_PsFdhV9  (Whole cell 폴더)
  2. pacyc_xr.dna              (Genes 폴더)

Mutations:
  E223A: GAA -> GCA  (position 223, Glu -> Ala)
  E223D: GAA -> GAC  (position 223, Glu -> Asp)
  S271A: TCC -> GCC  (position 271, Ser -> Ala)
"""

import os
from datetime import datetime
from Bio.Seq import Seq

from ..subst_primer_mode import iPCRSubstDesigner
from ..snapgene_parser import parse_snapgene


# ══════════════════════════════════════════════════════════════════════════
#  PsXyl1 Mutation Finder
# ══════════════════════════════════════════════════════════════════════════

def find_xr_cds(features, seq):
    """Find PsXyl1/XR CDS feature and return gene info."""
    # Try multiple name patterns
    xr_names = ["PsXyl1", "XR", "xyl1", "Xyl1", "xylose reductase"]

    candidates = []
    for f in features:
        for name_pattern in xr_names:
            if name_pattern.lower() in f["name"].lower():
                if f["type"] == "CDS" or f["type"] == "gene":
                    candidates.append(f)
                    break

    if not candidates:
        # Fallback: scan all CDS features for PsXR-sized proteins (~310-320 aa)
        for f in features:
            if f["type"] == "CDS":
                cds_len = f["end"] - f["start"]
                if 900 < cds_len < 1100:
                    candidates.append(f)

    if not candidates:
        # Last resort: scan the whole plasmid for PsXR motif
        return find_xr_by_motif(seq)

    # Pick longest
    return max(candidates, key=lambda f: f["end"] - f["start"])


def find_xr_by_motif(seq):
    """
    Find PsXR CDS by searching for conserved protein motifs.
    PsXR has distinctive motifs near E223 and S271 regions.
    Returns a synthetic feature dict or None.
    """
    # PsXR conserved motif around E223: IPKS + NIEL + E + NQGR
    # Search all 3 reading frames on both strands
    for strand_name, search_seq in [("fwd", seq), ("rev", str(Seq(seq).reverse_complement()))]:
        for frame in range(3):
            protein = str(Seq(search_seq[frame:]).translate())
            # Look for NIELE motif (around E223 in PsXR)
            idx = protein.find("NIELE")
            if idx == -1:
                idx = protein.find("IELE")
            if idx >= 0 and len(protein) > idx + 100:
                # Found! Now find the full ORF
                # Search backwards for Met
                met_pos = protein.rfind("M", 0, idx)
                if met_pos >= 0 and (idx - met_pos) > 180:
                    # Check total length
                    stop_pos = protein.find("*", met_pos)
                    if stop_pos > 0 and 280 < (stop_pos - met_pos) < 350:
                        # Calculate DNA positions
                        if strand_name == "fwd":
                            dna_start = frame + met_pos * 3
                            dna_end = frame + stop_pos * 3 + 3
                        else:
                            dna_end = len(seq) - frame - met_pos * 3
                            dna_start = len(seq) - frame - stop_pos * 3 - 3
                        return {
                            "name": "PsXyl1 (motif-found)",
                            "type": "CDS",
                            "start": dna_start,
                            "end": dna_end,
                            "strand": 1 if strand_name == "fwd" else 2,
                        }
    return None


def find_mutation_site_by_motif(protein, target_aa, motif_patterns):
    """
    Find the correct position of a target amino acid using conserved motifs.
    Returns 1-indexed position in the protein, or None.

    PsXR E223 context: ...IPKS-NIEL-E-NQGR... (E is the target)
    PsXR S271 context: ...ASNF-SNTV-S-VPRI... (S is the target, but check NTVS or SNTV)
    """
    for motif in motif_patterns:
        idx = protein.find(motif)
        if idx >= 0:
            # Find the target_aa within or adjacent to the motif
            # The motif should contain or be immediately before/after the target
            for offset in range(len(motif)):
                if protein[idx + offset] == target_aa:
                    return idx + offset + 1  # 1-indexed
            # Check positions just after the motif
            if idx + len(motif) < len(protein) and protein[idx + len(motif)] == target_aa:
                return idx + len(motif) + 1
    return None


def find_codon_position(cds_seq, aa_pos, expected_codon=None):
    """
    Find DNA position of amino acid at aa_pos (1-indexed).
    Returns (dna_offset_within_cds, codon_str).
    """
    dna_offset = (aa_pos - 1) * 3
    codon = cds_seq[dna_offset:dna_offset + 3].upper()

    if expected_codon and codon != expected_codon.upper():
        print(f"  WARNING: Expected codon {expected_codon} at aa {aa_pos}, "
              f"found {codon} (translates to {str(Seq(codon).translate())})")

    return dna_offset, codon


# ══════════════════════════════════════════════════════════════════════════
#  Main Analysis
# ══════════════════════════════════════════════════════════════════════════

def analyze_construct(filepath, construct_name):
    """Analyze a single construct and design primers for all mutations."""
    sep = "=" * 70
    results = {}
    designer = iPCRSubstDesigner()

    print(f"\n{'#' * 70}")
    print(f"  Construct: {construct_name}")
    print(f"  File: {os.path.basename(filepath)}")
    print(f"{'#' * 70}")

    # 1. Parse .dna file
    seq, circular, features = parse_snapgene(filepath)
    print(f"\n  Plasmid: {len(seq)} bp ({'circular' if circular else 'linear'})")
    print(f"  Features found: {len(features)}")

    # List all features
    print(f"\n  All features:")
    for f in features:
        print(f"    {f['name']:35s} {f['type']:12s} {f['start']:>5d}-{f['end']:<5d} "
              f"({'fwd' if f['strand'] == 1 else 'rev'})")

    # 2. Find PsXyl1 CDS
    xr = find_xr_cds(features, seq)
    if xr is None:
        print(f"\n  ERROR: PsXyl1/XR CDS not found!")
        print(f"  Trying manual search by ATG + translation...")
        # Try to find XR by searching for known motifs
        return None

    xr_start = xr["start"]
    xr_end = xr["end"]
    cds_seq = seq[xr_start:xr_end]
    protein = str(Seq(cds_seq).translate())

    print(f"\n  PsXyl1 CDS: {xr['name']}")
    print(f"    Position: {xr_start}-{xr_end} (0-indexed)")
    print(f"    Length: {len(cds_seq)} bp -> {len(protein)} aa")
    print(f"    Start: {cds_seq[:12]}... -> {protein[:4]}...")
    print(f"    Stop:  ...{cds_seq[-12:]} -> ...{protein[-4:]}")

    # 3. Find E223 and S271 by MOTIF matching (paper numbering vs construct numbering)
    print(f"\n  Motif-based residue identification:")

    # E223 motifs: conserved region around Glu223 in PsXR
    e223_motifs = ["NIELE", "IELE", "NIELEN", "IPKSNIELE"]
    e223_pos = find_mutation_site_by_motif(protein, "E", e223_motifs)

    if e223_pos is None:
        # Fallback: scan for E near expected position
        for i in range(max(0, 218), min(len(protein), 228)):
            if protein[i] == "E":
                context = protein[max(0, i-3):i+4]
                if "N" in context and "L" in context:  # NIELE pattern
                    e223_pos = i + 1
                    break

    if e223_pos:
        e223_offset = (e223_pos - 1) * 3
        e223_codon = cds_seq[e223_offset:e223_offset + 3].upper()
        e223_abs = xr_start + e223_offset
        e223_aa = str(Seq(e223_codon).translate())
        ctx_start = max(0, e223_pos - 5)
        ctx_end = min(len(protein), e223_pos + 4)
        print(f"    E223 (paper) -> construct position {e223_pos}: "
              f"codon {e223_codon} -> {e223_aa}")
        print(f"      Context: ...{protein[ctx_start:e223_pos-1]}"
              f"[{protein[e223_pos-1]}]"
              f"{protein[e223_pos:ctx_end]}...")
        assert e223_aa == "E", f"Expected E, found {e223_aa} at position {e223_pos}"
    else:
        print(f"    ERROR: Could not locate E223 equivalent!")
        return None

    # S271 motifs: conserved region around Ser271 in PsXR
    s271_motifs = ["SNTVS", "NTVS", "TVSVP", "SNTVSV"]
    s271_pos = find_mutation_site_by_motif(protein, "S", s271_motifs)

    if s271_pos is None:
        # Fallback: scan for S near expected position with context
        for i in range(max(0, 266), min(len(protein), 276)):
            if protein[i] == "S":
                context = protein[max(0, i-3):i+4]
                if "V" in context and "P" in context:
                    s271_pos = i + 1
                    break

    if s271_pos:
        s271_offset = (s271_pos - 1) * 3
        s271_codon = cds_seq[s271_offset:s271_offset + 3].upper()
        s271_abs = xr_start + s271_offset
        s271_aa = str(Seq(s271_codon).translate())
        ctx_start = max(0, s271_pos - 5)
        ctx_end = min(len(protein), s271_pos + 4)
        print(f"    S271 (paper) -> construct position {s271_pos}: "
              f"codon {s271_codon} -> {s271_aa}")
        print(f"      Context: ...{protein[ctx_start:s271_pos-1]}"
              f"[{protein[s271_pos-1]}]"
              f"{protein[s271_pos:ctx_end]}...")
        assert s271_aa == "S", f"Expected S, found {s271_aa} at position {s271_pos}"
    else:
        print(f"    ERROR: Could not locate S271 equivalent!")
        return None

    # Paper numbering offset
    offset = 223 - e223_pos
    print(f"\n    Paper numbering offset: {offset} "
          f"(paper pos = construct pos + {offset})")
    print(f"    Verification: S271 paper -> construct {271 - offset} "
          f"(actual found: {s271_pos}) "
          f"{'MATCH' if (271 - offset) == s271_pos else 'MISMATCH!'}")

    # 4. Define mutations using CORRECT absolute positions
    mutations = {
        "E223A": {
            "aa_pos": e223_pos,
            "paper_pos": 223,
            "old_codon": e223_codon,
            "new_codon": "GCA",
            "description": "Glu223 -> Ala (NADH binding elimination)",
            "abs_pos": e223_abs,
        },
        "E223D": {
            "aa_pos": e223_pos,
            "paper_pos": 223,
            "old_codon": e223_codon,
            "new_codon": "GAC",
            "description": "Glu223 -> Asp (NADH binding reduction)",
            "abs_pos": e223_abs,
        },
        "S271A": {
            "aa_pos": s271_pos,
            "paper_pos": 271,
            "old_codon": s271_codon,
            "new_codon": "GCC",
            "description": "Ser271 -> Ala (NADPH preference enhancement)",
            "abs_pos": s271_abs,
        },
    }

    # 5. Design primers for each single mutation
    print(f"\n{'=' * 70}")
    print(f"  PRIMER DESIGN RESULTS")
    print(f"{'=' * 70}")

    for mut_name, mut_info in mutations.items():
        print(f"\n{'-' * 70}")
        print(f"  Mutation: {mut_name} ({mut_info['description']})")
        print(f"  Codon change: {mut_info['old_codon']} -> {mut_info['new_codon']}")
        print(f"  Plasmid position: {mut_info['abs_pos']}")
        print(f"{'-' * 70}")

        try:
            result = designer.design(
                seq=seq,
                subst_pos=mut_info['abs_pos'],
                old_seq=mut_info['old_codon'],
                new_seq=mut_info['new_codon'],
                target_tm=61.0,
                overlap_len=18,
                min_len=18,
                max_len=35,
            )

            # Primer names
            f_name = f"iPCR_{mut_name}_F"
            r_name = f"iPCR_{mut_name}_R"

            # Quality results (integrated in designer.design)
            f_qc = result['f_qc']
            r_qc = result['r_qc']
            het = result['het']

            print(f"\n  Forward: {f_name}")
            print(f"    5'-{result['f_full']}-3'")
            print(f"    Length: {result['f_len']} nt")
            print(f"    Annealing (non-overlap): {result['f_ann']}  ({len(result['f_ann'])} bp)")
            print(f"    Effective binding (dn_tail+ann): {result['f_eff_bind']}  ({len(result['f_eff_bind'])} bp)")
            print(f"    Eff. binding Tm={result['f_tm']}C  GC={result['f_gc']}%")
            print(f"    Overlap tail: [{result['up_tail']}][{mut_info['new_codon']}][{result['dn_tail']}]")

            print(f"\n  Reverse: {r_name}")
            print(f"    5'-{result['r_full']}-3'")
            print(f"    Length: {result['r_len']} nt")
            print(f"    Annealing (non-overlap): {result['r_ann']}  ({len(result['r_ann'])} bp)")
            print(f"    Effective binding (up_tail_rc+ann): {result['r_eff_bind']}  ({len(result['r_eff_bind'])} bp)")
            print(f"    Eff. binding Tm={result['r_tm']}C  GC={result['r_gc']}%")

            print(f"\n  Overlap: {result['overlap_region']} ({len(result['overlap_region'])} bp)")
            print(f"    Verified: {'YES' if result['overlap_verified'] else 'NO'}")
            print(f"  Anneal temp: {result['anneal_temp']}C (Q5: min Tm + 1C)")

            print(f"\n  Quality (eff. binding, Tm-based at {result['anneal_temp']}C):")
            print(f"    F: {f_qc['verdict']}  "
                  f"hairpin Tm={f_qc.get('hairpin_tm','N/A')}C (dG={f_qc.get('hairpin_dg','N/A')})  "
                  f"homodimer Tm={f_qc.get('homodimer_tm','N/A')}C (dG={f_qc.get('homodimer_dg','N/A')})")
            print(f"    R: {r_qc['verdict']}  "
                  f"hairpin Tm={r_qc.get('hairpin_tm','N/A')}C (dG={r_qc.get('hairpin_dg','N/A')})  "
                  f"homodimer Tm={r_qc.get('homodimer_tm','N/A')}C (dG={r_qc.get('homodimer_dg','N/A')})")
            print(f"    Heterodimer: dG={het.get('dg','N/A')} kcal/mol (expected: 18bp overlap)")

            if result['warnings']:
                print(f"\n  Warnings:")
                for w in result['warnings']:
                    print(f"    - {w}")

            if f_qc.get('issues'):
                for iss in f_qc['issues']:
                    print(f"    - F: {iss}")
            if r_qc.get('issues'):
                for iss in r_qc['issues']:
                    print(f"    - R: {iss}")

            results[mut_name] = {
                "result": result,
                "f_name": f_name,
                "r_name": r_name,
                "f_qc": f_qc,
                "r_qc": r_qc,
                "het": het,
                "mut_info": mut_info,
            }

        except Exception as e:
            print(f"\n  ERROR: {e}")
            results[mut_name] = {"error": str(e)}

    # 6. Double mutant strategy summary
    print(f"\n{'=' * 70}")
    print(f"  DOUBLE MUTANT STRATEGY")
    print(f"{'=' * 70}")

    print(f"\n  AA (E223A/S271A) - Strictly NADPH-dependent, 106% WT activity:")
    print(f"    Round 1: iPCR with E223A primers -> DpnI -> KLD -> Transform")
    print(f"    Round 2: iPCR with S271A primers on E223A template -> DpnI -> KLD -> Transform")
    print(f"    (or reverse order: S271A first, then E223A)")

    print(f"\n  DA (E223D/S271A) - NADPH 1.27x, kcat/Km 1.26x:")
    print(f"    Round 1: iPCR with E223D primers -> DpnI -> KLD -> Transform")
    print(f"    Round 2: iPCR with S271A primers on E223D template -> DpnI -> KLD -> Transform")

    # PCR conditions
    product_size = len(seq)
    ext_time = product_size // 1000 + 1

    print(f"\n  PCR Conditions (both rounds):")
    print(f"    Polymerase: Q5 High-Fidelity (NEB #M0491)")
    print(f"    Template: {len(seq)} bp (circular) -> product ~{product_size} bp")
    print(f"    Extension: {ext_time} min (Q5: 1 min/kb)")
    print(f"    Post-PCR: KLD (Kinase-Ligase-DpnI, NEB #M0554)")

    return {
        "construct_name": construct_name,
        "filepath": filepath,
        "plasmid_size": len(seq),
        "circular": circular,
        "xr_start": xr_start,
        "xr_end": xr_end,
        "cds_length": len(cds_seq),
        "protein_length": len(protein),
        "mutations": results,
        "product_size": product_size,
        "ext_time": ext_time,
    }


# ══════════════════════════════════════════════════════════════════════════
#  Report Generator
# ══════════════════════════════════════════════════════════════════════════

def generate_report(all_results, output_path):
    """Generate markdown report."""
    designer = iPCRSubstDesigner()
    lines = []
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    lines.append("# PsXR (PsXyl1) Cofactor Specificity Mutagenesis - iPCR Primer Report")
    lines.append("")
    lines.append(f"**Date**: {now}")
    lines.append(f"**Reference**: Khattab et al. (2011) BBRC 404:634-637")
    lines.append(f"**Method**: Inverse PCR (Q5 High-Fidelity) + KLD")
    lines.append("")

    lines.append("---")
    lines.append("")
    lines.append("## 1. Background")
    lines.append("")
    lines.append("PsXR (Pichia stipitis xylose reductase, PsXyl1)의 cofactor specificity를 ")
    lines.append("변경하여 strictly NADPH-dependent 변이체를 제작한다.")
    lines.append("")
    lines.append("| Mutant | Mutations | Cofactor | NADPH Activity | NADH Activity |")
    lines.append("|--------|-----------|----------|----------------|---------------|")
    lines.append("| WT | - | Dual (NADPH preferred) | 100% | 100% |")
    lines.append("| **AA** | **E223A/S271A** | **Strictly NADPH** | **106%** | **ND** |")
    lines.append("| DA | E223D/S271A | NADPH (nearly strict) | 127% | ~15% WT |")
    lines.append("")
    lines.append("### Mutation Rationale")
    lines.append("")
    lines.append("- **E223A**: Glu223는 NAD+ adenosine의 2'/3'-OH와 bidentate 수소결합 형성. ")
    lines.append("  Ala 치환으로 비극성 잔기가 되어 NADH 결합 완전 차단.")
    lines.append("- **E223D**: Asp는 Glu보다 짧은 산성 잔기. NADH 결합을 줄이면서 부분 활성 유지.")
    lines.append("- **S271A**: Ser271은 NADP+의 2'-phosphate와 상호작용. ")
    lines.append("  Ala 치환이 오히려 NADPH 결합을 개선 (E223A에 의한 활성 감소 보상).")
    lines.append("")

    lines.append("---")
    lines.append("")
    lines.append("## 2. Target Constructs")
    lines.append("")

    for res in all_results:
        if res is None:
            continue
        lines.append(f"### {res['construct_name']}")
        lines.append("")
        lines.append(f"- **File**: `{os.path.basename(res['filepath'])}`")
        lines.append(f"- **Plasmid size**: {res['plasmid_size']} bp ({'circular' if res['circular'] else 'linear'})")
        lines.append(f"- **PsXyl1 CDS**: pos {res['xr_start']}-{res['xr_end']} "
                     f"({res['cds_length']} bp, {res['protein_length']} aa)")
        lines.append("")

    lines.append("---")
    lines.append("")
    lines.append("## 3. Primer Design")
    lines.append("")
    lines.append("### Design Parameters")
    lines.append("")
    lines.append("| Parameter | Value |")
    lines.append("|-----------|-------|")
    lines.append("| Tm calculation | Nearest-neighbor (Owczarzy 2008, saltcorr=7) |")
    lines.append("| Target Tm | 61C (annealing region) |")
    lines.append(f"| Salt conditions | Na={designer.na} mM, Mg={designer.mg} mM |")
    lines.append(f"| Primer conc. | {designer.dnac1} nM |")
    lines.append(f"| dNTPs | {designer.dntps} mM |")
    lines.append("| Overlap length | 18 bp |")
    lines.append("| Min/Max annealing | 18-35 nt |")
    lines.append("")

    for res in all_results:
        if res is None:
            continue

        lines.append(f"### {res['construct_name']}")
        lines.append("")

        for mut_name in ["E223A", "E223D", "S271A"]:
            mut_data = res['mutations'].get(mut_name, {})
            if "error" in mut_data:
                lines.append(f"#### {mut_name} - ERROR: {mut_data['error']}")
                lines.append("")
                continue

            if "result" not in mut_data:
                continue

            r = mut_data['result']
            mi = mut_data['mut_info']

            lines.append(f"#### {mut_name}: {mi['old_codon']} -> {mi['new_codon']} ({mi['description']})")
            lines.append("")
            lines.append(f"| Primer | Name | Sequence (5'->3') | Length | Tm (ann) | GC% (ann) |")
            lines.append(f"|--------|------|-------------------|--------|----------|-----------|")
            lines.append(f"| Forward | {mut_data['f_name']} | `{r['f_full']}` | {r['f_len']} nt | {r['f_tm']}C | {r['f_gc']}% |")
            lines.append(f"| Reverse | {mut_data['r_name']} | `{r['r_full']}` | {r['r_len']} nt | {r['r_tm']}C | {r['r_gc']}% |")
            lines.append("")
            lines.append(f"- **Overlap**: `{r['overlap_region']}` ({len(r['overlap_region'])} bp)")
            lines.append(f"- **Overlap verified**: {'YES' if r['overlap_verified'] else 'NO'}")
            lines.append(f"- **Recommended annealing temp**: {r['anneal_temp']}C")

            # Quality
            f_qc = mut_data.get('f_qc', {})
            r_qc = mut_data.get('r_qc', {})
            het = mut_data.get('het', {})

            lines.append(f"- **F quality**: {f_qc.get('verdict', 'N/A')} "
                        f"(hairpin dG={f_qc.get('hairpin_dg', 'N/A')} kcal/mol)")
            lines.append(f"- **R quality**: {r_qc.get('verdict', 'N/A')} "
                        f"(hairpin dG={r_qc.get('hairpin_dg', 'N/A')} kcal/mol)")
            lines.append(f"- **Heterodimer**: dG={het.get('dg', 'N/A')} kcal/mol")

            if r['warnings']:
                lines.append(f"- **Warnings**: {', '.join(r['warnings'])}")

            lines.append("")

    # 4. Protocol
    lines.append("---")
    lines.append("")
    lines.append("## 4. Double Mutant Construction Protocol")
    lines.append("")
    lines.append("### AA 변이체 (E223A/S271A) - **권장**")
    lines.append("")
    lines.append("Strictly NADPH-dependent, WT 대비 106% NADPH 활성, NADH 활성 ND")
    lines.append("")
    lines.append("```")
    lines.append("Round 1: E223A 도입")
    lines.append("  Template: WT plasmid (pETduet-1 or pACYCduet-1)")
    lines.append("  Primers: iPCR_E223A_F + iPCR_E223A_R")
    lines.append("  -> Q5 iPCR -> KLD (NEB #M0554) -> Transform -> Colony PCR -> Sequencing")
    lines.append("")
    lines.append("Round 2: S271A 도입")
    lines.append("  Template: E223A mutant plasmid (Round 1 결과)")
    lines.append("  Primers: iPCR_S271A_F + iPCR_S271A_R")
    lines.append("  -> Q5 iPCR -> KLD -> Transform -> Colony PCR -> Sequencing")
    lines.append("```")
    lines.append("")
    lines.append("### DA 변이체 (E223D/S271A)")
    lines.append("")
    lines.append("NADPH 1.27x, kcat/Km 1.26x (대안)")
    lines.append("")
    lines.append("```")
    lines.append("Round 1: E223D 도입")
    lines.append("  Primers: iPCR_E223D_F + iPCR_E223D_R")
    lines.append("")
    lines.append("Round 2: S271A 도입")
    lines.append("  Primers: iPCR_S271A_F + iPCR_S271A_R")
    lines.append("```")
    lines.append("")

    # 5. PCR conditions table
    lines.append("### PCR Conditions")
    lines.append("")
    lines.append("| Step | Temperature | Time | Cycles |")
    lines.append("|------|-------------|------|--------|")
    lines.append("| Initial denaturation | 98C | 30 sec | 1 |")

    # Get anneal temps from results
    anneal_temps = set()
    for res in all_results:
        if res is None:
            continue
        for mut_name, mut_data in res['mutations'].items():
            if 'result' in mut_data:
                anneal_temps.add(mut_data['result']['anneal_temp'])

    anneal_str = "/".join(f"{t}C" for t in sorted(anneal_temps)) if anneal_temps else "62C"

    lines.append(f"| Denaturation | 98C | 10 sec | 25 |")
    lines.append(f"| Annealing | {anneal_str} | 20 sec | 25 |")

    ext_times = set()
    for res in all_results:
        if res is None:
            continue
        ext_times.add(res['ext_time'])
    ext_str = "/".join(f"{t} min" for t in sorted(ext_times)) if ext_times else "4 min"

    lines.append(f"| Extension | 72C | {ext_str} | 25 |")
    lines.append("| Final extension | 72C | 2 min | 1 |")
    lines.append("| Hold | 4C | - | - |")
    lines.append("")
    lines.append("### Post-PCR Treatment")
    lines.append("")
    lines.append("```")
    lines.append("KLD Reaction (NEB #M0554):")
    lines.append("  1 uL PCR product")
    lines.append("  1 uL 10x KLD Enzyme Mix")
    lines.append("  5 uL 2x KLD Reaction Buffer")
    lines.append("  3 uL nuclease-free water")
    lines.append("  -> 25C, 5 min -> Transform 5 uL into competent cells")
    lines.append("```")
    lines.append("")

    # 6. Sequencing primers
    lines.append("---")
    lines.append("")
    lines.append("## 5. Confirmation Strategy")
    lines.append("")
    lines.append("### Sequencing Primers")
    lines.append("")
    lines.append("| Primer | Use | Target |")
    lines.append("|--------|-----|--------|")
    lines.append("| pET-upstream | pETduet-1 MCS1 | E223 region |")
    lines.append("| DuetDown1 | pETduet-1 MCS1 reverse | S271 region |")
    lines.append("| T7 promoter | pET/pACYC forward | Upstream confirmation |")
    lines.append("| T7 terminator | pET/pACYC reverse | Downstream confirmation |")
    lines.append("")
    lines.append("### Expected Sequencing Results")
    lines.append("")
    lines.append("| Mutation | WT Codon | Mutant Codon | AA Change |")
    lines.append("|----------|----------|--------------|-----------|")

    for res in all_results:
        if res is None:
            continue
        for mut_name in ["E223A", "E223D", "S271A"]:
            mut_data = res['mutations'].get(mut_name, {})
            if 'mut_info' in mut_data:
                mi = mut_data['mut_info']
                old_aa = str(Seq(mi['old_codon']).translate())
                new_aa = str(Seq(mi['new_codon']).translate())
                lines.append(f"| {mut_name} | {mi['old_codon']} ({old_aa}) | "
                           f"{mi['new_codon']} ({new_aa}) | {old_aa}{mi['aa_pos']}{new_aa} |")
        break  # Only need one table

    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 6. Primer Order Summary")
    lines.append("")
    lines.append("### pETduet-1 construct")
    lines.append("")
    lines.append("| # | Primer Name | Sequence (5'->3') | Length | Purpose |")
    lines.append("|---|-------------|-------------------|--------|---------|")

    idx = 1
    for res in all_results:
        if res is None:
            continue
        for mut_name in ["E223A", "E223D", "S271A"]:
            mut_data = res['mutations'].get(mut_name, {})
            if 'result' in mut_data:
                r = mut_data['result']
                lines.append(f"| {idx} | {mut_data['f_name']} | `{r['f_full']}` | {r['f_len']} nt | {mut_name} |")
                idx += 1
                lines.append(f"| {idx} | {mut_data['r_name']} | `{r['r_full']}` | {r['r_len']} nt | {mut_name} |")
                idx += 1
        break  # First construct table

    lines.append("")

    # Second construct
    if len([r for r in all_results if r is not None]) > 1:
        lines.append("### pACYCduet-1 construct")
        lines.append("")
        lines.append("| # | Primer Name | Sequence (5'->3') | Length | Purpose |")
        lines.append("|---|-------------|-------------------|--------|---------|")
        idx = 1
        for i, res in enumerate(all_results):
            if res is None or i == 0:
                continue
            for mut_name in ["E223A", "E223D", "S271A"]:
                mut_data = res['mutations'].get(mut_name, {})
                if 'result' in mut_data:
                    r = mut_data['result']
                    lines.append(f"| {idx} | {mut_data['f_name']} | `{r['f_full']}` | {r['f_len']} nt | {mut_name} |")
                    idx += 1
                    lines.append(f"| {idx} | {mut_data['r_name']} | `{r['r_full']}` | {r['r_len']} nt | {mut_name} |")
                    idx += 1
        lines.append("")

    lines.append("---")
    lines.append("")
    lines.append(f"*Generated by psxr_cofactor_mutagenesis.py on {now}*")
    lines.append(f"*Tm: Nearest-neighbor (Owczarzy 2008, saltcorr=7), "
                f"Na={designer.na} mM, Mg={designer.mg} mM, "
                f"[primer]={designer.dnac1} nM, dNTPs={designer.dntps} mM*")

    report = "\n".join(lines)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(report)

    print(f"\nReport saved: {output_path}")
    return report


# ══════════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════════

def main():
    base_dir = (
        r"C:\Users\Jahyun\OneDrive - 고려대학교\저장소"
        r"\8. D-Gal to MA and D-tagatose\Genes"
    )

    constructs = [
        {
            "path": os.path.join(base_dir, "Whole cell",
                                 "Cloned_pETduet-1_PsXyl1_PsFdhV9.dna"),
            "name": "pETduet-1_PsXyl1_PsFdhV9",
        },
        {
            "path": os.path.join(base_dir, "pacyc_xr.dna"),
            "name": "pACYCduet-1_PsXyl1 (pacyc_xr)",
        },
    ]

    all_results = []
    for construct in constructs:
        if not os.path.exists(construct["path"]):
            print(f"\nWARNING: File not found: {construct['path']}")
            all_results.append(None)
            continue
        result = analyze_construct(construct["path"], construct["name"])
        all_results.append(result)

    # Generate report
    output_dir = r"C:\Users\Jahyun\PycharmProjects\pythonProject1"
    report_path = os.path.join(output_dir, "PsXR_cofactor_mutagenesis_report.md")
    generate_report(all_results, report_path)

    print("\n" + "=" * 70)
    print("  DONE")
    print("=" * 70)


if __name__ == "__main__":
    main()
