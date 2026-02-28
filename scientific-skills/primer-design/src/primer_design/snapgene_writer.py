"""SnapGene .dna binary file writer — in-silico RE cloning.

Base vector .dna 파일에 insert CDS를 RE cloning 시뮬레이션으로 삽입하여
재조합 플라스미드 .dna 파일을 생성.

Binary format:
  Block = [type: 1 byte][length: 4 bytes big-endian][data: N bytes]
  Type 9  : Cookie/header (magic "SnapGene")
  Type 0  : DNA sequence + flags (circular bit)
  Type 10 : Features (XML)
"""

from __future__ import annotations

import re
import struct
import xml.etree.ElementTree as ET
from pathlib import Path

from Bio.Seq import Seq

from .snapgene_parser import parse_snapgene

# ── Genetic code for translation ─────────────────────────────────────────────

_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def _translate(dna: str) -> str:
    """DNA → amino acid string (stop='*')."""
    aa = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3].upper()
        aa.append(_CODON_TABLE.get(codon, "X"))
    return "".join(aa)


def _mol_weight(protein: str) -> float:
    """Approximate molecular weight (Da) from amino acid string."""
    avg_aa_mw = 110.0  # average amino acid MW
    water = 18.015
    # Remove stop codons
    protein = protein.replace("*", "")
    return len(protein) * avg_aa_mw - (len(protein) - 1) * water + water


# ── Low-level block writers ──────────────────────────────────────────────────

def _write_block(block_type: int, data: bytes) -> bytes:
    """Encode a single SnapGene block."""
    return struct.pack("B", block_type) + struct.pack(">I", len(data)) + data


def _build_cookie() -> bytes:
    """Cookie/header block (type 9)."""
    cookie = b"SnapGene\x00"
    cookie += b"\x01"                     # file type: DNA
    cookie += struct.pack(">H", 1)        # export version
    cookie += struct.pack(">H", 14)       # import version
    return _write_block(9, cookie)


def _build_sequence_block(sequence: str, circular: bool = True) -> bytes:
    """DNA sequence block (type 0)."""
    flags = 0x01 if circular else 0x00
    data = struct.pack("B", flags) + sequence.upper().encode("ascii")
    return _write_block(0, data)


def _build_features_block(features: list[dict]) -> bytes:
    """Features block (type 10) with full SnapGene-compatible XML.

    Each feature dict:
        name, type, start (0-based), end (0-based, inclusive),
        strand (1=fwd, 2=rev), color, translated (bool),
        translation (str), qualifiers (list of (name, value) tuples)
    """
    root = ET.Element("Features")
    for i, feat in enumerate(features):
        attrs = {
            "recentID": str(i),
            "name": feat["name"],
            "type": feat.get("type", "misc_feature"),
            "allowSegmentOverlaps": "0",
            "consecutiveTranslationNumbering": "1",
        }
        strand = feat.get("strand", 1)
        if strand in (1, 2):
            attrs["directionality"] = str(strand)

        if feat.get("translated"):
            prot = feat.get("translation", "")
            mw = _mol_weight(prot)
            attrs["translationMW"] = f"{mw:.2f}"
            # Reading frame direction (required by SnapGene for correct
            # codon display): negative = reverse strand, positive = forward
            if strand == 2:
                attrs["readingFrame"] = "-1"
                attrs["swappedSegmentNumbering"] = "1"
            else:
                attrs["readingFrame"] = "1"

        f_elem = ET.SubElement(root, "Feature", attrs)

        seg_attrs = {
            "range": f"{feat['start'] + 1}-{feat['end'] + 1}",  # 0-based → 1-based XML
            "type": "standard",
        }
        if "color" in feat:
            seg_attrs["color"] = feat["color"]
        if feat.get("translated"):
            seg_attrs["translated"] = "1"
        ET.SubElement(f_elem, "Segment", seg_attrs)

        # Qualifiers
        def _add_q(name: str, text: str | None = None, intval: int | None = None):
            q = ET.SubElement(f_elem, "Q", {"name": name})
            if text is not None:
                ET.SubElement(q, "V", {"text": text})
            elif intval is not None:
                ET.SubElement(q, "V", {"int": str(intval)})

        _add_q("label", text=feat["name"])

        if feat.get("translated"):
            _add_q("codon_start", intval=1)
            _add_q("transl_table", intval=1)
            prot = feat.get("translation", "")
            if prot:
                _add_q("translation", text=prot.rstrip("*"))

        for qname, qval in feat.get("qualifiers", []):
            _add_q(qname, text=qval)

    xml_str = ET.tostring(root, encoding="unicode")
    return _write_block(10, xml_str.encode("utf-8"))


# ── In-silico cloning ────────────────────────────────────────────────────────

def _find_site(seq: str, recognition: str) -> list[int]:
    """Find all occurrences of a recognition site in a sequence."""
    positions = []
    start = 0
    rec_upper = recognition.upper()
    seq_upper = seq.upper()
    while True:
        idx = seq_upper.find(rec_upper, start)
        if idx == -1:
            break
        positions.append(idx)
        start = idx + 1
    return positions


def _parse_vector_features(filepath: str | Path) -> list[dict]:
    """Parse feature XML from a .dna file, returning feature dicts."""
    with open(filepath, "rb") as fh:
        data = fh.read()

    features = []
    i = 0
    while i < len(data) - 5:
        btype = data[i]
        blen = int.from_bytes(data[i + 1:i + 5], "big")
        if i + 5 + blen > len(data):
            break

        if btype == 10:
            xml_str = data[i + 5:i + 5 + blen].decode("utf-8", errors="replace")
            root = ET.fromstring(xml_str)
            for feat in root.findall(".//Feature"):
                name = feat.get("name", "")
                ftype = feat.get("type", "")
                direc = int(feat.get("directionality", "0"))
                translated = False

                for seg in feat.findall("Segment"):
                    rng = seg.get("range", "")
                    if seg.get("translated") == "1":
                        translated = True
                    if "-" in rng:
                        s, e = rng.split("-", 1)
                        color = seg.get("color", "")
                        feat_dict = {
                            "name": name,
                            "type": ftype,
                            "start": int(s) - 1,   # 1-based XML → 0-based
                            "end": int(e) - 1,      # 1-based XML → 0-based
                            "strand": direc,
                            "translated": translated,
                        }
                        if color:
                            feat_dict["color"] = color
                        features.append(feat_dict)
            break

        i += 5 + blen

    return features


# ── Public API ───────────────────────────────────────────────────────────────

def write_cloning_construct(
    design_result: dict,
    insert_seq: str,
    output_path: Path | str,
    gene_name: str = "Insert",
    vector_dna_path: str | Path | None = None,
) -> Path:
    """Write a SnapGene .dna file for the RE cloning result.

    When ``vector_dna_path`` is given, performs in-silico cloning:
    parses the base vector, replaces the MCS region between the two RE
    sites with the insert, and writes a circular recombinant plasmid
    with all original vector features preserved + new CDS annotation.

    Parameters
    ----------
    design_result : dict
        Output from ``RestrictionCloningDesigner.design()``.
    insert_seq : str
        Original insert CDS sequence (same as passed to design()).
    output_path : Path or str
        Output ``.dna`` file path.
    gene_name : str
        Gene name for CDS annotation.
    vector_dna_path : str, Path, or None
        Path to the base vector ``.dna`` file. When provided, generates
        a full recombinant plasmid. When None, generates PCR product only.

    Returns
    -------
    Path : Written file path.
    """
    output_path = Path(output_path)
    insert_seq = insert_seq.upper().replace(" ", "")

    if vector_dna_path and Path(vector_dna_path).exists():
        return _write_recombinant_plasmid(
            design_result, insert_seq, output_path, gene_name,
            Path(vector_dna_path),
        )
    else:
        return _write_pcr_product(
            design_result, insert_seq, output_path, gene_name,
        )


def _write_recombinant_plasmid(
    dr: dict,
    insert_seq: str,
    output_path: Path,
    gene_name: str,
    vector_path: Path,
) -> Path:
    """In-silico cloning: insert CDS into base vector."""

    # ── 1. Parse base vector ──────────────────────────────────────────────
    vec_seq, vec_circular, _ = parse_snapgene(str(vector_path))
    vec_features = _parse_vector_features(vector_path)

    re5_site = dr["re_5prime_site"]
    re3_site = dr["re_3prime_site"]
    sp5 = dr.get("spacer_5prime", "")
    sp3 = dr.get("spacer_3prime", "")
    include_stop = dr.get("include_stop_codon", False)
    stop_codon = "TAA"  # default

    # ── 2. Find RE sites in vector ────────────────────────────────────────
    re5_positions = _find_site(vec_seq, re5_site)
    re3_positions = _find_site(vec_seq, re3_site)

    if not re5_positions:
        raise ValueError(f"{dr['re_5prime']} site ({re5_site}) not found in vector")
    if not re3_positions:
        raise ValueError(f"{dr['re_3prime']} site ({re3_site}) not found in vector")

    # Use first occurrence (should be in MCS)
    re5_pos = re5_positions[0]
    re3_pos = re3_positions[0]

    # ── 3. Determine insert orientation and gap ───────────────────────────
    # In pET vectors: 5' RE (promoter-proximal) has HIGHER position
    # because the gene is on the reverse complement strand.
    # The insert's forward CDS strand becomes the bottom strand.

    if re5_pos > re3_pos:
        # Gene on reverse strand (typical pET)
        gap_start = re3_pos + len(re3_site)   # after RE3 site
        gap_end = re5_pos                      # before RE5 site
        is_reverse = True
    else:
        # Gene on forward strand
        gap_start = re5_pos + len(re5_site)   # after RE5 site
        gap_end = re3_pos                      # before RE3 site
        is_reverse = False

    # Build the insert region between the RE sites
    stop_extra = stop_codon if include_stop else ""
    rc_sp3 = str(Seq(sp3).reverse_complement()) if sp3 else ""
    forward_insert = sp5 + insert_seq + stop_extra + rc_sp3

    if is_reverse:
        insert_on_top = str(Seq(forward_insert).reverse_complement())
    else:
        insert_on_top = forward_insert

    # ── 4. Build recombinant sequence ─────────────────────────────────────
    new_seq = vec_seq[:gap_start] + insert_on_top + vec_seq[gap_end:]
    delta = len(insert_on_top) - (gap_end - gap_start)

    # ── 5. Adjust vector features ─────────────────────────────────────────
    new_features = []
    for feat in vec_features:
        fs, fe = feat["start"], feat["end"]

        # Feature entirely before gap: keep as-is
        if fe < gap_start:
            new_features.append(feat)
        # Feature entirely after gap: shift
        elif fs >= gap_end:
            shifted = dict(feat)
            shifted["start"] = fs + delta
            shifted["end"] = fe + delta
            new_features.append(shifted)
        # Feature overlaps gap: skip (part of replaced MCS)
        # But keep features that only touch the RE sites themselves
        elif fe < gap_start or fs >= gap_end:
            new_features.append(feat)
        # else: overlaps the replaced region → remove

    # ── 6. Build full ORF CDS (vector ATG → tags → insert → stop) ──────
    if is_reverse:
        cds_offset = len(sp3) + (3 if include_stop else 0)
        cds_start = gap_start + cds_offset
        cds_end = cds_start + len(insert_seq) - 1
        cds_strand = 2
    else:
        cds_offset = len(sp5)
        cds_start = gap_start + cds_offset
        cds_end = cds_start + len(insert_seq) - 1
        cds_strand = 1

    insert_protein = _translate(insert_seq).rstrip("*")

    insert_has_stop = (
        insert_seq[-3:].upper() in ("TAA", "TAG", "TGA")
        or include_stop
    )

    # ── 6.1. Remove non-functional downstream CDS features ────────────
    if insert_has_stop:
        _downstream_range = 500
        if is_reverse:
            new_features = [
                f for f in new_features
                if not (
                    f.get("type") == "CDS"
                    and f.get("strand") == cds_strand
                    and f["end"] < cds_start
                    and f["start"] >= cds_start - _downstream_range
                )
            ]
        else:
            new_features = [
                f for f in new_features
                if not (
                    f.get("type") == "CDS"
                    and f.get("strand") == cds_strand
                    and f["start"] > cds_end
                    and f["end"] <= cds_end + _downstream_range
                )
            ]

    # ── 6.2. Keep insert CDS as-is; preserve vector N-terminal tag features ──
    # Do NOT extend the CDS back to the vector's ATG. Instead, annotate the
    # insert gene as its own CDS (with its own start codon). The vector's
    # N-terminal tag CDS features (His6, Thrombin, T7-tag, ATG) remain as
    # separate annotations, giving a cleaner SnapGene display.
    frame_check = dr.get("frame_check", {})
    protein = insert_protein

    # ── 6.3. Extend C-terminal if no stop codon ──────────────────────
    if not insert_has_stop:
        in_frame_3 = frame_check.get("in_frame_3prime", False)
        linker_3_aa = frame_check.get("linker_3prime_aa", "")
        linker_3_nt = frame_check.get("linker_3prime_nt", "")
        if in_frame_3 and linker_3_aa and linker_3_nt:
            protein = protein + linker_3_aa.rstrip("*")
            if is_reverse:
                cds_start = cds_start - len(sp3) - len(linker_3_nt)
            else:
                cds_end = cds_end + len(sp3) + len(linker_3_nt)

    protein = protein.rstrip("*")

    new_features.append({
        "name": gene_name,
        "type": "CDS",
        "start": cds_start,
        "end": cds_end,
        "strand": cds_strand,
        "color": "#66bb6a",
        "translated": True,
        "translation": protein,
    })

    # ── 7. Write .dna file ────────────────────────────────────────────────
    file_data = _build_cookie()
    file_data += _build_sequence_block(new_seq, circular=True)
    file_data += _build_features_block(new_features)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_bytes(file_data)

    return output_path


def _write_pcr_product(
    dr: dict,
    insert_seq: str,
    output_path: Path,
    gene_name: str,
) -> Path:
    """Fallback: write PCR product (linear) when no vector .dna is available."""

    f_tail = dr["f_tail"]
    r_tail = dr["r_tail"]
    r_tail_rc = str(Seq(r_tail).reverse_complement())

    pcr_product = f_tail + insert_seq + r_tail_rc
    total_len = len(pcr_product)

    f_tail_len = len(f_tail)
    prot5_len = len(dr["protection_5"])
    re5_site_len = len(dr["re_5prime_site"])
    prot3_len = len(dr["protection_3"])
    re3_site_len = len(dr["re_3prime_site"])

    features = []

    # 5' protection
    if prot5_len > 0:
        features.append({
            "name": "5' protection",
            "type": "misc_feature",
            "start": 0,
            "end": prot5_len - 1,
            "strand": 1,
            "color": "#b0bec5",
        })

    # 5' RE site
    re5_start = prot5_len
    features.append({
        "name": dr["re_5prime"],
        "type": "misc_binding",
        "start": re5_start,
        "end": re5_start + re5_site_len - 1,
        "strand": 1,
        "color": "#ef5350",
    })

    # Insert CDS
    insert_start = f_tail_len
    insert_end = f_tail_len + len(insert_seq) - 1
    protein = _translate(insert_seq).rstrip("*")
    features.append({
        "name": gene_name,
        "type": "CDS",
        "start": insert_start,
        "end": insert_end,
        "strand": 1,
        "color": "#66bb6a",
        "translated": True,
        "translation": protein,
    })

    # 3' RE site (reverse strand)
    re3_end_pos = total_len - prot3_len - 1
    re3_start_pos = re3_end_pos - re3_site_len + 1
    features.append({
        "name": dr["re_3prime"],
        "type": "misc_binding",
        "start": re3_start_pos,
        "end": re3_end_pos,
        "strand": 2,
        "color": "#ef5350",
    })

    # 3' protection
    if prot3_len > 0:
        features.append({
            "name": "3' protection",
            "type": "misc_feature",
            "start": total_len - prot3_len,
            "end": total_len - 1,
            "strand": 2,
            "color": "#b0bec5",
        })

    # Primer binding sites
    f_ann_len = dr["f_ann_len"]
    features.append({
        "name": f"F_{gene_name}",
        "type": "primer_bind",
        "start": 0,
        "end": f_tail_len + f_ann_len - 1,
        "strand": 1,
        "color": "#42a5f5",
    })

    r_ann_len = dr["r_ann_len"]
    features.append({
        "name": f"R_{gene_name}",
        "type": "primer_bind",
        "start": total_len - len(dr["r_tail"]) - r_ann_len,
        "end": total_len - 1,
        "strand": 2,
        "color": "#42a5f5",
    })

    file_data = _build_cookie()
    file_data += _build_sequence_block(pcr_product, circular=False)
    file_data += _build_features_block(features)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_bytes(file_data)

    return output_path
