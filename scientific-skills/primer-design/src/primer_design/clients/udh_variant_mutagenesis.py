#!/usr/bin/env python3
"""
UDH Variant Mutagenesis — iPCR Primer Designer
================================================
Issue #57 Phase 2 variant primer design.

핵심 기능: 6-key position 번호 체계와 실제 단백질 번호 간 자동 변환 + 검증.

6-key position 체계:
  Active site subclustering 정렬 기반 reference positions [75, 165, 166, 170, 171, 174].
  각 parent UDH는 실제 단백질 번호와 일정한 offset이 존재.

Offset (actual_pos = sixkey_pos + offset):
  AtUdh  : -1   (예: 6-key 170 → actual 169)
  PsUdh  : +10  (예: 6-key 170 → actual 180)
  PmUdh  : +10  (예: 6-key 170 → actual 180)
  m3Udh  : +4   (예: 6-key 170 → actual 174)

검증 과정:
  1. Construct .dna 파일에서 UDH CDS를 찾고 단백질로 번역
  2. 6-key WT 잔기가 offset 적용 후 실제 위치에서 일치하는지 확인
  3. 불일치 시 ValueError → numbering 오류 방지
"""

import os
from pathlib import Path
from dataclasses import dataclass
from datetime import datetime

from Bio.Data.CodonTable import standard_dna_table
from Bio.Seq import Seq

from ..subst_primer_mode import iPCRSubstDesigner
from ..del_primer_mode import iPCRDelDesigner
from ..snapgene_parser import parse_snapgene


# ══════════════════════════════════════════════════════════════════════════
#  Registry — Parent UDH 정의
# ══════════════════════════════════════════════════════════════════════════

SIXKEY_POSITIONS = (75, 165, 166, 170, 171, 174)


@dataclass
class UdhParent:
    """Parent UDH construct 정보."""
    name: str
    organism: str
    construct_path: str   # base_dir 기준 상대 경로
    offset: int           # actual_pos = sixkey_pos + offset
    sixkey_wt: dict       # {sixkey_pos: expected_wt_aa}


UDH_REGISTRY = {
    "AtUdh": UdhParent(
        name="AtUdh",
        organism="Agrobacterium tumefaciens C58",
        construct_path="AtUDH primers/pET28a_BamHI_AtUdh_XhoI.dna",
        offset=-1,
        sixkey_wt={75: "V", 165: "C", 166: "T", 170: "N", 171: "N", 174: "M"},
    ),
    "PsUdh": UdhParent(
        name="PsUdh",
        organism="Pseudomonas syringae DC3000",
        construct_path="PsUdh primers/pET28a_BamHI_PsUdh_NotI.dna",
        offset=10,
        sixkey_wt={75: "V", 165: "S", 166: "F", 170: "Q", 171: "N", 174: "M"},
    ),
    "PmUdh": UdhParent(
        name="PmUdh",
        organism="Pseudomonas mendocina",
        construct_path="PmUdh primers/pET28a_PmUDH_BamHI_XhoI.dna",
        offset=10,
        sixkey_wt={75: "V", 165: "S", 166: "F", 170: "A", 171: "N", 174: "M"},
    ),
    "m3Udh": UdhParent(
        name="m3Udh",
        organism="Metagenome-derived #3",
        construct_path="m3Udh primers/pET28a_m3UDH_BamHI_XhoI.dna",
        offset=4,
        sixkey_wt={75: "T", 165: "S", 166: "F", 170: "K", 171: "D", 174: "M"},
    ),
}


# ══════════════════════════════════════════════════════════════════════════
#  SixKeyMapper — 번호 변환 + 검증
# ══════════════════════════════════════════════════════════════════════════

class SixKeyMapper:
    """6-key alignment position <-> actual protein position 변환기.

    모든 변환 시 WT 잔기 검증을 수행하여 numbering 오류를 방지한다.
    """

    def __init__(self, parent_name):
        if parent_name not in UDH_REGISTRY:
            raise ValueError(
                f"Unknown parent: {parent_name}. "
                f"Available: {list(UDH_REGISTRY.keys())}"
            )
        self.parent = UDH_REGISTRY[parent_name]

    def to_actual(self, sixkey_pos):
        """6-key position → actual protein position (1-indexed)."""
        return sixkey_pos + self.parent.offset

    def to_sixkey(self, actual_pos):
        """Actual protein position → 6-key position."""
        return actual_pos - self.parent.offset

    def get_wt_aa(self, sixkey_pos):
        """6-key position의 WT amino acid."""
        wt = self.parent.sixkey_wt.get(sixkey_pos)
        if wt is None:
            raise ValueError(
                f"6-key position {sixkey_pos} not defined for {self.parent.name}. "
                f"Valid: {sorted(self.parent.sixkey_wt.keys())}"
            )
        return wt

    def verify_protein(self, protein_seq):
        """전체 6-key WT 잔기가 단백질 서열에서 일치하는지 검증.

        Returns
        -------
        list[str] : 에러 메시지 리스트 (비어있으면 검증 통과)
        """
        errors = []
        for sixkey_pos, expected_aa in sorted(self.parent.sixkey_wt.items()):
            actual_pos = self.to_actual(sixkey_pos)
            if actual_pos < 1 or actual_pos > len(protein_seq):
                errors.append(
                    f"6-key {sixkey_pos} -> actual {actual_pos}: "
                    f"out of range (protein length: {len(protein_seq)})"
                )
                continue
            actual_aa = protein_seq[actual_pos - 1]
            if actual_aa != expected_aa:
                errors.append(
                    f"6-key {sixkey_pos} -> actual {actual_pos}: "
                    f"found {actual_aa}, expected {expected_aa}"
                )
        return errors

    def format_mutation(self, sixkey_mutations):
        """변이 표기 생성.

        Parameters
        ----------
        sixkey_mutations : dict  {sixkey_pos: new_aa}

        Returns
        -------
        str : e.g., "N170L/N171D" or "T166D"
        """
        parts = []
        for pos in sorted(sixkey_mutations):
            wt = self.get_wt_aa(pos)
            new = sixkey_mutations[pos]
            parts.append(f"{wt}{pos}{new}")
        return "/".join(parts)

    def format_actual_mutation(self, sixkey_mutations):
        """실제 단백질 번호 기반 변이 표기.

        Returns
        -------
        str : e.g., "N169L/N170D" (AtUdh) or "Q180L/N181D" (PsUdh)
        """
        parts = []
        for pos in sorted(sixkey_mutations):
            wt = self.get_wt_aa(pos)
            actual = self.to_actual(pos)
            new = sixkey_mutations[pos]
            parts.append(f"{wt}{actual}{new}")
        return "/".join(parts)


# ══════════════════════════════════════════════════════════════════════════
#  Codon utilities
# ══════════════════════════════════════════════════════════════════════════

# E. coli BL21(DE3) preferred codons
_ECOLI_PREFERRED = {
    'A': 'GCG', 'R': 'CGT', 'N': 'AAC', 'D': 'GAT', 'C': 'TGC',
    'E': 'GAA', 'Q': 'CAG', 'G': 'GGC', 'H': 'CAT', 'I': 'ATT',
    'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTC', 'P': 'CCG',
    'S': 'AGC', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAT', 'V': 'GTG',
}

# Build reverse table: amino acid -> list of codons
_AA_TO_CODONS = {}
for _codon, _aa in standard_dna_table.forward_table.items():
    _AA_TO_CODONS.setdefault(_aa, []).append(_codon)


def pick_min_change_codon(new_aa, old_codon):
    """new_aa를 인코딩하면서 old_codon과 최소 변경되는 코돈 선택.

    Tie-breaking: E. coli preferred codon 우선.
    """
    new_aa = new_aa.upper()
    old_upper = old_codon.upper()

    possible = _AA_TO_CODONS.get(new_aa)
    if not possible:
        raise ValueError(f"No codons for amino acid: {new_aa}")

    preferred = _ECOLI_PREFERRED.get(new_aa, "")

    def score(codon):
        changes = sum(a != b for a, b in zip(codon, old_upper))
        return (changes, codon != preferred)

    return min(possible, key=score)


# ══════════════════════════════════════════════════════════════════════════
#  UdhMutagenesisDesigner
# ══════════════════════════════════════════════════════════════════════════

class UdhMutagenesisDesigner:
    """UDH variant iPCR primer designer.

    6-key 번호 체계로 변이를 지정하면
    자동으로 실제 단백질 번호로 변환 + 검증 후 primer를 설계한다.

    Parameters
    ----------
    base_dir : str
        UDH construct .dna 파일이 위치한 기본 디렉토리.
        예: "C:/Users/.../Uronate dehydrogenase (UDH)/"
    """

    def __init__(self, base_dir):
        self.base_dir = base_dir
        self.subst_designer = iPCRSubstDesigner()
        self.del_designer = iPCRDelDesigner()
        self._cache = {}  # parent_name -> (plasmid_seq, cds_start, cds_end, protein)

    def _load_and_verify(self, parent_name):
        """Construct 로드 + CDS 탐색 + 6-key 검증. 결과 캐싱."""
        if parent_name in self._cache:
            return self._cache[parent_name]

        parent = UDH_REGISTRY[parent_name]
        mapper = SixKeyMapper(parent_name)

        # 1. Parse construct .dna
        filepath = os.path.join(self.base_dir, parent.construct_path)
        if not os.path.exists(filepath):
            raise FileNotFoundError(
                f"Construct file not found: {filepath}\n"
                f"Expected: {parent.construct_path} in {self.base_dir}"
            )

        plasmid_seq, is_circular, features = parse_snapgene(filepath)
        if plasmid_seq is None:
            raise RuntimeError(f"Failed to parse sequence from {filepath}")

        # 2. Find UDH CDS by ORF scanning + 6-key validation
        cds_start, cds_end, protein = self._find_verified_cds(
            plasmid_seq, mapper
        )

        result = (plasmid_seq, cds_start, cds_end, protein)
        self._cache[parent_name] = result

        return result

    @staticmethod
    def _find_verified_cds(plasmid_seq, mapper):
        """ORF scanning으로 UDH CDS를 찾고 6-key 검증 통과하는 것을 반환.

        pET28a construct에서 UDH는 forward strand, ~750-900 bp 범위.
        모든 6-key WT 잔기가 일치해야 반환 → numbering 오류 불가.

        Returns
        -------
        tuple : (cds_start, cds_end, protein_str)
            cds_start, cds_end: 0-indexed, plasmid_seq[start:end] = CDS
        """
        seq_upper = plasmid_seq.upper()
        pos = 0

        while True:
            idx = seq_upper.find("ATG", pos)
            if idx == -1:
                break

            remaining = seq_upper[idx:]
            protein_full = str(Seq(remaining).translate())
            stop = protein_full.find("*")

            if stop > 0:
                cds_len = (stop + 1) * 3
                if 700 < cds_len < 1000:
                    protein = protein_full[:stop]
                    errors = mapper.verify_protein(protein)
                    if not errors:
                        return idx, idx + cds_len, protein

            pos = idx + 3

        raise RuntimeError(
            f"Could not find verified UDH CDS for {mapper.parent.name}. "
            f"Checked all forward-strand ORFs (700-1000 bp) but none passed "
            f"6-key WT residue verification."
        )

    def design_substitution(self, parent_name, sixkey_mutations,
                            target_tm=61.0, overlap_len=18,
                            min_len=18, max_len=35):
        """치환 변이 primer 설계.

        Parameters
        ----------
        parent_name : str
            "AtUdh", "PsUdh", "PmUdh", "m3Udh"
        sixkey_mutations : dict
            {sixkey_pos: new_aa}.
            단일: {170: "E"}, 이중: {170: "L", 171: "D"}
        target_tm : float
        overlap_len : int
        min_len : int
        max_len : int

        Returns
        -------
        dict : primer design result + mapping metadata
        """
        mapper = SixKeyMapper(parent_name)
        plasmid_seq, cds_start, cds_end, protein = self._load_and_verify(
            parent_name
        )

        # Validate mutations
        sorted_sixkeys = sorted(sixkey_mutations.keys())
        actual_positions = []
        for sk in sorted_sixkeys:
            wt_aa = mapper.get_wt_aa(sk)
            new_aa = sixkey_mutations[sk]
            actual = mapper.to_actual(sk)
            actual_positions.append(actual)

            # 실제 위치의 잔기 확인 (이중 검증)
            found = protein[actual - 1]
            if found != wt_aa:
                raise ValueError(
                    f"{parent_name} 6-key {sk} (actual {actual}): "
                    f"expected {wt_aa}, found {found}. Mapping error!"
                )
            if new_aa == wt_aa:
                raise ValueError(
                    f"{parent_name} 6-key {sk}: new_aa '{new_aa}' == WT '{wt_aa}'. "
                    f"Not a mutation."
                )

        # Check adjacency for multi-position mutations
        if len(actual_positions) > 1:
            for i in range(1, len(actual_positions)):
                if actual_positions[i] != actual_positions[i - 1] + 1:
                    raise ValueError(
                        f"Non-adjacent mutations not supported in single primer pair: "
                        f"actual positions {actual_positions}. "
                        f"Design separately or implement sequential mutagenesis."
                    )

        # Build old/new codon sequences
        old_codons = []
        new_codons = []
        for sk in sorted_sixkeys:
            actual = mapper.to_actual(sk)
            codon_offset_in_cds = (actual - 1) * 3
            old_codon = plasmid_seq[cds_start + codon_offset_in_cds:
                                    cds_start + codon_offset_in_cds + 3].upper()
            old_codons.append(old_codon)

            # Verify codon translates to WT
            old_aa = str(Seq(old_codon).translate())
            expected_wt = mapper.get_wt_aa(sk)
            if old_aa != expected_wt:
                raise ValueError(
                    f"{parent_name} actual {actual}: codon {old_codon} -> {old_aa}, "
                    f"expected {expected_wt}. CDS alignment error!"
                )

            new_codon = pick_min_change_codon(sixkey_mutations[sk], old_codon)
            new_codons.append(new_codon)

        old_seq = "".join(old_codons)
        new_seq = "".join(new_codons)

        # Absolute position in plasmid
        first_actual = actual_positions[0]
        abs_pos = cds_start + (first_actual - 1) * 3

        # Design primers
        result = self.subst_designer.design(
            seq=plasmid_seq,
            subst_pos=abs_pos,
            old_seq=old_seq,
            new_seq=new_seq,
            target_tm=target_tm,
            overlap_len=overlap_len,
            min_len=min_len,
            max_len=max_len,
        )

        # Add mapping metadata
        sixkey_label = mapper.format_mutation(sixkey_mutations)
        actual_label = mapper.format_actual_mutation(sixkey_mutations)

        result["parent"] = parent_name
        result["mutation_sixkey"] = sixkey_label
        result["mutation_actual"] = actual_label
        result["sixkey_positions"] = sorted_sixkeys
        result["actual_positions"] = actual_positions
        result["old_codons"] = old_codons
        result["new_codons"] = new_codons
        result["abs_pos"] = abs_pos
        result["cds_start"] = cds_start
        result["cds_end"] = cds_end
        result["plasmid_size"] = len(plasmid_seq)

        return result

    def design_deletion(self, parent_name, sixkey_position,
                        target_tm=61.0, overlap_len=18,
                        min_len=18, max_len=35):
        """아미노산 결실 primer 설계.

        Parameters
        ----------
        parent_name : str
        sixkey_position : int
            결실할 6-key position (단일 아미노산 = 3 bp in-frame deletion)

        Returns
        -------
        dict : primer design result + mapping metadata
        """
        mapper = SixKeyMapper(parent_name)
        plasmid_seq, cds_start, cds_end, protein = self._load_and_verify(
            parent_name
        )

        wt_aa = mapper.get_wt_aa(sixkey_position)
        actual = mapper.to_actual(sixkey_position)

        # Verify
        found = protein[actual - 1]
        if found != wt_aa:
            raise ValueError(
                f"{parent_name} 6-key {sixkey_position} (actual {actual}): "
                f"expected {wt_aa}, found {found}. Mapping error!"
            )

        codon_offset = (actual - 1) * 3
        del_start = cds_start + codon_offset
        del_end = del_start + 3

        result = self.del_designer.design(
            seq=plasmid_seq,
            del_start=del_start,
            del_end=del_end,
            target_tm=target_tm,
            overlap_len=overlap_len,
            min_len=min_len,
            max_len=max_len,
            cds_start=cds_start,
        )

        sixkey_label = f"{wt_aa}{sixkey_position}del"
        actual_label = f"{wt_aa}{actual}del"

        result["parent"] = parent_name
        result["mutation_sixkey"] = sixkey_label
        result["mutation_actual"] = actual_label
        result["sixkey_position"] = sixkey_position
        result["actual_position"] = actual
        result["deleted_codon"] = plasmid_seq[del_start:del_end].upper()
        result["cds_start"] = cds_start
        result["cds_end"] = cds_end
        result["plasmid_size"] = len(plasmid_seq)

        return result


# ══════════════════════════════════════════════════════════════════════════
#  Result Printer
# ══════════════════════════════════════════════════════════════════════════

def print_substitution_result(result):
    """치환 결과 출력."""
    sep = "=" * 70
    line = "-" * 70
    parent = result["parent"]
    sixkey = result["mutation_sixkey"]
    actual = result["mutation_actual"]

    print(sep)
    print(f"  {parent} {sixkey}")
    print(f"  (actual numbering: {actual})")
    print(sep)

    # Codon changes
    for i, (sk, ap) in enumerate(zip(
        result["sixkey_positions"], result["actual_positions"]
    )):
        old_c = result["old_codons"][i]
        new_c = result["new_codons"][i]
        old_aa = str(Seq(old_c).translate())
        new_aa = str(Seq(new_c).translate())
        changes = sum(a != b for a, b in zip(old_c, new_c))
        print(f"  6-key {sk} (actual {ap}): {old_aa} -> {new_aa}  "
              f"codon {old_c} -> {new_c}  ({changes} nt change)")

    print(f"  Plasmid position: {result['abs_pos']} "
          f"(CDS: {result['cds_start']}-{result['cds_end']})")

    # Primers
    f_qc = result["f_qc"]
    r_qc = result["r_qc"]

    print(f"\n{line}")
    f_name = f"iPCR_{parent}_{sixkey.replace('/', '_')}_F"
    r_name = f"iPCR_{parent}_{sixkey.replace('/', '_')}_R"

    print(f"  Forward: {f_name}")
    print(f"    5'-{result['f_full']}-3'")
    print(f"    {result['f_len']} nt, Tm={result['f_tm']}C, GC={result['f_gc']}%")
    print(f"    QC: {f_qc['verdict']}  "
          f"hairpin Tm={f_qc.get('hairpin_tm', 'N/A')}C  "
          f"homodimer Tm={f_qc.get('homodimer_tm', 'N/A')}C")

    print(f"\n  Reverse: {r_name}")
    print(f"    5'-{result['r_full']}-3'")
    print(f"    {result['r_len']} nt, Tm={result['r_tm']}C, GC={result['r_gc']}%")
    print(f"    QC: {r_qc['verdict']}  "
          f"hairpin Tm={r_qc.get('hairpin_tm', 'N/A')}C  "
          f"homodimer Tm={r_qc.get('homodimer_tm', 'N/A')}C")

    print(f"\n  Overlap: {result['overlap_region']} "
          f"({len(result['overlap_region'])} bp)")
    print(f"  Overlap verified: "
          f"{'YES' if result['overlap_verified'] else 'NO'}")
    print(f"  Anneal temp: {result['anneal_temp']}C")

    het = result["het"]
    print(f"  Heterodimer: dG={het.get('dg', 'N/A')} kcal/mol")

    if result["warnings"]:
        print(f"\n  Warnings:")
        for w in result["warnings"]:
            print(f"    - {w}")
    print()


def print_deletion_result(result):
    """결실 결과 출력."""
    sep = "=" * 70

    print(sep)
    print(f"  {result['parent']} {result['mutation_sixkey']}")
    print(f"  (actual numbering: {result['mutation_actual']})")
    print(sep)

    print(f"  Deleted codon: {result['deleted_codon']} "
          f"({str(Seq(result['deleted_codon']).translate())})")
    print(f"  Deletion: {result['del_len']} bp (in-frame)")

    iPCRDelDesigner.print_result(result)


# ══════════════════════════════════════════════════════════════════════════
#  Phase 2 Mutation Definitions (Issue #57)
# ══════════════════════════════════════════════════════════════════════════

PHASE2_SUBSTITUTIONS = [
    # AtUdh Phase 2
    ("AtUdh", {166: "D"}),     # T166D
    ("AtUdh", {166: "E"}),     # T166E
    ("AtUdh", {170: "L", 171: "D"}),  # N170L/N171D

    # PmUdh Phase 2
    ("PmUdh", {170: "K"}),     # A170K
    ("PmUdh", {170: "R"}),     # A170R
    ("PmUdh", {75: "T"}),      # V75T

    # m3Udh Phase 2
    ("m3Udh", {170: "R", 171: "E"}),  # K170R/D171E
    ("m3Udh", {75: "I"}),      # T75I
    ("m3Udh", {75: "V"}),      # T75V
]

PHASE2_DELETIONS = [
    # m3Udh Phase 2
    ("m3Udh", 75),             # T75del (in-frame codon deletion)
]


# ══════════════════════════════════════════════════════════════════════════
#  Report Generator
# ══════════════════════════════════════════════════════════════════════════

def generate_report(all_results, output_path):
    """Markdown report 생성."""
    lines = []
    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    lines.append("# UDH Variant Mutagenesis - iPCR Primer Report")
    lines.append("")
    lines.append(f"**Date**: {now}")
    lines.append(f"**Reference**: Issue #57 Phase 2")
    lines.append(f"**Method**: Inverse PCR (Q5 High-Fidelity) + KLD")
    lines.append("")

    lines.append("---")
    lines.append("")
    lines.append("## Numbering Convention")
    lines.append("")
    lines.append("| Parent | Offset | 6-key 170 | Actual 170 |")
    lines.append("|--------|--------|-----------|------------|")
    for name, parent in UDH_REGISTRY.items():
        sk170 = parent.sixkey_wt.get(170, "?")
        actual_170 = 170 + parent.offset
        lines.append(
            f"| {name} | {parent.offset:+d} | "
            f"{sk170}{170} | {sk170}{actual_170} |"
        )
    lines.append("")

    lines.append("---")
    lines.append("")
    lines.append("## Primer Design Results")
    lines.append("")

    current_parent = None
    for r in all_results:
        parent = r["parent"]
        if parent != current_parent:
            lines.append(f"### {parent} ({UDH_REGISTRY[parent].organism})")
            lines.append("")
            current_parent = parent

        sixkey = r["mutation_sixkey"]
        actual = r["mutation_actual"]
        is_del = "del" in sixkey

        lines.append(f"#### {sixkey} (actual: {actual})")
        lines.append("")

        if is_del:
            lines.append(f"- **Type**: In-frame deletion "
                         f"({r['del_len']} bp)")
            lines.append(f"- **Deleted codon**: "
                         f"`{r.get('deleted_codon', 'N/A')}`")
        else:
            codon_info = []
            for i, sk in enumerate(r.get("sixkey_positions", [])):
                old_c = r["old_codons"][i]
                new_c = r["new_codons"][i]
                changes = sum(a != b for a, b in zip(old_c, new_c))
                codon_info.append(
                    f"`{old_c}` -> `{new_c}` ({changes} nt)")
            lines.append(f"- **Codon change**: {', '.join(codon_info)}")

        f_name = f"iPCR_{parent}_{sixkey.replace('/', '_')}_F"
        r_name = f"iPCR_{parent}_{sixkey.replace('/', '_')}_R"

        lines.append("")
        lines.append("| Primer | Name | Sequence (5'->3') | Length | Tm | QC |")
        lines.append("|--------|------|-------------------|--------|----|----|")
        lines.append(
            f"| F | {f_name} | `{r['f_full']}` | "
            f"{r['f_len']} nt | {r['f_tm']}C | "
            f"{r['f_qc']['verdict']} |"
        )
        lines.append(
            f"| R | {r_name} | `{r['r_full']}` | "
            f"{r['r_len']} nt | {r['r_tm']}C | "
            f"{r['r_qc']['verdict']} |"
        )
        lines.append("")
        lines.append(
            f"- **Anneal temp**: {r['anneal_temp']}C")
        lines.append(
            f"- **Overlap verified**: "
            f"{'YES' if r['overlap_verified'] else 'NO'}")

        if r.get("warnings"):
            lines.append(f"- **Warnings**: {', '.join(r['warnings'])}")
        lines.append("")

    # Order summary
    lines.append("---")
    lines.append("")
    lines.append("## Primer Order Summary")
    lines.append("")
    lines.append("| # | Parent | Mutation | Primer | Sequence (5'->3') | "
                 "Length |")
    lines.append("|---|--------|----------|--------|-------------------|"
                 "--------|")

    idx = 1
    for r in all_results:
        parent = r["parent"]
        sixkey = r["mutation_sixkey"]
        f_name = f"iPCR_{parent}_{sixkey.replace('/', '_')}_F"
        r_name = f"iPCR_{parent}_{sixkey.replace('/', '_')}_R"
        lines.append(
            f"| {idx} | {parent} | {sixkey} | {f_name} | "
            f"`{r['f_full']}` | {r['f_len']} nt |"
        )
        idx += 1
        lines.append(
            f"| {idx} | {parent} | {sixkey} | {r_name} | "
            f"`{r['r_full']}` | {r['r_len']} nt |"
        )
        idx += 1
    lines.append("")

    # Protocol
    lines.append("---")
    lines.append("")
    lines.append("## Protocol")
    lines.append("")
    lines.append("### PCR Conditions (Q5 High-Fidelity)")
    lines.append("")
    lines.append("| Step | Temperature | Time | Cycles |")
    lines.append("|------|-------------|------|--------|")
    lines.append("| Initial denaturation | 98C | 30 sec | 1 |")
    lines.append("| Denaturation | 98C | 10 sec | 25 |")

    anneal_temps = {r["anneal_temp"] for r in all_results}
    anneal_str = "/".join(f"{t}C" for t in sorted(anneal_temps))
    lines.append(f"| Annealing | {anneal_str} | 20 sec | 25 |")

    plasmid_sizes = {r.get("plasmid_size", 5000) for r in all_results}
    ext_time = max(s // 1000 + 1 for s in plasmid_sizes)
    lines.append(f"| Extension | 72C | {ext_time} min | 25 |")
    lines.append("| Final extension | 72C | 2 min | 1 |")
    lines.append("| Hold | 4C | - | - |")
    lines.append("")
    lines.append("### Post-PCR: KLD (NEB #M0554)")
    lines.append("")
    lines.append("```")
    lines.append("1 uL PCR product")
    lines.append("1 uL 10x KLD Enzyme Mix")
    lines.append("5 uL 2x KLD Reaction Buffer")
    lines.append("3 uL nuclease-free water")
    lines.append("-> 25C, 5 min -> Transform 5 uL")
    lines.append("```")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append(f"*Generated by udh_variant_mutagenesis.py on {now}*")

    report = "\n".join(lines)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(report)
    print(f"\nReport saved: {output_path}")
    return report


# ══════════════════════════════════════════════════════════════════════════
#  Main
# ══════════════════════════════════════════════════════════════════════════

def main():
    base_dir = str(
        Path.home() / "OneDrive - 고려대학교" / "저장소"
        / "9. 한미해조류_해수부" / "Uronate dehydrogenase (UDH)"
    )

    designer = UdhMutagenesisDesigner(base_dir)
    all_results = []
    sep = "#" * 70

    # Substitution mutations
    for parent_name, sixkey_muts in PHASE2_SUBSTITUTIONS:
        mapper = SixKeyMapper(parent_name)
        label = mapper.format_mutation(sixkey_muts)

        print(f"\n{sep}")
        print(f"  {parent_name} {label}")
        print(sep)

        try:
            result = designer.design_substitution(parent_name, sixkey_muts)
            print_substitution_result(result)
            all_results.append(result)
        except Exception as e:
            print(f"  ERROR: {e}")

    # Deletion mutations
    for parent_name, sixkey_pos in PHASE2_DELETIONS:
        mapper = SixKeyMapper(parent_name)
        wt = mapper.get_wt_aa(sixkey_pos)

        print(f"\n{sep}")
        print(f"  {parent_name} {wt}{sixkey_pos}del")
        print(sep)

        try:
            result = designer.design_deletion(parent_name, sixkey_pos)
            print_deletion_result(result)
            all_results.append(result)
        except Exception as e:
            print(f"  ERROR: {e}")

    # Generate report
    if all_results:
        output_dir = os.path.dirname(os.path.abspath(__file__))
        report_path = os.path.join(
            output_dir, "..", "..", "..",
            "results", "primer_design",
            "UDH_Phase2_primer_report.md"
        )
        os.makedirs(os.path.dirname(report_path), exist_ok=True)
        generate_report(all_results, report_path)

    print(f"\n{'=' * 70}")
    print(f"  DONE - {len(all_results)} primer pairs designed")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
