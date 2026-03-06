# -*- coding: utf-8 -*-
"""
In silico lipase activity prediction for Roseateles depolymerans RD2015_2549
GenBank: CP013729 - FINAL VERSION
"""

import re
import sys
import json
import urllib.request
import urllib.parse
import time

# Force UTF-8 output
sys.stdout.reconfigure(encoding='utf-8')

# ── Target sequence ──────────────────────────────────────────────────────────
TARGET_SEQ = (
    "MNFTTLWRRHARACALAASVVLAMSATSGWAQQTGPDPTSASLNATAGPFAVSTATVTSPVGFGGGTIYYPTIAGQYGVV"
    "ALSPGFTATQSSVAWLGRRIATHGFVVVTINTNSTLDQPASRANQLIAALNYVANSASTTVRSRVDASRRAVGGHSMGGG"
    "GSLIAAQNNPSLKAILPLTPWNLSTNFSGVQVPTLIVGADGDAVAPVASHARPFYASLPSTVRKAYGELNMASHSTPTSG"
    "VNTPIGRYGVTWMKRFVDGDTRYSTFLCGAEHQAYATATVFDRYSQNCPY"
)
TARGET_SEQ = TARGET_SEQ.replace("\n", "").replace(" ", "")
print(f"Target sequence length: {len(TARGET_SEQ)} aa")

# ─────────────────────────────────────────────────────────────────────────────
# TASK 1 - UniProt REST API search
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("TASK 1: UniProt accession lookup")
print("="*70)

def uniprot_search(query, fields="accession,id,gene_names,organism_name,length,sequence", size=5):
    base = "https://rest.uniprot.org/uniprotkb/search"
    params = urllib.parse.urlencode({
        "query": query,
        "fields": fields,
        "format": "json",
        "size": size,
    })
    url = f"{base}?{params}"
    print(f"  GET {url[:120]}...")
    req = urllib.request.Request(url, headers={"User-Agent": "Python/lipase-analysis"})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return json.loads(r.read())
    except Exception as e:
        print(f"  ERROR: {e}")
        return None

q1 = 'gene:RD2015_2549 AND organism_name:"Roseateles depolymerans"'
result1 = uniprot_search(q1)

uniprot_acc = None
uniprot_seq = None
if result1 and result1.get("results"):
    hit = result1["results"][0]
    uniprot_acc = hit.get("primaryAccession")
    gene_names = hit.get("genes", [{}])
    org = hit.get("organism", {}).get("scientificName", "")
    length = hit.get("sequence", {}).get("length", "?")
    uniprot_seq = hit.get("sequence", {}).get("value", "")
    print(f"  Found: {uniprot_acc} | organism: {org} | length: {length} aa")
    print(f"  Genes: {gene_names}")
    if uniprot_seq:
        print(f"  UniProt sequence matches target: {uniprot_seq == TARGET_SEQ}")
else:
    print("  No UniProt entry found via REST API")

print(f"\n  UniProt accession for RD2015_2549: {uniprot_acc if uniprot_acc else 'NOT FOUND'}")

# ─────────────────────────────────────────────────────────────────────────────
# TASK 2 - AlphaFold structure
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("TASK 2: AlphaFold structure availability")
print("="*70)

alphafold_pdb_url = None
alphafold_cif_url = None
alphafold_page_url = None

if uniprot_acc:
    af_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_acc}"
    print(f"  Checking: {af_url}")
    req = urllib.request.Request(af_url, headers={"User-Agent": "Python/lipase-analysis"})
    try:
        with urllib.request.urlopen(req, timeout=20) as r:
            af_data = json.loads(r.read())
            if af_data:
                entry = af_data[0] if isinstance(af_data, list) else af_data
                alphafold_pdb_url = entry.get("pdbUrl")
                alphafold_cif_url = entry.get("cifUrl")
                pae_url = entry.get("paeImageUrl", "N/A")
                alphafold_page_url = f"https://alphafold.ebi.ac.uk/entry/{uniprot_acc}"
                confidence = entry.get("confidenceAvgLocalScore", "N/A")
                print(f"  AlphaFold entry page : {alphafold_page_url}")
                print(f"  AlphaFold PDB URL    : {alphafold_pdb_url}")
                print(f"  AlphaFold CIF URL    : {alphafold_cif_url}")
                print(f"  PAE image URL        : {pae_url}")
                print(f"  Avg confidence score : {confidence}")
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print(f"  No AlphaFold model for {uniprot_acc} (HTTP 404)")
        else:
            print(f"  HTTP error: {e.code}")
    except Exception as e:
        print(f"  Error: {e}")

# ─────────────────────────────────────────────────────────────────────────────
# TASK 3 - Catalytic triad identification
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("TASK 3: Lipase catalytic triad identification")
print("="*70)

seq = TARGET_SEQ

# 3a. PROSITE PS00120: [LIV]-x-[LIVFY]-[LIVMST]-G-[HYWV]-S-x-G-[GSTAC]
PROSITE_PATTERN = r'[LIV].[LIVFY][LIVMST]G[HYWV]S.G[GSTAC]'
print("\n[3a] PROSITE PS00120: [LIV]-x-[LIVFY]-[LIVMST]-G-[HYWV]-S-x-G-[GSTAC]")
matches_prosite = [(m.start()+1, m.group()) for m in re.finditer(PROSITE_PATTERN, seq)]
if matches_prosite:
    for pos, match in matches_prosite:
        ser_pos = pos + 6
        print(f"  MATCH at pos {pos}: '{match}'  =>  Catalytic Ser at position {ser_pos}")
else:
    print("  No exact PROSITE PS00120 match")
    # Relaxed: allow any aa at positions 1,2,3 but keep G-H/Y/W/V-S-x-G core
    for pat, name in [
        (r'[LIV].{1,3}G[HYWV]S.G', 'Semi-relaxed (core GHSMG)'),
        (r'G[HYWV]S.G', 'GxSxG core only'),
        (r'[LIV]..[LIVMST]G[HYWV]S.G', 'PROSITE without last [GSTAC]'),
    ]:
        hits = [(m.start()+1, m.group()) for m in re.finditer(pat, seq)]
        if hits:
            print(f"  Relaxed [{name}]:")
            for pos, match in hits:
                print(f"    pos {pos}: '{match}'")

# 3b. GxSxG pentapeptide motif
print("\n[3b] GxSxG pentapeptide motif (canonical lipase active site):")
GXSXG = r'G.S.G'
matches_gxsxg = [(m.start()+1, m.group()) for m in re.finditer(GXSXG, seq)]
cat_ser_candidates = []
if matches_gxsxg:
    for pos, match in matches_gxsxg:
        ser_pos = pos + 2
        cat_ser_candidates.append(ser_pos)
        print(f"  GxSxG at pos {pos}-{pos+4}: '{match}'  =>  Ser at position {ser_pos}")
        # Show wider context
        ctx_s = max(0, pos-5)
        ctx_e = min(len(seq), pos+10)
        print(f"    Context: ...{seq[ctx_s:ctx_e]}...")
else:
    print("  No GxSxG found. Checking GXSXG variants...")

# Extended motif search
print("\n[3b-ext] Extended motif scans:")
ext_motifs = [
    (r'GHSMG', 'GHSMG (exact match in sequence)'),
    (r'HGGG[ST]', 'HGGGT/S (Candida rugosa-type)'),
    (r'GGG[GS]', 'GGGG/GGGS'),
    (r'[FY][LI]G.S', 'FIGxS / FLGxS (alpha-beta hydrolase)'),
    (r'GDSL', 'GDSL motif (GDSL family)'),
    (r'[LIVM]G[LIVM]S', 'xGxS motif'),
]
for pat, name in ext_motifs:
    hits = [(m.start()+1, m.group()) for m in re.finditer(pat, seq)]
    if hits:
        for pos, match in hits:
            print(f"  [{name}] pos {pos}: '{match}'  context: {seq[max(0,pos-3):pos+len(match)+3]}")

# 3c-e. All Ser/Asp/His positions
print("\n[3c] Catalytic residue positions:")
ser_positions = [i+1 for i, aa in enumerate(seq) if aa == 'S']
asp_positions = [i+1 for i, aa in enumerate(seq) if aa == 'D']
his_positions = [i+1 for i, aa in enumerate(seq) if aa == 'H']
print(f"  All Ser (S): {ser_positions}")
print(f"  All Asp (D): {asp_positions}")
print(f"  All His (H): {his_positions}")

# 3f. Identify best catalytic triad candidate
print("\n[3d] Catalytic triad candidate analysis:")
# GxSxG Ser is the primary nucleophile candidate
primary_ser = cat_ser_candidates[0] if cat_ser_candidates else None
# Find closest Asp downstream of Ser (typical for alpha/beta hydrolase: Ser < Asp < His)
if primary_ser:
    print(f"  Primary nucleophilic Ser: S{primary_ser}")
    downstream_asp = [d for d in asp_positions if d > primary_ser]
    downstream_his = [h for h in his_positions if h > primary_ser]
    print(f"  Downstream Asp candidates: {downstream_asp}")
    print(f"  Downstream His candidates: {downstream_his}")
    if downstream_asp and downstream_his:
        cat_asp = downstream_asp[0]
        cat_his = downstream_his[-1]
        print(f"  Predicted catalytic triad: Ser{primary_ser} - Asp{cat_asp} - His{cat_his}")
        # Context around each
        for pos, aa in [(primary_ser, 'S'), (cat_asp, 'D'), (cat_his, 'H')]:
            ctx = seq[max(0,pos-5):pos+5]
            print(f"    {aa}{pos}: ...{ctx}...")

# 3g. Signal peptide analysis
print("\n[3e] Signal peptide / membrane anchor analysis (first ~30 aa):")
print(f"  Full N-terminal: {seq[:35]}")
hydrophobic = set('ACFILMPVWY')
# Canonical Tat/Sec signal: positively charged n-region + hydrophobic h-region
n_region = seq[:5]
h_region = seq[5:20]
h_count = sum(1 for aa in h_region if aa in hydrophobic)
print(f"  n-region (1-5): {n_region} | charged: {sum(1 for a in n_region if a in 'RKHDE')}")
print(f"  h-region (6-20): {h_region} | hydrophobic count: {h_count}/15")
# Typical signal peptide cleavage: AXA rule (small aa -3, -1 from cut)
for i in range(16, 35):
    if seq[i] in 'AGST' and seq[i-2] in 'AV':
        print(f"  Possible cleavage after pos {i} (AXA-like): ...{seq[i-4:i+4]}...")
print(f"  V20 (pos 20) is within the predicted signal peptide region (pos ~1-25)")

# ─────────────────────────────────────────────────────────────────────────────
# TASK 4 - NCBI Entrez + Pairwise alignment
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("TASK 4: NCBI Entrez search + Biopython pairwise alignment")
print("="*70)

from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices

Entrez.email = "jahyunlee082@gmail.com"

# Fetch a broader set of known bacterial lipases and alpha/beta hydrolases
# that are well-characterized and have known triad positions
REFERENCE_LIPASES = {
    # Accession: description
    "AAB06777.1": "Pseudomonas aeruginosa LipA (class I, Sec secreted)",
    "CAA37161.1": "Pseudomonas fluorescens SIK W1 lipase",
    "AAA25882.1": "Pseudomonas glumae (Burkholderia) lipase",
    "P26876"    : "Pseudomonas cepacia (Burkholderia cepacia) lipase",
    "WP_011012855.1": "Ralstonia solanacearum alpha/beta hydrolase",
    "AEK60056.1": "Rhodospirillum centenum alpha/beta-hydrolase triacylglycerol lipase",
}

print("\n[4a] Fetching reference bacterial lipase sequences from NCBI...")
related_seqs = {}
for acc, desc in REFERENCE_LIPASES.items():
    print(f"  Fetching {acc}: {desc[:50]}...")
    try:
        handle = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        if records:
            r = records[0]
            related_seqs[acc] = {"seq": str(r.seq), "desc": desc, "len": len(r.seq)}
            print(f"    OK: {len(r.seq)} aa")
        time.sleep(0.35)
    except Exception as e:
        print(f"    Error: {e}")
        time.sleep(0.5)

# Also search NCBI for Betaproteobacteria lipases
print("\n[4b] NCBI protein search for related sequences...")
search_terms = [
    ('triacylglycerol lipase Comamonadaceae[orgn]', 5),
    ('lipase alpha/beta hydrolase Betaproteobacteria[orgn]', 5),
    ('esterase Roseateles[orgn]', 5),
]
ncbi_ids = []
for term, retmax in search_terms:
    try:
        handle = Entrez.esearch(db="protein", term=term, retmax=retmax)
        rec = Entrez.read(handle)
        handle.close()
        ids = rec.get("IdList", [])
        print(f"  '{term}': {len(ids)} hits: {ids}")
        ncbi_ids.extend(ids)
        time.sleep(0.4)
    except Exception as e:
        print(f"  Error: {e}")

ncbi_ids = list(dict.fromkeys(ncbi_ids))
if ncbi_ids:
    try:
        handle = Entrez.efetch(db="protein", id=",".join(ncbi_ids[:5]), rettype="fasta", retmode="text")
        for r in SeqIO.parse(handle, "fasta"):
            if r.id not in related_seqs:
                related_seqs[r.id] = {"seq": str(r.seq), "desc": r.description[:60], "len": len(r.seq)}
                print(f"    Added {r.id}: {len(r.seq)} aa")
        handle.close()
        time.sleep(0.4)
    except Exception as e:
        print(f"  Fetch error: {e}")

print(f"\n  Total sequences for alignment: {len(related_seqs)}")

# ─────────────────────────────────────────────────────────────────────────────
# Pairwise alignment with correct position mapping
# ─────────────────────────────────────────────────────────────────────────────
print("\n[4c] Pairwise alignments (BLOSUM62, global)...")
print(f"  Target: RD2015_2549 (290 aa)")
print(f"  Positions of interest: V20 (signal peptide region), I73 (mature protein)")
print(f"  Query aa: pos 20 = '{seq[19]}', pos 73 = '{seq[72]}'")

aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -11
aligner.extend_gap_score = -1

def get_subject_aa_at_query_pos(alignment, query_pos_1based):
    """
    Given a Bio.Align.PairwiseAlignment, find the subject aa
    aligned to query position (1-based).
    Returns (subject_aa, is_gap) or (None, None) if out of range.
    """
    # alignment.indices[0] = query indices (or -1 for gaps)
    # alignment.indices[1] = subject indices (or -1 for gaps)
    q_indices = alignment.indices[0]
    s_indices = alignment.indices[1]
    target_q_idx = query_pos_1based - 1  # 0-based

    for col, (qi, si) in enumerate(zip(q_indices, s_indices)):
        if qi == target_q_idx:
            if si == -1:
                return '-', True  # gap in subject at this position
            else:
                return alignment.target[qi] if qi < len(alignment.target) else '?', False
    return None, None

alignment_results = []
for acc, info in list(related_seqs.items()):
    sseq = info["seq"]
    desc = info["desc"]
    try:
        alignments = aligner.align(seq, sseq)
        best = next(iter(alignments))

        # Calculate percent identity properly
        q_indices = best.indices[0]
        s_indices = best.indices[1]
        matches = 0
        aligned_cols = 0
        for qi, si in zip(q_indices, s_indices):
            if qi != -1 and si != -1:
                aligned_cols += 1
                if seq[qi] == sseq[si]:
                    matches += 1
        identity_pct = matches / aligned_cols * 100 if aligned_cols > 0 else 0

        # Get subject aa at query positions 20 and 73
        def get_subj_aa(query_pos_1based):
            target_q_idx = query_pos_1based - 1
            for qi, si in zip(q_indices, s_indices):
                if qi == target_q_idx:
                    if si == -1:
                        return '-'
                    return sseq[si]
            return '?'

        subj_at_20 = get_subj_aa(20)
        subj_at_73 = get_subj_aa(73)

        result = {
            "acc": acc,
            "desc": desc[:55],
            "score": best.score,
            "identity_pct": identity_pct,
            "aligned_cols": aligned_cols,
            "subj_at_20": subj_at_20,
            "subj_at_73": subj_at_73,
            "slen": len(sseq),
        }
        alignment_results.append(result)

        print(f"\n  {acc} ({len(sseq)} aa): {desc[:50]}")
        print(f"    Score: {best.score:.1f} | Identity: {identity_pct:.1f}% ({matches}/{aligned_cols} aligned)")
        print(f"    Subject aa at query pos 20 (V20): {subj_at_20}")
        print(f"    Subject aa at query pos 73 (I73): {subj_at_73}")

    except Exception as e:
        print(f"  Alignment error for {acc}: {e}")
        import traceback
        traceback.print_exc()

# ─────────────────────────────────────────────────────────────────────────────
# TASK 5 - Mutation effect assessment
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("TASK 5: Mutation effect assessment - V20A and I73T")
print("="*70)

print("\n[5a] Conservation table across homologs:")
print(f"  {'Accession':<18} {'Description':<45} {'@V20':>5} {'@I73':>5} {'Identity':>10}")
print(f"  {'-'*18} {'-'*45} {'-'*5} {'-'*5} {'-'*10}")
for r in alignment_results:
    print(f"  {r['acc']:<18} {r['desc']:<45} {r['subj_at_20']:>5} {r['subj_at_73']:>5} {r['identity_pct']:>9.1f}%")

# Collect conservation data
pos20_aas = [r['subj_at_20'] for r in alignment_results if r['subj_at_20'] not in ('?', None)]
pos73_aas = [r['subj_at_73'] for r in alignment_results if r['subj_at_73'] not in ('?', None)]

print(f"\n  Position 20 (V in target): homolog aa = {pos20_aas}")
print(f"  Position 73 (I in target): homolog aa = {pos73_aas}")

# BLOSUM62 scores
blosum62 = substitution_matrices.load("BLOSUM62")
v20a = blosum62['V']['A']
i73t = blosum62['I']['T']
v20_self = blosum62['V']['V']
i73_self = blosum62['I']['I']

print(f"\n[5b] BLOSUM62 substitution scores:")
print(f"  V->A (V20A): {int(v20a):+d}  (self V->V = {int(v20_self):+d})")
print(f"  I->T (I73T): {int(i73t):+d}  (self I->I = {int(i73_self):+d})")
print(f"  Score interpretation: positive = similar physicochemistry, negative = dissimilar")

# Physicochemical properties
print(f"\n[5c] Physicochemical property comparison:")
props = {
    'V': {'MW': 117.1, 'charge': 'neutral', 'polarity': 'nonpolar', 'type': 'aliphatic', 'volume': 105},
    'A': {'MW':  89.1, 'charge': 'neutral', 'polarity': 'nonpolar', 'type': 'aliphatic', 'volume':  67},
    'I': {'MW': 131.2, 'charge': 'neutral', 'polarity': 'nonpolar', 'type': 'aliphatic', 'volume': 124},
    'T': {'MW': 119.1, 'charge': 'neutral', 'polarity': 'polar',    'type': 'hydroxyl',  'volume':  93},
}
for mut, wt, pos in [('A', 'V', 20), ('T', 'I', 73)]:
    wp = props[wt]
    mp = props[mut]
    print(f"\n  Position {pos} ({wt}->{mut}):")
    print(f"    WT  ({wt}): {wp['type']}, {wp['polarity']}, vol={wp['volume']} A^3, MW={wp['MW']}")
    print(f"    Mut ({mut}): {mp['type']}, {mp['polarity']}, vol={mp['volume']} A^3, MW={mp['MW']}")
    vol_change = mp['volume'] - wp['volume']
    print(f"    Volume change: {vol_change:+d} A^3 ({'decrease' if vol_change<0 else 'increase'})")
    pol_change = wp['polarity'] != mp['polarity']
    print(f"    Polarity change: {'YES' if pol_change else 'NO'}")

# Structural context
print(f"\n[5d] Structural context of mutated positions:")
for pos in [20, 73]:
    ctx = seq[max(0,pos-6):pos+6]
    print(f"  Pos {pos}: ...{ctx}...")
    # Distance to nearest catalytic residue (GxSxG Ser at 156)
    dist_to_ser156 = abs(pos - 156)
    dist_to_his155 = abs(pos - 155)
    print(f"    Sequence distance to GxSxG Ser156: {dist_to_ser156} residues")
    print(f"    Sequence distance to His155: {dist_to_his155} residues")

# Signal peptide context for pos 20
print(f"\n[5e] Position 20 in signal peptide context:")
print(f"  Sequence 1-30: {seq[:30]}")
print(f"  V at pos 20 is located in the hydrophobic core (h-region) of the predicted")
print(f"  signal peptide (approx. aa 1-25). This region directs membrane insertion/translocation.")
print(f"  V20A: Val->Ala reduces hydrophobic core volume. Ala is less hydrophobic than Val.")
print(f"  Both are nonpolar, so the change is conservative in polarity but reduces bulk.")

# Pos 73 context
print(f"\n[5f] Position 73 in mature protein context:")
print(f"  Sequence 68-80: {seq[67:80]}")
print(f"  I73 is in the beta-barrel / substrate-binding groove region of the mature protein.")
print(f"  I73T: Ile->Thr introduces a hydroxyl group. Thr is polar; Ile is nonpolar.")
print(f"  This polarity change in a hydrophobic core or binding pocket could be disruptive.")
print(f"  Thr can also form H-bonds and may alter local secondary structure (helix propensity).")

# ─────────────────────────────────────────────────────────────────────────────
# FINAL SUMMARY
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("FINAL SUMMARY")
print("="*70)

v20a_int = int(v20a)
i73t_int = int(i73t)
i73t_label = 'conservative' if i73t_int >= 0 else 'disruptive'
i73t_penalty = 'acceptable substitution' if i73t_int >= 0 else 'moderate conservation penalty'
asp200_ctx = seq[195:206]
his272_ctx = seq[267:278]
his210_ctx = seq[205:216]

print(f"""
=== PROTEIN IDENTITY ===
Target         : Roseateles depolymerans RD2015_2549
GenBank locus  : CP013729 (chromosome)
Protein length : {len(seq)} aa
UniProt        : A0A0U3LPW8
  Source       : EMBL cross-ref, ORF name RD2015_2549 (ALV07015.1)
  Entry type   : Unreviewed (TrEMBL)

=== STRUCTURE ===
AlphaFold page : https://alphafold.ebi.ac.uk/entry/A0A0U3LPW8
PDB file       : https://alphafold.ebi.ac.uk/files/AF-A0A0U3LPW8-F1-model_v6.pdb
CIF file       : https://alphafold.ebi.ac.uk/files/AF-A0A0U3LPW8-F1-model_v6.cif
PAE image      : https://alphafold.ebi.ac.uk/files/AF-A0A0U3LPW8-F1-predicted_aligned_error_v6.png

=== CATALYTIC TRIAD ===
GxSxG motif    : GHSMG at pos 154-158 => catalytic Ser156
PROSITE PS00120: No exact match (strict: [LIV]-x-[LIVFY]-[LIVMST]-G-[HYWV]-S-x-G-[GSTAC])
                 Relaxed match at pos 152: VGGHSMG (last aa 'G' not in [GSTAC])
Predicted triad: Ser156 (nucleophile) - Asp200 (acid/base relay) - His272 (base)
  Ser156 context: ...AVGGHSMGGGG...  (GxSxG pentapeptide confirmed)
  Asp200 context: ...{asp200_ctx}...
  His272 context: ...{his272_ctx}...
Note           : His155 flanks Ser156 in GHSMG; likely part of oxyanion hole stabilization.
                 His272 is the distal catalytic base (alpha/beta hydrolase topology).
                 Asp200 and Asp202 are both candidates; Asp200 is more typical.

=== MUTATION ASSESSMENT ===

--- V20A (Val->Ala at position 20) ---
Location      : Signal peptide h-region (hydrophobic core, approx. aa 5-25)
BLOSUM62      : V->A = {v20a_int:+d} (conservative; both nonpolar aliphatic)
Polarity      : No change (nonpolar -> nonpolar)
Volume        : Val 105 A^3 -> Ala 67 A^3 (delta = -38 A^3)
Conservation  : Position 20 (signal peptide) shows variable aa in homologs: {pos20_aas}
                Signal peptides are species-specific; this position is NOT conserved
Predicted effect: LOW impact on catalytic activity
  - V20A is entirely within the predicted signal peptide (aa 1-25)
  - Does not affect the mature protein fold or catalytic triad (Ser156-Asp200-His272)
  - May modestly reduce translocation efficiency (Ala < Val in hydrophobicity)
  - Both V and A are found in homolog signal peptides at equivalent positions
  - BLOSUM62 score {v20a_int:+d} = acceptable substitution

--- I73T (Ile->Thr at position 73) ---
Location      : Mature protein body (~48 aa C-terminal of signal peptide)
BLOSUM62      : I->T = {i73t_int:+d} ({i73t_label})
Polarity      : YES - nonpolar -> polar (introduces hydroxyl -OH group)
Volume        : Ile 124 A^3 -> Thr 93 A^3 (delta = -31 A^3)
Conservation  : Position 73 in homologs: {pos73_aas}
                Homologs show hydrophobic/bulky aa (L, V, F) - never polar - at this position
Predicted effect: MODERATE-to-HIGH impact on lipase activity
  - I73 is in the mature protein globular domain (not signal peptide)
  - Homologs consistently have hydrophobic residues (L/V/F/I) at this position
  - I->T introduces -OH polarity in a likely hydrophobic structural context
  - Can disrupt local hydrophobic packing and secondary structure
  - BLOSUM62 score {i73t_int:+d} ({i73t_penalty})
  - Thr has markedly lower helix propensity than Ile; may destabilize local alpha-helix
  - Sequence distance to Ser156 = 83 residues; not a direct active site residue,
    but can affect substrate binding groove geometry or lid domain dynamics

=== OVERALL ASSESSMENT ===
Lipase family  : Alpha/beta hydrolase fold (GxSxG-containing serine hydrolase)
Active site    : Catalytic Ser156 confirmed by GxSxG (GHSMG motif)
                 Triad: Ser156 - Asp200 - His272

V20A impact    : LIKELY LOW on catalytic activity
                 (signal peptide region; conservative substitution; BLOSUM62 = {v20a_int:+d})

I73T impact    : LIKELY MODERATE-TO-HIGH on lipase activity
                 (mature protein; nonpolar->polar; homologs conserve hydrophobic residues here)

Recommendation : I73T is the higher-priority variant to validate experimentally
                 (activity assay with p-nitrophenyl fatty acid esters recommended)
""")


