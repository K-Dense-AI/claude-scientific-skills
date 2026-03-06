"""
In silico lipase activity prediction for Roseateles depolymerans RD2015_2549
GenBank: CP013729
"""

import re
import sys
import json
import urllib.request
import urllib.parse
import time

# ── Target sequence ──────────────────────────────────────────────────────────
TARGET_SEQ = (
    "MNFTTLWRRHARACALAASVVLAMSATSGWAQQTGPDPTSASLNATAGPFAVSTATVTSPVGFGGGTIYYPTIAGQYGVV"
    "ALSPGFTATQSSVAWLGRRIATHGFVVVTINTNSTLDQPASRANQLIAALNYVANSASTTVRSRVDASRRAVGGHSMGGG"
    "GSLIAAQNNPSLKAILPLTPWNLSTNFSGVQVPTLIVGADGDAVAPVASHARPFYASLPSTVRKAYGELNMASHSTPTSG"
    "VNTPIGRYGVTWMKRFVDGDTRYSTFLCGAEHQAYATATVFDRYSQNCPY"
)
# Flatten to single string (remove line breaks if any)
TARGET_SEQ = TARGET_SEQ.replace("\n", "").replace(" ", "")
print(f"Target sequence length: {len(TARGET_SEQ)} aa")

# ─────────────────────────────────────────────────────────────────────────────
# TASK 1 – UniProt REST API search
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("TASK 1: UniProt accession lookup")
print("="*70)

def uniprot_search(query, fields="accession,id,gene_names,organism_name,length", size=5):
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

# Primary search: gene name + organism
q1 = 'gene:RD2015_2549 AND organism_name:"Roseateles depolymerans"'
result1 = uniprot_search(q1)

uniprot_acc = None
if result1 and result1.get("results"):
    hit = result1["results"][0]
    uniprot_acc = hit.get("primaryAccession")
    gene_names = hit.get("genes", [{}])
    org = hit.get("organism", {}).get("scientificName", "")
    length = hit.get("sequence", {}).get("length", "?")
    print(f"  Found: {uniprot_acc} | organism: {org} | length: {length} aa")
    print(f"  Genes: {gene_names}")
else:
    print("  No result for primary query. Trying broader search...")
    # Fallback: locus tag search
    q2 = 'xref:embl-RD2015_2549 OR gene:RD2015_2549'
    result2 = uniprot_search(q2)
    if result2 and result2.get("results"):
        hit = result2["results"][0]
        uniprot_acc = hit.get("primaryAccession")
        print(f"  Fallback found: {uniprot_acc}")
    else:
        print("  No UniProt entry found via REST API (protein may not be in UniProt)")

# Try also searching by organism taxon directly
q3 = 'organism_id:415244 AND (lipase OR esterase)'
result3 = uniprot_search(q3, size=10)
if result3 and result3.get("results"):
    print(f"\n  Roseateles depolymerans lipase/esterase entries in UniProt:")
    for h in result3["results"]:
        acc = h.get("primaryAccession")
        genes = [g.get("geneName", {}).get("value", "") for g in h.get("genes", [])]
        pname = h.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "")
        print(f"    {acc}: {pname} | genes: {genes}")
else:
    # Try without taxon ID
    q3b = 'organism_name:"Roseateles depolymerans"'
    result3b = uniprot_search(q3b, size=10)
    if result3b and result3b.get("results"):
        print(f"\n  All Roseateles depolymerans UniProt entries:")
        for h in result3b["results"]:
            acc = h.get("primaryAccession")
            genes = [g.get("geneName", {}).get("value", "") for g in h.get("genes", [])]
            pname = (h.get("proteinDescription", {}).get("recommendedName", {})
                     .get("fullName", {}).get("value", ""))
            print(f"    {acc}: {pname} | genes: {genes}")
    else:
        print("  No Roseateles depolymerans entries found in UniProt")

print(f"\n  UniProt accession for RD2015_2549: {uniprot_acc if uniprot_acc else 'NOT FOUND'}")

# ─────────────────────────────────────────────────────────────────────────────
# TASK 2 – AlphaFold structure
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("TASK 2: AlphaFold structure availability")
print("="*70)

alphafold_url = None
if uniprot_acc:
    af_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_acc}"
    print(f"  Checking: {af_url}")
    req = urllib.request.Request(af_url, headers={"User-Agent": "Python/lipase-analysis"})
    try:
        with urllib.request.urlopen(req, timeout=20) as r:
            af_data = json.loads(r.read())
            if af_data:
                entry = af_data[0] if isinstance(af_data, list) else af_data
                alphafold_url = entry.get("pdbUrl") or entry.get("cifUrl")
                pae_url = entry.get("paeImageUrl", "N/A")
                print(f"  AlphaFold PDB URL : {entry.get('pdbUrl', 'N/A')}")
                print(f"  AlphaFold CIF URL : {entry.get('cifUrl', 'N/A')}")
                print(f"  PAE image URL     : {pae_url}")
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print(f"  No AlphaFold model for {uniprot_acc} (HTTP 404)")
        else:
            print(f"  HTTP error: {e.code}")
    except Exception as e:
        print(f"  Error: {e}")
else:
    print("  Skipped (no UniProt accession). AlphaFold requires UniProt accession.")
    print("  Note: AlphaFold DB URL pattern: https://alphafold.ebi.ac.uk/entry/<UniProt_ACC>")

# ─────────────────────────────────────────────────────────────────────────────
# TASK 3 – Catalytic triad identification
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("TASK 3: Lipase catalytic triad identification")
print("="*70)

seq = TARGET_SEQ

# 3a. PROSITE PS00120: [LIV]-x-[LIVFY]-[LIVMST]-G-[HYWV]-S-x-G-[GSTAC]
#     Pattern definition (1-indexed positions relative to match)
PROSITE_PATTERN = r'[LIV].[LIVFY][LIVMST]G[HYWV]S.G[GSTAC]'
print("\n[3a] PROSITE PS00120 pattern: [LIV]-x-[LIVFY]-[LIVMST]-G-[HYWV]-S-x-G-[GSTAC]")
matches_prosite = [(m.start()+1, m.group()) for m in re.finditer(PROSITE_PATTERN, seq)]
if matches_prosite:
    for pos, match in matches_prosite:
        ser_pos = pos + 6  # S is 7th residue in pattern (1-indexed)
        print(f"  Match at pos {pos}: '{match}'  →  Catalytic Ser at position {ser_pos}")
else:
    print("  No PROSITE PS00120 match found")
    # Try relaxed version
    PROSITE_RELAXED = r'[LIV].{0,2}G[HYWV]S.G'
    matches_relaxed = [(m.start()+1, m.group()) for m in re.finditer(PROSITE_RELAXED, seq)]
    if matches_relaxed:
        print("  Relaxed pattern matches:")
        for pos, match in matches_relaxed:
            print(f"    pos {pos}: '{match}'")

# 3b. GxSxG motif (pentapeptide lipase motif)
print("\n[3b] GxSxG pentapeptide motif search:")
GXSXG = r'G.S.G'
matches_gxsxg = [(m.start()+1, m.group()) for m in re.finditer(GXSXG, seq)]
if matches_gxsxg:
    for pos, match in matches_gxsxg:
        ser_pos = pos + 2  # S is middle (3rd)
        print(f"  GxSxG at pos {pos}-{pos+4}: '{match}'  →  Ser at position {ser_pos}")
else:
    print("  No GxSxG motif found")

# Also check GxSxAG and variants
print("\n[3b-ext] Extended motifs (GxSxG variants, HGGG, etc.):")
EXT_MOTIFS = [
    (r'G[A-Z]S[A-Z]G', 'GxSxG'),
    (r'HGGG', 'HGGG (oxyanion hole)'),
    (r'GGG[ST]', 'GGGx'),
    (r'[FY][LI]G.S', 'FIGxS alpha/beta hydrolase'),
]
for pat, name in EXT_MOTIFS:
    m_list = [(m.start()+1, m.group()) for m in re.finditer(pat, seq)]
    if m_list:
        for pos, match in m_list:
            print(f"  [{name}] pos {pos}: '{match}'")

# 3c. Serine nucleophile candidates
print("\n[3c] All Ser (S) positions in sequence:")
ser_positions = [i+1 for i, aa in enumerate(seq) if aa == 'S']
print(f"  S positions: {ser_positions}")

# 3d. Conserved Asp candidates (potential catalytic Asp)
print("\n[3d] All Asp (D) positions in sequence:")
asp_positions = [i+1 for i, aa in enumerate(seq) if aa == 'D']
print(f"  D positions: {asp_positions}")

# 3e. Conserved His candidates (potential catalytic His)
print("\n[3e] All His (H) positions in sequence:")
his_positions = [i+1 for i, aa in enumerate(seq) if aa == 'H']
print(f"  H positions: {his_positions}")

# 3f. Alpha/beta hydrolase fold signature: look for GXSXG in context
print("\n[3f] Contextual analysis around Ser residues (± 5 aa):")
gxsxg_sers = [pos+2 for pos, match in matches_gxsxg] if matches_gxsxg else []
for s in ser_positions:
    ctx_start = max(0, s-6)
    ctx_end = min(len(seq), s+5)
    ctx = seq[ctx_start:ctx_end]
    is_gxsxg = s in gxsxg_sers
    flag = " <-- GxSxG motif (CATALYTIC SER CANDIDATE)" if is_gxsxg else ""
    print(f"  S{s}: ...{ctx}...{flag}")

# 3g. Signal peptide / mature protein boundary
print("\n[3g] Signal peptide analysis (first 30 aa):")
print(f"  First 30 aa: {seq[:30]}")
# SignalP-like heuristic: hydrophobic core
hydrophobic = set('ACFILMPVWY')
h_count = sum(1 for aa in seq[1:20] if aa in hydrophobic)
print(f"  Hydrophobic aa in pos 2-20: {h_count}/19")
# Find potential cleavage site (after hydrophobic core)
for i in range(15, 35):
    window = seq[i-3:i]
    if seq[i] in 'AGST' and all(aa in hydrophobic for aa in seq[i-5:i-2]):
        print(f"  Potential SP cleavage near pos {i+1}: ...{seq[i-5:i+3]}...")

print("\n[3g-summary] Predicted mature sequence start:")
# Simple: look for AXA motif pattern near pos 20-30
for i in range(18, 35):
    if seq[i-1] == 'A' and seq[i+1] == 'A':
        print(f"  AXA motif at pos {i}: {seq[i-3:i+4]}")

# ─────────────────────────────────────────────────────────────────────────────
# TASK 4 – NCBI Entrez + Pairwise alignment
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("TASK 4: NCBI Entrez search + Biopython pairwise alignment")
print("="*70)

from Bio import Entrez, SeqIO, pairwise2
from Bio.Align import PairwiseAligner

Entrez.email = "jahyunlee082@gmail.com"

# Search for related lipases/esterases in bacteria
print("\n[4a] Searching NCBI protein database for related lipases...")
search_queries = [
    "lipase[title] Betaproteobacteria[orgn] AND alpha/beta hydrolase[title]",
    "esterase lipase Roseateles[orgn]",
    "alpha/beta hydrolase lipase Comamonadaceae[orgn]",
]

hit_ids = []
for q in search_queries:
    print(f"  Query: {q}")
    try:
        handle = Entrez.esearch(db="protein", term=q, retmax=5)
        record = Entrez.read(handle)
        handle.close()
        ids = record["IdList"]
        print(f"    Found {len(ids)} hits: {ids}")
        hit_ids.extend(ids)
        time.sleep(0.5)
    except Exception as e:
        print(f"    Error: {e}")

# Deduplicate
hit_ids = list(dict.fromkeys(hit_ids))[:8]
print(f"\n  Total unique hits to fetch: {len(hit_ids)}: {hit_ids}")

# Fetch sequences
related_seqs = {}
if hit_ids:
    print("\n[4b] Fetching sequences from NCBI...")
    try:
        handle = Entrez.efetch(db="protein", id=",".join(hit_ids), rettype="fasta", retmode="text")
        for record in SeqIO.parse(handle, "fasta"):
            related_seqs[record.id] = str(record.seq)
            org_hint = record.description[:80]
            print(f"  {record.id}: {len(record.seq)} aa | {org_hint}")
        handle.close()
        time.sleep(0.5)
    except Exception as e:
        print(f"  Fetch error: {e}")

# If not enough, fetch known lipase accessions manually
KNOWN_LIPASES = {
    "AAB06777.1": "Pseudomonas aeruginosa LipA (LipA, class I lipase)",
    "P26876"    : "Pseudomonas cepacia lipase (BC lipase)",
    "CAA40558.1": "Burkholderia glumae lipase",
}
print("\n[4b-ext] Fetching known reference lipases...")
for acc, desc in KNOWN_LIPASES.items():
    print(f"  Fetching {acc}: {desc}")
    try:
        handle = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
        record = list(SeqIO.parse(handle, "fasta"))[0]
        handle.close()
        related_seqs[f"{acc}_{desc.split()[0]}_{desc.split()[1]}"] = str(record.seq)
        print(f"    OK: {len(record.seq)} aa")
        time.sleep(0.3)
    except Exception as e:
        print(f"    Error: {e}")

# ─────────────────────────────────────────────────────────────────────────────
# Pairwise alignment and conservation analysis
# ─────────────────────────────────────────────────────────────────────────────
print("\n[4c] Pairwise alignment vs target (RD2015_2549)...")
print(f"  Query mutations of interest: V20A (position 20), I73T (position 73)")
print(f"  Target seq[19] = {seq[19]} (pos 20), seq[72] = {seq[72]} (pos 73)")

aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.substitution_matrix = None  # use default
# Use built-in BLOSUM62
from Bio.Align import substitution_matrices
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

results = []
for sid, sseq in list(related_seqs.items())[:6]:  # max 6 alignments
    try:
        alignments = aligner.align(seq, sseq)
        best = next(iter(alignments))
        score = best.score
        # Calculate identity
        aligned_pairs = list(best.aligned)
        # Get aligned sequences
        aln_target = str(best).split("\n")[0]
        aln_subject = str(best).split("\n")[2]
        matches = sum(a == b for a, b in zip(aln_target, aln_subject) if a != '-' and b != '-')
        aln_len = len(aln_target.replace('-', '').replace(' ',''))
        identity = matches / max(len(seq), 1) * 100

        # Find equivalent positions for pos 20 and 73 in the alignment
        # Map target position → alignment column
        target_aln = str(best).split("\n")[0]
        subject_aln = str(best).split("\n")[2]

        def get_aligned_aa(target_aln, subject_aln, query_pos):
            """Given 1-based query position, find the subject aa at that column."""
            count = 0
            for col, (ta, sa) in enumerate(zip(target_aln, subject_aln)):
                if ta != '-':
                    count += 1
                    if count == query_pos:
                        return sa, col
            return None, None

        subj_at_20, col_20 = get_aligned_aa(target_aln, subject_aln, 20)
        subj_at_73, col_73 = get_aligned_aa(target_aln, subject_aln, 73)

        results.append({
            "id": sid[:50],
            "score": score,
            "identity_pct": identity,
            "aa_at_20": subj_at_20,
            "aa_at_73": subj_at_73,
        })
        print(f"\n  Alignment: {sid[:50]}")
        print(f"    Score: {score:.1f}, Identity: {identity:.1f}%")
        print(f"    Subject aa at pos 20: {subj_at_20} (target: V)")
        print(f"    Subject aa at pos 73: {subj_at_73} (target: I)")
    except Exception as e:
        print(f"  Alignment error for {sid}: {e}")

# ─────────────────────────────────────────────────────────────────────────────
# TASK 5 – Mutation effect assessment
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*70)
print("TASK 5: Mutation effect assessment – V20A and I73T")
print("="*70)

# Conservation summary
print("\n[5a] Conservation of positions 20 and 73 across homologs:")
print(f"  {'Homolog':<50} {'aa@20':>6} {'aa@73':>6} {'Identity':>10}")
print(f"  {'-'*50} {'-'*6} {'-'*6} {'-'*10}")
for r in results:
    print(f"  {r['id']:<50} {str(r['aa_at_20']):>6} {str(r['aa_at_73']):>6} {r['identity_pct']:>9.1f}%")

# BLOSUM62 scores for mutations
from Bio.Align.substitution_matrices import load as load_matrix
blosum62 = load_matrix("BLOSUM62")

def blosum_score(aa1, aa2):
    try:
        return blosum62[aa1][aa2]
    except Exception:
        try:
            return blosum62[aa2][aa1]
        except:
            return None

v20a_score = blosum_score('V', 'A')
i73t_score = blosum_score('I', 'T')
print(f"\n[5b] BLOSUM62 substitution scores:")
print(f"  V→A (V20A): {v20a_score}  (>0 = conservative, <0 = disruptive)")
print(f"  I→T (I73T): {i73t_score}  (>0 = conservative, <0 = disruptive)")

# Structural context
print(f"\n[5c] Structural context of positions 20 and 73:")
window = 5
for pos, mut, wt in [(20, 'A', 'V'), (73, 'T', 'I')]:
    ctx = seq[max(0, pos-1-window): pos+window]
    print(f"  Position {pos} ({wt}→{mut}): ...{ctx}...")
    print(f"    WT residue: {wt} (hydrophobic, {'small' if wt=='V' else 'branched'})")
    print(f"    Mut residue: {mut} ({'small hydrophobic' if mut=='A' else 'polar hydroxyl'})")

# Check if pos 20/73 are near catalytic residues
print(f"\n[5d] Distance to catalytic residues:")
cat_ser = [pos+2 for pos, match in matches_gxsxg] if matches_gxsxg else []
# Also include PROSITE match Ser positions
cat_ser += [pos + 6 for pos, match in matches_prosite]
cat_ser = sorted(set(cat_ser))
print(f"  Catalytic Ser candidates: {cat_ser}")
print(f"  Asp positions: {asp_positions}")
print(f"  His positions: {his_positions}")

for mut_pos in [20, 73]:
    for cat_pos in cat_ser:
        dist = abs(mut_pos - cat_pos)
        print(f"  |pos {mut_pos} - Ser{cat_pos}| = {dist} residues (sequence distance)")

# Signal peptide context
print(f"\n[5e] Signal peptide context:")
print(f"  Position 20 (V20) is in the predicted signal peptide / membrane anchor region")
print(f"  First 30 aa: {seq[:30]}")
print(f"  Pos 20 context: {seq[15:25]}")

print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"""
Target protein : Roseateles depolymerans RD2015_2549 ({len(seq)} aa)
GenBank        : CP013729

UniProt        : {uniprot_acc if uniprot_acc else 'Not found in UniProt (uncommon organism)'}
AlphaFold      : {'Available at ' + str(alphafold_url) if alphafold_url else 'Not available (requires UniProt accession)'}

Catalytic triad identification:
  PROSITE PS00120 hits : {matches_prosite if matches_prosite else 'None (may need relaxed pattern)'}
  GxSxG motif hits     : {matches_gxsxg if matches_gxsxg else 'None'}
  GxSxG Ser candidates : {[pos+2 for pos, m in matches_gxsxg] if matches_gxsxg else 'None'}
  All Ser positions    : {ser_positions}
  All Asp positions    : {asp_positions}
  All His positions    : {his_positions}

Mutation assessment:
  V20A: BLOSUM62 = {v20a_score} | Position in signal peptide region
  I73T: BLOSUM62 = {i73t_score} | Position in mature protein / substrate binding
""")
