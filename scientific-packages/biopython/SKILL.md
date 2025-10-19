---
name: biopython
description: Comprehensive toolkit for computational molecular biology using BioPython. Use this skill when working with biological sequences (DNA, RNA, protein), parsing sequence files (FASTA, GenBank, FASTQ), accessing NCBI databases (Entrez, BLAST), performing sequence alignments, building phylogenetic trees, analyzing protein structures (PDB), or any bioinformatics task requiring BioPython modules.
---

# BioPython

## Overview

BioPython is a comprehensive Python library for computational molecular biology and bioinformatics. This skill provides guidance on using BioPython's extensive modules for sequence manipulation, file I/O, database access, sequence similarity searches, alignments, phylogenetics, structural biology, and population genetics.

## When to Use This Skill

Use this skill when:
- Working with biological sequences (DNA, RNA, protein)
- Reading or writing sequence files (FASTA, GenBank, FASTQ, etc.)
- Accessing NCBI databases (GenBank, PubMed, Protein, Nucleotide)
- Running or parsing BLAST searches
- Performing sequence alignments (pairwise or multiple)
- Building or analyzing phylogenetic trees
- Analyzing protein structures (PDB files)
- Calculating sequence properties (GC content, melting temp, molecular weight)
- Converting between sequence file formats
- Performing population genetics analysis
- Any bioinformatics task requiring BioPython

## Core Capabilities

### 1. Sequence Manipulation

Create and manipulate biological sequences using `Bio.Seq`:

```python
from Bio.Seq import Seq

dna_seq = Seq("ATGGTGCATCTGACT")
rna_seq = dna_seq.transcribe()           # DNA → RNA
protein = dna_seq.translate()             # DNA → Protein
rev_comp = dna_seq.reverse_complement()   # Reverse complement
```

**Common operations:**
- Transcription and back-transcription
- Translation with custom genetic codes
- Complement and reverse complement
- Sequence slicing and concatenation
- Pattern searching and counting

**Reference:** See `references/core_modules.md` (section: Bio.Seq) for detailed operations and examples.

### 2. File Input/Output

Read and write sequence files in multiple formats using `Bio.SeqIO`:

```python
from Bio import SeqIO

# Read sequences
for record in SeqIO.parse("sequences.fasta", "fasta"):
    print(record.id, len(record.seq))

# Write sequences
SeqIO.write(records, "output.gb", "genbank")

# Convert formats
SeqIO.convert("input.fasta", "fasta", "output.gb", "genbank")
```

**Supported formats:** FASTA, FASTQ, GenBank, EMBL, Swiss-Prot, PDB, Clustal, PHYLIP, NEXUS, Stockholm, and many more.

**Common workflows:**
- Format conversion (FASTA ↔ GenBank ↔ FASTQ)
- Filtering sequences by length, ID, or content
- Batch processing large files with iterators
- Random access with `SeqIO.index()` for large files

**Script:** Use `scripts/file_io.py` for file I/O examples and patterns.

**Reference:** See `references/core_modules.md` (section: Bio.SeqIO) for comprehensive format details and workflows.

### 3. NCBI Database Access

Access NCBI databases (GenBank, PubMed, Protein, etc.) using `Bio.Entrez`:

```python
from Bio import Entrez

Entrez.email = "your.email@example.com"  # Required!

# Search database
handle = Entrez.esearch(db="nucleotide", term="human kinase", retmax=100)
record = Entrez.read(handle)
id_list = record["IdList"]

# Fetch sequences
handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta", retmode="text")
records = SeqIO.parse(handle, "fasta")
```

**Key Entrez functions:**
- `esearch()`: Search databases, retrieve IDs
- `efetch()`: Download full records
- `esummary()`: Get document summaries
- `elink()`: Find related records across databases
- `einfo()`: Get database information
- `epost()`: Upload ID lists for large queries

**Important:** Always set `Entrez.email` before using Entrez functions.

**Script:** Use `scripts/ncbi_entrez.py` for complete Entrez workflows including batch downloads and WebEnv usage.

**Reference:** See `references/database_tools.md` (section: Bio.Entrez) for detailed function documentation and parameters.

### 4. BLAST Searches

Run BLAST searches and parse results using `Bio.Blast`:

```python
from Bio.Blast import NCBIWWW, NCBIXML

# Run BLAST online
result_handle = NCBIWWW.qblast("blastn", "nt", sequence)

# Save results
with open("blast_results.xml", "w") as out:
    out.write(result_handle.read())

# Parse results
with open("blast_results.xml") as result_handle:
    blast_record = NCBIXML.read(result_handle)

    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < 0.001:
                print(f"Hit: {alignment.title}")
                print(f"E-value: {hsp.expect}")
                print(f"Identity: {hsp.identities}/{hsp.align_length}")
```

**BLAST programs:** blastn, blastp, blastx, tblastn, tblastx

**Key result attributes:**
- `alignment.title`: Hit description
- `hsp.expect`: E-value
- `hsp.identities`: Number of identical residues
- `hsp.query`, `hsp.match`, `hsp.sbjct`: Aligned sequences

**Script:** Use `scripts/blast_search.py` for complete BLAST workflows including result filtering and extraction.

**Reference:** See `references/database_tools.md` (section: Bio.Blast) for detailed parsing and filtering strategies.

### 5. Sequence Alignment

Perform pairwise and multiple sequence alignments using `Bio.Align`:

**Pairwise alignment:**
```python
from Bio import Align

aligner = Align.PairwiseAligner()
aligner.mode = 'global'  # or 'local'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.gap_score = -2

alignments = aligner.align(seq1, seq2)
print(alignments[0])
print(f"Score: {alignments.score}")
```

**Multiple sequence alignment I/O:**
```python
from Bio import AlignIO

# Read alignment
alignment = AlignIO.read("alignment.clustal", "clustal")

# Write alignment
AlignIO.write(alignment, "output.phylip", "phylip")

# Convert formats
AlignIO.convert("input.clustal", "clustal", "output.fasta", "fasta")
```

**Supported formats:** Clustal, PHYLIP, Stockholm, NEXUS, FASTA, MAF

**Script:** Use `scripts/alignment_phylogeny.py` for alignment examples and workflows.

**Reference:** See `references/core_modules.md` (sections: Bio.Align, Bio.AlignIO) for detailed alignment capabilities.

### 6. Phylogenetic Analysis

Build and analyze phylogenetic trees using `Bio.Phylo`:

```python
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# Read alignment
alignment = AlignIO.read("sequences.fasta", "fasta")

# Calculate distance matrix
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(alignment)

# Build tree (UPGMA or Neighbor-Joining)
constructor = DistanceTreeConstructor(calculator)
tree = constructor.upgma(dm)  # or constructor.nj(dm)

# Visualize tree
Phylo.draw_ascii(tree)
Phylo.draw(tree)  # matplotlib visualization

# Save tree
Phylo.write(tree, "tree.nwk", "newick")
```

**Tree manipulation:**
- `tree.ladderize()`: Sort branches
- `tree.root_at_midpoint()`: Root at midpoint
- `tree.prune()`: Remove taxa
- `tree.collapse_all()`: Collapse short branches
- `tree.distance()`: Calculate distances between clades

**Supported formats:** Newick, NEXUS, PhyloXML, NeXML

**Script:** Use `scripts/alignment_phylogeny.py` for tree construction and manipulation examples.

**Reference:** See `references/specialized_modules.md` (section: Bio.Phylo) for comprehensive tree analysis capabilities.

### 7. Structural Bioinformatics

Analyze protein structures using `Bio.PDB`:

```python
from Bio.PDB import PDBParser, PDBList

# Download structure
pdbl = PDBList()
pdbl.retrieve_pdb_file("1ABC", file_format="pdb", pdir=".")

# Parse structure
parser = PDBParser()
structure = parser.get_structure("protein", "1abc.pdb")

# Navigate hierarchy: Structure → Model → Chain → Residue → Atom
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                print(atom.name, atom.coord)

# Secondary structure with DSSP
from Bio.PDB import DSSP
dssp = DSSP(model, "structure.pdb")

# Structural alignment
from Bio.PDB import Superimposer
sup = Superimposer()
sup.set_atoms(ref_atoms, alt_atoms)
print(f"RMSD: {sup.rms}")
```

**Key capabilities:**
- Parse PDB, mmCIF, MMTF formats
- Secondary structure analysis (DSSP)
- Solvent accessibility calculations
- Structural superimposition
- Distance and angle calculations
- Structure quality validation

**Reference:** See `references/specialized_modules.md` (section: Bio.PDB) for complete structural analysis capabilities.

### 8. Sequence Analysis Utilities

Calculate sequence properties using `Bio.SeqUtils`:

```python
from Bio.SeqUtils import gc_fraction, MeltingTemp as mt
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# DNA analysis
gc = gc_fraction(dna_seq) * 100
tm = mt.Tm_NN(dna_seq)  # Melting temperature

# Protein analysis
protein_analysis = ProteinAnalysis(str(protein_seq))
mw = protein_analysis.molecular_weight()
pi = protein_analysis.isoelectric_point()
aromaticity = protein_analysis.aromaticity()
instability = protein_analysis.instability_index()
```

**Available analyses:**
- GC content and GC skew
- Melting temperature (multiple methods)
- Molecular weight
- Isoelectric point
- Aromaticity
- Instability index
- Secondary structure prediction
- Sequence checksums

**Script:** Use `scripts/sequence_operations.py` for sequence analysis examples.

**Reference:** See `references/core_modules.md` (section: Bio.SeqUtils) for all available utilities.

### 9. Specialized Modules

**Restriction enzymes:**
```python
from Bio import Restriction
enzyme = Restriction.EcoRI
sites = enzyme.search(seq)
```

**Motif analysis:**
```python
from Bio import motifs
m = motifs.create([seq1, seq2, seq3])
pwm = m.counts.normalize(pseudocounts=0.5)
```

**Population genetics:**
Use `Bio.PopGen` for allele frequencies, Hardy-Weinberg equilibrium, FST calculations.

**Clustering:**
Use `Bio.Cluster` for hierarchical clustering, k-means, PCA on biological data.

**Reference:** See `references/core_modules.md` and `references/specialized_modules.md` for specialized module documentation.

## Common Workflows

### Workflow 1: Download and Analyze NCBI Sequences

1. Search NCBI database with `Entrez.esearch()`
2. Fetch sequences with `Entrez.efetch()`
3. Parse with `SeqIO.parse()`
4. Analyze sequences (GC content, translation, etc.)
5. Save results to file

**Script:** Use `scripts/ncbi_entrez.py` for complete implementation.

### Workflow 2: Sequence Similarity Search

1. Run BLAST with `NCBIWWW.qblast()` or parse existing results
2. Parse XML results with `NCBIXML.read()`
3. Filter hits by E-value, identity, coverage
4. Extract and save significant hits
5. Perform downstream analysis

**Script:** Use `scripts/blast_search.py` for complete implementation.

### Workflow 3: Phylogenetic Tree Construction

1. Read multiple sequence alignment with `AlignIO.read()`
2. Calculate distance matrix with `DistanceCalculator`
3. Build tree with `DistanceTreeConstructor` (UPGMA or NJ)
4. Manipulate tree (ladderize, root, prune)
5. Visualize with `Phylo.draw()` or `Phylo.draw_ascii()`
6. Save tree with `Phylo.write()`

**Script:** Use `scripts/alignment_phylogeny.py` for complete implementation.

### Workflow 4: Format Conversion Pipeline

1. Read sequences in original format with `SeqIO.parse()`
2. Filter or modify sequences as needed
3. Write to new format with `SeqIO.write()`
4. Or use `SeqIO.convert()` for direct conversion

**Script:** Use `scripts/file_io.py` for format conversion examples.

## Best Practices

### Email Configuration
Always set `Entrez.email` before using NCBI services:
```python
Entrez.email = "your.email@example.com"
```

### Rate Limiting
Be polite to NCBI servers:
- Use `time.sleep()` between requests
- Use WebEnv for large queries
- Batch downloads in reasonable chunks (100-500 sequences)

### Memory Management
For large files:
- Use iterators (`SeqIO.parse()`) instead of lists
- Use `SeqIO.index()` for random access without loading entire file
- Process in batches when possible

### Error Handling
Always handle potential errors:
```python
try:
    record = SeqIO.read(handle, format)
except Exception as e:
    print(f"Error: {e}")
```

### File Format Selection
Choose appropriate formats:
- FASTA: Simple sequences, no annotations
- GenBank: Rich annotations, features, references
- FASTQ: Sequences with quality scores
- PDB: 3D structural data

## Resources

### scripts/
Executable Python scripts demonstrating common BioPython workflows:

- `sequence_operations.py`: Basic sequence manipulation (transcription, translation, complement, GC content, melting temp)
- `file_io.py`: Reading, writing, and converting sequence files; filtering; indexing large files
- `ncbi_entrez.py`: Searching and downloading from NCBI databases; batch processing with WebEnv
- `blast_search.py`: Running BLAST searches online; parsing and filtering results
- `alignment_phylogeny.py`: Pairwise and multiple sequence alignment; phylogenetic tree construction and manipulation

Run any script with `python3 scripts/<script_name>.py` to see examples.

### references/
Comprehensive reference documentation for BioPython modules:

- `core_modules.md`: Core sequence handling (Seq, SeqRecord, SeqIO, AlignIO, Align, SeqUtils, CodonTable, motifs, Restriction)
- `database_tools.md`: Database access and searches (Entrez, BLAST, SearchIO, BioSQL)
- `specialized_modules.md`: Advanced analyses (PDB, Phylo, PAML, PopGen, Cluster, Graphics)

Reference these files when:
- Learning about specific module capabilities
- Looking up function parameters and options
- Understanding supported file formats
- Finding example code patterns

Use `grep` to search references for specific topics:
```bash
grep -n "secondary structure" references/specialized_modules.md
grep -n "efetch" references/database_tools.md
```

## Additional Resources

**Official Documentation:** https://biopython.org/docs/latest/

**Tutorial:** https://biopython.org/docs/latest/Tutorial/index.html

**API Reference:** https://biopython.org/docs/latest/api/index.html

**Cookbook:** https://biopython.org/wiki/Category:Cookbook
