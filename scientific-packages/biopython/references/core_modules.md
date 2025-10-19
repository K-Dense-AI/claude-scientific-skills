# BioPython Core Modules Reference

This document provides detailed information about BioPython's core modules and their capabilities.

## Sequence Handling

### Bio.Seq - Sequence Objects

Seq objects are BioPython's fundamental data structure for biological sequences, providing biological methods on top of string-like behavior.

**Creation:**
```python
from Bio.Seq import Seq
my_seq = Seq("AGTACACTGGT")
```

**Key Operations:**
- String methods: `find()`, `count()`, `count_overlap()` (for overlapping patterns)
- Complement/Reverse complement: Returns complementary sequences
- Transcription: DNA → RNA (T → U)
- Back transcription: RNA → DNA
- Translation: DNA/RNA → protein with customizable genetic codes and stop codon handling

**Use Cases:**
- DNA/RNA sequence manipulation
- Converting between nucleic acid types
- Protein translation from coding sequences
- Sequence searching and pattern counting

### Bio.SeqRecord - Sequence Metadata

SeqRecord wraps Seq objects with metadata like ID, description, and features.

**Attributes:**
- `seq`: The sequence itself (Seq object)
- `id`: Unique identifier
- `name`: Short name
- `description`: Longer description
- `features`: List of SeqFeature objects
- `annotations`: Dictionary of additional information
- `letter_annotations`: Per-letter annotations (e.g., quality scores)

### Bio.SeqFeature - Sequence Annotations

Manages sequence annotations and features such as genes, promoters, and coding regions.

**Common Features:**
- Gene locations
- CDS (coding sequences)
- Promoters and regulatory elements
- Exons and introns
- Protein domains

## File Input/Output

### Bio.SeqIO - Sequence File I/O

Unified interface for reading and writing sequence files in multiple formats.

**Supported Formats:**
- FASTA/FASTQ: Standard sequence formats
- GenBank/EMBL: Feature-rich annotation formats
- Clustal/Stockholm/PHYLIP: Alignment formats
- ABI/SFF: Trace and flowgram data
- Swiss-Prot/PIR: Protein databases
- PDB: Protein structure files

**Key Functions:**

**SeqIO.parse()** - Iterator for reading multiple records:
```python
from Bio import SeqIO
for record in SeqIO.parse("file.fasta", "fasta"):
    print(record.id, len(record.seq))
```

**SeqIO.read()** - Read single record:
```python
record = SeqIO.read("file.fasta", "fasta")
```

**SeqIO.write()** - Write sequences:
```python
SeqIO.write(sequences, "output.fasta", "fasta")
```

**SeqIO.convert()** - Direct format conversion:
```python
count = SeqIO.convert("input.gb", "genbank", "output.fasta", "fasta")
```

**SeqIO.index()** - Memory-efficient random access for large files:
```python
record_dict = SeqIO.index("large_file.fasta", "fasta")
sequence = record_dict["seq_id"]
```

**SeqIO.to_dict()** - Load all records into dictionary (memory-based):
```python
record_dict = SeqIO.to_dict(SeqIO.parse("file.fasta", "fasta"))
```

**Common Patterns:**
- Format conversion between FASTA, GenBank, FASTQ
- Filtering sequences by length, ID, or content
- Extracting subsequences
- Batch processing large files with iterators

### Bio.AlignIO - Multiple Sequence Alignment I/O

Handles multiple sequence alignment files.

**Key Functions:**
- `write()`: Save alignments
- `parse()`: Read multiple alignments
- `read()`: Read single alignment
- `convert()`: Convert between formats

**Supported Formats:**
- Clustal
- PHYLIP (sequential and interleaved)
- Stockholm
- NEXUS
- FASTA (aligned)
- MAF (Multiple Alignment Format)

## Sequence Alignment

### Bio.Align - Alignment Tools

**PairwiseAligner** - High-performance pairwise alignment:
```python
from Bio import Align
aligner = Align.PairwiseAligner()
aligner.mode = 'global'  # or 'local'
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.gap_score = -2.5
alignments = aligner.align(seq1, seq2)
```

**CodonAligner** - Codon-aware alignment

**MultipleSeqAlignment** - Container for MSA with column access

### Bio.pairwise2 (Legacy)

Legacy pairwise alignment module with functions like `align.globalxx()`, `align.localxx()`.

## Sequence Analysis Utilities

### Bio.SeqUtils - Sequence Analysis

Collection of utility functions:

**CheckSum** - Calculate sequence checksums (CRC32, CRC64, GCG)

**MeltingTemp** - DNA melting temperature calculations:
- Nearest-neighbor method
- Wallace rule
- GC content method

**IsoelectricPoint** - Protein pI calculation

**ProtParam** - Protein analysis:
- Molecular weight
- Aromaticity
- Instability index
- Secondary structure fractions

**GC/GC_skew** - Calculate GC content and GC skew for sequence windows

### Bio.Data.CodonTable - Genetic Codes

Access to NCBI genetic code tables:
```python
from Bio.Data import CodonTable
standard_table = CodonTable.unambiguous_dna_by_id[1]
print(standard_table.forward_table)  # codon to amino acid
print(standard_table.back_table)     # amino acid to codons
print(standard_table.start_codons)
print(standard_table.stop_codons)
```

**Available codes:**
- Standard code (1)
- Vertebrate mitochondrial (2)
- Yeast mitochondrial (3)
- And many more organism-specific codes

## Sequence Motifs and Patterns

### Bio.motifs - Sequence Motif Analysis

Tools for working with sequence motifs:

**Position Weight Matrices (PWM):**
- Create PWM from aligned sequences
- Calculate information content
- Search sequences for motif matches
- Generate consensus sequences

**Position Specific Scoring Matrices (PSSM):**
- Convert PWM to PSSM
- Score sequences against motifs
- Determine significance thresholds

**Supported Formats:**
- JASPAR
- TRANSFAC
- MEME
- AlignAce

### Bio.Restriction - Restriction Enzymes

Comprehensive restriction enzyme database and analysis:

**Capabilities:**
- Search for restriction sites
- Predict digestion products
- Analyze restriction maps
- Access enzyme properties (recognition site, cut positions, isoschizomers)

**Example usage:**
```python
from Bio import Restriction
from Bio.Seq import Seq

seq = Seq("GAATTC...")
enzyme = Restriction.EcoRI
results = enzyme.search(seq)
```
