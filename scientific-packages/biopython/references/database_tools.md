# BioPython Database Access and Search Tools

This document covers BioPython's capabilities for accessing biological databases and performing sequence searches.

## NCBI Database Access

### Bio.Entrez - NCBI E-utilities Interface

Provides programmatic access to NCBI databases including PubMed, GenBank, Protein, Nucleotide, and more.

**Important:** Always set your email before using Entrez:
```python
from Bio import Entrez
Entrez.email = "your.email@example.com"
```

#### Core Query Functions

**esearch** - Search databases and retrieve IDs:
```python
handle = Entrez.esearch(db="nucleotide", term="Homo sapiens[Organism] AND COX1")
record = Entrez.read(handle)
id_list = record["IdList"]
```

Parameters:
- `db`: Database to search (nucleotide, protein, pubmed, etc.)
- `term`: Search query
- `retmax`: Maximum number of IDs to return
- `sort`: Sort order (relevance, pub_date, etc.)
- `usehistory`: Store results on server (useful for large queries)

**efetch** - Retrieve full records:
```python
handle = Entrez.efetch(db="nucleotide", id="123456", rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
```

Parameters:
- `db`: Database name
- `id`: Single ID or comma-separated list
- `rettype`: Return type (gb, fasta, gp, xml, etc.)
- `retmode`: Return mode (text, xml, asn.1)
- Automatically uses POST for >200 IDs

**elink** - Find related records across databases:
```python
handle = Entrez.elink(dbfrom="protein", db="gene", id="15718680")
result = Entrez.read(handle)
```

Parameters:
- `dbfrom`: Source database
- `db`: Target database
- `id`: ID(s) to link from
- Returns LinkOut providers and relevancy scores

**esummary** - Get document summaries:
```python
handle = Entrez.esummary(db="protein", id="15718680")
summary = Entrez.read(handle)
print(summary[0]['Title'])
```

Returns quick overviews without full records.

**einfo** - Get database statistics:
```python
handle = Entrez.einfo(db="nucleotide")
info = Entrez.read(handle)
```

Provides field indices, term counts, update dates, and available links.

**epost** - Upload ID lists to server:
```python
handle = Entrez.epost("nucleotide", id="123456,789012")
result = Entrez.read(handle)
webenv = result["WebEnv"]
query_key = result["QueryKey"]
```

Useful for large queries split across multiple requests.

**espell** - Get spelling suggestions:
```python
handle = Entrez.espell(term="brest cancer")
result = Entrez.read(handle)
print(result["CorrectedQuery"])  # "breast cancer"
```

**ecitmatch** - Convert citations to PubMed IDs:
```python
citation = "proc natl acad sci u s a|1991|88|3248|mann bj|"
handle = Entrez.ecitmatch(db="pubmed", bdata=citation)
```

#### Data Processing Functions

**Entrez.read()** - Parse XML to Python dictionary:
```python
handle = Entrez.esearch(db="protein", term="insulin")
record = Entrez.read(handle)
```

**Entrez.parse()** - Generator for large XML results:
```python
handle = Entrez.efetch(db="protein", id=id_list, rettype="gp", retmode="xml")
for record in Entrez.parse(handle):
    process(record)
```

#### Common Workflows

**Download sequences by accession:**
```python
handle = Entrez.efetch(db="nucleotide", id="NM_001301717", rettype="fasta", retmode="text")
record = SeqIO.read(handle, "fasta")
```

**Search and download multiple sequences:**
```python
# Search
search_handle = Entrez.esearch(db="nucleotide", term="human kinase", retmax="100")
search_results = Entrez.read(search_handle)

# Download
fetch_handle = Entrez.efetch(db="nucleotide", id=search_results["IdList"], rettype="gb", retmode="text")
for record in SeqIO.parse(fetch_handle, "genbank"):
    print(record.id)
```

**Use WebEnv for large queries:**
```python
# Post IDs
post_handle = Entrez.epost(db="nucleotide", id=",".join(large_id_list))
post_result = Entrez.read(post_handle)

# Fetch in batches
batch_size = 500
for start in range(0, count, batch_size):
    fetch_handle = Entrez.efetch(
        db="nucleotide",
        rettype="fasta",
        retmode="text",
        retstart=start,
        retmax=batch_size,
        webenv=post_result["WebEnv"],
        query_key=post_result["QueryKey"]
    )
    # Process batch
```

### Bio.GenBank - GenBank Format Parsing

Low-level GenBank file parser (SeqIO is usually preferred).

### Bio.SwissProt - Swiss-Prot/UniProt Parsing

Parse Swiss-Prot and UniProtKB flat file format:
```python
from Bio import SwissProt
with open("uniprot.dat") as handle:
    for record in SwissProt.parse(handle):
        print(record.entry_name, record.organism)
```

## Sequence Similarity Searches

### Bio.Blast - BLAST Interface

Tools for running BLAST searches and parsing results.

#### Running BLAST

**NCBI QBLAST (online):**
```python
from Bio.Blast import NCBIWWW
result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
```

Parameters:
- Program: blastn, blastp, blastx, tblastn, tblastx
- Database: nt, nr, refseq_rna, pdb, etc.
- Sequence: string or Seq object
- Additional parameters: `expect`, `word_size`, `hitlist_size`, `format_type`

**Local BLAST:**
Run standalone BLAST from command line, then parse results.

#### Parsing BLAST Results

**XML format (recommended):**
```python
from Bio.Blast import NCBIXML

result_handle = open("blast_results.xml")
blast_records = NCBIXML.parse(result_handle)

for blast_record in blast_records:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < 0.001:
                print(f"Hit: {alignment.title}")
                print(f"Length: {alignment.length}")
                print(f"E-value: {hsp.expect}")
                print(f"Identities: {hsp.identities}/{hsp.align_length}")
```

**Functions:**
- `NCBIXML.read()`: Single query
- `NCBIXML.parse()`: Multiple queries (generator)

**Key Record Attributes:**
- `alignments`: List of matching sequences
- `query`: Query sequence ID
- `query_length`: Length of query

**Alignment Attributes:**
- `title`: Description of hit
- `length`: Length of hit sequence
- `hsps`: High-scoring segment pairs

**HSP Attributes:**
- `expect`: E-value
- `score`: Bit score
- `identities`: Number of identical residues
- `positives`: Number of positive scoring matches
- `gaps`: Number of gaps
- `align_length`: Length of alignment
- `query`: Aligned query sequence
- `match`: Match indicators
- `sbjct`: Aligned subject sequence
- `query_start`, `query_end`: Query coordinates
- `sbjct_start`, `sbjct_end`: Subject coordinates

#### Common BLAST Workflows

**Find homologs:**
```python
result = NCBIWWW.qblast("blastp", "nr", protein_sequence, expect=1e-10)
with open("results.xml", "w") as out:
    out.write(result.read())
```

**Filter results by criteria:**
```python
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < 1e-5 and hsp.identities/hsp.align_length > 0.5:
            # Process high-quality hits
            pass
```

### Bio.SearchIO - Unified Search Results Parser

Modern interface for parsing various search tool outputs (BLAST, HMMER, BLAT, etc.).

**Key Functions:**
- `read()`: Parse single query
- `parse()`: Parse multiple queries (generator)
- `write()`: Write results to file
- `convert()`: Convert between formats

**Supported Tools:**
- BLAST (XML, tabular, plain text)
- HMMER (hmmscan, hmmsearch, phmmer)
- BLAT
- FASTA
- InterProScan
- Exonerate

**Example:**
```python
from Bio import SearchIO
results = SearchIO.parse("blast_output.xml", "blast-xml")
for result in results:
    for hit in result:
        if hit.hsps[0].evalue < 0.001:
            print(hit.id, hit.hsps[0].evalue)
```

## Local Database Management

### BioSQL - SQL Database Interface

Store and manage biological sequences in SQL databases (PostgreSQL, MySQL, SQLite).

**Features:**
- Store SeqRecord objects with annotations
- Efficient querying and retrieval
- Cross-reference sequences
- Track relationships between sequences

**Example:**
```python
from BioSQL import BioSeqDatabase
server = BioSeqDatabase.open_database(driver="MySQLdb", user="user", passwd="pass", host="localhost", db="bioseqdb")
db = server["my_db"]

# Store sequences
db.load(SeqIO.parse("sequences.gb", "genbank"))

# Query
seq = db.lookup(accession="NC_005816")
```
