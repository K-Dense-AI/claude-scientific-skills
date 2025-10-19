#!/usr/bin/env python3
"""
BLAST searches and result parsing using BioPython.

This script demonstrates:
- Running BLAST searches via NCBI (qblast)
- Parsing BLAST XML output
- Filtering and analyzing results
- Working with alignments and HSPs
"""

from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO


def run_blast_online(sequence, program="blastn", database="nt", expect=0.001):
    """
    Run BLAST search via NCBI's qblast.

    Parameters:
    - sequence: Sequence string or Seq object
    - program: blastn, blastp, blastx, tblastn, tblastx
    - database: nt (nucleotide), nr (protein), refseq_rna, etc.
    - expect: E-value threshold
    """

    print(f"Running {program} search against {database} database...")
    print(f"E-value threshold: {expect}")
    print("-" * 60)

    # Run BLAST
    result_handle = NCBIWWW.qblast(
        program=program,
        database=database,
        sequence=sequence,
        expect=expect,
        hitlist_size=50,  # Number of sequences to show alignments for
    )

    # Save results
    output_file = "blast_results.xml"
    with open(output_file, "w") as out:
        out.write(result_handle.read())

    result_handle.close()

    print(f"BLAST search complete. Results saved to {output_file}")
    print()

    return output_file


def parse_blast_results(xml_file, max_hits=10, evalue_threshold=0.001):
    """Parse BLAST XML results."""

    print(f"Parsing BLAST results from: {xml_file}")
    print(f"E-value threshold: {evalue_threshold}")
    print("=" * 60)

    with open(xml_file) as result_handle:
        blast_record = NCBIXML.read(result_handle)

    print(f"Query: {blast_record.query}")
    print(f"Query length: {blast_record.query_length} residues")
    print(f"Database: {blast_record.database}")
    print(f"Number of alignments: {len(blast_record.alignments)}")
    print()

    hit_count = 0

    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect <= evalue_threshold:
                hit_count += 1

                if hit_count <= max_hits:
                    print(f"Hit {hit_count}:")
                    print(f"  Sequence: {alignment.title}")
                    print(f"  Length: {alignment.length}")
                    print(f"  E-value: {hsp.expect:.2e}")
                    print(f"  Score: {hsp.score}")
                    print(f"  Identities: {hsp.identities}/{hsp.align_length} ({hsp.identities / hsp.align_length * 100:.1f}%)")
                    print(f"  Positives: {hsp.positives}/{hsp.align_length} ({hsp.positives / hsp.align_length * 100:.1f}%)")
                    print(f"  Gaps: {hsp.gaps}/{hsp.align_length}")
                    print(f"  Query range: {hsp.query_start} - {hsp.query_end}")
                    print(f"  Subject range: {hsp.sbjct_start} - {hsp.sbjct_end}")
                    print()

                    # Show alignment (first 100 characters)
                    print("  Alignment preview:")
                    print(f"  Query:  {hsp.query[:100]}")
                    print(f"  Match:  {hsp.match[:100]}")
                    print(f"  Sbjct:  {hsp.sbjct[:100]}")
                    print()

    print(f"Total significant hits (E-value <= {evalue_threshold}): {hit_count}")
    print()

    return blast_record


def parse_multiple_queries(xml_file):
    """Parse BLAST results with multiple queries."""

    print(f"Parsing multiple queries from: {xml_file}")
    print("=" * 60)

    with open(xml_file) as result_handle:
        blast_records = NCBIXML.parse(result_handle)

        for i, blast_record in enumerate(blast_records, 1):
            print(f"\nQuery {i}: {blast_record.query}")
            print(f"  Number of hits: {len(blast_record.alignments)}")

            if blast_record.alignments:
                best_hit = blast_record.alignments[0]
                best_hsp = best_hit.hsps[0]
                print(f"  Best hit: {best_hit.title[:80]}...")
                print(f"  Best E-value: {best_hsp.expect:.2e}")


def filter_blast_results(blast_record, min_identity=0.7, min_coverage=0.5):
    """Filter BLAST results by identity and coverage."""

    print(f"Filtering results:")
    print(f"  Minimum identity: {min_identity * 100}%")
    print(f"  Minimum coverage: {min_coverage * 100}%")
    print("-" * 60)

    filtered_hits = []

    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            identity_fraction = hsp.identities / hsp.align_length
            coverage = hsp.align_length / blast_record.query_length

            if identity_fraction >= min_identity and coverage >= min_coverage:
                filtered_hits.append(
                    {
                        "title": alignment.title,
                        "length": alignment.length,
                        "evalue": hsp.expect,
                        "identity": identity_fraction,
                        "coverage": coverage,
                        "alignment": alignment,
                        "hsp": hsp,
                    }
                )

    print(f"Found {len(filtered_hits)} hits matching criteria")
    print()

    # Sort by E-value
    filtered_hits.sort(key=lambda x: x["evalue"])

    # Display top hits
    for i, hit in enumerate(filtered_hits[:5], 1):
        print(f"{i}. {hit['title'][:80]}")
        print(f"   Identity: {hit['identity']*100:.1f}%, Coverage: {hit['coverage']*100:.1f}%, E-value: {hit['evalue']:.2e}")
        print()

    return filtered_hits


def extract_hit_sequences(blast_record, output_file="blast_hits.fasta"):
    """Extract aligned sequences from BLAST results."""

    print(f"Extracting hit sequences to {output_file}...")

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    records = []

    for i, alignment in enumerate(blast_record.alignments[:10]):  # Top 10 hits
        hsp = alignment.hsps[0]  # Best HSP for this alignment

        # Extract accession from title
        accession = alignment.title.split()[0]

        # Create SeqRecord from aligned subject sequence
        record = SeqRecord(
            Seq(hsp.sbjct.replace("-", "")),  # Remove gaps
            id=accession,
            description=f"E-value: {hsp.expect:.2e}, Identity: {hsp.identities}/{hsp.align_length}",
        )

        records.append(record)

    # Write to FASTA
    SeqIO.write(records, output_file, "fasta")

    print(f"Extracted {len(records)} sequences")
    print()


def analyze_blast_statistics(blast_record):
    """Compute statistics from BLAST results."""

    print("BLAST Result Statistics:")
    print("-" * 60)

    if not blast_record.alignments:
        print("No hits found")
        return

    evalues = []
    identities = []
    scores = []

    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            evalues.append(hsp.expect)
            identities.append(hsp.identities / hsp.align_length)
            scores.append(hsp.score)

    import statistics

    print(f"Total HSPs: {len(evalues)}")
    print(f"\nE-values:")
    print(f"  Min: {min(evalues):.2e}")
    print(f"  Max: {max(evalues):.2e}")
    print(f"  Median: {statistics.median(evalues):.2e}")
    print(f"\nIdentity percentages:")
    print(f"  Min: {min(identities)*100:.1f}%")
    print(f"  Max: {max(identities)*100:.1f}%")
    print(f"  Mean: {statistics.mean(identities)*100:.1f}%")
    print(f"\nBit scores:")
    print(f"  Min: {min(scores):.1f}")
    print(f"  Max: {max(scores):.1f}")
    print(f"  Mean: {statistics.mean(scores):.1f}")
    print()


def example_workflow():
    """Demonstrate BLAST workflow."""

    print("=" * 60)
    print("BioPython BLAST Example Workflow")
    print("=" * 60)
    print()

    # Example sequence (human beta-globin)
    example_sequence = """
    ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAAGCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGC
    """.replace("\n", "").replace(" ", "")

    print("Example: Human beta-globin sequence")
    print(f"Length: {len(example_sequence)} bp")
    print()

    # Note: Uncomment to run actual BLAST search (takes time)
    # xml_file = run_blast_online(example_sequence, program="blastn", database="nt", expect=0.001)

    # For demonstration, use a pre-existing results file
    print("To run a real BLAST search, uncomment the run_blast_online() line")
    print("For now, demonstrating parsing with example results file")
    print()

    # If you have results, parse them:
    # blast_record = parse_blast_results("blast_results.xml", max_hits=5)
    # filtered = filter_blast_results(blast_record, min_identity=0.9)
    # analyze_blast_statistics(blast_record)
    # extract_hit_sequences(blast_record)


if __name__ == "__main__":
    example_workflow()

    print()
    print("Note: BLAST searches can take several minutes.")
    print("For production use, consider running local BLAST instead.")
