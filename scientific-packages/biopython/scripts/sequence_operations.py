#!/usr/bin/env python3
"""
Common sequence operations using BioPython.

This script demonstrates basic sequence manipulation tasks like:
- Creating and manipulating Seq objects
- Transcription and translation
- Complement and reverse complement
- Calculating GC content and melting temperature
"""

from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, MeltingTemp as mt


def demonstrate_seq_operations():
    """Show common Seq object operations."""

    # Create DNA sequence
    dna_seq = Seq("ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTG")

    print("Original DNA sequence:")
    print(dna_seq)
    print()

    # Transcription (DNA -> RNA)
    rna_seq = dna_seq.transcribe()
    print("Transcribed to RNA:")
    print(rna_seq)
    print()

    # Translation (DNA -> Protein)
    protein_seq = dna_seq.translate()
    print("Translated to protein:")
    print(protein_seq)
    print()

    # Translation with stop codon handling
    protein_to_stop = dna_seq.translate(to_stop=True)
    print("Translated to first stop codon:")
    print(protein_to_stop)
    print()

    # Complement
    complement = dna_seq.complement()
    print("Complement:")
    print(complement)
    print()

    # Reverse complement
    reverse_complement = dna_seq.reverse_complement()
    print("Reverse complement:")
    print(reverse_complement)
    print()

    # GC content
    gc = gc_fraction(dna_seq) * 100
    print(f"GC content: {gc:.2f}%")
    print()

    # Melting temperature
    tm = mt.Tm_NN(dna_seq)
    print(f"Melting temperature (nearest-neighbor): {tm:.2f}°C")
    print()

    # Sequence searching
    codon_start = dna_seq.find("ATG")
    print(f"Start codon (ATG) position: {codon_start}")

    # Count occurrences
    g_count = dna_seq.count("G")
    print(f"Number of G nucleotides: {g_count}")
    print()


def translate_with_genetic_code():
    """Demonstrate translation with different genetic codes."""

    dna_seq = Seq("ATGGTGCATCTGACTCCTGAGGAGAAGTCT")

    # Standard genetic code (table 1)
    standard = dna_seq.translate(table=1)
    print("Standard genetic code translation:")
    print(standard)

    # Vertebrate mitochondrial code (table 2)
    mito = dna_seq.translate(table=2)
    print("Vertebrate mitochondrial code translation:")
    print(mito)
    print()


def working_with_codons():
    """Access genetic code tables."""
    from Bio.Data import CodonTable

    # Get standard genetic code
    standard_table = CodonTable.unambiguous_dna_by_id[1]

    print("Standard genetic code:")
    print(f"Start codons: {standard_table.start_codons}")
    print(f"Stop codons: {standard_table.stop_codons}")
    print()

    # Show some codon translations
    print("Example codons:")
    for codon in ["ATG", "TGG", "TAA", "TAG", "TGA"]:
        if codon in standard_table.stop_codons:
            print(f"{codon} -> STOP")
        else:
            aa = standard_table.forward_table.get(codon, "Unknown")
            print(f"{codon} -> {aa}")


if __name__ == "__main__":
    print("=" * 60)
    print("BioPython Sequence Operations Demo")
    print("=" * 60)
    print()

    demonstrate_seq_operations()
    print("-" * 60)
    translate_with_genetic_code()
    print("-" * 60)
    working_with_codons()
