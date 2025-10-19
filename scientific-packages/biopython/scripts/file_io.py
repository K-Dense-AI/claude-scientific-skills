#!/usr/bin/env python3
"""
File I/O operations using BioPython SeqIO.

This script demonstrates:
- Reading sequences from various formats
- Writing sequences to files
- Converting between formats
- Filtering and processing sequences
- Working with large files efficiently
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def read_sequences(filename, format_type):
    """Read and display sequences from a file."""

    print(f"Reading {format_type} file: {filename}")
    print("-" * 60)

    count = 0
    for record in SeqIO.parse(filename, format_type):
        count += 1
        print(f"ID: {record.id}")
        print(f"Name: {record.name}")
        print(f"Description: {record.description}")
        print(f"Sequence length: {len(record.seq)}")
        print(f"Sequence: {record.seq[:50]}...")
        print()

        # Only show first 3 sequences
        if count >= 3:
            break

    # Count total sequences
    total = len(list(SeqIO.parse(filename, format_type)))
    print(f"Total sequences in file: {total}")
    print()


def read_single_sequence(filename, format_type):
    """Read a single sequence from a file."""

    record = SeqIO.read(filename, format_type)

    print("Single sequence record:")
    print(f"ID: {record.id}")
    print(f"Sequence: {record.seq}")
    print()


def write_sequences(records, output_filename, format_type):
    """Write sequences to a file."""

    count = SeqIO.write(records, output_filename, format_type)
    print(f"Wrote {count} sequences to {output_filename} in {format_type} format")
    print()


def convert_format(input_file, input_format, output_file, output_format):
    """Convert sequences from one format to another."""

    count = SeqIO.convert(input_file, input_format, output_file, output_format)
    print(f"Converted {count} sequences from {input_format} to {output_format}")
    print()


def filter_sequences(input_file, format_type, min_length=100, max_length=1000):
    """Filter sequences by length."""

    filtered = []

    for record in SeqIO.parse(input_file, format_type):
        if min_length <= len(record.seq) <= max_length:
            filtered.append(record)

    print(f"Found {len(filtered)} sequences between {min_length} and {max_length} bp")
    return filtered


def extract_subsequence(input_file, format_type, seq_id, start, end):
    """Extract a subsequence from a specific record."""

    # Index for efficient access
    record_dict = SeqIO.index(input_file, format_type)

    if seq_id in record_dict:
        record = record_dict[seq_id]
        subseq = record.seq[start:end]
        print(f"Extracted subsequence from {seq_id} ({start}:{end}):")
        print(subseq)
        return subseq
    else:
        print(f"Sequence {seq_id} not found")
        return None


def create_sequence_records():
    """Create SeqRecord objects from scratch."""

    # Simple record
    simple_record = SeqRecord(
        Seq("ATGCATGCATGC"),
        id="seq001",
        name="MySequence",
        description="Example sequence"
    )

    # Record with annotations
    annotated_record = SeqRecord(
        Seq("ATGGTGCATCTGACTCCTGAGGAG"),
        id="seq002",
        name="GeneX",
        description="Important gene"
    )
    annotated_record.annotations["molecule_type"] = "DNA"
    annotated_record.annotations["organism"] = "Homo sapiens"

    return [simple_record, annotated_record]


def index_large_file(filename, format_type):
    """Index a large file for random access without loading into memory."""

    # Create index
    record_index = SeqIO.index(filename, format_type)

    print(f"Indexed {len(record_index)} sequences")
    print(f"Available IDs: {list(record_index.keys())[:10]}...")
    print()

    # Access specific record by ID
    if len(record_index) > 0:
        first_id = list(record_index.keys())[0]
        record = record_index[first_id]
        print(f"Accessed record: {record.id}")
        print()

    # Close index
    record_index.close()


def parse_with_quality_scores(fastq_file):
    """Parse FASTQ files with quality scores."""

    print("Parsing FASTQ with quality scores:")
    print("-" * 60)

    for record in SeqIO.parse(fastq_file, "fastq"):
        print(f"ID: {record.id}")
        print(f"Sequence: {record.seq[:50]}...")
        print(f"Quality scores (first 10): {record.letter_annotations['phred_quality'][:10]}")

        # Calculate average quality
        avg_quality = sum(record.letter_annotations["phred_quality"]) / len(record)
        print(f"Average quality: {avg_quality:.2f}")
        print()
        break  # Just show first record


def batch_process_large_file(input_file, format_type, batch_size=100):
    """Process large files in batches to manage memory."""

    batch = []
    count = 0

    for record in SeqIO.parse(input_file, format_type):
        batch.append(record)
        count += 1

        if len(batch) == batch_size:
            # Process batch
            print(f"Processing batch of {len(batch)} sequences...")
            # Do something with batch
            batch = []  # Clear for next batch

    # Process remaining records
    if batch:
        print(f"Processing final batch of {len(batch)} sequences...")

    print(f"Total sequences processed: {count}")


def example_workflow():
    """Demonstrate a complete workflow."""

    print("=" * 60)
    print("BioPython SeqIO Workflow Example")
    print("=" * 60)
    print()

    # Create example sequences
    records = create_sequence_records()

    # Write as FASTA
    write_sequences(records, "example_output.fasta", "fasta")

    # Write as GenBank
    write_sequences(records, "example_output.gb", "genbank")

    # Convert FASTA to GenBank (would work if file exists)
    # convert_format("input.fasta", "fasta", "output.gb", "genbank")

    print("Example workflow completed!")


if __name__ == "__main__":
    example_workflow()

    print()
    print("Note: This script demonstrates BioPython SeqIO operations.")
    print("Uncomment and adapt the functions for your specific files.")
