#!/usr/bin/env python3
"""
NCBI Entrez database access using BioPython.

This script demonstrates:
- Searching NCBI databases
- Downloading sequences by accession
- Retrieving PubMed articles
- Batch downloading with WebEnv
- Proper error handling and rate limiting
"""

import time
from Bio import Entrez, SeqIO

# IMPORTANT: Always set your email
Entrez.email = "your.email@example.com"  # Change this!


def search_nucleotide(query, max_results=10):
    """Search NCBI nucleotide database."""

    print(f"Searching nucleotide database for: {query}")
    print("-" * 60)

    handle = Entrez.esearch(db="nucleotide", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()

    print(f"Found {record['Count']} total matches")
    print(f"Returning top {len(record['IdList'])} IDs:")
    print(record["IdList"])
    print()

    return record["IdList"]


def fetch_sequence_by_accession(accession):
    """Download a sequence by accession number."""

    print(f"Fetching sequence: {accession}")

    try:
        handle = Entrez.efetch(
            db="nucleotide", id=accession, rettype="gb", retmode="text"
        )
        record = SeqIO.read(handle, "genbank")
        handle.close()

        print(f"Successfully retrieved: {record.id}")
        print(f"Description: {record.description}")
        print(f"Length: {len(record.seq)} bp")
        print(f"Organism: {record.annotations.get('organism', 'Unknown')}")
        print()

        return record

    except Exception as e:
        print(f"Error fetching {accession}: {e}")
        return None


def fetch_multiple_sequences(id_list, output_file="downloaded_sequences.fasta"):
    """Download multiple sequences and save to file."""

    print(f"Fetching {len(id_list)} sequences...")

    try:
        # For >200 IDs, efetch automatically uses POST
        handle = Entrez.efetch(
            db="nucleotide", id=id_list, rettype="fasta", retmode="text"
        )

        # Parse and save
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()

        SeqIO.write(records, output_file, "fasta")

        print(f"Successfully downloaded {len(records)} sequences to {output_file}")
        print()

        return records

    except Exception as e:
        print(f"Error fetching sequences: {e}")
        return []


def search_and_download(query, output_file, max_results=100):
    """Complete workflow: search and download sequences."""

    print(f"Searching and downloading: {query}")
    print("=" * 60)

    # Search
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()

    id_list = record["IdList"]
    print(f"Found {len(id_list)} sequences")

    if not id_list:
        print("No results found")
        return

    # Download in batches to be polite
    batch_size = 100
    all_records = []

    for start in range(0, len(id_list), batch_size):
        end = min(start + batch_size, len(id_list))
        batch_ids = id_list[start:end]

        print(f"Downloading batch {start // batch_size + 1} ({len(batch_ids)} sequences)...")

        handle = Entrez.efetch(
            db="nucleotide", id=batch_ids, rettype="fasta", retmode="text"
        )
        batch_records = list(SeqIO.parse(handle, "fasta"))
        handle.close()

        all_records.extend(batch_records)

        # Be polite - wait between requests
        time.sleep(0.5)

    # Save all records
    SeqIO.write(all_records, output_file, "fasta")
    print(f"Downloaded {len(all_records)} sequences to {output_file}")
    print()


def use_history_for_large_queries(query, max_results=1000):
    """Use NCBI History server for large queries."""

    print("Using NCBI History server for large query")
    print("-" * 60)

    # Search with history
    search_handle = Entrez.esearch(
        db="nucleotide", term=query, retmax=max_results, usehistory="y"
    )
    search_results = Entrez.read(search_handle)
    search_handle.close()

    count = int(search_results["Count"])
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]

    print(f"Found {count} total sequences")
    print(f"WebEnv: {webenv[:20]}...")
    print(f"QueryKey: {query_key}")
    print()

    # Fetch in batches using history
    batch_size = 500
    all_records = []

    for start in range(0, min(count, max_results), batch_size):
        end = min(start + batch_size, max_results)

        print(f"Downloading records {start + 1} to {end}...")

        fetch_handle = Entrez.efetch(
            db="nucleotide",
            rettype="fasta",
            retmode="text",
            retstart=start,
            retmax=batch_size,
            webenv=webenv,
            query_key=query_key,
        )

        batch_records = list(SeqIO.parse(fetch_handle, "fasta"))
        fetch_handle.close()

        all_records.extend(batch_records)

        # Be polite
        time.sleep(0.5)

    print(f"Downloaded {len(all_records)} sequences total")
    return all_records


def search_pubmed(query, max_results=10):
    """Search PubMed for articles."""

    print(f"Searching PubMed for: {query}")
    print("-" * 60)

    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()

    id_list = record["IdList"]
    print(f"Found {record['Count']} total articles")
    print(f"Returning {len(id_list)} PMIDs:")
    print(id_list)
    print()

    return id_list


def fetch_pubmed_abstracts(pmid_list):
    """Fetch PubMed article summaries."""

    print(f"Fetching summaries for {len(pmid_list)} articles...")

    handle = Entrez.efetch(db="pubmed", id=pmid_list, rettype="abstract", retmode="text")
    abstracts = handle.read()
    handle.close()

    print(abstracts[:500])  # Show first 500 characters
    print("...")
    print()


def get_database_info(database="nucleotide"):
    """Get information about an NCBI database."""

    print(f"Getting info for database: {database}")
    print("-" * 60)

    handle = Entrez.einfo(db=database)
    record = Entrez.read(handle)
    handle.close()

    db_info = record["DbInfo"]
    print(f"Name: {db_info['DbName']}")
    print(f"Description: {db_info['Description']}")
    print(f"Record count: {db_info['Count']}")
    print(f"Last update: {db_info['LastUpdate']}")
    print()


def link_databases(db_from, db_to, id_):
    """Find related records in other databases."""

    print(f"Finding links from {db_from} ID {id_} to {db_to}")
    print("-" * 60)

    handle = Entrez.elink(dbfrom=db_from, db=db_to, id=id_)
    record = Entrez.read(handle)
    handle.close()

    if record[0]["LinkSetDb"]:
        linked_ids = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]]
        print(f"Found {len(linked_ids)} linked records")
        print(f"IDs: {linked_ids[:10]}")
    else:
        print("No linked records found")

    print()


def example_workflow():
    """Demonstrate complete Entrez workflow."""

    print("=" * 60)
    print("BioPython Entrez Example Workflow")
    print("=" * 60)
    print()

    # Note: These are examples - uncomment to run with your email set

    # # Example 1: Search and get IDs
    # ids = search_nucleotide("Homo sapiens[Organism] AND COX1[Gene]", max_results=5)
    #
    # # Example 2: Fetch a specific sequence
    # fetch_sequence_by_accession("NM_001301717")
    #
    # # Example 3: Complete search and download
    # search_and_download("Escherichia coli[Organism] AND 16S", "ecoli_16s.fasta", max_results=50)
    #
    # # Example 4: PubMed search
    # pmids = search_pubmed("CRISPR[Title] AND 2023[PDAT]", max_results=5)
    # fetch_pubmed_abstracts(pmids[:2])
    #
    # # Example 5: Get database info
    # get_database_info("nucleotide")

    print("Examples are commented out. Uncomment and set your email to run.")


if __name__ == "__main__":
    example_workflow()

    print()
    print("IMPORTANT: Always set Entrez.email before using these functions!")
    print("NCBI requires an email address for their E-utilities.")
