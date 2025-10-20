#!/usr/bin/env python3
"""
UniProt REST API Client

A Python client for interacting with the UniProt REST API.
Provides helper functions for common operations including search,
retrieval, ID mapping, and streaming.

Usage examples:
    # Search for proteins
    results = search_proteins("insulin AND organism_name:human", format="json")

    # Get a single protein
    protein = get_protein("P12345", format="fasta")

    # Map IDs
    mapped = map_ids(["P12345", "P04637"], from_db="UniProtKB_AC-ID", to_db="PDB")

    # Stream large results
    for batch in stream_results("taxonomy_id:9606 AND reviewed:true", format="fasta"):
        process(batch)
"""

import requests
import time
import json
from typing import List, Dict, Optional, Generator
from urllib.parse import urlencode

BASE_URL = "https://rest.uniprot.org"
POLLING_INTERVAL = 3  # seconds


def search_proteins(query: str, format: str = "json",
                   fields: Optional[List[str]] = None,
                   size: int = 25) -> Dict:
    """
    Search UniProt database with a query.

    Args:
        query: Search query (e.g., "insulin AND organism_name:human")
        format: Response format (json, tsv, xlsx, xml, fasta, txt, rdf)
        fields: List of fields to return (e.g., ["accession", "gene_names", "organism_name"])
        size: Number of results per page (default 25, max 500)

    Returns:
        Response data in requested format
    """
    endpoint = f"{BASE_URL}/uniprotkb/search"

    params = {
        "query": query,
        "format": format,
        "size": size
    }

    if fields:
        params["fields"] = ",".join(fields)

    response = requests.get(endpoint, params=params)
    response.raise_for_status()

    if format == "json":
        return response.json()
    else:
        return response.text


def get_protein(accession: str, format: str = "json") -> str:
    """
    Retrieve a single protein entry by accession number.

    Args:
        accession: UniProt accession number (e.g., "P12345")
        format: Response format (json, txt, xml, fasta, gff, rdf)

    Returns:
        Protein data in requested format
    """
    endpoint = f"{BASE_URL}/uniprotkb/{accession}.{format}"

    response = requests.get(endpoint)
    response.raise_for_status()

    if format == "json":
        return response.json()
    else:
        return response.text


def batch_retrieve(accessions: List[str], format: str = "json",
                  fields: Optional[List[str]] = None) -> str:
    """
    Retrieve multiple protein entries efficiently.

    Args:
        accessions: List of UniProt accession numbers
        format: Response format
        fields: List of fields to return

    Returns:
        Combined results in requested format
    """
    query = " OR ".join([f"accession:{acc}" for acc in accessions])
    return search_proteins(query, format=format, fields=fields, size=len(accessions))


def stream_results(query: str, format: str = "fasta",
                  fields: Optional[List[str]] = None,
                  chunk_size: int = 8192) -> Generator[str, None, None]:
    """
    Stream large result sets without pagination.

    Args:
        query: Search query
        format: Response format
        fields: List of fields to return
        chunk_size: Size of chunks to yield

    Yields:
        Chunks of response data
    """
    endpoint = f"{BASE_URL}/uniprotkb/stream"

    params = {
        "query": query,
        "format": format
    }

    if fields:
        params["fields"] = ",".join(fields)

    response = requests.get(endpoint, params=params, stream=True)
    response.raise_for_status()

    for chunk in response.iter_content(chunk_size=chunk_size, decode_unicode=True):
        if chunk:
            yield chunk


def map_ids(ids: List[str], from_db: str, to_db: str,
           format: str = "json") -> Dict:
    """
    Map protein identifiers between different database systems.

    Args:
        ids: List of identifiers to map (max 100,000)
        from_db: Source database (e.g., "UniProtKB_AC-ID", "Gene_Name")
        to_db: Target database (e.g., "PDB", "Ensembl", "RefSeq_Protein")
        format: Response format

    Returns:
        Mapping results

    Note:
        - Maximum 100,000 IDs per job
        - Results stored for 7 days
        - See id_mapping_databases.md for all supported databases
    """
    if len(ids) > 100000:
        raise ValueError("Maximum 100,000 IDs allowed per mapping job")

    # Step 1: Submit job
    submit_endpoint = f"{BASE_URL}/idmapping/run"

    data = {
        "from": from_db,
        "to": to_db,
        "ids": ",".join(ids)
    }

    response = requests.post(submit_endpoint, data=data)
    response.raise_for_status()
    job_id = response.json()["jobId"]

    # Step 2: Poll for completion
    status_endpoint = f"{BASE_URL}/idmapping/status/{job_id}"

    while True:
        response = requests.get(status_endpoint)
        response.raise_for_status()
        status = response.json()

        if "results" in status or "failedIds" in status:
            break

        time.sleep(POLLING_INTERVAL)

    # Step 3: Retrieve results
    results_endpoint = f"{BASE_URL}/idmapping/results/{job_id}"

    params = {"format": format}
    response = requests.get(results_endpoint, params=params)
    response.raise_for_status()

    if format == "json":
        return response.json()
    else:
        return response.text


def get_available_fields() -> List[Dict]:
    """
    Get list of all available fields for queries.

    Returns:
        List of field definitions with names and descriptions
    """
    endpoint = f"{BASE_URL}/configure/uniprotkb/result-fields"

    response = requests.get(endpoint)
    response.raise_for_status()

    return response.json()


def get_id_mapping_databases() -> Dict:
    """
    Get list of all supported databases for ID mapping.

    Returns:
        Dictionary of database groups and their supported databases
    """
    endpoint = f"{BASE_URL}/configure/idmapping/fields"

    response = requests.get(endpoint)
    response.raise_for_status()

    return response.json()


# Example usage
if __name__ == "__main__":
    # Example 1: Search for human insulin proteins
    print("Searching for human insulin proteins...")
    results = search_proteins(
        "insulin AND organism_name:human AND reviewed:true",
        format="json",
        fields=["accession", "id", "gene_names", "protein_name"],
        size=5
    )
    print(json.dumps(results, indent=2))

    # Example 2: Get a specific protein in FASTA format
    print("\nRetrieving protein P01308 (human insulin)...")
    protein = get_protein("P01308", format="fasta")
    print(protein)

    # Example 3: Map UniProt IDs to PDB IDs
    print("\nMapping UniProt IDs to PDB...")
    mapping = map_ids(
        ["P01308", "P04637"],
        from_db="UniProtKB_AC-ID",
        to_db="PDB"
    )
    print(json.dumps(mapping, indent=2))
