"""Literature retrieval from PubMed and import from other databases."""

import json
import logging
import time
from pathlib import Path
from typing import Any

import requests
from Bio import Entrez, Medline

from .config import Settings
from .models import Citation

logger = logging.getLogger(__name__)


class PubMedRetriever:
    """Retrieve literature from PubMed using NCBI E-utilities."""

    def __init__(self, settings: Settings) -> None:
        """
        Initialize PubMed retriever.

        Args:
            settings: Application settings with API credentials
        """
        self.settings = settings
        Entrez.email = settings.ncbi_email
        if settings.ncbi_api_key:
            Entrez.api_key = settings.ncbi_api_key
            self.requests_per_second = 10
        else:
            self.requests_per_second = 3
        self.delay = 1.0 / self.requests_per_second

    def search(self, query: str, max_results: int = 10000) -> list[str]:
        """
        Search PubMed and return list of PMIDs.

        Args:
            query: PubMed search query
            max_results: Maximum number of results to retrieve

        Returns:
            List of PMIDs
        """
        logger.info(f"Searching PubMed with query: {query[:100]}...")

        try:
            handle = Entrez.esearch(
                db="pubmed", term=query, retmax=max_results, usehistory="y"
            )
            search_results = Entrez.read(handle)
            handle.close()

            count = int(search_results["Count"])
            logger.info(f"Found {count} results")

            # Get PMIDs
            pmids = search_results["IdList"]

            # If more results than retmax, use history
            if count > len(pmids):
                web_env = search_results["WebEnv"]
                query_key = search_results["QueryKey"]
                pmids = self._fetch_pmids_from_history(web_env, query_key, count)

            return pmids

        except Exception as e:
            logger.error(f"PubMed search failed: {e}")
            raise

    def _fetch_pmids_from_history(
        self, web_env: str, query_key: str, count: int
    ) -> list[str]:
        """Fetch PMIDs using history server."""
        batch_size = 500
        pmids = []

        for start in range(0, count, batch_size):
            try:
                handle = Entrez.esearch(
                    db="pubmed",
                    WebEnv=web_env,
                    query_key=query_key,
                    retstart=start,
                    retmax=batch_size,
                )
                batch_results = Entrez.read(handle)
                handle.close()
                pmids.extend(batch_results["IdList"])
                time.sleep(self.delay)
            except Exception as e:
                logger.warning(f"Error fetching batch at {start}: {e}")

        return pmids

    def fetch_records(self, pmids: list[str]) -> list[Citation]:
        """
        Fetch full records for given PMIDs.

        Args:
            pmids: List of PubMed IDs

        Returns:
            List of Citation objects
        """
        logger.info(f"Fetching {len(pmids)} PubMed records...")

        citations = []
        batch_size = 200

        for i in range(0, len(pmids), batch_size):
            batch = pmids[i : i + batch_size]
            try:
                handle = Entrez.efetch(
                    db="pubmed", id=batch, rettype="medline", retmode="text"
                )
                records = Medline.parse(handle)

                for record in records:
                    citation = self._parse_medline_record(record)
                    citations.append(citation)

                handle.close()
                time.sleep(self.delay)

                logger.info(f"Fetched {len(citations)}/{len(pmids)} records")

            except Exception as e:
                logger.error(f"Error fetching batch {i}: {e}")
                continue

        return citations

    def _parse_medline_record(self, record: dict[str, Any]) -> Citation:
        """Parse Medline record to Citation object."""
        # Extract authors
        authors = record.get("AU", [])

        # Extract year
        year = None
        if "DP" in record:
            try:
                year = int(record["DP"].split()[0])
            except (ValueError, IndexError):
                pass

        # Extract DOI
        doi = ""
        if "AID" in record:
            for aid in record["AID"]:
                if "[doi]" in aid.lower():
                    doi = aid.replace("[doi]", "").strip()
                    break

        return Citation(
            id=record.get("PMID", ""),
            pmid=record.get("PMID", ""),
            title=record.get("TI", ""),
            authors=authors,
            year=year,
            journal=record.get("JT", ""),
            volume=record.get("VI", ""),
            issue=record.get("IP", ""),
            pages=record.get("PG", ""),
            abstract=record.get("AB", ""),
            doi=doi,
            keywords=record.get("OT", []),
            mesh_terms=record.get("MH", []),
            language=record.get("LA", [""])[0] if "LA" in record else "",
            source_database="PubMed",
        )

    def search_and_fetch(self, query: str, max_results: int = 10000) -> list[Citation]:
        """
        Search PubMed and fetch full records in one call.

        Args:
            query: PubMed search query
            max_results: Maximum results to retrieve

        Returns:
            List of Citations
        """
        pmids = self.search(query, max_results)
        if not pmids:
            logger.warning("No results found")
            return []

        return self.fetch_records(pmids)


def save_citations_jsonl(citations: list[Citation], output_path: Path) -> None:
    """
    Save citations to JSONL file.

    Args:
        citations: List of citations
        output_path: Output file path
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as f:
        for citation in citations:
            f.write(citation.model_dump_json() + "\n")

    logger.info(f"Saved {len(citations)} citations to {output_path}")


def load_citations_jsonl(input_path: Path) -> list[Citation]:
    """
    Load citations from JSONL file.

    Args:
        input_path: Input file path

    Returns:
        List of citations
    """
    citations = []
    with open(input_path, "r", encoding="utf-8") as f:
        for line in f:
            data = json.loads(line)
            citations.append(Citation(**data))

    logger.info(f"Loaded {len(citations)} citations from {input_path}")
    return citations
