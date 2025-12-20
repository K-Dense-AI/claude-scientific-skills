#!/usr/bin/env python3
"""
Search PubMed and OpenAlex for recent meta-analyses and systematic reviews
on smartphone addiction and sleep disorders (2020-2025).
"""

import requests
import json
import time
import sys
import os
from typing import List, Dict, Any

# Add OpenAlex scripts to path
sys.path.insert(0, 'scientific-skills/openalex-database/scripts')
from openalex_client import OpenAlexClient


def search_pubmed(query: str, max_results: int = 100) -> List[Dict[str, Any]]:
    """
    Search PubMed using E-utilities API.

    Args:
        query: PubMed search query
        max_results: Maximum number of results to retrieve

    Returns:
        List of article information dictionaries
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

    # Step 1: Search for PMIDs
    print(f"Searching PubMed with query: {query}")
    search_url = f"{base_url}esearch.fcgi"
    search_params = {
        "db": "pubmed",
        "term": query,
        "retmax": max_results,
        "retmode": "json",
        "sort": "relevance"
    }

    response = requests.get(search_url, params=search_params)
    response.raise_for_status()
    search_results = response.json()

    pmids = search_results.get("esearchresult", {}).get("idlist", [])
    total_count = search_results.get("esearchresult", {}).get("count", 0)

    print(f"Found {total_count} total results, retrieving {len(pmids)} articles")

    if not pmids:
        return []

    # Step 2: Fetch article details in batches
    articles = []
    batch_size = 20

    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i+batch_size]
        time.sleep(0.4)  # Rate limiting: ~3 requests/second

        fetch_url = f"{base_url}esummary.fcgi"
        fetch_params = {
            "db": "pubmed",
            "id": ",".join(batch),
            "retmode": "json"
        }

        response = requests.get(fetch_url, params=fetch_params)
        response.raise_for_status()
        fetch_results = response.json()

        results = fetch_results.get("result", {})
        for pmid in batch:
            if pmid in results:
                article = results[pmid]
                articles.append({
                    "pmid": pmid,
                    "title": article.get("title", ""),
                    "authors": [author.get("name", "") for author in article.get("authors", [])],
                    "journal": article.get("fulljournalname", ""),
                    "pub_date": article.get("pubdate", ""),
                    "source": "PubMed",
                    "doi": article.get("elocationid", "").replace("doi: ", ""),
                    "pub_type": article.get("pubtype", [])
                })

    return articles


def search_openalex(search_terms: str, years: str = "2020-2025", max_results: int = 100) -> List[Dict[str, Any]]:
    """
    Search OpenAlex for reviews on a topic.

    Args:
        search_terms: Search query
        years: Year range filter
        max_results: Maximum number of results to retrieve

    Returns:
        List of work information dictionaries
    """
    print(f"Searching OpenAlex for: {search_terms}")

    client = OpenAlexClient(email="research@example.com")

    # Search for works with type filter for reviews
    results = client.search_works(
        search=search_terms,
        filter_params={
            "publication_year": years,
            "type": "review|meta-analysis"
        },
        per_page=min(max_results, 200),
        sort="cited_by_count:desc"
    )

    works = results.get("results", [])
    total_count = results.get("meta", {}).get("count", 0)

    print(f"Found {total_count} total results, retrieving {len(works)} works")

    # Extract relevant information
    articles = []
    for work in works:
        # Get author names
        authorships = work.get("authorships", [])
        authors = [a.get("author", {}).get("display_name", "") for a in authorships[:5]]

        # Get publication info
        primary_location = work.get("primary_location", {})
        source = primary_location.get("source", {})

        articles.append({
            "openalex_id": work.get("id", ""),
            "doi": work.get("doi", "").replace("https://doi.org/", ""),
            "title": work.get("title", ""),
            "authors": authors,
            "journal": source.get("display_name", ""),
            "pub_date": str(work.get("publication_year", "")),
            "cited_by_count": work.get("cited_by_count", 0),
            "source": "OpenAlex",
            "type": work.get("type", ""),
            "is_open_access": work.get("open_access", {}).get("is_oa", False),
            "oa_url": work.get("open_access", {}).get("oa_url", "")
        })

    return articles


def deduplicate_by_doi(articles: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Remove duplicate articles based on DOI."""
    seen_dois = set()
    unique_articles = []

    for article in articles:
        doi = article.get("doi", "")
        if doi and doi in seen_dois:
            continue
        if doi:
            seen_dois.add(doi)
        unique_articles.append(article)

    return unique_articles


def main():
    # PubMed search query for meta-analyses and systematic reviews on smartphone addiction and sleep
    pubmed_query = (
        "(smartphone[tiab] OR phone[tiab] OR mobile[tiab] OR screen time[tiab] OR "
        "phone addiction[tiab] OR smartphone addiction[tiab]) AND "
        "(sleep[tiab] OR insomnia[tiab] OR sleep disorder[tiab] OR sleep quality[tiab] OR "
        "sleep disturbance[tiab]) AND "
        "(Meta-Analysis[pt] OR Systematic Review[pt]) AND "
        "2020:2025[dp]"
    )

    # OpenAlex search terms
    openalex_query = "smartphone addiction sleep insomnia"

    print("=" * 80)
    print("SEARCHING FOR SMARTPHONE ADDICTION AND SLEEP DISORDER LITERATURE")
    print("Focus: Meta-analyses and Systematic Reviews (2020-2025)")
    print("=" * 80)
    print()

    # Search PubMed
    print("\n--- PUBMED SEARCH ---")
    pubmed_articles = search_pubmed(pubmed_query, max_results=50)
    print(f"Retrieved {len(pubmed_articles)} articles from PubMed\n")

    # Search OpenAlex
    print("\n--- OPENALEX SEARCH ---")
    openalex_articles = search_openalex(openalex_query, years="2020-2025", max_results=50)
    print(f"Retrieved {len(openalex_articles)} articles from OpenAlex\n")

    # Combine and deduplicate
    all_articles = pubmed_articles + openalex_articles
    unique_articles = deduplicate_by_doi(all_articles)

    # Sort by year (most recent first) and citation count
    unique_articles.sort(
        key=lambda x: (
            -int(x.get("pub_date", "0")[:4]) if x.get("pub_date") else 0,
            -x.get("cited_by_count", 0)
        )
    )

    print("\n" + "=" * 80)
    print(f"TOTAL UNIQUE ARTICLES FOUND: {len(unique_articles)}")
    print("=" * 80)

    # Save results
    output_file = "smartphone_sleep_reviews.json"
    with open(output_file, "w") as f:
        json.dump(unique_articles, f, indent=2)

    print(f"\nResults saved to: {output_file}")

    # Display summary
    print("\n" + "=" * 80)
    print("TOP 20 RESULTS (SORTED BY YEAR AND CITATIONS)")
    print("=" * 80)

    for i, article in enumerate(unique_articles[:20], 1):
        print(f"\n{i}. {article.get('title', 'No title')}")

        authors = article.get('authors', [])
        if authors:
            author_str = ", ".join(authors[:3])
            if len(authors) > 3:
                author_str += f" et al. ({len(authors)} authors)"
            print(f"   Authors: {author_str}")

        print(f"   Journal: {article.get('journal', 'Unknown')}")
        print(f"   Year: {article.get('pub_date', 'Unknown')}")

        if article.get('cited_by_count'):
            print(f"   Citations: {article.get('cited_by_count')}")

        if article.get('doi'):
            print(f"   DOI: {article.get('doi')}")

        if article.get('pmid'):
            print(f"   PMID: {article.get('pmid')}")

        if article.get('is_open_access'):
            print(f"   Open Access: Yes")
            if article.get('oa_url'):
                print(f"   OA URL: {article.get('oa_url')}")

        print(f"   Source: {article.get('source')}")

    return unique_articles


if __name__ == "__main__":
    articles = main()
