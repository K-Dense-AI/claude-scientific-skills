#!/usr/bin/env python3
"""
Module to convert paper queries (titles or title+authors) into BibTeX format
using OpenAlex for metadata retrieval.
"""

import sys
import os
from typing import List, Dict, Optional, Tuple

# Add the parent directory to path so we can import from scripts
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.openalex_client import OpenAlexClient


import unicodedata

def normalize_text(text: str) -> str:
    """Normalize text to removing accents and lowercasing."""
    if not text:
        return ""
    return unicodedata.normalize('NFKD', text).encode('ASCII', 'ignore').decode('utf-8').lower()


def _format_authors(authorships: List[Dict]) -> str:
    """Format authors list for BibTeX."""
    names = []
    for authorship in authorships:
        author = authorship.get('author', {})
        name = author.get('display_name', '')
        if name:
            names.append(name)
    return ' and '.join(names)


def _generate_citation_key(work: Dict) -> str:
    """Generate a citation key (FirstAuthorYear)."""
    try:
        if not work.get('authorships'):
            return f"unknown{work.get('publication_year', '0000')}"
        
        first_author = work['authorships'][0]['author']['display_name'].split()[-1].lower()
        year = work.get('publication_year', '0000')
        # Simple sanitization
        first_author = "".join([c for c in first_author if c.isalnum()])
        return f"{first_author}{year}"
    except Exception:
        return f"openalex-{work.get('id', '').split('/')[-1]}"


def work_to_bibtex(work: Dict) -> str:
    """Convert an OpenAlex work object to a BibTeX entry."""
    
    work_type = work.get('type', 'article')
    # Map OpenAlex types to BibTeX types
    bib_type = 'misc'
    if work_type == 'article':
        bib_type = 'article'
    elif work_type == 'book':
        bib_type = 'book'
    elif work_type in ['book-chapter', 'proceedings-article']:
        bib_type = 'inproceedings'
    elif work_type == 'dissertation':
        bib_type = 'phdthesis'
    
    citation_key = _generate_citation_key(work)
    
    fields = [
        f"title = {{{work.get('title', '')}}}",
        f"author = {{{_format_authors(work.get('authorships', []))}}}",
        f"year = {{{work.get('publication_year', '')}}}",
    ]
    
    # Optional fields
    if work.get('ids', {}).get('doi'):
        fields.append(f"doi = {{{work['ids']['doi'].replace('https://doi.org/', '')}}}")
    
    if work.get('primary_location', {}).get('source', {}):
        source = work['primary_location']['source']
        if source:
             fields.append(f"journal = {{{source.get('display_name', '')}}}")
             
    # Volume/Issue
    loc = work.get('biblio', {})
    if loc.get('volume'):
        fields.append(f"volume = {{{loc['volume']}}}")
    if loc.get('issue'):
        fields.append(f"number = {{{loc['issue']}}}")
    if loc.get('first_page') and loc.get('last_page'):
        fields.append(f"pages = {{{loc['first_page']}--{loc['last_page']}}}")
    elif loc.get('first_page'):
         fields.append(f"pages = {{{loc['first_page']}}}")

    # URL
    if work.get('id'):
         fields.append(f"url = {{{work['id']}}}")

    inner = ",\n  ".join(fields)
    return f"@{bib_type}{{{citation_key},\n  {inner}\n}}\n"


def titles_to_bib(titles: List[str], client: Optional[OpenAlexClient] = None) -> str:
    """
    Convert a list of paper titles to BibTeX entries.
    
    Args:
        titles: List of paper titles
        client: Optional OpenAlexClient instance
        
    Returns:
        String containing BibTeX entries
    """
    if not client:
        client = OpenAlexClient()
        
    bib_entries = []
    
    for title in titles:
        # Search for the work
        results = client.search_works(search=title, per_page=1)
        if results.get('results'):
            best_match = results['results'][0]
            # Verify if it's a good match? For now, assume top result is correct
            bib_entries.append(work_to_bibtex(best_match))
        else:
            bib_entries.append(f"% No match found for title: {title}\n")
            
    return "\n".join(bib_entries)


def titles_authors_to_bib(papers: List[Tuple[str, str]], client: Optional[OpenAlexClient] = None) -> str:
    """
    Convert a list of (title, author) tuples to BibTeX entries.
    
    Args:
        papers: List of (title, author) tuples/lists
        client: Optional OpenAlexClient instance
        
    Returns:
        String containing BibTeX entries
    """
    if not client:
        client = OpenAlexClient()
        
    bib_entries = []
    
    for title, target_author in papers:
        # Search primarily by title
        results = client.search_works(search=title, per_page=10)
        
        found = False
        work_found = None
        
        if results.get('results'):
            # Try to find the specific author in the results
            for work in results['results']:
                authors = [normalize_text(a.get('author', {}).get('display_name', '')) 
                          for a in work.get('authorships', [])]
                
                # Check if target author is in this work's authors
                # We do a loose check: is 'zhou' in 'jing zhou'?
                target_norm = normalize_text(target_author)
                if any(target_norm in a for a in authors):
                    work_found = work
                    found = True
                    break
            
            # If not found by author filter, maybe just take the first one 
            # if the title is an exact match? 
            # For now, let's rely on the filter, or fallback to first result 
            # if we are fairly confident (e.g. strict title match).
            # But safe bet is fallback to first if we can't match author, 
            # with a warning note?
            if not found:
                 # Fallback: Search with "Title Author" string
                 fallback_query = f"{title} {target_author}"
                 fallback_res = client.search_works(search=fallback_query, per_page=1)
                 if fallback_res.get('results'):
                     work_found = fallback_res['results'][0]
                     found = True
        
        if found and work_found:
            bib_entries.append(work_to_bibtex(work_found))
        else:
            bib_entries.append(f"% No match found for: {title} by {target_author}\n")
            
    return "\n".join(bib_entries)


if __name__ == "__main__":
    # Example usage
    client = OpenAlexClient()
    
    # Test case 1: Titles only
    print("Testing titles_to_bib...")
    titles = [
        "Attention Is All You Need",
        "Deep Residual Learning for Image Recognition"
    ]
    print(titles_to_bib(titles, client))
    
    # Test case 2: Title + Author
    print("\nTesting titles_authors_to_bib...")
    papers = [
        ("Recurrent Convolutional Neural Network Regression for Continuous Pain Intensity Estimation in Video", "Zhou"),
        ("Deep Pain: Exploiting Long Short-Term Memory Networks for Facial Expression Classification", "Rodriguez")
    ]
    print(titles_authors_to_bib(papers, client))
