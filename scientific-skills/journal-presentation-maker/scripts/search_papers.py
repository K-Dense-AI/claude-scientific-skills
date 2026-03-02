#!/usr/bin/env python3
"""
Search and filter academic papers from multiple sources including preprints.
Supports keyword-based search with citation and date filtering.
"""

import sys
import json
from typing import List, Dict, Optional
from datetime import datetime

def search_papers(
    keywords: str,
    min_papers: int = 5,
    max_papers: int = 100,
    high_citation_threshold: int = 200,
    medium_citation_threshold: int = 50
) -> List[Dict]:
    """
    Search for papers using keywords and apply smart filtering.
    
    Args:
        keywords: Search keywords
        min_papers: Minimum number of papers to return (default: 5)
        max_papers: Maximum number of papers to return (default: 100)
        high_citation_threshold: Citation count for "classic" papers (default: 200)
        medium_citation_threshold: Citation count for medium-cited papers (default: 50)
    
    Returns:
        List of paper dictionaries with metadata
    """
    print(f"Searching for papers with keywords: {keywords}", file=sys.stderr)
    print(f"Paper range: {min_papers}-{max_papers}", file=sys.stderr)
    
    # This is a placeholder - actual implementation would use:
    # - PubMed API
    # - Google Scholar scraping
    # - Crossref API
    # - bioRxiv/medRxiv/arXiv APIs
    
    # Return format example
    example_results = []
    
    return example_results


def filter_papers(
    papers: List[Dict],
    high_citation_threshold: int = 200,
    medium_citation_threshold: int = 50
) -> List[Dict]:
    """
    Apply smart filtering to papers based on citation and date.
    
    Filtering strategy:
    - 2023-2025: Include regardless of citations
    - 2020-2022: Include if citations >= medium_threshold
    - Before 2020: Include only if citations >= high_citation_threshold
    - Preprints: Include regardless of date/citations
    
    Args:
        papers: List of paper dictionaries
        high_citation_threshold: Minimum citations for pre-2020 papers
        medium_citation_threshold: Minimum citations for 2020-2022 papers
    
    Returns:
        Filtered and sorted list of papers
    """
    current_year = datetime.now().year
    filtered = []
    
    for paper in papers:
        year = paper.get('year', 0)
        citations = paper.get('citations', 0)
        is_preprint = paper.get('is_preprint', False)
        
        # Always include preprints
        if is_preprint:
            paper['category'] = 'preprint'
            filtered.append(paper)
            continue
        
        # Recent papers (2023-2025)
        if year >= 2023:
            paper['category'] = 'recent'
            filtered.append(paper)
        
        # Medium recent (2020-2022)
        elif 2020 <= year < 2023:
            if citations >= medium_citation_threshold:
                paper['category'] = 'medium-recent'
                filtered.append(paper)
        
        # Classic papers (pre-2020)
        elif year < 2020:
            if citations >= high_citation_threshold:
                paper['category'] = 'classic'
                filtered.append(paper)
    
    # Sort: recent first, then by citations
    filtered.sort(key=lambda x: (
        -x.get('year', 0),
        -x.get('citations', 0)
    ))
    
    return filtered


def format_paper_list(papers: List[Dict]) -> str:
    """
    Format papers for user selection.
    
    Args:
        papers: List of paper dictionaries
    
    Returns:
        Formatted string for display
    """
    output = []
    
    for i, paper in enumerate(papers, 1):
        title = paper.get('title', 'Unknown')
        year = paper.get('year', 'N/A')
        authors = paper.get('authors', 'Unknown')
        journal = paper.get('journal', 'Unknown')
        citations = paper.get('citations', 0)
        doi = paper.get('doi', 'N/A')
        is_preprint = paper.get('is_preprint', False)
        category = paper.get('category', '')
        
        # Format entry
        if is_preprint:
            entry = f"{i}. [{title}] ({year}, {journal} - PREPRINT)\n"
        elif category == 'classic':
            entry = f"{i}. [{title}] ({year}, {journal}, cited: {citations}) ★ High-impact classic\n"
        else:
            entry = f"{i}. [{title}] ({year}, {journal}, cited: {citations})\n"
        
        entry += f"   {authors}\n"
        entry += f"   DOI: {doi}\n"
        
        output.append(entry)
    
    return "\n".join(output)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: search_papers.py <keywords> [min_papers] [max_papers]")
        sys.exit(1)
    
    keywords = sys.argv[1]
    min_papers = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    max_papers = int(sys.argv[3]) if len(sys.argv) > 3 else 100
    
    papers = search_papers(keywords, min_papers, max_papers)
    filtered = filter_papers(papers)
    
    # Limit to requested range
    result = filtered[0:max_papers]
    
    if len(result) < min_papers:
        print(f"Warning: Only found {len(result)} papers (minimum requested: {min_papers})", 
              file=sys.stderr)
    
    print(format_paper_list(result))
