#!/usr/bin/env python3
"""
Search and filter academic papers - Class-based refactored version.

This module provides the PaperSearcher class for searching academic databases
with smart filtering based on citation counts and publication dates.
"""

import sys
import json
from typing import List, Dict, Optional
from datetime import datetime
import logging

logger = logging.getLogger(__name__)


class PaperSearcher:
    """Search and filter academic papers from multiple sources.

    Example:
        >>> searcher = PaperSearcher("enzyme cascade", min_papers=5, max_papers=50)
        >>> searcher.search()
        >>> filtered = searcher.filter_papers()
        >>> print(searcher.format_paper_list(filtered))
    """

    def __init__(
        self,
        keywords: str,
        min_papers: int = 5,
        max_papers: int = 100,
        high_citation_threshold: int = 200,
        medium_citation_threshold: int = 50
    ):
        """Initialize paper searcher.

        Args:
            keywords: Search keywords
            min_papers: Minimum number of papers (default: 5)
            max_papers: Maximum number of papers (default: 100)
            high_citation_threshold: Citations for classics (default: 200)
            medium_citation_threshold: Citations for medium (default: 50)

        Raises:
            ValueError: If parameters are invalid
        """
        if not keywords or not keywords.strip():
            raise ValueError("Keywords cannot be empty")

        if min_papers < 1:
            raise ValueError(f"min_papers must be >= 1, got {min_papers}")

        if max_papers < min_papers:
            raise ValueError(f"max_papers must be >= min_papers")

        self.keywords = keywords.strip()
        self.min_papers = min_papers
        self.max_papers = max_papers
        self.high_citation_threshold = high_citation_threshold
        self.medium_citation_threshold = medium_citation_threshold
        self.results = []

        logger.info(
            f"Initialized: keywords='{self.keywords}', "
            f"range={self.min_papers}-{self.max_papers}"
        )

    def search(self) -> List[Dict]:
        """Execute search across databases.

        Returns:
            List of paper dictionaries
        """
        logger.info(f"Searching: {self.keywords}")

        # Placeholder - actual uses web_search/web_fetch
        self.results = []

        logger.info(f"Found {len(self.results)} papers")
        return self.results

    def filter_papers(self, papers: Optional[List[Dict]] = None) -> List[Dict]:
        """Apply smart filtering.

        Strategy:
        - 2023+: Include all
        - 2020-2022: Include if citations >= medium threshold
        - Pre-2020: Include if citations >= high threshold
        - Preprints: Always include

        Args:
            papers: Papers to filter (uses self.results if None)

        Returns:
            Filtered and sorted papers
        """
        if papers is None:
            if not self.results:
                raise ValueError("No papers. Run search() first.")
            papers = self.results

        logger.info(f"Filtering {len(papers)} papers")

        filtered = []

        for paper in papers:
            year = paper.get('year', 0)
            citations = paper.get('citations', 0)
            is_preprint = paper.get('is_preprint', False)

            if is_preprint:
                paper['category'] = 'preprint'
                filtered.append(paper)
            elif year >= 2023:
                paper['category'] = 'recent'
                filtered.append(paper)
            elif 2020 <= year < 2023 and citations >= self.medium_citation_threshold:
                paper['category'] = 'medium-recent'
                filtered.append(paper)
            elif year < 2020 and citations >= self.high_citation_threshold:
                paper['category'] = 'classic'
                filtered.append(paper)

        filtered.sort(key=lambda x: (-x.get('year', 0), -x.get('citations', 0)))

        logger.info(f"Filtered to {len(filtered)} papers")
        return filtered

    def format_paper_list(self, papers: List[Dict]) -> str:
        """Format papers for display.

        Args:
            papers: Papers to format

        Returns:
            Formatted string
        """
        if not papers:
            return "No papers to display."

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

            if is_preprint:
                entry = f"{i}. [{title}] ({year}, {journal} - PREPRINT)\n"
            elif category == 'classic':
                entry = f"{i}. [{title}] ({year}, {journal}, cited: {citations}) ★ High-impact\n"
            else:
                entry = f"{i}. [{title}] ({year}, {journal}, cited: {citations})\n"

            entry += f"   {authors}\n"
            entry += f"   DOI: {doi}\n"

            output.append(entry)

        return "\n".join(output)


def main():
    """CLI entry point."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    if len(sys.argv) < 2:
        print("Usage: search_papers_v2.py <keywords> [min] [max]")
        sys.exit(1)

    keywords = sys.argv[1]
    min_papers = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    max_papers = int(sys.argv[3]) if len(sys.argv) > 3 else 100

    try:
        searcher = PaperSearcher(keywords, min_papers, max_papers)
        searcher.search()
        filtered = searcher.filter_papers()
        result = filtered[:max_papers]

        if len(result) < min_papers:
            logger.warning(f"Only {len(result)} papers (min: {min_papers})")

        print(searcher.format_paper_list(result))

    except Exception as e:
        logger.error(f"Failed: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
