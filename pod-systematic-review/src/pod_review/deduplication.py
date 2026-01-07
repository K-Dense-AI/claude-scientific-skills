"""De-duplication of citations using exact and fuzzy matching."""

import logging
from dataclasses import dataclass
from typing import Any

from thefuzz import fuzz

from .models import Citation

logger = logging.getLogger(__name__)


@dataclass
class DuplicateMatch:
    """Represents a duplicate match between two citations."""

    original_id: str
    duplicate_id: str
    match_type: str  # 'doi', 'pmid', 'fuzzy_title', 'fuzzy_author_title'
    similarity_score: float  # For fuzzy matches


class Deduplicator:
    """De-duplicate citations using multiple strategies."""

    def __init__(
        self, title_threshold: float = 0.90, author_threshold: float = 0.85
    ) -> None:
        """
        Initialize deduplicator.

        Args:
            title_threshold: Minimum similarity score for title matching (0-1)
            author_threshold: Minimum similarity score for author matching (0-1)
        """
        self.title_threshold = title_threshold
        self.author_threshold = author_threshold

    def deduplicate(self, citations: list[Citation]) -> tuple[list[Citation], list[DuplicateMatch]]:
        """
        Deduplicate citations.

        Args:
            citations: List of citations to deduplicate

        Returns:
            Tuple of (unique_citations, duplicate_matches)
        """
        logger.info(f"Deduplicating {len(citations)} citations...")

        # Track which citations are duplicates
        duplicate_ids: set[str] = set()
        matches: list[DuplicateMatch] = []

        # Build indices for efficient lookup
        doi_index: dict[str, Citation] = {}
        pmid_index: dict[str, Citation] = {}
        title_index: list[tuple[str, Citation]] = []

        for citation in citations:
            if citation.id in duplicate_ids:
                continue

            # Check DOI match
            if citation.doi:
                doi_key = citation.doi.lower().strip()
                if doi_key in doi_index:
                    matches.append(
                        DuplicateMatch(
                            original_id=doi_index[doi_key].id,
                            duplicate_id=citation.id,
                            match_type="doi",
                            similarity_score=1.0,
                        )
                    )
                    duplicate_ids.add(citation.id)
                    continue
                else:
                    doi_index[doi_key] = citation

            # Check PMID match
            if citation.pmid:
                pmid_key = citation.pmid.strip()
                if pmid_key in pmid_index:
                    matches.append(
                        DuplicateMatch(
                            original_id=pmid_index[pmid_key].id,
                            duplicate_id=citation.id,
                            match_type="pmid",
                            similarity_score=1.0,
                        )
                    )
                    duplicate_ids.add(citation.id)
                    continue
                else:
                    pmid_index[pmid_key] = citation

            # Fuzzy matching on title and authors
            is_duplicate = False
            for existing_title, existing_citation in title_index:
                if self._is_fuzzy_duplicate(citation, existing_citation):
                    similarity = self._calculate_similarity(citation, existing_citation)
                    matches.append(
                        DuplicateMatch(
                            original_id=existing_citation.id,
                            duplicate_id=citation.id,
                            match_type="fuzzy_title",
                            similarity_score=similarity,
                        )
                    )
                    duplicate_ids.add(citation.id)
                    is_duplicate = True
                    break

            if not is_duplicate:
                title_index.append((citation.title.lower(), citation))

        # Filter out duplicates
        unique_citations = [c for c in citations if c.id not in duplicate_ids]

        logger.info(
            f"Found {len(matches)} duplicates. {len(unique_citations)} unique citations remain."
        )

        return unique_citations, matches

    def _is_fuzzy_duplicate(self, citation1: Citation, citation2: Citation) -> bool:
        """
        Check if two citations are fuzzy duplicates.

        Args:
            citation1: First citation
            citation2: Second citation

        Returns:
            True if likely duplicates
        """
        # Title similarity
        title1 = citation1.title.lower().strip()
        title2 = citation2.title.lower().strip()

        if not title1 or not title2:
            return False

        title_similarity = fuzz.ratio(title1, title2) / 100.0

        if title_similarity < self.title_threshold:
            return False

        # If titles are very similar, check authors and year
        # Authors
        authors1 = [a.lower() for a in citation1.authors]
        authors2 = [a.lower() for a in citation2.authors]

        if authors1 and authors2:
            # Check if first authors match
            first_author_similarity = fuzz.ratio(authors1[0], authors2[0]) / 100.0
            if first_author_similarity < self.author_threshold:
                return False

        # Year check (must match if both available)
        if citation1.year and citation2.year:
            if citation1.year != citation2.year:
                return False

        return True

    def _calculate_similarity(self, citation1: Citation, citation2: Citation) -> float:
        """Calculate overall similarity score between two citations."""
        title_sim = fuzz.ratio(citation1.title.lower(), citation2.title.lower()) / 100.0

        if citation1.authors and citation2.authors:
            author_sim = (
                fuzz.ratio(citation1.authors[0].lower(), citation2.authors[0].lower()) / 100.0
            )
            return (title_sim + author_sim) / 2
        else:
            return title_sim


def generate_dedup_report(
    original_count: int, unique_count: int, matches: list[DuplicateMatch]
) -> dict[str, Any]:
    """
    Generate deduplication report.

    Args:
        original_count: Original number of citations
        unique_count: Number of unique citations after dedup
        matches: List of duplicate matches

    Returns:
        Report dictionary
    """
    match_types: dict[str, int] = {}
    for match in matches:
        match_types[match.match_type] = match_types.get(match.match_type, 0) + 1

    report = {
        "original_count": original_count,
        "unique_count": unique_count,
        "duplicates_removed": original_count - unique_count,
        "duplicate_rate": (original_count - unique_count) / original_count if original_count > 0 else 0,
        "matches_by_type": match_types,
        "total_matches": len(matches),
    }

    return report
