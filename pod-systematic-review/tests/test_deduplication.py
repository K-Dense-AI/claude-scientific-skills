"""Tests for deduplication module."""

import pytest
from pod_review.deduplication import Deduplicator
from pod_review.models import Citation


def test_exact_doi_match() -> None:
    """Test exact DOI matching."""
    dedup = Deduplicator()

    c1 = Citation(id="1", title="Test Study", doi="10.1234/test", authors=["Smith J"])
    c2 = Citation(id="2", title="Test Study", doi="10.1234/test", authors=["Smith J"])
    c3 = Citation(id="3", title="Different Study", doi="10.5678/other", authors=["Jones A"])

    unique, matches = dedup.deduplicate([c1, c2, c3])

    assert len(unique) == 2
    assert len(matches) == 1
    assert matches[0].match_type == "doi"


def test_exact_pmid_match() -> None:
    """Test exact PMID matching."""
    dedup = Deduplicator()

    c1 = Citation(id="1", title="Test Study", pmid="12345678", authors=["Smith J"])
    c2 = Citation(id="2", title="Test Study", pmid="12345678", authors=["Smith J"])

    unique, matches = dedup.deduplicate([c1, c2])

    assert len(unique) == 1
    assert len(matches) == 1
    assert matches[0].match_type == "pmid"


def test_fuzzy_title_match() -> None:
    """Test fuzzy title matching."""
    dedup = Deduplicator(title_threshold=0.90)

    c1 = Citation(
        id="1",
        title="Postoperative delirium in elderly patients after hip arthroplasty",
        authors=["Smith J"],
        year=2020,
    )
    c2 = Citation(
        id="2",
        title="Postoperative delirium in elderly patients after hip arthroplasty.",
        authors=["Smith J"],
        year=2020,
    )

    unique, matches = dedup.deduplicate([c1, c2])

    assert len(unique) == 1
    assert len(matches) == 1
    assert matches[0].match_type == "fuzzy_title"


def test_no_match_different_year() -> None:
    """Test that different years prevent fuzzy match."""
    dedup = Deduplicator()

    c1 = Citation(
        id="1",
        title="Postoperative delirium study",
        authors=["Smith J"],
        year=2020,
    )
    c2 = Citation(
        id="2",
        title="Postoperative delirium study",
        authors=["Smith J"],
        year=2021,
    )

    unique, matches = dedup.deduplicate([c1, c2])

    assert len(unique) == 2
    assert len(matches) == 0


def test_empty_list() -> None:
    """Test deduplication of empty list."""
    dedup = Deduplicator()
    unique, matches = dedup.deduplicate([])

    assert len(unique) == 0
    assert len(matches) == 0
