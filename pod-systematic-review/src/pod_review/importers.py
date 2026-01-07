"""Import citations from various file formats (RIS, BibTeX, CSV)."""

import csv
import logging
from pathlib import Path
from typing import Any

import bibtexparser
import rispy

from .models import Citation

logger = logging.getLogger(__name__)


class RISImporter:
    """Import citations from RIS format."""

    def import_file(self, file_path: Path) -> list[Citation]:
        """
        Import citations from RIS file.

        Args:
            file_path: Path to RIS file

        Returns:
            List of Citations
        """
        logger.info(f"Importing RIS file: {file_path}")

        with open(file_path, "r", encoding="utf-8") as f:
            entries = rispy.load(f)

        citations = []
        for entry in entries:
            citation = self._parse_ris_entry(entry)
            citations.append(citation)

        logger.info(f"Imported {len(citations)} citations from RIS")
        return citations

    def _parse_ris_entry(self, entry: dict[str, Any]) -> Citation:
        """Parse RIS entry to Citation."""
        # RIS field mappings
        title = entry.get("title", entry.get("primary_title", ""))
        authors = entry.get("authors", [])
        year_str = entry.get("year", entry.get("publication_year", ""))
        year = None
        if year_str:
            try:
                year = int(year_str)
            except (ValueError, TypeError):
                pass

        # Generate ID from DOI or first author + year
        doi = entry.get("doi", "")
        pmid = entry.get("pmid", "")
        if doi:
            citation_id = doi
        elif pmid:
            citation_id = f"PMID{pmid}"
        elif authors and year:
            citation_id = f"{authors[0].split()[-1]}{year}"
        else:
            citation_id = str(hash(title))[:16]

        return Citation(
            id=citation_id,
            title=title,
            authors=authors,
            year=year,
            journal=entry.get("journal_name", entry.get("secondary_title", "")),
            volume=entry.get("volume", ""),
            issue=entry.get("number", ""),
            pages=entry.get("start_page", ""),
            abstract=entry.get("abstract", ""),
            doi=doi,
            pmid=pmid,
            keywords=entry.get("keywords", []),
            url=entry.get("url", ""),
            language=entry.get("language", ""),
            source_database="RIS_Import",
        )


class BibTeXImporter:
    """Import citations from BibTeX format."""

    def import_file(self, file_path: Path) -> list[Citation]:
        """
        Import citations from BibTeX file.

        Args:
            file_path: Path to BibTeX file

        Returns:
            List of Citations
        """
        logger.info(f"Importing BibTeX file: {file_path}")

        with open(file_path, "r", encoding="utf-8") as f:
            bib_database = bibtexparser.load(f)

        citations = []
        for entry in bib_database.entries:
            citation = self._parse_bibtex_entry(entry)
            citations.append(citation)

        logger.info(f"Imported {len(citations)} citations from BibTeX")
        return citations

    def _parse_bibtex_entry(self, entry: dict[str, Any]) -> Citation:
        """Parse BibTeX entry to Citation."""
        # Parse authors
        authors_str = entry.get("author", "")
        authors = [a.strip() for a in authors_str.split(" and ")] if authors_str else []

        # Parse year
        year = None
        year_str = entry.get("year", "")
        if year_str:
            try:
                year = int(year_str)
            except ValueError:
                pass

        # Citation ID
        citation_id = entry.get("ID", entry.get("doi", str(hash(entry.get("title", "")))[:16]))

        return Citation(
            id=citation_id,
            title=entry.get("title", ""),
            authors=authors,
            year=year,
            journal=entry.get("journal", ""),
            volume=entry.get("volume", ""),
            issue=entry.get("number", ""),
            pages=entry.get("pages", ""),
            abstract=entry.get("abstract", ""),
            doi=entry.get("doi", ""),
            pmid=entry.get("pmid", ""),
            url=entry.get("url", ""),
            source_database="BibTeX_Import",
        )


class CSVImporter:
    """Import citations from CSV format."""

    def __init__(self, column_mapping: dict[str, str] | None = None) -> None:
        """
        Initialize CSV importer.

        Args:
            column_mapping: Mapping from CSV column names to Citation fields
                           If None, uses default mapping
        """
        self.column_mapping = column_mapping or {
            "Title": "title",
            "Authors": "authors",
            "Year": "year",
            "Journal": "journal",
            "Abstract": "abstract",
            "DOI": "doi",
            "PMID": "pmid",
            "Volume": "volume",
            "Issue": "issue",
            "Pages": "pages",
            "URL": "url",
        }

    def import_file(self, file_path: Path) -> list[Citation]:
        """
        Import citations from CSV file.

        Args:
            file_path: Path to CSV file

        Returns:
            List of Citations
        """
        logger.info(f"Importing CSV file: {file_path}")

        citations = []
        with open(file_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                citation = self._parse_csv_row(row)
                citations.append(citation)

        logger.info(f"Imported {len(citations)} citations from CSV")
        return citations

    def _parse_csv_row(self, row: dict[str, str]) -> Citation:
        """Parse CSV row to Citation."""
        # Map CSV columns to Citation fields
        data: dict[str, Any] = {}

        for csv_col, citation_field in self.column_mapping.items():
            if csv_col in row and row[csv_col]:
                value = row[csv_col].strip()

                # Special handling for certain fields
                if citation_field == "authors":
                    # Split authors by semicolon or comma
                    data[citation_field] = [
                        a.strip() for a in value.replace(";", ",").split(",")
                    ]
                elif citation_field == "year":
                    try:
                        data[citation_field] = int(value)
                    except ValueError:
                        data[citation_field] = None
                else:
                    data[citation_field] = value

        # Generate ID
        doi = data.get("doi", "")
        pmid = data.get("pmid", "")
        title = data.get("title", "")

        if doi:
            citation_id = doi
        elif pmid:
            citation_id = f"PMID{pmid}"
        else:
            citation_id = str(hash(title))[:16]

        data["id"] = citation_id
        data["source_database"] = "CSV_Import"

        return Citation(**data)


def import_citations(file_path: Path, file_format: str | None = None) -> list[Citation]:
    """
    Import citations from file, auto-detecting format if needed.

    Args:
        file_path: Path to import file
        file_format: File format (ris, bibtex, csv). If None, detect from extension

    Returns:
        List of Citations
    """
    if file_format is None:
        # Auto-detect from extension
        suffix = file_path.suffix.lower()
        if suffix in [".ris", ".txt"]:
            file_format = "ris"
        elif suffix in [".bib", ".bibtex"]:
            file_format = "bibtex"
        elif suffix == ".csv":
            file_format = "csv"
        else:
            raise ValueError(f"Cannot determine format for file: {file_path}")

    if file_format == "ris":
        importer = RISImporter()
    elif file_format == "bibtex":
        importer = BibTeXImporter()
    elif file_format == "csv":
        importer = CSVImporter()
    else:
        raise ValueError(f"Unsupported format: {file_format}")

    return importer.import_file(file_path)
