"""Export citations to various formats (RIS, BibTeX, CSV)."""

import csv
import logging
from pathlib import Path
from typing import Any

from .models import Citation

logger = logging.getLogger(__name__)


def export_to_ris(citations: list[Citation], output_path: Path) -> None:
    """
    Export citations to RIS format.

    Args:
        citations: List of citations
        output_path: Output file path
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as f:
        for citation in citations:
            # RIS format
            f.write("TY  - JOUR\n")  # Journal article

            # Title
            if citation.title:
                f.write(f"TI  - {citation.title}\n")

            # Authors
            for author in citation.authors:
                f.write(f"AU  - {author}\n")

            # Year
            if citation.year:
                f.write(f"PY  - {citation.year}\n")

            # Journal
            if citation.journal:
                f.write(f"JO  - {citation.journal}\n")
                f.write(f"T2  - {citation.journal}\n")

            # Volume
            if citation.volume:
                f.write(f"VL  - {citation.volume}\n")

            # Issue
            if citation.issue:
                f.write(f"IS  - {citation.issue}\n")

            # Pages
            if citation.pages:
                f.write(f"SP  - {citation.pages}\n")

            # Abstract
            if citation.abstract:
                f.write(f"AB  - {citation.abstract}\n")

            # DOI
            if citation.doi:
                f.write(f"DO  - {citation.doi}\n")

            # PMID
            if citation.pmid:
                f.write(f"AN  - {citation.pmid}\n")

            # URL
            if citation.url:
                f.write(f"UR  - {citation.url}\n")

            # Keywords
            for keyword in citation.keywords:
                f.write(f"KW  - {keyword}\n")

            # Language
            if citation.language:
                f.write(f"LA  - {citation.language}\n")

            # Database
            f.write(f"DB  - {citation.source_database}\n")

            # End of record
            f.write("ER  - \n\n")

    logger.info(f"Exported {len(citations)} citations to RIS: {output_path}")


def export_to_csv(citations: list[Citation], output_path: Path) -> None:
    """
    Export citations to CSV format.

    Args:
        citations: List of citations
        output_path: Output file path
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "ID",
        "Title",
        "Authors",
        "Year",
        "Journal",
        "Volume",
        "Issue",
        "Pages",
        "Abstract",
        "DOI",
        "PMID",
        "URL",
        "Keywords",
        "Language",
        "Source_Database",
    ]

    with open(output_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for citation in citations:
            writer.writerow(
                {
                    "ID": citation.id,
                    "Title": citation.title,
                    "Authors": "; ".join(citation.authors),
                    "Year": citation.year or "",
                    "Journal": citation.journal,
                    "Volume": citation.volume,
                    "Issue": citation.issue,
                    "Pages": citation.pages,
                    "Abstract": citation.abstract,
                    "DOI": citation.doi,
                    "PMID": citation.pmid,
                    "URL": citation.url,
                    "Keywords": "; ".join(citation.keywords),
                    "Language": citation.language,
                    "Source_Database": citation.source_database,
                }
            )

    logger.info(f"Exported {len(citations)} citations to CSV: {output_path}")


def export_to_bibtex(citations: list[Citation], output_path: Path) -> None:
    """
    Export citations to BibTeX format.

    Args:
        citations: List of citations
        output_path: Output file path
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as f:
        for citation in citations:
            # Create citation key
            if citation.authors and citation.year:
                first_author = citation.authors[0].split()[-1]
                cite_key = f"{first_author}{citation.year}"
            elif citation.doi:
                cite_key = citation.doi.replace("/", "_").replace(".", "_")
            else:
                cite_key = citation.id

            f.write(f"@article{{{cite_key},\n")

            # Required and optional fields
            if citation.title:
                f.write(f'  title = {{{citation.title}}},\n')

            if citation.authors:
                authors_str = " and ".join(citation.authors)
                f.write(f'  author = {{{authors_str}}},\n')

            if citation.year:
                f.write(f"  year = {{{citation.year}}},\n")

            if citation.journal:
                f.write(f'  journal = {{{citation.journal}}},\n')

            if citation.volume:
                f.write(f"  volume = {{{citation.volume}}},\n")

            if citation.issue:
                f.write(f"  number = {{{citation.issue}}},\n")

            if citation.pages:
                f.write(f"  pages = {{{citation.pages}}},\n")

            if citation.abstract:
                # Escape special characters in abstract
                abstract = citation.abstract.replace("{", "\\{").replace("}", "\\}")
                f.write(f'  abstract = {{{abstract}}},\n')

            if citation.doi:
                f.write(f"  doi = {{{citation.doi}}},\n")

            if citation.pmid:
                f.write(f"  pmid = {{{citation.pmid}}},\n")

            if citation.url:
                f.write(f'  url = {{{citation.url}}},\n')

            f.write("}\n\n")

    logger.info(f"Exported {len(citations)} citations to BibTeX: {output_path}")
