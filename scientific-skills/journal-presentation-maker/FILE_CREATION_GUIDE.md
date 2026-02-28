# File Creation Best Practices

This guide outlines best practices for creating new Python files and documentation in this skill.

## Core Principles

### 1. Check Existing Code First

**ALWAYS** before creating any new file:

1. **Search for existing implementations**
   ```bash
   # Find similar files
   glob "**/*.py"

   # Search for related classes/functions
   grep "class.*Parser" --type py
   grep "def search" --type py
   ```

2. **Read and understand existing patterns**
   - Use `Read` tool to examine similar files
   - Identify coding conventions (naming, structure, patterns)
   - Look for reusable components

3. **Check project structure**
   - Review `STRUCTURE.md` if it exists
   - Understand where new files should be placed
   - Follow existing directory organization

### 2. Use Class-Based Architecture

**ALWAYS** prefer class-based design over procedural code.

#### ✅ Good Example: Class-Based

```python
"""
Module for searching academic papers across multiple databases.
"""

from typing import List, Dict, Optional
import logging

logger = logging.getLogger(__name__)


class PaperSearcher:
    """Search and filter academic papers from multiple sources.

    This class handles searching PubMed, Google Scholar, bioRxiv,
    arXiv, and other academic databases with smart filtering.

    Attributes:
        keywords (str): Search keywords
        config (Dict): Search configuration including filters
        results (List[Dict]): Cached search results
    """

    def __init__(self, keywords: str, config: Optional[Dict] = None):
        """Initialize paper searcher.

        Args:
            keywords: Keywords to search for
            config: Optional configuration with filters and limits

        Raises:
            ValueError: If keywords are empty
        """
        if not keywords:
            raise ValueError("Keywords cannot be empty")

        self.keywords = keywords
        self.config = config or {}
        self.results = []
        logger.info(f"Initialized searcher with keywords: {keywords}")

    def search(self) -> List[Dict]:
        """Execute search across all configured databases.

        Returns:
            List of paper dictionaries with metadata
        """
        self.results = []

        # Search each database
        pubmed_results = self._search_pubmed()
        scholar_results = self._search_scholar()
        preprint_results = self._search_preprints()

        # Combine and deduplicate
        self.results = self._deduplicate_results(
            pubmed_results + scholar_results + preprint_results
        )

        return self.results

    def filter_by_date(self, papers: List[Dict]) -> List[Dict]:
        """Filter papers by publication date.

        Args:
            papers: List of papers to filter

        Returns:
            Filtered list of papers
        """
        min_year = self.config.get('min_year', 2020)
        return [p for p in papers if p.get('year', 0) >= min_year]

    def filter_by_citations(self, papers: List[Dict]) -> List[Dict]:
        """Filter papers by citation count.

        Args:
            papers: List of papers to filter

        Returns:
            Filtered list of papers
        """
        min_citations = self.config.get('min_citations', 50)
        return [p for p in papers if p.get('citations', 0) >= min_citations]

    def _search_pubmed(self) -> List[Dict]:
        """Search PubMed database."""
        # Implementation
        pass

    def _search_scholar(self) -> List[Dict]:
        """Search Google Scholar."""
        # Implementation
        pass

    def _search_preprints(self) -> List[Dict]:
        """Search preprint servers."""
        # Implementation
        pass

    def _deduplicate_results(self, papers: List[Dict]) -> List[Dict]:
        """Remove duplicate papers based on DOI."""
        seen_dois = set()
        unique = []

        for paper in papers:
            doi = paper.get('doi')
            if doi and doi not in seen_dois:
                seen_dois.add(doi)
                unique.append(paper)

        return unique
```

#### ❌ Bad Example: Procedural

```python
"""Paper search functions."""

def search_papers(keywords):
    """Search for papers."""
    results = []
    # Global state, hard to test, no type hints
    # ...
    return results

def filter_papers(papers, min_year):
    """Filter papers."""
    # Separated from search logic
    # Hard to maintain state
    return [p for p in papers if p['year'] >= min_year]
```

### 3. Maintain STRUCTURE.md

Keep project structure documentation **always up to date**.

#### When to Update

- ✅ After creating any new file
- ✅ After adding new classes/modules
- ✅ After reorganizing directories
- ✅ After changing file responsibilities

#### STRUCTURE.md Template

```markdown
# Project Structure

Last updated: 2025-11-19

## Directory Layout

```
journal-presentation-maker/
├── SKILL.md                      # Main skill documentation
├── FILE_CREATION_GUIDE.md        # This file
├── STRUCTURE.md                  # Project structure (keep updated!)
├── references/                   # Reference documentation
│   ├── api_reference.md          # Database API docs
│   └── styling_guide.md          # Presentation styling guide
└── scripts/                      # Implementation templates
    ├── search_papers.py          # PaperSearcher class
    └── generate_presentation.py  # PresentationGenerator class
```

## Key Components

### Paper Search Module (`scripts/search_papers.py`)
- **PaperSearcher**: Main class for searching academic databases
  - Methods: `search()`, `filter_by_date()`, `filter_by_citations()`
  - Searches: PubMed, Google Scholar, bioRxiv, arXiv
  - Handles deduplication and smart filtering

### Presentation Generation (`scripts/generate_presentation.py`)
- **PresentationGenerator**: Creates slides with citation tracking
  - Methods: `extract_content()`, `create_slide_structure()`, `generate_html_slides()`
  - Manages reference numbering and speaker notes
  - Outputs HTML files for conversion to PPTX

## Design Patterns

- **Class-based architecture**: All major functionality in classes
- **Type hints**: All public methods have complete type annotations
- **Docstrings**: Google style for all classes and public methods
- **Error handling**: Specific exceptions with context
- **Logging**: Using Python logging module throughout

## Dependencies

- Python 3.8+
- typing module (built-in)
- logging module (built-in)
- No external dependencies in templates (actual implementation uses web_search/web_fetch tools)
```

### 4. Python Code Quality Standards

Every Python file must include:

#### ✅ Complete Type Hints

```python
from typing import List, Dict, Optional, Union, Tuple

def process_papers(
    papers: List[Dict],
    filters: Optional[Dict] = None
) -> Tuple[List[Dict], int]:
    """Process papers with optional filters.

    Args:
        papers: List of paper metadata dictionaries
        filters: Optional filter configuration

    Returns:
        Tuple of (filtered papers, count of removed papers)
    """
    pass
```

#### ✅ Google-Style Docstrings

```python
class CitationManager:
    """Manage citation numbering and reference formatting.

    This class tracks citations across slides, assigns reference numbers
    in order of first appearance, and formats reference lists.

    Attributes:
        papers (List[Dict]): List of papers being cited
        citation_map (Dict[str, int]): Maps DOI to reference number
        first_appearance (Dict[int, str]): Maps ref number to slide ID

    Example:
        >>> manager = CitationManager(papers)
        >>> ref_num = manager.cite(doi="10.1038/xxxxx", slide_id="slide_03")
        >>> print(ref_num)  # 1 (first citation)
    """

    def cite(self, doi: str, slide_id: str) -> int:
        """Add citation and return reference number.

        Assigns a new reference number on first citation, or returns
        existing number for subsequent citations.

        Args:
            doi: Paper DOI being cited
            slide_id: ID of slide containing citation

        Returns:
            Reference number for this paper

        Raises:
            ValueError: If DOI not found in papers list
        """
        pass
```

#### ✅ Proper Error Handling

```python
class PaperFetcher:
    """Fetch paper content from online sources."""

    def fetch(self, doi: str) -> Dict:
        """Fetch paper metadata and content.

        Args:
            doi: Paper DOI

        Returns:
            Paper metadata dictionary

        Raises:
            ValueError: If DOI format is invalid
            ConnectionError: If unable to reach paper source
            NotFoundError: If paper not found at DOI
        """
        if not self._validate_doi(doi):
            raise ValueError(f"Invalid DOI format: {doi}")

        try:
            response = self._make_request(doi)
        except Exception as e:
            raise ConnectionError(f"Failed to fetch {doi}: {e}")

        if response.status_code == 404:
            raise NotFoundError(f"Paper not found: {doi}")

        return self._parse_response(response)
```

#### ✅ Logging Instead of Print

```python
import logging

logger = logging.getLogger(__name__)

class ContentExtractor:
    """Extract content from papers."""

    def extract(self, paper: Dict) -> Dict:
        """Extract structured content from paper.

        Args:
            paper: Paper metadata and content

        Returns:
            Extracted content dictionary
        """
        logger.info(f"Extracting content from: {paper.get('title', 'Unknown')}")

        try:
            content = self._parse_sections(paper)
            logger.debug(f"Extracted {len(content)} sections")
            return content
        except Exception as e:
            logger.error(f"Extraction failed: {e}", exc_info=True)
            raise
```

### 5. Modular Design

Break complex functionality into **small, focused classes** following Single Responsibility Principle.

#### ✅ Good: Focused Classes

```python
# Each class has ONE clear responsibility

class PaperFetcher:
    """Fetch paper content from online sources."""

    def fetch(self, doi: str) -> Dict:
        """Fetch paper by DOI."""
        pass

class ContentParser:
    """Parse paper content into structured sections."""

    def parse(self, raw_content: str) -> Dict:
        """Parse raw content."""
        pass

class FigureExtractor:
    """Extract figures and captions from papers."""

    def extract_figures(self, paper: Dict) -> List[Dict]:
        """Extract all figures."""
        pass

class CitationManager:
    """Manage citation numbering and references."""

    def cite(self, doi: str) -> int:
        """Get reference number for DOI."""
        pass

class SlideGenerator:
    """Generate presentation slides from content."""

    def generate(self, content: Dict) -> List[str]:
        """Generate HTML slides."""
        pass
```

#### ❌ Bad: Monolithic Class

```python
class EverythingProcessor:
    """Does everything related to presentations."""

    def process(self):
        """Does all processing - 500 lines of code."""
        # Fetches papers
        # Parses content
        # Extracts figures
        # Manages citations
        # Generates slides
        # Too many responsibilities!
        pass
```

## Workflow for Creating New Files

### Step-by-Step Process

1. **Search existing code**
   ```
   - glob "**/*.py"
   - grep "class.*[RelatedName]"
   - Read similar files
   ```

2. **Plan the structure**
   - Identify classes needed (keep focused!)
   - Define clear responsibilities
   - Plan method interfaces with types

3. **Create the file with proper structure**
   - Module docstring at top
   - Imports organized (stdlib, third-party, local)
   - Classes with complete docstrings
   - Type hints on all signatures
   - Error handling throughout

4. **Update STRUCTURE.md immediately**
   - Add file to directory layout
   - Document new classes/components
   - Update design patterns if changed

5. **Verify quality**
   - All functions have type hints ✓
   - All classes have docstrings ✓
   - Error handling in place ✓
   - Logging instead of print ✓
   - Follows existing patterns ✓

## Example: Creating a New Module

Let's say we need to add figure downloading capability.

### 1. Check Existing Code

```bash
# Search for related functionality
grep "download\|fetch\|retrieve" --type py
grep "class.*Figure" --type py

# Read existing fetcher
read scripts/search_papers.py
```

### 2. Plan the Module

```
Module: figure_downloader.py
Class: FigureDownloader
Methods:
  - download(url: str, output_path: Path) -> Path
  - download_batch(urls: List[str], output_dir: Path) -> List[Path]
  - _validate_url(url: str) -> bool
  - _get_extension(url: str) -> str
```

### 3. Create the File

```python
"""
Download and save figures from academic papers.

This module handles downloading figures from various sources including
publisher websites, preprint servers, and image hosting services.
"""

from typing import List
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


class FigureDownloader:
    """Download figures from online sources.

    Handles downloading images from paper sources with proper error
    handling, retry logic, and file type validation.

    Attributes:
        output_dir (Path): Default directory for downloads
        timeout (int): Request timeout in seconds
    """

    def __init__(self, output_dir: Path, timeout: int = 30):
        """Initialize figure downloader.

        Args:
            output_dir: Directory to save downloaded figures
            timeout: Download timeout in seconds

        Raises:
            ValueError: If output_dir doesn't exist
        """
        if not output_dir.exists():
            raise ValueError(f"Output directory doesn't exist: {output_dir}")

        self.output_dir = output_dir
        self.timeout = timeout
        logger.info(f"Initialized downloader, output: {output_dir}")

    def download(self, url: str, output_path: Path) -> Path:
        """Download single figure from URL.

        Args:
            url: Figure URL
            output_path: Path to save figure

        Returns:
            Path to downloaded file

        Raises:
            ValueError: If URL is invalid
            ConnectionError: If download fails
        """
        if not self._validate_url(url):
            raise ValueError(f"Invalid URL: {url}")

        logger.info(f"Downloading figure: {url}")

        try:
            # Download implementation
            pass
        except Exception as e:
            logger.error(f"Download failed: {e}")
            raise ConnectionError(f"Failed to download {url}: {e}")

        return output_path

    def download_batch(
        self,
        urls: List[str],
        output_dir: Path
    ) -> List[Path]:
        """Download multiple figures.

        Args:
            urls: List of figure URLs
            output_dir: Directory to save figures

        Returns:
            List of paths to downloaded files
        """
        downloaded = []

        for i, url in enumerate(urls):
            try:
                output_path = output_dir / f"figure_{i:03d}{self._get_extension(url)}"
                path = self.download(url, output_path)
                downloaded.append(path)
            except Exception as e:
                logger.warning(f"Skipped {url}: {e}")

        logger.info(f"Downloaded {len(downloaded)}/{len(urls)} figures")
        return downloaded

    def _validate_url(self, url: str) -> bool:
        """Validate URL format."""
        # Implementation
        return True

    def _get_extension(self, url: str) -> str:
        """Extract file extension from URL."""
        # Implementation
        return ".png"
```

### 4. Update STRUCTURE.md

```markdown
## Key Components

### Figure Download Module (`scripts/figure_downloader.py`)
- **FigureDownloader**: Download figures from online sources
  - Methods: `download()`, `download_batch()`
  - Features: Retry logic, timeout handling, format validation
  - Output: Saves figures to specified directory
```

### 5. Verify Quality Checklist

- ✅ Module docstring
- ✅ All imports organized
- ✅ Complete type hints
- ✅ Google-style docstrings
- ✅ Error handling with specific exceptions
- ✅ Logging throughout
- ✅ STRUCTURE.md updated
- ✅ Follows existing patterns

## Summary

**Always remember:**

1. 🔍 **Check existing code first** - Search before creating
2. 📦 **Use classes** - Class-based, not procedural
3. 📝 **Update STRUCTURE.md** - Keep documentation current
4. 🎯 **Follow standards** - Type hints, docstrings, logging, errors
5. 🧩 **Stay modular** - Small, focused classes

**The goal:** Maintainable, testable, well-documented code that follows consistent patterns.
