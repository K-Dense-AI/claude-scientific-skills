# Project Structure

Last updated: 2025-11-19

## Directory Layout

```
journal-presentation-maker/
├── SKILL.md                      # Main skill documentation and workflow
├── FILE_CREATION_GUIDE.md        # Best practices for creating new files
├── STRUCTURE.md                  # This file - project structure overview
├── references/                   # Reference documentation
│   ├── api_reference.md          # Database API docs, search strategies
│   └── styling_guide.md          # Presentation styling guide
└── scripts/                      # Implementation templates
    ├── search_papers.py          # Paper search (procedural - legacy)
    ├── search_papers_v2.py       # Paper search (class-based - NEW)
    └── generate_presentation.py  # Slide generation with citations
```

## Key Components

### Paper Search Module v2 (`scripts/search_papers_v2.py`) ⭐ RECOMMENDED

**PaperSearcher Class** (Class-based):
- `__init__()`: Initialize with keywords and thresholds
- `search()`: Search databases and cache results
- `filter_papers()`: Apply smart filtering
- `format_paper_list()`: Format for display

**Features**:
- Complete type hints and docstrings
- Error handling with specific exceptions
- Logging (no print statements)
- State management

**Filtering Strategy**:
- 2023-2025: Include all
- 2020-2022: Include if citations ≥ 50
- Pre-2020: Include if citations ≥ 200
- Preprints: Always include

### Paper Search Module (Legacy) (`scripts/search_papers.py`)

**Functions** (Procedural):
- `search_papers()`, `filter_papers()`, `format_paper_list()`

**Note**: Use v2 for new code.

### Presentation Generation (`scripts/generate_presentation.py`)

**PresentationGenerator Class**:
- `extract_content()`: Extract and organize content from papers
- `create_slide_structure()`: Build slide layout with sections
- `generate_html_slides()`: Create HTML files for slides
- `_generate_notes()`: Create speaker notes with source tracking
- `_format_references()`: Format reference list in standard format

**Features**:
- Citation tracking with reference numbers [1,2,3]
- Figure extraction with original captions
- Speaker notes with detailed source attribution
- Support for max 30 slides per presentation

## Design Patterns

### Current State
- **generate_presentation.py**: ✅ Class-based architecture
- **search_papers.py**: ⚠️ Procedural (needs refactoring to class-based)

### Standards
- **Architecture**: Class-based for all major functionality
- **Type Hints**: Complete type annotations on all public methods
- **Docstrings**: Google style for classes and public methods
- **Error Handling**: Specific exceptions with context
- **Logging**: Python logging module (not print)
- **Configuration**: Separate config from code

## File Organization

### Documentation Files
- **SKILL.md**: Main workflow and usage instructions
- **FILE_CREATION_GUIDE.md**: Coding standards and best practices
- **STRUCTURE.md**: This file - current project structure
- **references/*.md**: API and styling reference materials

### Implementation Files
- **scripts/*.py**: Python templates showing code structure
- Note: Actual implementation uses Claude's web_search and web_fetch tools

## Dependencies

### Python Requirements
- Python 3.8+
- typing module (built-in)
- logging module (built-in)
- json module (built-in)

### External Dependencies (for actual implementation)
- None in templates (placeholders only)
- Actual skill uses Claude's built-in tools:
  - `web_search`: Search academic databases
  - `web_fetch`: Fetch paper content
  - `pptx skill`: Generate PowerPoint presentations

## Planned Improvements

### Near Term
1. Refactor `search_papers.py` to class-based architecture
2. Add `PaperSearcher` class with proper methods
3. Add comprehensive error handling throughout
4. Add logging to all modules

### Future Enhancements
1. Add `FigureDownloader` class for figure retrieval
2. Add `CitationManager` class for reference tracking
3. Add `ContentParser` class for paper content extraction
4. Create unit tests for all components

## Update Protocol

**When creating/modifying files**, update this file immediately:

1. Add new files to directory layout
2. Document new classes and their methods
3. Update design patterns section
4. Add to dependencies if needed
5. Note any architectural changes

**Example update**:
```markdown
### New Module (`scripts/figure_downloader.py`)
- **FigureDownloader**: Download figures from online sources
  - Methods: `download()`, `download_batch()`
  - Features: Retry logic, timeout handling
```

## Contact & Maintenance

This skill is maintained as part of Claude Code's skill system.

For updates or questions, refer to:
- SKILL.md for usage instructions
- FILE_CREATION_GUIDE.md for coding standards
