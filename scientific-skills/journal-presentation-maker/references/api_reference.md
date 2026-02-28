# Academic Database APIs and Search Sources

This document provides reference information for accessing academic papers from various sources.

## Data Sources

### 1. PubMed / PMC
- **API**: NCBI E-utilities
- **Base URL**: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/`
- **Coverage**: Biomedical and life sciences
- **Endpoints**:
  - `esearch.fcgi`: Search for papers
  - `efetch.fcgi`: Fetch full paper metadata
  - `elink.fcgi`: Find related papers

### 2. Crossref
- **API**: Crossref REST API
- **Base URL**: `https://api.crossref.org/`
- **Coverage**: General academic literature with DOIs
- **Endpoints**:
  - `/works`: Search and retrieve metadata
  - `/works/{doi}`: Get specific paper by DOI

### 3. Google Scholar
- **Method**: Web scraping (no official API)
- **Coverage**: Comprehensive across all disciplines
- **Note**: Use serpapi.com or scholarly Python package

### 4. Preprint Servers

#### bioRxiv / medRxiv
- **API**: bioRxiv/medRxiv API
- **Base URL**: `https://api.biorxiv.org/`
- **Coverage**: Biology and medicine preprints
- **Endpoints**:
  - `/details/{server}/{interval}`: Get papers by date range
  - `/publisher/{server}/{interval}`: Get published status

#### arXiv
- **API**: arXiv API
- **Base URL**: `http://export.arxiv.org/api/`
- **Coverage**: Physics, mathematics, CS, biology
- **Method**: Atom/RSS feeds

#### ChemRxiv
- **Platform**: Figshare-based
- **Coverage**: Chemistry preprints
- **Access**: Through Figshare API

## Search Strategy

### Keyword Search
```
1. Construct query with field-specific search
   - Title: TI[title]
   - Abstract: AB[abstract]
   - Author: AU[author]
   
2. Combine with Boolean operators
   - AND, OR, NOT
   
3. Apply filters
   - Publication date
   - Article type
```

### Citation Metrics
```
Sources for citation counts:
- Google Scholar: Most comprehensive
- Crossref: DOI-based citations
- Semantic Scholar API: Free citation data
```

### Filtering Logic

#### Date-based filtering
```python
def filter_by_date_and_citations(paper):
    year = paper['year']
    citations = paper['citations']
    
    if year >= 2023:
        return True  # Always include recent
    elif 2020 <= year < 2023:
        return citations >= 50  # Medium threshold
    elif year < 2020:
        return citations >= 200  # High threshold
    
    return False
```

#### Preprint handling
```python
def is_preprint(paper):
    preprint_servers = [
        'biorxiv', 'medrxiv', 'arxiv', 
        'chemrxiv', 'ssrn'
    ]
    
    source = paper.get('source', '').lower()
    return any(server in source for server in preprint_servers)
```

## Data Structure

### Paper Metadata Format
```json
{
  "title": "Paper Title",
  "authors": "First Author, Second Author, et al.",
  "first_author": "First Author",
  "year": 2024,
  "journal": "Nature",
  "doi": "10.1038/xxxxx",
  "pmid": "12345678",
  "abstract": "Full abstract text...",
  "citations": 150,
  "is_preprint": false,
  "source": "pubmed",
  "keywords": ["keyword1", "keyword2"],
  "figures": [
    {
      "number": "Figure 1",
      "caption": "Figure caption..."
    }
  ],
  "tables": [
    {
      "number": "Table 1",
      "caption": "Table caption..."
    }
  ]
}
```

## Citation Formats

### APA Style
```
Author, A. A., Author, B. B., & Author, C. C. (Year). Title of article. 
Journal Name, volume(issue), pages. https://doi.org/xxxxx
```

### Nature Style
```
Author, A. A. et al. Title of article. Journal Name volume, pages (Year).
```

### Vancouver Style
```
1. Author AA, Author BB. Title of article. Journal Name. Year;volume(issue):pages.
```

## Rate Limits and Best Practices

### NCBI E-utilities
- Maximum 3 requests/second without API key
- Maximum 10 requests/second with API key
- Use `&email=` parameter in all requests

### Crossref
- Polite pool: 50 requests/second
- Public pool: 1 request/second
- Include email in User-Agent header

### Best Practices
1. Cache results to avoid repeated requests
2. Implement exponential backoff for rate limits
3. Respect robots.txt for web scraping
4. Use batch requests when available
5. Store DOIs for future reference
