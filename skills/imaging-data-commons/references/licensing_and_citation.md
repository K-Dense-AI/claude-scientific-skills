# Licensing and Citation Guide for IDC

## When to Use This Guide

Load this guide when:
- A user asks whether IDC data can be used commercially, redistributed, or included in a product
- You are assembling a cohort that mixes collections and need to know which terms govern it
- A user is publishing results and needs formatted citations (APA, BibTeX, CSL JSON, RDF Turtle)
- You need the parameters or output formats for citation generation

The obligation summary — check the license, generate citations — lives in `SKILL.md`. This
guide holds the detail behind it.

**These are the two IDC tasks least tied to any one access path.** Licenses and citations are
available identically from `idc-index`, the REST API, and the hosted MCP server. Use whichever
route the session is already on rather than installing Python to answer a licensing question,
or dropping out of an MCP session to run a script.

| Task | `idc-index` (Python) | REST API | MCP server |
|------|----------------------|----------|------------|
| License breakdown for a selection | `sql_query` on `license_short_name` | `POST /v3/licenses` | `get_licenses` |
| Citations for a selection | `citations_from_selection()` | `POST /v3/citations` | `get_citations` |

Route-specific detail lives in `references/rest_api_guide.md` (endpoint reference, filter
syntax, and the body-shape pitfall that makes a mis-shaped filter return all of IDC) and
`references/mcp_guide.md` (tool inventory). The license semantics below apply to all three.

## Licenses in IDC

Every DICOM file in IDC is tagged with its license in the file metadata, and every row in the
`index` table carries a `license_short_name` column. There is no single IDC-wide license.

| License | Share of data | Commercial use | Attribution required |
|---------|---------------|----------------|----------------------|
| CC BY 4.0 | 74.7% | Yes | Yes |
| CC BY 3.0 | 22.1% | Yes | Yes |
| CC BY-NC 4.0 | 2.1% | **No** | Yes |
| CC BY-NC 3.0 | 0.8% | **No** | Yes |
| NLM Terms and Conditions | 0.3% | Read the terms | Yes |

About 97% of IDC data by size permits commercial reuse; just under 3% is non-commercial. Treat
any `license_short_name` that is not a recognizable Creative Commons string as custom, and
report the exact value to the user rather than assuming it permits reuse.

**Licenses attach to individual series, not to whole collections.** 39 of IDC's 179 collections
carry more than one license — analysis results and original images within one collection can
differ, as can series from different sources. Never conclude that a collection is
commercially usable from one series, or from the collection's headline license: group by
`license_short_name` over the exact selection you intend to use.

**When a cohort mixes licenses, the most restrictive term governs the combined dataset.** If a
selection contains any CC BY-NC series, either drop those series or tell the user the whole
derived dataset is non-commercial.

Commercially restricted data is also physically separated in cloud storage: the
`idc-open-data-cr` (AWS) / `idc-open-cr` (GCS) buckets hold the CC BY-NC collections. See
`references/cloud_storage_guide.md` for bucket details.

## Checking licenses

### Via `idc-index`

```python
from idc_index import IDCClient
client = IDCClient()

# Licenses across all collections
licenses = client.sql_query("""
    SELECT DISTINCT
      collection_id,
      license_short_name,
      COUNT(DISTINCT SeriesInstanceUID) as series_count
    FROM index
    GROUP BY collection_id, license_short_name
    ORDER BY collection_id
""")
print(licenses)
```

```python
# Licenses present in one specific cohort — run this before handing a dataset to a user
cohort_licenses = client.sql_query("""
    SELECT license_short_name, COUNT(DISTINCT SeriesInstanceUID) as series_count
    FROM index
    WHERE Modality = 'MR' AND BodyPartExamined = 'BREAST'
    GROUP BY license_short_name
""")
print(cohort_licenses)
```

```python
# Commercial-safe subset: exclude non-commercial collections outright
commercial_ok = client.sql_query("""
    SELECT collection_id, SeriesInstanceUID
    FROM index
    WHERE Modality = 'CT'
      AND license_short_name NOT LIKE '%NC%'
    LIMIT 20
""")
```

### Via the REST API

`POST /v3/licenses` takes the filter object **directly** (not wrapped in a `filters` key) and
returns the per-license breakdown with series counts and sizes:

```bash
B=https://api.imaging.datacommons.cancer.gov/v3
curl -s $B/licenses \
  -H 'content-type: application/json' \
  -d '{"terms": {"Modality": ["MR"], "BodyPartExamined": ["BREAST"]}}'
```

Response shape: `licenses[{license_short_name, series, size_TB}]`. A collection's licenses are
also included in `GET /v3/collections/{id}`.

### Via the MCP server

Call `get_licenses` with the same selection you built with `build_cohort`. The result carries
the same per-license breakdown; the CC BY vs CC BY-NC distinction above applies unchanged.

## Citations and attribution

The `source_DOI` column links to the publications describing how each dataset was generated.
All three routes turn a selection into formatted citations that satisfy the attribution
requirement common to every IDC license.

Generate citations from the *same* selection you downloaded, not from the collection as a
whole — a five-series subset of a collection that spans several source publications should
cite only the publications it actually draws on.

### Via `idc-index`

```python
# Citations for a collection (APA is the default format)
citations = client.citations_from_selection(collection_id="rider_pilot")
for citation in citations:
    print(citation)
```

```python
# Citations for a specific set of series — matches what you actually downloaded
results = client.sql_query("""
    SELECT SeriesInstanceUID FROM index
    WHERE collection_id = 'tcga_luad' LIMIT 5
""")
citations = client.citations_from_selection(
    seriesInstanceUID=list(results['SeriesInstanceUID'].values)
)
```

```python
# BibTeX, for LaTeX manuscripts
bibtex_citations = client.citations_from_selection(
    collection_id="tcga_luad",
    citation_format=IDCClient.CITATION_FORMAT_BIBTEX
)
```

`citations_from_selection()` takes the same selection filters as the download methods —
`collection_id`, `patientId`, `studyInstanceUID`, `seriesInstanceUID` — plus `citation_format`.

### Via the REST API

`POST /v3/citations` **wraps** the filter in a `filters` key (unlike `/v3/licenses` — this
asymmetry is the single most common REST mistake; see `references/rest_api_guide.md`):

```bash
curl -s $B/citations \
  -H 'content-type: application/json' \
  -d '{"filters": {"terms": {"collection_id": ["rider_pilot"]}}, "citation_format": "bibtex"}'
```

The response separates the per-dataset `citations[]` from `idc_acknowledgment` (the IDC paper)
and `recommendation`. Include both parts — see *What to include when publishing* below.

### Via the MCP server

Call `get_citations` for a selection. It returns the per-dataset citations plus the IDC paper
to acknowledge IDC itself, matching the REST response.

### Citation formats

| `idc-index` constant | REST / MCP `citation_format` | Output |
|----------------------|------------------------------|--------|
| `IDCClient.CITATION_FORMAT_APA` | `apa` (default) | APA string |
| `IDCClient.CITATION_FORMAT_BIBTEX` | `bibtex` | BibTeX entry, for LaTeX |
| `IDCClient.CITATION_FORMAT_JSON` | `csl-json` | CSL JSON |
| `IDCClient.CITATION_FORMAT_TURTLE` | `turtle` | RDF Turtle |

## What to include when publishing

1. **The dataset citations** for every collection or series set used.
2. **The IDC data version** — `client.get_idc_version()`, `GET /v3/version`, or the MCP
   `get_idc_version` tool. IDC releases are versioned and series are added and revised between
   them, so the version is what makes the selection reproducible.
3. **The IDC platform citation**, to acknowledge IDC itself. The REST and MCP routes return
   this as `idc_acknowledgment`; when using `idc-index`, add it yourself:

   > Fedorov, A., et al. "National Cancer Institute Imaging Data Commons: Toward Transparency,
   > Reproducibility, and Scalability in Imaging Artificial Intelligence." *RadioGraphics* 43.12
   > (2023). https://doi.org/10.1148/rg.230180

4. **The series manifest** — save the `SeriesInstanceUID` list alongside the analysis so the
   exact cohort can be rebuilt.

## Troubleshooting

### Issue: Fewer citations returned than collections selected

- **Cause:** Citations are derived from `source_DOI`, and several collections can share one
  DOI, so a multi-collection selection may legitimately produce a shorter list.
- **Solution:** Query `SELECT DISTINCT collection_id, source_DOI FROM index WHERE ...` to see
  the mapping directly.

### Issue: `POST /v3/citations` returns citations for all of IDC

- **Cause:** The filter was passed directly instead of wrapped in `filters`. `/v3/licenses`
  takes the filter directly; `/v3/citations` wraps it. A mis-shaped body is not an error — it
  is treated as an empty filter.
- **Solution:** Check the response counts against a `POST /v3/cohort/counts` for the same
  selection. See `references/rest_api_guide.md`.

## Resources

- **IDC Portal** — https://portal.imaging.datacommons.cancer.gov/
- **IDC data licensing documentation** — https://learn.canceridc.dev/data/licensing
- **`references/rest_api_guide.md`** — `/v3/licenses` and `/v3/citations` endpoint reference
- **`references/mcp_guide.md`** — `get_licenses` and `get_citations` tool inventory
- **`references/cloud_storage_guide.md`** — bucket separation for commercially restricted data
