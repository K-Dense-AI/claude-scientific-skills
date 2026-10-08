# Index Tables Guide for IDC

**Tested with:** idc-index 0.13.0 (IDC data version v25)

This guide covers the structure and access patterns for IDC index tables: programmatic schema discovery, DataFrame access, and join column references. For the overview of available tables and their purposes, see the "Index Tables" section in the main SKILL.md.

**Complete index table documentation:** https://idc-index.readthedocs.io/en/latest/indices_reference.html

## When to Use This Guide

Load this guide when you need to:
- Discover table schemas and column types programmatically
- Access index tables as pandas DataFrames (not via SQL)
- Understand key columns and join relationships between tables

For SQL query examples (filter discovery, finding annotations, size estimation), see `references/sql_patterns.md`.

## Prerequisites

Needs `idc-index` installed — run `scripts/check_version.py`, which reports the installed
version and prints the install command for the interpreter you are running.

## Available Tables

`SKILL.md` carries a compact map of the table families. This is the full inventory with row
granularity and contents. Always call `client.fetch_index("table_name")` before querying any
of them — it is safe and idempotent for all tables, including those loaded automatically at
startup.

| Table | Row Granularity | Description |
|-------|-----------------|-------------|
| `index` | 1 row = 1 DICOM series | Primary metadata for all current IDC data |
| `version_metadata_index` | 1 row = 1 IDC release version | IDC version release timestamps; join on `idc_version` to correlate series with their release date |
| `collections_index` | 1 row = 1 collection | Collection-level metadata and descriptions |
| `analysis_results_index` | 1 row = 1 analysis result collection | Metadata about derived datasets (annotations, segmentations) |
| `clinical_index` | 1 row = 1 (collection, table, column) triple | Dictionary mapping clinical data table columns to collections |
| `sm_index` | 1 row = 1 slide microscopy series | Slide Microscopy (pathology) series metadata |
| `sm_instance_index` | 1 row = 1 slide microscopy instance | Instance-level (SOPInstanceUID) metadata for slide microscopy |
| `seg_index` | 1 row = 1 DICOM Segmentation series | Segmentation metadata: algorithm, segment count, reference to source image series |
| `ann_index` | 1 row = 1 DICOM ANN series | Microscopy Bulk Simple Annotations series metadata; references annotated image series |
| `ann_group_index` | 1 row = 1 annotation group | Detailed annotation group metadata: graphic type, annotation count, property codes, algorithm |
| `contrast_index` | 1 row = 1 series with contrast info | Contrast agent metadata: agent name, ingredient, administration route (CT, MR, PT, XA, RF) |
| `volume_geometry_index` | 1 row = 1 CT/MR/PT series | 3D volume geometry validation for single-frame CT, MR, and PT series; boolean checks for orientation, spacing, dimensions, and slice positions; composite `regularly_spaced_3d_volume` flag |
| `rtstruct_index` | 1 row = 1 RTSTRUCT series | RT Structure Set metadata: total ROI count, ROI names, generation algorithms, interpreted types, and the referenced image series UID |
| `ct_index` | 1 row = 1 CT series | CT acquisition/reconstruction parameters: pixel spacing, slice thickness, kVp, convolution kernel, tube current (min/max for dose-modulated), exposure, spiral pitch, scan options |
| `mr_index` | 1 row = 1 MR series | MR acquisition/sequence parameters: field strength, scanning sequence, TE (array for multi-echo), TR, flip angle, DiffusionBValue (array for DWI), pixel bandwidth, receive coil, number of temporal positions |
| `pt_index` | 1 row = 1 PET series | PET acquisition/reconstruction/radiopharmaceutical parameters: series type, units, decay/scatter/attenuation correction, reconstruction method, radionuclide, injected dose, frame duration (array for dynamic PET) |
| `prior_versions_index` | 1 row = 1 DICOM series | **Reproducibility only.** Contains series permanently removed from IDC (all `max_idc_version` < current version; zero overlap with `index`). Use ONLY when a user explicitly needs to reproduce work from a prior IDC version using data no longer in the current release. Do NOT use for version history or "what's new" questions — those use `series_init_idc_version`/`series_revised_idc_version` in the main `index` table. Column names `min_idc_version`/`max_idc_version` here are NOT equivalent to `series_init_idc_version`/`series_revised_idc_version` in `index`. |

## Accessing Index Tables

### Via SQL (recommended for filtering/aggregation)

```python
from idc_index import IDCClient
client = IDCClient()

# Query the primary index (always available)
results = client.sql_query("SELECT * FROM index WHERE Modality = 'CT' LIMIT 10")

# Fetch and query additional indices
client.fetch_index("collections_index")
collections = client.sql_query("SELECT collection_id, cancer_types, tumor_locations FROM collections_index")

client.fetch_index("analysis_results_index")
analysis = client.sql_query("SELECT * FROM analysis_results_index LIMIT 5")
```

### As pandas DataFrames (direct access)

```python
# Primary index (always available after client initialization)
df = client.index

# Fetch and access on-demand indices
client.fetch_index("sm_index")
sm_df = client.sm_index
```

## Discovering Table Schemas

The `indices_overview` dictionary contains complete schema information for all tables. **Always consult this when writing queries or exploring data structure.**

**DICOM attribute mapping:** Many columns are populated directly from DICOM attributes in the source files. The column description in the schema indicates when a column corresponds to a DICOM attribute (e.g., "DICOM Modality attribute" or references a DICOM tag). This allows leveraging DICOM knowledge when querying — standard DICOM attribute names like `PatientID`, `StudyInstanceUID`, `Modality`, `BodyPartExamined` work as expected.

```python
from idc_index import IDCClient
client = IDCClient()

# List all available indices with descriptions
for name, info in client.indices_overview.items():
    print(f"\n{name}:")
    print(f"  Installed: {info['installed']}")
    print(f"  Description: {info['description']}")

# Get complete schema for a specific index (columns, types, descriptions)
schema = client.indices_overview["index"]["schema"]
print(f"\nTable: {schema['table_description']}")
print("\nColumns:")
for col in schema['columns']:
    desc = col.get('description', 'No description')
    # Description indicates if column is from DICOM attribute
    print(f"  {col['name']} ({col['type']}): {desc}")

# Find columns that are DICOM attributes (check description for "DICOM" reference)
dicom_cols = [c['name'] for c in schema['columns'] if 'DICOM' in c.get('description', '').upper()]
print(f"\nDICOM-sourced columns: {dicom_cols}")
```

**Alternative: use `get_index_schema()` method:**
```python
schema = client.get_index_schema("index")
# Returns same schema dict: {'table_description': ..., 'columns': [...]}
```

### Finding which table contains a column

The most common schema question is "where does `SliceThickness` live?" — the primary `index`
holds series-level metadata only, so modality-specific acquisition parameters are in dedicated
tables. Search the overview rather than guessing; neither call fetches anything:

```python
# Find which table(s) contain a specific column (no fetch required)
target = "SliceThickness"
for table_name, info in client.indices_overview.items():
    if any(c["name"] == target for c in info["schema"]["columns"]):
        print(f"'{target}' is in: {table_name}")
# → 'SliceThickness' is in: ct_index

# List all columns in a table from the schema (no fetch required)
ct_cols = [c["name"] for c in client.indices_overview["ct_index"]["schema"]["columns"]]
print("ct_index columns:", ct_cols)
# → ['SeriesInstanceUID', 'PixelSpacing_row_mm', 'PixelSpacing_col_mm', 'Rows',
#    'Columns', 'SliceThickness', 'KVP', 'ConvolutionKernel', ...]
```

Then `client.fetch_index("ct_index")` and join to `index` on `SeriesInstanceUID`.

## Key Columns Reference

Most common columns in the primary `index` table (use `indices_overview` for complete list and descriptions):

| Column | Type | DICOM | Description |
|--------|------|-------|-------------|
| `collection_id` | STRING | No | IDC collection identifier |
| `analysis_result_id` | STRING | No | If applicable, indicates what analysis results collection given series is part of |
| `source_DOI` | STRING | No | DOI linking to dataset details; use for learning more about the content and for attribution (see citations below) |
| `PatientID` | STRING | Yes | Patient identifier |
| `StudyInstanceUID` | STRING | Yes | DICOM Study UID |
| `SeriesInstanceUID` | STRING | Yes | DICOM Series UID — use for downloads/viewing |
| `Modality` | STRING | Yes | Imaging modality (CT, MR, PT, SM, SEG, ANN, RTSTRUCT, etc.) |
| `BodyPartExamined` | STRING | Yes | Anatomical region |
| `SeriesDescription` | STRING | Yes | Description of the series |
| `Manufacturer` | STRING | Yes | Equipment manufacturer |
| `StudyDate` | STRING | Yes | Date study was performed |
| `PatientSex` | STRING | Yes | Patient sex |
| `PatientAge` | STRING | Yes | Patient age at time of study |
| `license_short_name` | STRING | No | License type (CC BY 4.0, CC BY-NC 4.0, etc.) |
| `series_size_MB` | FLOAT | No | Size of series in megabytes |
| `instanceCount` | INTEGER | No | Number of DICOM instances in series |
| `SOPClassUID` | STRING | Yes | DICOM SOP Class UID (identifies the object/service class, e.g., CT Image Storage) |
| `TransferSyntaxUID` | STRING | Yes | DICOM Transfer Syntax UID (encoding/compression method) |
| `sop_class_name` | STRING | No | Human-readable `SOPClassUID` (e.g., "CT Image Storage", "Segmentation Storage"); the filterable form, and more specific than `Modality` for distinguishing object types |
| `transfer_syntax_name` | STRING | No | Human-readable `TransferSyntaxUID` (e.g., "JPEG 2000", "Explicit VR Little Endian") |
| `PhotometricInterpretation` | STRING | Yes | Intended interpretation of pixel data — `MONOCHROME2` (grayscale), `RGB` / `YBR_FULL_422` (color) |
| `PixelRepresentation` | STRING | Yes | Pixel sample representation: `"0"` unsigned, `"1"` signed |

**Pixel-encoding columns are per-series summaries, not per-instance values.** `TransferSyntaxUID`,
`transfer_syntax_name`, `PhotometricInterpretation`, and `PixelRepresentation` are
**comma-separated** when instances within one series differ — about 72 000 series in v25, mostly
SM — and the latter two are `NULL` for non-image objects (SEG, SR, RTSTRUCT). Match them with
`LIKE '%...%'` or split before comparing; an `=` test silently drops every mixed-encoding series.
Filter on `sop_class_name` rather than parsing `SOPClassUID`: it is one of the 19 attributes the
cohort and REST filter APIs accept, while the raw UID is not.

**DICOM = Yes**: Column value extracted from the DICOM attribute with the same name. Refer to the [DICOM standard](https://dicom.nema.org/medical/dicom/current/output/chtml/part06/chapter_6.html) for numeric tag mappings. Use standard DICOM knowledge for expected values and formats.

## Data Provenance (new in v25)

IDC v25 added a structured, queryable record of **where each piece of data came from and who
handled it** — previously answerable only by heuristics over `collection_id` and DOI strings.
It lives in two places, both as a `provenance` STRUCT with the same four fields:

| Field | Meaning |
|-------|---------|
| `data_contributor` | Who contributed the data to IDC (`TCIA`, `IDC`) |
| `source_data_provider` | Who held the upstream source material (`GDC`, `TCIA`, `HTAN`, `NLM`, …) |
| `deidentification_party` | Who performed de-identification (`TCIA`, `GDC`, `HTAN`, `NCH`, …) |
| `dicom_conversion_by` | Who produced the DICOM representation |

**Provenance is per-source, not per-collection.** A collection routinely mixes components with
different lineage — TCGA radiology arrives as DICOM from TCIA while its pathology originates as
vendor SVS from GDC and is converted by IDC under a separate Zenodo DOI. The unit is therefore
`source_DOI`:

- `collections_index.sources` is a **list of structs**, one per source, each carrying its own
  `license`, `citation`, and `provenance`. Unnest before reaching into it.
- `analysis_results_index.provenance` is a **plain struct** on the row (new column in v25).

### Reading it

```python
client.fetch_index("collections_index")

# Who de-identified each component of a collection, and did IDC convert it?
client.sql_query("""
    SELECT collection_id,
           src.source_doi,
           src.provenance.data_contributor,
           src.provenance.deidentification_party,
           src.provenance.dicom_conversion_by
    FROM (SELECT collection_id, unnest(sources) AS src FROM collections_index)
    WHERE collection_id = 'tcga_brca'
""")
```

Attach provenance to actual series by joining `index.source_DOI` to the unnested `source_doi`
(the DOIs match exactly; lowercase both sides to be safe):

```python
client.sql_query("""
    SELECT i.collection_id, i.Modality, COUNT(*) AS series,
           s.src.provenance.data_contributor AS contributed_by,
           s.src.provenance.dicom_conversion_by AS converted_by
    FROM index i
    JOIN (SELECT collection_id, unnest(sources) AS src FROM collections_index) s
      ON s.collection_id = i.collection_id
     AND lower(s.src.source_doi) = lower(i.source_DOI)
    GROUP BY ALL
""")
```

### Two traps

**`NULL` in `collections_index` does not mean "unknown".** 92 of the 308 sources in v25 have a
`NULL` provenance struct, and all 92 are analysis-result DOIs whose provenance lives in
`analysis_results_index.provenance` instead — which is populated for all 26 rows. Across the two
tables coverage is complete, so a `NULL` here is a **pointer to the other table**, not missing
data. Look the `source_DOI` up there before reporting provenance as unavailable.

**Absent and unknown are different values.** The vocabulary uses explicit sentinels, so test for
them rather than treating any non-empty string as an answer:

| Sentinel | Means |
|----------|-------|
| `NOT_APPLICABLE` | The step did not happen — e.g. `dicom_conversion_by` for the 141 sources that arrived already as DICOM |
| `UNKNOWN` | The step happened but the party is not recorded |
| `NOT_DOCUMENTED` | Seen in `analysis_results_index.dicom_conversion_by` |

In v25, `data_contributor` is `TCIA` for 141 sources and `IDC` for 75; the 75 IDC-contributed
sources are exactly those with `dicom_conversion_by = 'IDC'`, since IDC contributing a component
implies IDC converted it from a non-DICOM upstream source.

## Join Column Reference

Use this table to identify join columns between index tables. Always call `client.fetch_index("table_name")` before using a table in SQL.

| Table A | Table B | Join Condition |
|---------|---------|----------------|
| `index` | `collections_index` | `index.collection_id = collections_index.collection_id` |
| `index` | `collections_index.sources` (provenance) | `index.collection_id = collections_index.collection_id` AND `lower(index.source_DOI) = lower(src.source_doi)` after `unnest(sources) AS src` |
| `index` | `analysis_results_index` | `lower(index.source_DOI) = lower(analysis_results_index.source_DOI)` |
| `index` | `sm_index` | `index.SeriesInstanceUID = sm_index.SeriesInstanceUID` |
| `index` | `seg_index` | `index.SeriesInstanceUID = seg_index.segmented_SeriesInstanceUID` |
| `index` | `ann_index` | `index.SeriesInstanceUID = ann_index.SeriesInstanceUID` |
| `ann_index` | `ann_group_index` | `ann_index.SeriesInstanceUID = ann_group_index.SeriesInstanceUID` |
| `index` | `clinical_index` | `index.collection_id = clinical_index.collection_id` (then filter by patient) |
| `index` | `contrast_index` | `index.SeriesInstanceUID = contrast_index.SeriesInstanceUID` |
| `index` | `volume_geometry_index` | `index.SeriesInstanceUID = volume_geometry_index.SeriesInstanceUID` |
| `index` | `rtstruct_index` | `index.SeriesInstanceUID = rtstruct_index.SeriesInstanceUID` |
| `rtstruct_index` | `index` (source images) | `rtstruct_index.referenced_SeriesInstanceUID = index.SeriesInstanceUID` |
| `index` | `ct_index` | `index.SeriesInstanceUID = ct_index.SeriesInstanceUID` |
| `index` | `mr_index` | `index.SeriesInstanceUID = mr_index.SeriesInstanceUID` |
| `index` | `pt_index` | `index.SeriesInstanceUID = pt_index.SeriesInstanceUID` |

For complete query examples using these joins, see `references/sql_patterns.md`.

## Troubleshooting

**Issue:** Column not found in table
- **Cause:** Column name misspelled or doesn't exist in that table
- **Solution:** Use `client.indices_overview["table_name"]["schema"]["columns"]` to list available columns

**Issue:** DataFrame access returns None
- **Cause:** Index not fetched or property name incorrect
- **Solution:** Fetch first with `client.fetch_index()`, then access via property matching the index name

## Resources

- Complete index table documentation: https://idc-index.readthedocs.io/en/latest/indices_reference.html
- `references/sql_patterns.md` for query examples using these tables
- `references/clinical_data_guide.md` for clinical data workflows
- `references/digital_pathology_guide.md` for pathology-specific indices
