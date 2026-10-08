# IDC REST API Guide

**Tested with:** API `3.0.0b5` (build `a7abd9e`), IDC data version v25, `idc_index_data_version` 25.0.0

IDC operates a hosted REST API that exposes discovery, cohort building, metadata SQL, and
download manifests over plain HTTP. No authentication, account, or credentials are required —
every example on this page can be run from any terminal with `curl`.

The API and the [MCP server](mcp_guide.md) are the same service behind two transports: the MCP
tools wrap these endpoints. Both are hosted against a current IDC release, independently of the
`idc-index` version installed locally.

## When to Use This Guide

Load this guide when you need to:
- Query IDC from a language or environment with no `idc-index` install (shell, R, Java, JS, a
  notebook without pip access)
- Build cohorts and manifests over HTTP for a pipeline or web application
- Run metadata SQL against a current IDC release without downloading local index files
- Hand a user copy-pasteable `curl` commands they can run anywhere

For downloads, pandas/notebook analysis, reading pixel data, DICOMweb, BigQuery, or digital
pathology tiling, use `idc-index` as documented in `SKILL.md`. The API never moves image bytes:
it returns public `s3://` URLs and manifests, and the transfer happens directly from cloud
storage to the client.

**Choosing among the three interfaces:**

| Situation | Use |
|-----------|-----|
| Session already has the hosted MCP server | MCP tools (`mcp_guide.md`) — same data, no HTTP plumbing |
| `idc-index` already installed, results feed pandas / downloads | `idc-index` (`SKILL.md`) |
| `idc-index` **not** installed and the task is read-only metadata | This REST API — do not install to answer a metadata question |
| No Python, or a non-Python client, or shell commands the user re-runs | This REST API |
| Installed `idc-index` is a whole IDC data release behind and cannot be upgraded here | REST API for the query, **direct bucket transfer** may be needed for the download — `idc-index` cannot fetch what its index does not list |

Being hosted, the API always serves a current IDC release — but so does an up-to-date
`idc-index`. "I want the newest data" is not on its own a reason to prefer the API: the normal
fix for a stale local index is to upgrade it. Reach for the API on version grounds only when
upgrading is not an option in that environment. Avoiding the install *is* a reason, though: the
packaged index data is ~77 MB before pandas, pyarrow, and duckdb, which a metadata question does
not need.

## Endpoint and Versioning

| Property | Value |
|----------|-------|
| Base URL | `https://api.imaging.datacommons.cancer.gov/v3` |
| Authentication | None |
| Content type | `application/json` (request and response), except `cohort/manifest.txt` → `text/plain` |
| Interactive docs | https://api.imaging.datacommons.cancer.gov/v3/docs (Swagger UI) |
| OpenAPI spec | https://api.imaging.datacommons.cancer.gov/v3/openapi.json |
| Source | https://github.com/ImagingDataCommons/IDC-REST-MCP |

**v3 is in beta.** The contract may still change before the final `3.0.0` release, so verify the
running build rather than assuming the values in this guide:

```bash
curl -s https://api.imaging.datacommons.cancer.gov/v3/version
# {"idc_version":"v25","idc_index_data_version":"25.0.0","api_version":"3.0.0b5","build":"a7abd9e"}
```

`idc_version` is the IDC data release the API serves and is the authority when the API is in
use — prefer it over the `idc-data-version` pinned in the `SKILL.md` frontmatter, which records
the release the skill was last verified against.

### Checking the API against a local idc-index

The API and `idc-index` are built on the same `idc-index-data` package, and both report its
version, so consistency is an exact check rather than a guess:

| Side | Data release (coarse) | `idc-index-data` version (exact) |
|------|-----------------------|----------------------------------|
| API | `GET /v3/version` → `idc_version` | `GET /v3/version` → `idc_index_data_version` |
| Local | `IDCClient().get_idc_version()` | `idc_index_data.__version__` |

```python
import idc_index_data
import requests

api = requests.get("https://api.imaging.datacommons.cancer.gov/v3/version", timeout=30).json()
api_version, local_version = api["idc_index_data_version"], idc_index_data.__version__

if api_version.split(".")[0] != local_version.split(".")[0]:
    print(f"Different IDC data release: API {api_version}, local {local_version} — upgrade")
elif api_version != local_version:
    print(f"Same data release, different index build: API {api_version}, local {local_version}")
    # → Same data release, different index build: API 24.2.2, local 24.2.0
```

**The major is the IDC data release; the rest is the index build.** `idc-index-data` `24.x.y`
serves IDC `v24` — 24.0.0 shipped with the v24 release, and 24.1.0 / 24.2.x are later builds of
the *same* release. What differs between them is the index itself (added tables and columns,
corrected metadata), never which series IDC contains.

So read a mismatch by its position:

| Difference | Means | Consequence |
|------------|-------|-------------|
| Major (24.x.y vs 25.x.y) | Different IDC data release | Series added, revised, or removed. Counts legitimately differ, and `idc-index` **cannot download** what its index does not list |
| Minor or patch (24.2.0 vs 24.2.2) | Same data release, different index build | Same series everywhere; downloads are unaffected. A metadata query can still differ if it touches a column that was added or corrected |

Comparing `idc_version` alone cannot make this distinction in the other direction either — the
`vNN` label is exactly the major, so matching `v25` on both sides tells you the release agrees
but says nothing about the index build.

When the two disagree, say so and name both versions, then reconcile rather than mixing
results: upgrading `idc-index` brings the local side to the newer `idc-index-data`, and
`scripts/check_version.py` reports whether an upgrade is available and prints the
command for the interpreter you are running. Do not present API-derived and locally-derived counts side by side as if they came
from one index.

**A major behind also breaks downloads.** `idc-index` resolves every `s3://` URL it is given
against its *own* index, so a manifest built from a newer IDC data release can name series it
has never heard of — see *When the local index is a data release behind the API* under
**Getting the Data**.

### Use v3 only — V1 and V2 are being retired

**Do not write new code against the V1 or V2 IDC APIs, and do not follow V1/V2 examples found
in older tutorials, notebooks, blog posts, or forum answers.** Both are superseded by v3 and
are scheduled to be deprecated and shut down; code written against them will stop working. If
a user brings V1 or V2 code, say so and port it to v3 rather than extending it.

Two signals that a snippet is V1/V2 rather than v3, both of which will fail against `/v3`:

- a base URL other than `https://api.imaging.datacommons.cancer.gov/v3` (for example a
  `/v1/` or `/v2/` path segment)
- per-attribute filter suffixes such as `Modality_btw` or `_lt` / `_gt`, instead of v3's
  separate `terms` and `ranges` objects

V1/V2 documentation survives only in the [IDC docs archive](https://learn.canceridc.dev/archive/archive)
for historical reference. Treat it as read-only history, not as a source of working examples.

## The Query Surfaces

The endpoints group into five surfaces that build on each other:

| Surface | Answers | Endpoints |
|---------|---------|-----------|
| Discovery | What exists? What can I filter on? | `GET /version`, `/stats`, `/collections`, `/collections/{id}`, `/analysis_results`, `/attributes`, `/attributes/{attr}/values` |
| Cohort | How big is my selection, and what's in it? | `POST /cohort/counts`, `POST /cohort/manifest` |
| Retrieval | Give me the download links | `POST /cohort/manifest.txt` |
| SQL | Anything relational or aggregate | `GET /tables`, `GET /tables/{table}`, `POST /sql` |
| Side tools | View / cite / license-check a cohort | `GET /viewer-url`, `POST /citations`, `POST /licenses` |

Discovery supplies the vocabulary (attribute names and their valid values) that the cohort
filters consume. SQL is the escape hatch for questions structured filters cannot express.

Clinical data has its own discovery endpoints (`GET /clinical/tables`,
`/clinical/tables/{table}`, `/clinical/tables/{table}/rows`) and is filtered or joined through
SQL against the `clinical` schema.

## Endpoint Reference

All paths are relative to `https://api.imaging.datacommons.cancer.gov/v3`.

| Method & path | Purpose | Key response fields |
|---------------|---------|---------------------|
| `GET /version` | IDC data release + API build | `idc_version`, `idc_index_data_version`, `api_version`, `build` |
| `GET /stats` | Headline totals | `collections`, `patients`, `studies`, `series`, `instances`, `size_TB` |
| `GET /collections` | List collections (JSON array) | `collection_id`, `collection_name`, `cancer_types`, `tumor_locations`, `species`, `subjects`, `description` |
| `GET /collections/{id}` | Collection detail | the above plus `patients`, `studies`, `series`, `instances`, `size_TB`, `modalities[]`, `licenses[]` |
| `GET /analysis_results` | Derived datasets (JSON array) | `analysis_result_id`, `analysis_result_title`, `source_DOI`, `subjects`, `collections`, `modalities`, `license_short_name` |
| `GET /attributes` | Filterable attributes | `name`, `table`, `data_type`, `kind` (`term` \| `range`), `categorical`, `description` |
| `GET /attributes/{attr}/values?limit=` | Distinct values with counts | `attribute`, `values[{value,count}]`, `truncated`, `note` |
| `GET /tables` | Tables available to SQL | `tables[{name,description,column_count}]` |
| `GET /tables/{table}` | Column schema | `name`, `description`, `columns[{name,type,description}]` |
| `GET /clinical/tables?collection_id=` | Clinical tables, optionally one collection | `tables[{table_name,sql_path,collection_id,column_count}]` |
| `GET /clinical/tables/{table}` | Clinical columns + labels | `name`, `columns[{name,type,description}]` |
| `GET /clinical/tables/{table}/rows?max_rows=` | Clinical rows (capped) | `columns`, `rows`, `row_count`, `truncated`, `max_rows` |
| `POST /cohort/counts` | Distinct counts for a filter (cheap) | `patients`, `studies`, `series`, `instances`, `size_TB` |
| `POST /cohort/manifest` | Counts + page of series + download payload | `counts`, `page`, `page_size`, `returned`, `total_series`, `series[]`, `download` |
| `POST /cohort/manifest.txt` | Full manifest as `text/plain` | one `s3://…/*` URL per line |
| `POST /sql` | Guarded read-only SQL (DuckDB) | `columns`, `rows`, `row_count`, `truncated`, `max_rows` |
| `GET /viewer-url` | OHIF / Slim viewer link | `viewer_url`, `viewer`, `study_instance_uid`, `series_instance_uid` |
| `POST /citations` | Citations for a cohort | `format`, `citations[]`, `idc_acknowledgment`, `recommendation` |
| `POST /licenses` | License breakdown for a cohort | `licenses[{license_short_name,series,size_TB}]` |

`GET /health` and `GET /v3` (API root) also exist for liveness checks.

## Filter Syntax

Cohort filters are shared by `cohort/counts`, `cohort/manifest`, `cohort/manifest.txt`,
`licenses`, and `citations`. **The filter object always goes under a `filters` key**, on every one
of them. It has two parts:

- **`terms`** — `{attribute: [values]}` for equality/membership. Values are **OR**'d within an
  attribute and **AND**'d across attributes.
- **`ranges`** — `{attribute: {"gte": x, "lte": y}}` for numeric and date attributes. Either
  bound may be omitted for an open-ended range.

```json
{
  "filters": {
    "terms": {"Modality": ["CT"], "collection_id": ["nlst"]},
    "ranges": {"instanceCount": {"gte": 100, "lte": 200}}
  }
}
```

Filters operate only on the `index` table's filterable attributes — 19 of them as of `3.0.0b5`:

| Kind | Attributes |
|------|------------|
| `term` | `collection_id`, `analysis_result_id`, `PatientID`, `StudyInstanceUID`, `SeriesInstanceUID`, `Modality`, `BodyPartExamined`, `Manufacturer`, `ManufacturerModelName`, `PatientSex`, `sop_class_name`, `license_short_name`, `source_DOI` |
| `range` | `instanceCount`, `series_size_MB`, `series_init_idc_version`, `series_revised_idc_version`, `StudyDate`, `SeriesDate` |

`SeriesInstanceUID`, `StudyInstanceUID`, and `PatientID` being filterable is what makes the
side tools work at any granularity — the licenses or citations for a single series are just a
one-value filter.

Anything outside this list — clinical values, segmented anatomy, per-modality acquisition
parameters — is not filterable here; use the SQL surface instead.

**Ground values before filtering.** Call `GET /attributes` for what is filterable and whether
it is a term or a range, then `GET /attributes/{attr}/values` for real values and their casing —
values are matched case-sensitively.

### The server reports what it filtered on

Every filtered response echoes `filters_applied` and `warnings`, and misuse is refused rather
than ignored. Together these make a result self-describing, so you do not have to sanity-check a
count against a number you happen to remember.

| Mistake | Response |
|---------|----------|
| Bare filter object, no `filters` key | `422` naming the correct shape |
| Unrecognized key at any depth (`term` for `terms`, a range bound spelled `min`) | `422` pointing at the key |
| Unknown filter attribute, or a range attribute used as a term | `400` naming the discovery call |
| Unfiltered `cohort/manifest` or `manifest.txt` | `400` — it will not enumerate the whole archive |
| Unfiltered `cohort/counts` or `licenses` | `200` plus an explicit "ENTIRE IDC archive" warning |
| A predicate that constrains nothing (`{"collection_id": []}`, `{"instanceCount": {}}`) | `200`, and `warnings` names the ignored predicate |
| Miscased value (`mr` for `MR`) | `200`, zero counts, and a warning naming the casing that exists |

```bash
curl -s $B/cohort/counts -H 'content-type: application/json' \
  -d '{"filters": {"terms": {"collection_id": ["rider_pilot"]}}}'
# {"patients":8,"studies":154,"series":774,"instances":21111,"size_TB":0.011,
#  "filters_applied":{"terms":{"collection_id":["rider_pilot"]},"ranges":{}},"warnings":[]}
```

**Read `warnings` before reporting any count.** A zero count with `warnings: []` means the filter
applied and matched nothing — that is a real answer. A zero count *with* a casing warning means
the filter was wrong. And `filters_applied` covers the case no shape check can: a request
carrying one good predicate plus one that constrains nothing returns a perfectly plausible
number, and only `warnings` reveals the dropped half.

`cohort/manifest.txt` returns `text/plain` and so cannot carry these fields; there the
required-predicate `400` does the same job, appending the ignored-predicate reason to its message.

## Worked Examples

### Discovery

```bash
B=https://api.imaging.datacommons.cancer.gov/v3

curl -s $B/version                      # data release + API build
curl -s $B/stats                        # headline totals
curl -s $B/collections                  # all 179 collections
curl -s $B/collections/rider_pilot      # one collection: counts, modalities, licenses
curl -s $B/analysis_results             # derived datasets (segmentations, annotations)
curl -s $B/attributes                   # what can be filtered, and how
curl -s "$B/attributes/Modality/values?limit=10"
```

### Cohort building

Check size first — `counts` is cheap and answers "is this download sane?":

```bash
curl -s $B/cohort/counts \
  -H 'content-type: application/json' \
  -d '{"filters": {"terms": {"Modality": ["MR"], "BodyPartExamined": ["BREAST"]}}}'
# {"patients":3718,"studies":5689,"series":47986,"instances":6493262,"size_TB":2.421,
#  "filters_applied":{...},"warnings":[]}
```

Then request a page of series plus the download payload — same filter, same `filters` key.

```bash
curl -s $B/cohort/manifest \
  -H 'content-type: application/json' \
  -d '{"filters": {"terms": {"Modality": ["MR"], "BodyPartExamined": ["BREAST"]}},
       "page": 0, "page_size": 3}'
```

Each `series[]` row carries `collection_id`, `PatientID`, `StudyInstanceUID`,
`SeriesInstanceUID`, `Modality`, `SeriesDescription`, `instanceCount`, `series_size_MB`,
`aws_bucket`, `crdc_series_uuid`, and `series_aws_url`. Set `"include_rows": false` to get
counts and the download payload without the rows.

### SQL

```bash
curl -s $B/sql \
  -H 'content-type: application/json' \
  -d '{"sql": "SELECT Modality, count(*) n FROM index GROUP BY 1 ORDER BY n DESC", "max_rows": 20}'
```

### Viewer, licenses, citations

```bash
curl -s "$B/viewer-url?study_instance_uid=1.3.6.1.4.1.14519.5.2.1.7695.4164.129908397467389975396031099306"

curl -s $B/licenses \
  -H 'content-type: application/json' \
  -d '{"filters": {"terms": {"collection_id": ["rider_pilot"]}}}'

curl -s $B/citations \
  -H 'content-type: application/json' \
  -d '{"filters": {"terms": {"collection_id": ["rider_pilot"]}}, "citation_format": "bibtex"}'
```

`viewer-url` takes `series_instance_uid` or `study_instance_uid` (and an optional `viewer`
override); it picks OHIF v3 for radiology and Slim for slide microscopy automatically.

`citation_format` is one of `apa` (default), `bibtex`, `csl-json`, `turtle`. The response
separates the per-dataset `citations[]` from `idc_acknowledgment`, the IDC paper — include
both when publishing.

### Python client

```python
import requests

BASE = "https://api.imaging.datacommons.cancer.gov/v3"
session = requests.Session()

def get(path, **params):
    r = session.get(f"{BASE}{path}", params=params, timeout=60)
    r.raise_for_status()
    return r.json()

def post(path, payload):
    r = session.post(f"{BASE}{path}", json=payload, timeout=120)
    r.raise_for_status()
    return r.json()

# 1. Confirm which IDC release the API is serving
print(get("/version")["idc_version"])

# 2. Ground the filter values before using them
modalities = {v["value"] for v in get("/attributes/Modality/values", limit=1000)["values"]}
assert "MR" in modalities

# 3. Size the cohort, then page through it — the filter always goes under `filters`
filters = {"terms": {"collection_id": ["rider_pilot"], "Modality": ["CT"]}}
counts = post("/cohort/counts", {"filters": filters})
assert not counts["warnings"], counts["warnings"]   # nothing was silently dropped
print(counts)

manifest = post("/cohort/manifest", {"filters": filters, "page": 0, "page_size": 100})
uids = [row["SeriesInstanceUID"] for row in manifest["series"]]

# 4. Hand off to idc-index for the download (see "Handing Off to idc-index" below)
```

Error responses raise through `raise_for_status()`; read `r.json()["error"]["message"]` for the
reason before retrying.

## The SQL Surface

`POST /sql` runs read-only DuckDB SQL against the same tables `idc-index` exposes locally.
Ground the schema with `GET /tables` and `GET /tables/{table}` — do not guess table or column
names.

**Guardrails** (verified against `3.0.0b5`):

- Only single read-only `SELECT` / `WITH … SELECT` statements are accepted. Anything else is
  rejected with `{"error": {"code": "invalid_query", "message": "Only read-only SELECT (or WITH ... SELECT) statements are allowed."}}` and HTTP 400.
- A server row cap and per-query timeout apply. `max_rows` defaults to **5000** and is clamped
  to **10000**; results carry `truncated: true` when the cap was hit, and echo the `max_rows`
  actually applied.
- Invalid SQL returns DuckDB's own error text, including its "Candidate bindings" suggestions —
  useful for fixing a misspelled column without another schema round-trip.

**Tables reachable from SQL** are the `index` table plus the specialized indices documented in
`SKILL.md` (`collections_index`, `analysis_results_index`, `version_metadata_index`,
`prior_versions_index`, `seg_index`, `ann_index`, `ann_group_index`, `rtstruct_index`,
`ct_index`, `mr_index`, `pt_index`, `sm_index`, `sm_instance_index`, `contrast_index`,
`volume_geometry_index`, `clinical_index`). Join them to `index` on `SeriesInstanceUID`, except
where `SKILL.md` documents a different key (`segmented_SeriesInstanceUID`,
`referenced_SeriesInstanceUID`, `collection_id`, `analysis_result_id`).

The SQL patterns in `references/sql_patterns.md` are written for `client.sql_query()` but the
SQL itself transfers unchanged — send it as the `sql` field.

**Array columns:** columns typed `STRING[]` (e.g. `SegmentedPropertyType_CodeMeanings`) hold a
list per row. Match with `list_contains(col, 'value')`, not `=` or `LIKE`:

```bash
curl -s $B/sql -H 'content-type: application/json' -d '{
  "sql": "SELECT i.collection_id, count(DISTINCT i.SeriesInstanceUID) AS slides FROM index i JOIN seg_index seg ON seg.segmented_SeriesInstanceUID = i.SeriesInstanceUID WHERE i.Modality = '\''SM'\'' AND list_contains(seg.SegmentedPropertyType_CodeMeanings, '\''Nucleus'\'') GROUP BY 1 ORDER BY slides DESC",
  "max_rows": 20}'
```

**A SQL result can be a manifest.** Selecting `series_aws_url` gives download URLs directly,
and unlike `cohort/manifest.txt` you can carry extra columns alongside them — a per-row
`license_short_name`, for instance, which the plain manifest does not include:

```sql
SELECT SeriesInstanceUID, license_short_name, series_aws_url
FROM index WHERE collection_id = 'rider_pilot'
```

For bulk series, still prefer `cohort/manifest.txt` — it is not subject to the SQL row cap.

### Clinical data

Clinical data comes in two layers:

- **`clinical_index`** — a data dictionary, one row per (collection, table, column) with a
  human-readable `column_label` and coded `values`. It is an ordinary SQL table and joins to
  `index` on `collection_id`.
- **Per-collection clinical tables** (e.g. `nlst_canc`) — the actual rows. There are ~150, so
  they are kept out of `GET /tables` and discovered through the `/clinical/tables` endpoints
  instead. In SQL they live under a separate schema and are addressed as `clinical.<table>`.

Clinical tables join to imaging on **`dicom_patient_id = index.PatientID`**, not on a series
UID. Clinical data is not harmonized across collections — table and column names vary, so
always discover before querying.

```bash
curl -s "$B/clinical/tables?collection_id=nlst"            # which tables this collection has
curl -s $B/clinical/tables/nlst_canc                       # columns + labels
curl -s "$B/clinical/tables/nlst_canc/rows?max_rows=100"   # rows, capped

curl -s $B/sql -H 'content-type: application/json' -d '{
  "sql": "SELECT count(DISTINCT i.PatientID) AS patients FROM index i JOIN clinical.nlst_canc c ON c.dicom_patient_id = i.PatientID WHERE i.collection_id = '\''nlst'\'' AND i.Modality = '\''CT'\'' AND c.clinical_stag = '\''400'\''"}'
```

See `references/clinical_data_guide.md` for value mapping and the wider clinical data model.

## Getting the Data

The API returns manifests; a client moves the bytes. Every URL points at public AWS S3 or GCS
buckets and needs no credentials.

```bash
# save the full manifest
curl -s $B/cohort/manifest.txt \
  -H 'content-type: application/json' \
  -d '{"filters": {"terms": {"collection_id": ["rider_pilot"]}}}' > idc_manifest.txt

# download it (needs idc-index installed)
idc download-from-manifest idc_manifest.txt --download-dir ./idc-data
```

For a filter that is a single `collection_id`, the `download` payload of `cohort/manifest`
emits the simpler `idc download <collection_id> --download-dir ./idc-data` form.

### When the local index is a data release behind the API

This applies when the two `idc-index-data` **majors** differ — the API serving `25.x.y` against
a local `24.x.y`, say. A newer build of the same release (24.2.2 vs 24.2.0) covers the same
series, so manifests from it resolve locally and downloads are unaffected.

`idc download-from-manifest` does not simply hand the URLs to a transfer client: it extracts
each `crdc_series_uuid` from the manifest and joins it against the **local** index (then
against `prior_versions_index`) to compute sizes and build the output hierarchy. A manifest
produced by an API serving a newer IDC data release can therefore contain series the local
index has never heard of.

Those rows are **not** downloaded, and the command does not fail — it logs, then continues
with the rest:

```
The total of N copy commands are not recognized as referencing any associated series in the
main index. ... they may correspond to files available in a release of IDC different from v25
used in this version of idc-index.
...
The corresponding files could not be downloaded.
```

The result is a partial download that otherwise looks successful. `download_from_selection(seriesInstanceUID=…)`
has the same blind spot from the other direction: it filters the local index, so UIDs it does
not contain are silently dropped from the selection.

**Fix it one of two ways:**

1. **Upgrade** — upgrade `idc-index` (`scripts/check_version.py` prints the command),
   then re-run. This is the right answer
   whenever it is possible; it restores the hierarchy, size checks, and progress reporting.
2. **Bypass the index** — transfer directly from the bucket. The manifest URLs are
   self-contained (`s3://<bucket>/<crdc_series_uuid>/*`), so no index is needed at all:

```bash
# one directory per series, named by crdc_series_uuid
awk -F/ '{print "cp " $0 " ./idc-data/" $4 "/"}' idc_manifest.txt > s5cmd_commands.txt
s5cmd --no-sign-request run s5cmd_commands.txt

# for source=gcs manifests, add the GCS endpoint
s5cmd --no-sign-request --endpoint-url https://storage.googleapis.com run s5cmd_commands.txt
```

`s5cmd run` expects one *command* per line, which is why the bare URLs are rewritten as `cp`
commands; drop the `$4` segment to land every file flat in `./idc-data/`. `aws s3 cp
--no-sign-request --recursive` works the same way per URL.

What you give up by bypassing `idc-index` is convenience, not data: files are laid out under
CRDC UUIDs instead of `%collection_id/%PatientID/%Modality`, and there is no local size or
disk-space check — so read `size_TB` from `cohort/counts` first. See
`references/cloud_storage_guide.md` for bucket layout, `aws`/`gsutil` equivalents, and
UUID-to-DICOM-UID mapping.

**`source`: `aws` (default) vs `gcs`.** Both return `s3://` URLs — GCS is reached through its
S3-compatible endpoint, never a `gs://` URL. That is why `idc download-from-manifest` only
recognizes `s3://` lines. Driving `s5cmd` yourself, use `--no-sign-request`, and for
`source=gcs` add `--endpoint-url https://storage.googleapis.com`.

IDC is ~100 TB across 179 collections. Always report `series` and `size_TB` from
`cohort/counts` and confirm with the user before starting a broad download.

## Limits, Defaults, and Errors

Measured against `3.0.0b5`. Values above a cap are silently clamped — the response echoes the
value actually used (`max_rows`, `page_size`), so read it back rather than assuming the request
was honored.

| Endpoint | Parameter | Default | Cap |
|----------|-----------|---------|-----|
| `GET /attributes/{attr}/values` | `limit` | 100 | 10000 |
| `POST /sql` | `max_rows` | 5000 | 10000 |
| `GET /clinical/tables/{table}/rows` | `max_rows` | 5000 | 100000 |
| `POST /cohort/manifest` | `page_size` | 100 | 5000 |
| `POST /cohort/manifest` | `page` | 0 | — |
| `POST /cohort/manifest` | `include_rows` | `true` | — |
| `POST /cohort/manifest.txt` | `limit` | 100000 | — |
| `POST /cohort/manifest.txt` | `source` | `aws` | — |
| `POST /citations` | `citation_format` | `apa` | — |

`cohort/manifest.txt` is the surface that is not *row*-capped the way `/sql` is — it returned all
774 lines for `rider_pilot` with no `limit` set, and enumerates up to 100 000 series. That is why
bulk series belong there rather than in a `/sql` dump.

There is **no per-caller rate limit or quota** and no `429`. What is bounded is the individual
request: a 30 s SQL statement timeout, 4 GB query memory, and the caps above. A burst is absorbed
by autoscaling and surfaces as slower responses or a `503` — back off and retry rather than
treating it as permanent. For sustained heavy metadata access, query the `idc-index` Parquet files
(`references/parquet_access_guide.md`) or BigQuery instead of driving this API hard.

Size-capped responses carry a `truncated` boolean: `false` means the result is complete, `true`
means raise the limit or aggregate/narrow instead. Explore narrow, then widen.

Errors come in two shapes. Semantic problems are HTTP 400 with a uniform body:

```json
{"error": {"code": "invalid_query", "message": "Unknown or non-term filter attribute: 'NotAnAttribute'. Use list_attributes to see valid attributes."}}
```

Request-shape problems are HTTP 422 with FastAPI's `detail[]` array, which names the offending
key — a bare filter object, an unrecognized key, or a misspelled range bound all land here. Both
kinds are actionable: an unknown attribute names the discovery call to make, and a bad column name
carries DuckDB's candidate bindings. Read the message and fix the request rather than retrying it
unchanged.

What does **not** produce an error is a filter that is valid but empty or over-broad. Those return
HTTP 200 with a `warnings` entry saying so — see *The server reports what it filtered on*. Read
`warnings`; do not infer from the count alone.

## Handing Off to idc-index

The boundary artifact is a list of `SeriesInstanceUID` values (or a saved `manifest.txt`).

```python
# UIDs from a cohort/manifest or /sql response
series_uids = [row["SeriesInstanceUID"] for row in manifest["series"]]

from idc_index import IDCClient
client = IDCClient()

client.download_from_selection(
    downloadDir="./data",
    seriesInstanceUID=series_uids,        # a list, not a DataFrame
    dirTemplate="%collection_id/%PatientID/%Modality",
)
```

Run `scripts/check_version.py` before the first `idc-index` call in a session, even when
discovery happened over the API — the two version independently. Compare
`idc_index_data_version` on both sides first (see *Checking the API against a local
idc-index*).

**This handoff is only valid while the two are on the same IDC data release** — the same
`idc-index-data` major. `idc-index` can only download series its own index lists, so when the
API is a release ahead, UIDs and manifest URLs it returned may resolve to nothing locally:
`download_from_selection` silently drops them and `download-from-manifest` logs them as
unrecognized and skips them. Do not report that as "no data": name both versions, then either
upgrade `idc-index` or download straight from the bucket, as described in *When the local index
is a data release behind the API*.

## What the API Does Not Cover

- **Image bytes.** The API returns URLs and manifests only; files transfer from S3/GCS.
- **Pixel data access and DICOMweb.** Use `references/dicomweb_guide.md`.
- **Full DICOM metadata, per-segment detail, SR quantitative/qualitative measurements, private
  DICOM elements.** Still BigQuery-only — see `references/bigquery_guide.md`.
- **Writes.** The service is read-only by construction; there are no POST endpoints that mutate
  state, and the SQL connection rejects anything but `SELECT`.
- **Local analysis.** DataFrames, plotting, pydicom/SimpleITK, pathology tiling all stay with
  `idc-index`.

## Related Documentation

- IDC REST API docs: https://learn.canceridc.dev/rest-api/api
- Swagger UI: https://api.imaging.datacommons.cancer.gov/v3/docs
- API and MCP server source: https://github.com/ImagingDataCommons/IDC-REST-MCP
- `references/mcp_guide.md` — the same capabilities as agent tools
- `references/sql_patterns.md` — SQL that transfers unchanged to `POST /sql`
- `references/cli_guide.md` — `idc download-from-manifest` and the rest of the CLI
