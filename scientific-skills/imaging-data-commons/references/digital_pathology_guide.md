# Digital Pathology Guide for IDC

**Tested with:** IDC data version v23, idc-index 0.11.9

For general IDC queries and downloads, use `idc-index` (see main SKILL.md). This guide covers slide microscopy (SM) imaging, microscopy bulk simple annotations (ANN), and segmentations (SEG) in the context of digital pathology in IDC.

## Index Tables for Digital Pathology

Five specialized index tables provide curated metadata without needing BigQuery:

| Table | Row Granularity | Description |
|-------|-----------------|-------------|
| `sm_index` | 1 row = 1 SM series | Slide Microscopy series metadata: lens power, pixel spacing, image dimensions |
| `sm_instance_index` | 1 row = 1 SM instance | Instance-level (SOPInstanceUID) metadata for individual slide images |
| `seg_index` | 1 row = 1 SEG series | DICOM Segmentation metadata: algorithm, segment count, reference to source series. Used for both radiology and pathology — filter by source Modality to find pathology-specific segmentations |
| `ann_index` | 1 row = 1 ANN series | Microscopy Bulk Simple Annotations series metadata; includes `referenced_SeriesInstanceUID` linking to the annotated slide |
| `ann_group_index` | 1 row = 1 annotation group | Annotation group details: `AnnotationGroupLabel`, `GraphicType`, `NumberOfAnnotations`, `AlgorithmName`, property codes |

All require `client.fetch_index("table_name")` before querying. Use `client.indices_overview` to inspect column schemas programmatically.

## Slide Microscopy Queries

### Basic SM metadata

```python
from idc_index import IDCClient
client = IDCClient()

# sm_index has detailed metadata; join with index for collection_id
client.fetch_index("sm_index")
client.sql_query("""
    SELECT i.collection_id, COUNT(*) as slides,
           MIN(s.min_PixelSpacing_2sf) as min_resolution
    FROM sm_index s
    JOIN index i ON s.SeriesInstanceUID = i.SeriesInstanceUID
    GROUP BY i.collection_id
    ORDER BY slides DESC
""")
```

### Find SM series with specific properties

```python
# Find high-resolution slides with specific objective lens power
client.fetch_index("sm_index")
client.sql_query("""
    SELECT
        i.collection_id,
        i.PatientID,
        s.ObjectiveLensPower,
        s.min_PixelSpacing_2sf
    FROM sm_index s
    JOIN index i ON s.SeriesInstanceUID = i.SeriesInstanceUID
    WHERE s.ObjectiveLensPower >= 40
    ORDER BY s.min_PixelSpacing_2sf
    LIMIT 20
""")
```

## Annotation Queries (ANN)

DICOM Microscopy Bulk Simple Annotations (Modality = 'ANN') are annotations **on** slide microscopy images. They appear in `ann_index` (series-level) and `ann_group_index` (group-level detail). Each ANN series references the slide it annotates via `referenced_SeriesInstanceUID`.

### Basic annotation discovery

```python
# Find annotation series and their referenced images
client.fetch_index("ann_index")
client.fetch_index("ann_group_index")

client.sql_query("""
    SELECT
        a.SeriesInstanceUID as ann_series,
        a.AnnotationCoordinateType,
        a.referenced_SeriesInstanceUID as source_series
    FROM ann_index a
    LIMIT 10
""")
```

### Annotation group statistics

```python
# Get annotation group details (graphic types, counts, algorithms)
client.sql_query("""
    SELECT
        GraphicType,
        SUM(NumberOfAnnotations) as total_annotations,
        COUNT(*) as group_count
    FROM ann_group_index
    GROUP BY GraphicType
    ORDER BY total_annotations DESC
""")
```

### Find annotations with source slide context

```python
# Find annotations with their source slide microscopy context
client.sql_query("""
    SELECT
        i.collection_id,
        g.GraphicType,
        g.AnnotationPropertyType_CodeMeaning,
        g.AlgorithmName,
        g.NumberOfAnnotations
    FROM ann_group_index g
    JOIN ann_index a ON g.SeriesInstanceUID = a.SeriesInstanceUID
    JOIN index i ON a.referenced_SeriesInstanceUID = i.SeriesInstanceUID
    WHERE g.AlgorithmName IS NOT NULL
    LIMIT 10
""")
```

## Segmentations on Slide Microscopy

DICOM Segmentations (Modality = 'SEG') are used for both radiology (e.g., organ segmentations on CT) and pathology (e.g., tissue region segmentations on whole slide images). Use `seg_index.segmented_SeriesInstanceUID` to find the source series, then filter by source Modality to isolate pathology segmentations.

```python
# Find segmentations whose source is a slide microscopy image
client.fetch_index("seg_index")
client.fetch_index("sm_index")
client.sql_query("""
    SELECT
        seg.SeriesInstanceUID as seg_series,
        seg.AlgorithmName,
        seg.total_segments,
        src.collection_id,
        src.Modality as source_modality
    FROM seg_index seg
    JOIN index src ON seg.segmented_SeriesInstanceUID = src.SeriesInstanceUID
    WHERE src.Modality = 'SM'
    LIMIT 20
""")
```

## Filter by AnnotationGroupLabel

`AnnotationGroupLabel` is the most direct column for finding annotation groups by name or semantic content. Use `LIKE` with wildcards for text search.

### Simple label filtering

```python
# Find annotation groups by label (e.g., groups mentioning "blast")
client.fetch_index("ann_group_index")
client.sql_query("""
    SELECT
        g.SeriesInstanceUID,
        g.AnnotationGroupLabel,
        g.GraphicType,
        g.NumberOfAnnotations,
        g.AlgorithmName
    FROM ann_group_index g
    WHERE LOWER(g.AnnotationGroupLabel) LIKE '%blast%'
    ORDER BY g.NumberOfAnnotations DESC
""")
```

### Label filtering with collection context

```python
# Find annotation groups matching a label within a specific collection
client.fetch_index("ann_index")
client.fetch_index("ann_group_index")
client.sql_query("""
    SELECT
        i.collection_id,
        g.AnnotationGroupLabel,
        g.GraphicType,
        g.NumberOfAnnotations,
        g.AnnotationPropertyType_CodeMeaning
    FROM ann_group_index g
    JOIN ann_index a ON g.SeriesInstanceUID = a.SeriesInstanceUID
    JOIN index i ON a.SeriesInstanceUID = i.SeriesInstanceUID
    WHERE i.collection_id = 'your_collection_id'
      AND LOWER(g.AnnotationGroupLabel) LIKE '%keyword%'
    ORDER BY g.NumberOfAnnotations DESC
""")
```

## Annotations on Slide Microscopy (SM + ANN Cross-Reference)

When looking for annotations related to slide microscopy data, use both SM and ANN tables together. The `ann_index.referenced_SeriesInstanceUID` links each annotation series to its source slide.

```python
# Find slide microscopy images and their annotations in a collection
client.fetch_index("sm_index")
client.fetch_index("ann_index")
client.fetch_index("ann_group_index")
client.sql_query("""
    SELECT
        i.collection_id,
        s.ObjectiveLensPower,
        g.AnnotationGroupLabel,
        g.NumberOfAnnotations,
        g.GraphicType
    FROM ann_group_index g
    JOIN ann_index a ON g.SeriesInstanceUID = a.SeriesInstanceUID
    JOIN sm_index s ON a.referenced_SeriesInstanceUID = s.SeriesInstanceUID
    JOIN index i ON a.SeriesInstanceUID = i.SeriesInstanceUID
    WHERE i.collection_id = 'your_collection_id'
    ORDER BY g.NumberOfAnnotations DESC
""")
```

## Join Patterns

### SM join (slide microscopy details with collection context)

```python
client.fetch_index("sm_index")
result = client.sql_query("""
    SELECT i.collection_id, i.PatientID, s.ObjectiveLensPower, s.min_PixelSpacing_2sf
    FROM index i
    JOIN sm_index s ON i.SeriesInstanceUID = s.SeriesInstanceUID
    LIMIT 10
""")
```

### ANN join (annotation groups with collection context)

```python
client.fetch_index("ann_index")
client.fetch_index("ann_group_index")
result = client.sql_query("""
    SELECT
        i.collection_id,
        g.AnnotationGroupLabel,
        g.GraphicType,
        g.NumberOfAnnotations,
        a.referenced_SeriesInstanceUID as source_series
    FROM ann_group_index g
    JOIN ann_index a ON g.SeriesInstanceUID = a.SeriesInstanceUID
    JOIN index i ON a.SeriesInstanceUID = i.SeriesInstanceUID
    LIMIT 10
""")
```

## Related Tools

The following tools work with DICOM format for digital pathology workflows:

**Python Libraries:**
- [highdicom](https://github.com/ImagingDataCommons/highdicom) - High-level DICOM abstractions for Python. Create and read DICOM Segmentations (SEG), Structured Reports (SR), and parametric maps for pathology and radiology. Developed by IDC.
- [wsidicom](https://github.com/imi-bigpicture/wsidicom) - Python package for reading DICOM WSI datasets. Parses metadata into easy-to-use dataclasses for whole slide image analysis.
- [TIA-Toolbox](https://github.com/TissueImageAnalytics/tiatoolbox) - End-to-end computational pathology library with DICOM support via `DICOMWSIReader`. Provides tile extraction, feature extraction, and pretrained deep learning models.
- [EZ-WSI-DICOMweb](https://github.com/GoogleCloudPlatform/EZ-WSI-DICOMweb) - Extract image patches from DICOM whole slide images via DICOMweb. Designed for AI/ML workflows with cloud DICOM stores.

**Viewers:**
- [Slim](https://github.com/ImagingDataCommons/slim) - Web-based DICOM slide microscopy viewer and annotation tool. Supports brightfield and multiplexed immunofluorescence imaging via DICOMweb. Developed by IDC.
- [QuPath](https://qupath.github.io/) - Cross-platform open source software for whole slide image analysis. Supports DICOM WSI via Bio-Formats and OpenSlide (v0.4.0+).

**Conversion:**
- [dicom_wsi](https://github.com/Steven-N-Hart/dicom_wsi) - Python implementation for converting proprietary WSI formats to DICOM-compliant files.
