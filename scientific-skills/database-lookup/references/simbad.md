# SIMBAD Astronomical Database (CDS Strasbourg)

## Base URLs

**TAP endpoint (recommended):**
```
https://simbad.cds.unistra.fr/simbad/sim-tap/sync
```

**Legacy script interface:**
```
https://simbad.cds.unistra.fr/simbad/sim-script
```

**Simple query endpoints:**
```
https://simbad.cds.unistra.fr/simbad/sim-id
https://simbad.cds.unistra.fr/simbad/sim-coo
```

## Authentication

No API key required. All endpoints are public.

## Key Endpoints

### 1. TAP Queries (ADQL — recommended for programmatic use)

```
GET /simbad/sim-tap/sync?request=doQuery&lang=adql&format={format}&query={ADQL}
```

| Parameter | Type   | Description |
|-----------|--------|-------------|
| `request` | string | `doQuery` |
| `lang`    | string | `adql` |
| `format`  | string | `json`, `votable`, `csv`, `tsv`. |
| `query`   | string | **Required.** ADQL query. |

**Example — look up object by name:**
```
https://simbad.cds.unistra.fr/simbad/sim-tap/sync?request=doQuery&lang=adql&format=json&query=SELECT basic.OID, ra, dec, main_id, otype FROM basic JOIN ident ON oid = ident.oidref WHERE id = 'M31'
```

**Example — cone search (objects within 5 arcmin of coordinates):**
```
https://simbad.cds.unistra.fr/simbad/sim-tap/sync?request=doQuery&lang=adql&format=json&query=SELECT TOP 50 main_id, ra, dec, otype FROM basic WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', 10.684, 41.269, 0.083)) = 1
```
Note: radius in CIRCLE is in degrees (5 arcmin = 0.083 deg).

**Example — objects by type (e.g., all pulsars):**
```
https://simbad.cds.unistra.fr/simbad/sim-tap/sync?request=doQuery&lang=adql&format=json&query=SELECT TOP 100 main_id, ra, dec, otype FROM basic WHERE otype = 'Pulsar'
```

### 2. Identifier Query (simple lookup)

```
GET /simbad/sim-id?Ident={name}&output.format=votable
```

| Parameter        | Type   | Description |
|------------------|--------|-------------|
| `Ident`          | string | **Required.** Object name (e.g., `M31`, `Sirius`, `NGC 1275`). |
| `output.format`  | string | `votable`, `html`. |

**Example:**
```
https://simbad.cds.unistra.fr/simbad/sim-id?Ident=M31&output.format=votable
```

### 3. Coordinate Query

```
GET /simbad/sim-coo?Coord={coords}&Radius={radius}&Radius.unit={unit}&output.format=votable
```

| Parameter      | Type   | Description |
|----------------|--------|-------------|
| `Coord`        | string | **Required.** Coordinates, e.g., `10.684 +41.269` or `00 42 44 +41 16 09`. |
| `Radius`       | float  | Search radius. Default: 2. |
| `Radius.unit`  | string | `arcmin`, `arcsec`, `deg`. Default: `arcmin`. |
| `output.format`| string | `votable`, `html`. |

**Example:**
```
https://simbad.cds.unistra.fr/simbad/sim-coo?Coord=10.684+%2B41.269&Radius=5&Radius.unit=arcmin&output.format=votable
```

### 4. Script Interface (for multi-command queries)

```
POST /simbad/sim-script
Content-Type: application/x-www-form-urlencoded
script=format+object+"%MAIN_ID+|+%RA+|+%DEC+|+%OTYPE"\nquery+id+M31
```

## Key TAP Tables

| Table           | Description |
|-----------------|-------------|
| `basic`         | Core data: coordinates, main_id, object type. |
| `ident`         | All known identifiers for objects. |
| `flux`          | Flux/magnitude measurements. |
| `mesVelocities` | Radial velocity measurements. |
| `mesDistance`    | Distance measurements. |
| `otypedef`      | Object type definitions/labels. |
| `allfluxes`     | All flux data joined. |

## Common Object Types (otype)

`Star`, `Galaxy`, `Pulsar`, `QSO`, `Nebula`, `GlobCluster`, `RadioSource`, `X-raySource`, `SNRemnant`

## Response Format (TAP JSON)

```json
{
  "metadata": [
    {"name": "main_id", "datatype": "char"},
    {"name": "ra", "datatype": "double"},
    {"name": "dec", "datatype": "double"}
  ],
  "data": [
    ["M 31", 10.6847, 41.2687]
  ]
}
```

## Rate Limits

No formal rate limits documented. SIMBAD requests that automated scripts include reasonable delays between queries. Very large TAP queries may time out; use `TOP N` to limit results or switch to async TAP at `/simbad/sim-tap/async`.
