---
name: xarray
description: Labeled multi-dimensional arrays for scientific data. Use for NetCDF/HDF5/Zarr I/O, climate and weather data analysis, labeled dimension operations (groupby, resample, rolling), broadcasting with named dimensions, Dask-backed out-of-core computation, and time series analysis on gridded datasets.
license: Apache-2.0 license
metadata:
    skill-author: Adem Usta
---

# xarray: Labeled Multi-Dimensional Data

## Overview

xarray extends NumPy with labeled dimensions, coordinates, and attributes, making it natural to work with multi-dimensional scientific datasets. It is the standard tool for climate science, weather forecasting, oceanography, remote sensing, and any domain that works with gridded data stored in NetCDF, HDF5, or Zarr formats.

## Installation

```bash
# Install xarray using uv
uv pip install xarray

# With NetCDF support (recommended)
uv pip install "xarray[io]"
# or individually:
uv pip install netCDF4 h5netcdf scipy

# With Dask for parallel/out-of-core computation
uv pip install "xarray[parallel]"
# or:
uv pip install dask[array]

# Full installation
uv pip install xarray netCDF4 h5netcdf dask[array] zarr bottleneck flox
```

## When to Use This Skill

Use the xarray skill when:

- Working with labeled multi-dimensional arrays (time × lat × lon × level)
- Reading/writing NetCDF, HDF5, Zarr, or GRIB files
- Analyzing climate, weather, or Earth observation data
- Performing groupby, resampling, or rolling operations on gridded data
- Needing dimension-aware broadcasting and alignment
- Processing datasets too large for memory (with Dask)
- Working with time series on regular or irregular grids
- Computing spatial/temporal statistics over named dimensions

Consider alternatives when:
- You need general-purpose dataframes → Use **pandas** or **polars**
- You need geospatial vector data → Use **geopandas**
- You need image processing → Use OpenCV or scikit-image
- You need in-memory tabular data only → Use **pandas**

## Quick Start

### Creating Data

```python
import xarray as xr
import numpy as np
import pandas as pd

# DataArray: a single labeled array
temperature = xr.DataArray(
    data=np.random.randn(4, 3),
    dims=['time', 'location'],
    coords={
        'time': pd.date_range('2024-01-01', periods=4),
        'location': ['NYC', 'LA', 'Chicago'],
    },
    attrs={'units': 'degC', 'long_name': 'Temperature'}
)

# Dataset: a collection of DataArrays sharing dimensions
ds = xr.Dataset({
    'temperature': (['time', 'lat', 'lon'], np.random.randn(365, 180, 360)),
    'precipitation': (['time', 'lat', 'lon'], np.abs(np.random.randn(365, 180, 360))),
}, coords={
    'time': pd.date_range('2024-01-01', periods=365),
    'lat': np.linspace(-89.5, 89.5, 180),
    'lon': np.linspace(0.5, 359.5, 360),
})
ds.attrs['title'] = 'Sample Climate Dataset'
```

### Reading Data

```python
# NetCDF (single file)
ds = xr.open_dataset('data.nc')

# NetCDF (multiple files → single dataset)
ds = xr.open_mfdataset('data_*.nc', combine='by_coords')

# Zarr (cloud-optimized)
ds = xr.open_zarr('data.zarr')

# HDF5
ds = xr.open_dataset('data.h5', engine='h5netcdf')

# Remote data via OpenDAP
ds = xr.open_dataset('https://example.com/opendap/data.nc')

# Inspect
print(ds)
print(ds.dims)       # Dimension names and sizes
print(ds.coords)     # Coordinate variables
print(ds.data_vars)  # Data variables
print(ds.attrs)      # Global attributes
```

### Writing Data

```python
# NetCDF
ds.to_netcdf('output.nc')

# Zarr (supports appending)
ds.to_zarr('output.zarr', mode='w')

# Subset before saving
ds[['temperature']].sel(lat=slice(-30, 30)).to_netcdf('tropics_temp.nc')
```

## Core Operations

### 1. Indexing and Selection

```python
# By label (recommended)
ds.sel(time='2024-06-15', lat=45.0, lon=10.0)
ds.sel(time=slice('2024-01', '2024-06'))          # Range
ds.sel(lat=45.5, method='nearest')                 # Nearest neighbor

# By position
ds.isel(time=0, lat=slice(0, 10))

# By condition
ds.where(ds['temperature'] > 0, drop=True)

# Specific variable
temp = ds['temperature']  # Returns DataArray
temp = ds.temperature      # Attribute access (same)
```

### 2. Computation with Named Dimensions

```python
# Mean over time → spatial average
spatial_mean = ds['temperature'].mean(dim='time')

# Weighted mean (area-weighted for lat-lon grids)
weights = np.cos(np.deg2rad(ds.lat))
weighted_mean = ds['temperature'].weighted(weights).mean(dim=['lat', 'lon'])

# Standard deviation over space
temporal_std = ds['temperature'].std(dim=['lat', 'lon'])

# Cumulative sum along time
cumsum = ds['precipitation'].cumsum(dim='time')
```

### 3. GroupBy Operations

```python
# Monthly climatology
climatology = ds['temperature'].groupby('time.month').mean(dim='time')

# Anomalies from climatology
anomalies = ds['temperature'].groupby('time.month') - climatology

# Seasonal averages
seasonal = ds['temperature'].groupby('time.season').mean(dim='time')

# GroupBy with custom labels
ds['hemisphere'] = xr.where(ds.lat >= 0, 'NH', 'SH')
hemisphere_mean = ds.groupby('hemisphere').mean()
```

### 4. Resampling (Time Series)

```python
# Monthly mean from daily data
monthly = ds.resample(time='MS').mean()

# Annual maximum
annual_max = ds['temperature'].resample(time='YS').max()

# 7-day rolling average
rolling_mean = ds['temperature'].rolling(time=7, center=True).mean()

# Custom aggregation
quarterly = ds.resample(time='QS').apply(lambda x: x.max() - x.min())
```

### 5. Interpolation and Regridding

```python
# Interpolate to new coordinates
new_lats = np.linspace(-90, 90, 361)
new_lons = np.linspace(0, 360, 721)
ds_hires = ds.interp(lat=new_lats, lon=new_lons, method='linear')

# Interpolate to specific points
point_values = ds['temperature'].interp(
    lat=[40.7, 34.0, 41.8],
    lon=[286.0, 241.7, 272.3],
    method='nearest'
)

# Fill NaN values
filled = ds['temperature'].interpolate_na(dim='time', method='linear')
```

### 6. Merging and Combining

```python
# Merge datasets sharing coordinates
ds_merged = xr.merge([ds_temp, ds_precip])

# Concatenate along a dimension
ds_all = xr.concat([ds_2020, ds_2021, ds_2022], dim='time')

# Combine by coordinates (auto-align)
ds_combined = xr.combine_by_coords([ds_a, ds_b])

# Align datasets to common coordinates
ds_a_aligned, ds_b_aligned = xr.align(ds_a, ds_b, join='inner')
```

### 7. Plotting

```python
import matplotlib.pyplot as plt

# Quick 2D heatmap
ds['temperature'].isel(time=0).plot(cmap='RdBu_r', figsize=(12, 6))
plt.title('Temperature Map')
plt.savefig('temp_map.png', dpi=150, bbox_inches='tight')

# Line plot (1D)
ds['temperature'].sel(lat=45, lon=10, method='nearest').plot(figsize=(10, 4))

# Contour plot
ds['temperature'].isel(time=0).plot.contourf(levels=20, cmap='coolwarm')

# Faceted plots (multiple panels)
ds['temperature'].isel(time=slice(0, 12)).plot(
    col='time', col_wrap=4, cmap='RdBu_r'
)

# Histogram
ds['temperature'].plot.hist(bins=50)
```

## Common Scientific Workflows

### Climate Anomaly Analysis

```python
import xarray as xr

# Load data
ds = xr.open_dataset('era5_temperature.nc')

# Compute climatology (1991-2020 baseline)
baseline = ds.sel(time=slice('1991', '2020'))
climatology = baseline['t2m'].groupby('time.month').mean(dim='time')

# Compute anomalies
anomalies = ds['t2m'].groupby('time.month') - climatology

# Annual mean anomaly time series
annual_anomaly = anomalies.mean(dim=['lat', 'lon']).resample(time='YS').mean()

# Trend analysis
from scipy.stats import linregress
years = annual_anomaly.time.dt.year.values
slope, intercept, r, p, se = linregress(years, annual_anomaly.values)
print(f"Warming trend: {slope*10:.3f} °C/decade (p={p:.4f})")
```

### Multi-File Processing with Dask

```python
import xarray as xr

# Lazy-load many files (data stays on disk)
ds = xr.open_mfdataset(
    'data/model_output_*.nc',
    combine='by_coords',
    chunks={'time': 100, 'lat': 90, 'lon': 180},  # Dask chunks
    parallel=True
)

# Operations are lazy until .compute() or .load()
daily_max = ds['temperature'].resample(time='D').max()
monthly_mean = daily_max.resample(time='MS').mean()

# Trigger computation
result = monthly_mean.compute()  # or .load() for in-place loading
result.to_netcdf('monthly_means.nc')
```

### Spatial Subsetting and Analysis

```python
# Select geographic region (e.g., Europe)
europe = ds.sel(lat=slice(70, 35), lon=slice(-10, 40))

# Area-weighted spatial average
weights = np.cos(np.deg2rad(europe.lat))
europe_mean = europe['temperature'].weighted(weights).mean(dim=['lat', 'lon'])

# Zonal mean (average over longitudes)
zonal_mean = ds['temperature'].mean(dim='lon')
```

## Key Parameters to Adjust

### I/O
- `engine`: Backend ('netcdf4', 'h5netcdf', 'scipy', 'zarr', 'pydap')
- `chunks`: Dask chunk sizes for lazy loading (`{'time': 100}`)
- `decode_times`: Auto-decode time coordinates (default True)

### Computing
- `dim`: Dimension(s) to reduce over
- `skipna`: Skip NaN values (default True)
- `min_count`: Minimum non-NaN values required

### Interpolation
- `method`: 'linear', 'nearest', 'cubic', 'polynomial'
- `fill_value`: Value for extrapolation ('extrapolate' or a number)

## Common Pitfalls and Best Practices

1. **Use `sel()` over `isel()`**: Label-based indexing is safer and more readable
2. **Use `chunks` for large files**: Prevents loading entire datasets into memory
3. **Check `decode_times`**: Some datasets have non-standard calendars (360-day, no-leap)
4. **Use `xr.open_mfdataset`** for multi-file datasets instead of manual loops
5. **Area-weight spatial averages**: Use `cos(lat)` weights for lat-lon grids
6. **Save intermediates in Zarr**: Faster I/O than NetCDF for large datasets
7. **Use `.compute()` strategically**: Only call when you need the result in memory
8. **Avoid chained `.sel()` calls**: Combine selections into one call for efficiency

## Integration with Other Tools

- **pandas**: `ds.to_dataframe()` / `xr.Dataset.from_dataframe(df)`
- **Dask**: Parallel/out-of-core via `chunks` parameter
- **matplotlib**: Built-in `.plot()` methods
- **cartopy**: Add map projections for geospatial plotting
- **scipy**: Use with `xr.apply_ufunc()` for element-wise operations
- **Zarr**: Cloud-optimized storage backend
- **PANGEO**: Cloud-native geoscience stack built on xarray

## Additional Resources

- **Official Documentation**: https://docs.xarray.dev/
- **Tutorial**: https://tutorial.xarray.dev/
- **Gallery**: https://docs.xarray.dev/en/stable/gallery.html
- **PANGEO**: https://pangeo.io/ (cloud-native xarray workflows)
- **cf_xarray**: https://cf-xarray.readthedocs.io/ (CF conventions support)

## Reference Documentation

- **API Reference**: See [references/api_reference.md](references/api_reference.md) for DataArray/Dataset method listing
