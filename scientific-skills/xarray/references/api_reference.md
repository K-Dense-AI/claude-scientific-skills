# xarray API Quick Reference

## Creating Data

| Constructor | Purpose | Key Parameters |
|-------------|---------|----------------|
| `xr.DataArray(data, dims, coords)` | Single labeled array | `attrs`, `name` |
| `xr.Dataset(data_vars, coords)` | Collection of DataArrays | `attrs` |
| `xr.zeros_like(other)` | Array of zeros matching shape | `dtype` |
| `xr.ones_like(other)` | Array of ones matching shape | `dtype` |
| `xr.full_like(other, fill_value)` | Array filled with value | `dtype` |

## I/O Functions

| Function | Purpose | Key Parameters |
|----------|---------|----------------|
| `xr.open_dataset(path)` | Open single file | `engine`, `chunks`, `decode_times` |
| `xr.open_mfdataset(paths)` | Open multiple files | `combine`, `chunks`, `parallel` |
| `xr.open_zarr(store)` | Open Zarr store | `chunks`, `consolidated` |
| `xr.open_dataarray(path)` | Open single DataArray | `engine`, `chunks` |
| `ds.to_netcdf(path)` | Save as NetCDF | `format`, `engine`, `encoding` |
| `ds.to_zarr(store)` | Save as Zarr | `mode`, `append_dim`, `consolidated` |
| `ds.to_dataframe()` | Convert to pandas | `dim_order` |

### Engines
`'netcdf4'`, `'h5netcdf'`, `'scipy'`, `'zarr'`, `'pydap'`

## DataArray Methods

### Selection & Indexing

| Method | Purpose | Example |
|--------|---------|---------|
| `.sel(**kwargs)` | Select by label | `da.sel(time='2024-01', lat=45.0)` |
| `.isel(**kwargs)` | Select by position | `da.isel(time=0, lat=slice(0,10))` |
| `.loc[...]` | Label-based indexing | `da.loc['2024-01':'2024-06']` |
| `.where(cond)` | Conditional selection | `da.where(da > 0, drop=True)` |
| `.drop_sel(**kwargs)` | Drop by label | `da.drop_sel(time='2024-01-01')` |
| `.squeeze()` | Remove size-1 dims | `da.squeeze('level')` |

### Computation

| Method | Purpose | Key Parameters |
|--------|---------|----------------|
| `.mean(dim)` | Mean | `skipna`, `min_count` |
| `.sum(dim)` | Sum | `skipna`, `min_count` |
| `.std(dim)` | Standard deviation | `ddof` |
| `.var(dim)` | Variance | `ddof` |
| `.min(dim)` / `.max(dim)` | Min/Max | `skipna` |
| `.median(dim)` | Median | `skipna` |
| `.quantile(q, dim)` | Quantile | `interpolation` |
| `.cumsum(dim)` | Cumulative sum | `skipna` |
| `.diff(dim)` | Differences | `n` (order) |
| `.weighted(weights)` | Weighted operations | chain with `.mean()`, `.sum()` |
| `.rank(dim)` | Rank values | `pct` |
| `.coarsen(**kwargs)` | Block reduce | chain with `.mean()`, etc. |

### GroupBy & Resample

| Method | Purpose | Example |
|--------|---------|---------|
| `.groupby(group)` | GroupBy operation | `da.groupby('time.month').mean()` |
| `.groupby_bins(group, bins)` | Binned groupby | `da.groupby_bins('lat', bins=18).mean()` |
| `.resample(**kwargs)` | Time resampling | `da.resample(time='MS').mean()` |
| `.rolling(**kwargs)` | Rolling window | `da.rolling(time=7).mean()` |

### Reshaping

| Method | Purpose |
|--------|---------|
| `.stack(new_dim=dims)` | Stack dimensions |
| `.unstack(dim)` | Unstack dimension |
| `.transpose(*dims)` | Reorder dimensions |
| `.expand_dims(dim)` | Add new dimension |
| `.swap_dims({old: new})` | Swap dimension names |
| `.rename({old: new})` | Rename dims/coords |

### Interpolation

| Method | Purpose | Key Parameters |
|--------|---------|----------------|
| `.interp(**kwargs)` | Interpolate to coords | `method`, `kwargs` |
| `.interp_like(other)` | Interpolate to match | `method` |
| `.interpolate_na(dim)` | Fill NaN | `method`, `limit`, `fill_value` |
| `.ffill(dim)` | Forward fill | `limit` |
| `.bfill(dim)` | Backward fill | `limit` |

### Plotting

| Method | Purpose | Key Parameters |
|--------|---------|----------------|
| `.plot()` | Auto-detect plot type | `cmap`, `vmin`, `vmax`, `figsize` |
| `.plot.contourf()` | Filled contours | `levels`, `cmap` |
| `.plot.contour()` | Line contours | `levels`, `colors` |
| `.plot.pcolormesh()` | Pseudocolor mesh | `cmap`, `edgecolors` |
| `.plot.line()` | Line plot | `hue`, `col` |
| `.plot.hist()` | Histogram | `bins` |
| `.plot.scatter()` | Scatter plot | `hue`, `size` |

## Dataset Methods

| Method | Purpose |
|--------|---------|
| `ds[var]` | Access DataArray |
| `ds.drop_vars(names)` | Drop variables |
| `ds.rename_vars({old: new})` | Rename variables |
| `ds.assign(new_var=da)` | Add variable |
| `ds.info()` | Summary information |

## Combining Data

| Function | Purpose | Key Parameters |
|----------|---------|----------------|
| `xr.merge(datasets)` | Merge by coordinates | `join`, `compat` |
| `xr.concat(datasets, dim)` | Concatenate along dim | `data_vars`, `coords` |
| `xr.combine_by_coords(datasets)` | Auto-combine | `combine_attrs` |
| `xr.align(*datasets)` | Align coordinates | `join` ('inner', 'outer', 'left', 'right') |

## apply_ufunc

Apply custom functions element-wise with dimension awareness:

```python
result = xr.apply_ufunc(
    np.gradient,           # Function to apply
    da,                    # Input DataArray
    kwargs={'axis': -1},   # Function kwargs
    input_core_dims=[['lon']],   # Dims to pass as array axes
    output_core_dims=[['lon']],  # Dims of output
    dask='parallelized',   # Enable Dask support
    output_dtypes=[float],
)
```

## Time Accessor (`.dt`)

| Accessor | Returns |
|----------|---------|
| `da.time.dt.year` | Year |
| `da.time.dt.month` | Month (1-12) |
| `da.time.dt.day` | Day |
| `da.time.dt.hour` | Hour |
| `da.time.dt.dayofweek` | Day of week (0=Mon) |
| `da.time.dt.dayofyear` | Day of year |
| `da.time.dt.season` | Season ('DJF', 'MAM', 'JJA', 'SON') |
| `da.time.dt.is_leap_year` | Leap year flag |
| `da.time.dt.strftime(fmt)` | Formatted string |

## Resample Frequency Aliases

| Alias | Meaning |
|-------|---------|
| `'D'` | Calendar day |
| `'W'` | Weekly |
| `'MS'` | Month start |
| `'ME'` | Month end |
| `'QS'` | Quarter start |
| `'YS'` | Year start |
| `'h'` | Hourly |
| `'min'` | Minutely |
