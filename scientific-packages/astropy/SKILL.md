---
name: astropy
description: "Astronomy toolkit. FITS I/O, celestial coordinate transforms, cosmology calculations, time systems, WCS, units, astronomical tables, for astronomical data analysis and imaging."
---

# Astropy

## Overview

Astropy is the community standard Python library for astronomy, providing core functionality for astronomical data analysis and computation. This skill provides comprehensive guidance and tools for working with astropy's extensive capabilities across coordinate systems, file I/O, units and quantities, time systems, cosmology, modeling, and more.

## When to Use This Skill

This skill should be used when:
- Working with FITS files (reading, writing, inspecting, modifying)
- Performing coordinate transformations between astronomical reference frames
- Calculating cosmological distances, ages, or other quantities
- Handling astronomical time systems and conversions
- Working with physical units and dimensional analysis
- Processing astronomical data tables with specialized column types
- Fitting models to astronomical data
- Converting between pixel and world coordinates (WCS)
- Performing robust statistical analysis on astronomical data
- Visualizing astronomical images with proper scaling and stretching

## Core Capabilities

### 1. FITS File Operations

FITS (Flexible Image Transport System) is the standard file format in astronomy. Astropy provides comprehensive FITS support.

**Quick FITS Inspection**:
Use the included `scripts/fits_info.py` script for rapid file inspection:
```bash
python scripts/fits_info.py observation.fits
python scripts/fits_info.py observation.fits --detailed
python scripts/fits_info.py observation.fits --ext 1
```

**Common FITS workflows**:
```python
from astropy.io import fits

# Read FITS file
with fits.open('image.fits') as hdul:
    hdul.info()  # Display structure
    data = hdul[0].data
    header = hdul[0].header

# Write FITS file
fits.writeto('output.fits', data, header, overwrite=True)

# Quick access (less efficient for multiple operations)
data = fits.getdata('image.fits', ext=0)
header = fits.getheader('image.fits', ext=0)

# Update specific header keyword
fits.setval('image.fits', 'OBJECT', value='M31')
```

**Multi-extension FITS**:
```python
from astropy.io import fits

# Create multi-extension FITS
primary = fits.PrimaryHDU(primary_data)
image_ext = fits.ImageHDU(science_data, name='SCI')
error_ext = fits.ImageHDU(error_data, name='ERR')

hdul = fits.HDUList([primary, image_ext, error_ext])
hdul.writeto('multi_ext.fits', overwrite=True)
```

**Binary tables**:
```python
from astropy.io import fits

# Read binary table
with fits.open('catalog.fits') as hdul:
    table_data = hdul[1].data
    ra = table_data['RA']
    dec = table_data['DEC']

# Better: use astropy.table for table operations (see section 5)
```

### 2. Coordinate Systems and Transformations

Astropy supports ~25 coordinate frames with seamless transformations.

**Quick Coordinate Conversion**:
Use the included `scripts/coord_convert.py` script:
```bash
python scripts/coord_convert.py 10.68 41.27 --from icrs --to galactic
python scripts/coord_convert.py --file coords.txt --from icrs --to galactic --output sexagesimal
```

**Basic coordinate operations**:
```python
from astropy.coordinates import SkyCoord
import astropy.units as u

# Create coordinate (multiple input formats supported)
c = SkyCoord(ra=10.68*u.degree, dec=41.27*u.degree, frame='icrs')
c = SkyCoord('00:42:44.3 +41:16:09', unit=(u.hourangle, u.deg))
c = SkyCoord('00h42m44.3s +41d16m09s')

# Transform between frames
c_galactic = c.galactic
c_fk5 = c.fk5

print(f"Galactic: l={c_galactic.l.deg:.3f}, b={c_galactic.b.deg:.3f}")
```

**Working with coordinate arrays**:
```python
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

# Arrays of coordinates
ra = np.array([10.1, 10.2, 10.3]) * u.degree
dec = np.array([40.1, 40.2, 40.3]) * u.degree
coords = SkyCoord(ra=ra, dec=dec, frame='icrs')

# Calculate separations
sep = coords[0].separation(coords[1])
print(f"Separation: {sep.to(u.arcmin)}")

# Position angle
pa = coords[0].position_angle(coords[1])
```

**Catalog matching**:
```python
from astropy.coordinates import SkyCoord
import astropy.units as u

catalog1 = SkyCoord(ra=[10, 11, 12]*u.degree, dec=[40, 41, 42]*u.degree)
catalog2 = SkyCoord(ra=[10.01, 11.02, 13]*u.degree, dec=[40.01, 41.01, 43]*u.degree)

# Find nearest neighbors
idx, sep2d, dist3d = catalog1.match_to_catalog_sky(catalog2)

# Filter by separation threshold
max_sep = 1 * u.arcsec
matched = sep2d < max_sep
```

**Horizontal coordinates (Alt/Az)**:
```python
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u

location = EarthLocation(lat=40*u.deg, lon=-70*u.deg, height=300*u.m)
obstime = Time('2023-01-01 03:00:00')
target = SkyCoord(ra=10*u.degree, dec=40*u.degree, frame='icrs')

altaz_frame = AltAz(obstime=obstime, location=location)
target_altaz = target.transform_to(altaz_frame)

print(f"Alt: {target_altaz.alt.deg:.2f}°, Az: {target_altaz.az.deg:.2f}°")
```

**Available coordinate frames**:
- `icrs` - International Celestial Reference System (default, preferred)
- `fk5`, `fk4` - Fifth/Fourth Fundamental Katalog
- `galactic` - Galactic coordinates
- `supergalactic` - Supergalactic coordinates
- `altaz` - Horizontal (altitude-azimuth) coordinates
- `gcrs`, `cirs`, `itrs` - Earth-based systems
- Ecliptic frames: `BarycentricMeanEcliptic`, `HeliocentricMeanEcliptic`, `GeocentricMeanEcliptic`

### 3. Units and Quantities

Physical units are fundamental to astronomical calculations. Astropy's units system provides dimensional analysis and automatic conversions.

**Basic unit operations**:
```python
import astropy.units as u

# Create quantities
distance = 5.2 * u.parsec
velocity = 300 * u.km / u.s
time = 10 * u.year

# Convert units
distance_ly = distance.to(u.lightyear)
velocity_mps = velocity.to(u.m / u.s)

# Arithmetic with units
wavelength = 500 * u.nm
frequency = wavelength.to(u.Hz, equivalencies=u.spectral())
```

**Working with arrays**:
```python
import numpy as np
import astropy.units as u

wavelengths = np.array([400, 500, 600]) * u.nm
frequencies = wavelengths.to(u.THz, equivalencies=u.spectral())

fluxes = np.array([1.2, 2.3, 1.8]) * u.Jy
luminosities = 4 * np.pi * (10*u.pc)**2 * fluxes
```

**Important equivalencies**:
- `u.spectral()` - Convert wavelength ↔ frequency ↔ energy
- `u.doppler_optical(rest)` - Optical Doppler velocity
- `u.doppler_radio(rest)` - Radio Doppler velocity
- `u.doppler_relativistic(rest)` - Relativistic Doppler
- `u.temperature()` - Temperature unit conversions
- `u.brightness_temperature(freq)` - Brightness temperature

**Physical constants**:
```python
from astropy import constants as const

print(const.c)      # Speed of light
print(const.G)      # Gravitational constant
print(const.M_sun)  # Solar mass
print(const.R_sun)  # Solar radius
print(const.L_sun)  # Solar luminosity
```

**Performance tip**: Use the `<<` operator for fast unit assignment to arrays:
```python
# Fast
result = large_array << u.m

# Slower
result = large_array * u.m
```

### 4. Time Systems

Astronomical time systems require high precision and multiple time scales.

**Creating time objects**:
```python
from astropy.time import Time
import astropy.units as u

# Various input formats
t1 = Time('2023-01-01T00:00:00', format='isot', scale='utc')
t2 = Time(2459945.5, format='jd', scale='utc')
t3 = Time(['2023-01-01', '2023-06-01'], format='iso')

# Convert formats
print(t1.jd)    # Julian Date
print(t1.mjd)   # Modified Julian Date
print(t1.unix)  # Unix timestamp
print(t1.iso)   # ISO format

# Convert time scales
print(t1.tai)   # International Atomic Time
print(t1.tt)    # Terrestrial Time
print(t1.tdb)   # Barycentric Dynamical Time
```

**Time arithmetic**:
```python
from astropy.time import Time, TimeDelta
import astropy.units as u

t1 = Time('2023-01-01T00:00:00')
dt = TimeDelta(1*u.day)

t2 = t1 + dt
diff = t2 - t1
print(diff.to(u.hour))

# Array of times
times = t1 + np.arange(10) * u.day
```

**Astronomical time calculations**:
```python
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u

location = EarthLocation(lat=40*u.deg, lon=-70*u.deg)
t = Time('2023-01-01T00:00:00')

# Local sidereal time
lst = t.sidereal_time('apparent', longitude=location.lon)

# Barycentric correction
target = SkyCoord(ra=10*u.deg, dec=40*u.deg)
ltt = t.light_travel_time(target, location=location)
t_bary = t.tdb + ltt
```

**Available time scales**:
- `utc` - Coordinated Universal Time
- `tai` - International Atomic Time
- `tt` - Terrestrial Time
- `tcb`, `tcg` - Barycentric/Geocentric Coordinate Time
- `tdb` - Barycentric Dynamical Time
- `ut1` - Universal Time

### 5. Data Tables

Astropy tables provide astronomy-specific enhancements over pandas.

**Creating and manipulating tables**:
```python
from astropy.table import Table
import astropy.units as u

# Create table
t = Table()
t['name'] = ['Star1', 'Star2', 'Star3']
t['ra'] = [10.5, 11.2, 12.3] * u.degree
t['dec'] = [41.2, 42.1, 43.5] * u.degree
t['flux'] = [1.2, 2.3, 0.8] * u.Jy

# Column metadata
t['flux'].description = 'Flux at 1.4 GHz'
t['flux'].format = '.2f'

# Add calculated column
t['flux_mJy'] = t['flux'].to(u.mJy)

# Filter and sort
bright = t[t['flux'] > 1.0 * u.Jy]
t.sort('flux')
```

**Table I/O**:
```python
from astropy.table import Table

# Read (format auto-detected from extension)
t = Table.read('data.fits')
t = Table.read('data.csv', format='ascii.csv')
t = Table.read('data.ecsv', format='ascii.ecsv')  # Preserves units!
t = Table.read('data.votable', format='votable')

# Write
t.write('output.fits', overwrite=True)
t.write('output.ecsv', format='ascii.ecsv', overwrite=True)
```

**Advanced operations**:
```python
from astropy.table import Table, join, vstack, hstack

# Join tables (like SQL)
joined = join(table1, table2, keys='id')

# Stack tables
combined_rows = vstack([t1, t2])
combined_cols = hstack([t1, t2])

# Grouping and aggregation
t.group_by('category').groups.aggregate(np.mean)
```

**Tables with astronomical objects**:
```python
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

coords = SkyCoord(ra=[10, 11, 12]*u.deg, dec=[40, 41, 42]*u.deg)
times = Time(['2023-01-01', '2023-01-02', '2023-01-03'])

t = Table([coords, times], names=['coords', 'obstime'])
print(t['coords'][0].ra)  # Access coordinate properties
```

### 6. Cosmological Calculations

Quick cosmology calculations using standard models.

**Using the cosmology calculator**:
```bash
python scripts/cosmo_calc.py 0.5 1.0 1.5
python scripts/cosmo_calc.py --range 0 3 0.5 --cosmology Planck18
python scripts/cosmo_calc.py 0.5 --verbose
python scripts/cosmo_calc.py --convert 1000 --from luminosity_distance
```

**Programmatic usage**:
```python
from astropy.cosmology import Planck18
import astropy.units as u
import numpy as np

cosmo = Planck18

# Calculate distances
z = 1.5
d_L = cosmo.luminosity_distance(z)
d_A = cosmo.angular_diameter_distance(z)
d_C = cosmo.comoving_distance(z)

# Time calculations
age = cosmo.age(z)
lookback = cosmo.lookback_time(z)

# Hubble parameter
H_z = cosmo.H(z)

print(f"At z={z}:")
print(f"  Luminosity distance: {d_L:.2f}")
print(f"  Age of universe: {age:.2f}")
```

**Convert observables**:
```python
from astropy.cosmology import Planck18
import astropy.units as u

cosmo = Planck18
z = 1.5

# Angular size to physical size
d_A = cosmo.angular_diameter_distance(z)
angular_size = 1 * u.arcsec
physical_size = (angular_size.to(u.radian) * d_A).to(u.kpc)

# Flux to luminosity
flux = 1e-17 * u.erg / u.s / u.cm**2
d_L = cosmo.luminosity_distance(z)
luminosity = flux * 4 * np.pi * d_L**2

# Find redshift for given distance
from astropy.cosmology import z_at_value
z = z_at_value(cosmo.luminosity_distance, 1000*u.Mpc)
```

**Available cosmologies**:
- `Planck18`, `Planck15`, `Planck13` - Planck satellite parameters
- `WMAP9`, `WMAP7`, `WMAP5` - WMAP satellite parameters
- Custom: `FlatLambdaCDM(H0=70*u.km/u.s/u.Mpc, Om0=0.3)`

### 7. Model Fitting

Fit mathematical models to astronomical data.

**1D fitting example**:
```python
from astropy.modeling import models, fitting
import numpy as np

# Generate data
x = np.linspace(0, 10, 100)
y_data = 10 * np.exp(-0.5 * ((x - 5) / 1)**2) + np.random.normal(0, 0.5, x.shape)

# Create and fit model
g_init = models.Gaussian1D(amplitude=8, mean=4.5, stddev=0.8)
fitter = fitting.LevMarLSQFitter()
g_fit = fitter(g_init, x, y_data)

# Results
print(f"Amplitude: {g_fit.amplitude.value:.3f}")
print(f"Mean: {g_fit.mean.value:.3f}")
print(f"Stddev: {g_fit.stddev.value:.3f}")

# Evaluate fitted model
y_fit = g_fit(x)
```

**Common 1D models**:
- `Gaussian1D` - Gaussian profile
- `Lorentz1D` - Lorentzian profile
- `Voigt1D` - Voigt profile
- `Moffat1D` - Moffat profile (PSF modeling)
- `Polynomial1D` - Polynomial
- `PowerLaw1D` - Power law
- `BlackBody` - Blackbody spectrum

**Common 2D models**:
- `Gaussian2D` - 2D Gaussian
- `Moffat2D` - 2D Moffat (stellar PSF)
- `AiryDisk2D` - Airy disk (diffraction pattern)
- `Disk2D` - Circular disk

**Fitting with constraints**:
```python
from astropy.modeling import models, fitting

g = models.Gaussian1D(amplitude=10, mean=5, stddev=1)

# Set bounds
g.amplitude.bounds = (0, None)  # Positive only
g.mean.bounds = (4, 6)          # Constrain center

# Fix parameters
g.stddev.fixed = True

# Compound models
model = models.Gaussian1D() + models.Polynomial1D(degree=1)
```

**Available fitters**:
- `LinearLSQFitter` - Linear least squares (fast, for linear models)
- `LevMarLSQFitter` - Levenberg-Marquardt (most common)
- `SimplexLSQFitter` - Downhill simplex
- `SLSQPLSQFitter` - Sequential Least Squares with constraints

### 8. World Coordinate System (WCS)

Transform between pixel and world coordinates in images.

**Basic WCS usage**:
```python
from astropy.io import fits
from astropy.wcs import WCS

# Read FITS with WCS
hdu = fits.open('image.fits')[0]
wcs = WCS(hdu.header)

# Pixel to world
ra, dec = wcs.pixel_to_world_values(100, 200)

# World to pixel
x, y = wcs.world_to_pixel_values(ra, dec)

# Using SkyCoord (more powerful)
from astropy.coordinates import SkyCoord
import astropy.units as u

coord = SkyCoord(ra=150*u.deg, dec=-30*u.deg)
x, y = wcs.world_to_pixel(coord)
```

**Plotting with WCS**:
```python
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import ImageNormalize, LogStretch, PercentileInterval
import matplotlib.pyplot as plt

hdu = fits.open('image.fits')[0]
wcs = WCS(hdu.header)
data = hdu.data

# Create figure with WCS projection
fig = plt.figure()
ax = fig.add_subplot(111, projection=wcs)

# Plot with coordinate grid
norm = ImageNormalize(data, interval=PercentileInterval(99.5),
                     stretch=LogStretch())
ax.imshow(data, norm=norm, origin='lower', cmap='viridis')

# Coordinate labels and grid
ax.set_xlabel('RA')
ax.set_ylabel('Dec')
ax.coords.grid(color='white', alpha=0.5)
```

### 9. Statistics and Data Processing

Robust statistical tools for astronomical data.

**Sigma clipping** (remove outliers):
```python
from astropy.stats import sigma_clip, sigma_clipped_stats

# Remove outliers
clipped = sigma_clip(data, sigma=3, maxiters=5)

# Get statistics on cleaned data
mean, median, std = sigma_clipped_stats(data, sigma=3)

# Use clipped data
background = median
signal = data - background
snr = signal / std
```

**Other statistical functions**:
```python
from astropy.stats import mad_std, biweight_location, biweight_scale

# Robust standard deviation
std_robust = mad_std(data)

# Robust central location
center = biweight_location(data)

# Robust scale
scale = biweight_scale(data)
```

### 10. Visualization

Display astronomical images with proper scaling.

**Image normalization and stretching**:
```python
from astropy.visualization import (ImageNormalize, MinMaxInterval,
                                   PercentileInterval, ZScaleInterval,
                                   SqrtStretch, LogStretch, PowerStretch,
                                   AsinhStretch)
import matplotlib.pyplot as plt

# Common combination: percentile interval + sqrt stretch
norm = ImageNormalize(data,
                     interval=PercentileInterval(99),
                     stretch=SqrtStretch())

plt.imshow(data, norm=norm, origin='lower', cmap='gray')
plt.colorbar()
```

**Available intervals** (determine min/max):
- `MinMaxInterval()` - Use actual min/max
- `PercentileInterval(percentile)` - Clip to percentile (e.g., 99%)
- `ZScaleInterval()` - IRAF's zscale algorithm
- `ManualInterval(vmin, vmax)` - Specify manually

**Available stretches** (nonlinear scaling):
- `LinearStretch()` - Linear (default)
- `SqrtStretch()` - Square root (common for images)
- `LogStretch()` - Logarithmic (for high dynamic range)
- `PowerStretch(power)` - Power law
- `AsinhStretch()` - Arcsinh (good for wide range)

## Bundled Resources

### scripts/

**`fits_info.py`** - Comprehensive FITS file inspection tool
```bash
python scripts/fits_info.py observation.fits
python scripts/fits_info.py observation.fits --detailed
python scripts/fits_info.py observation.fits --ext 1
```

**`coord_convert.py`** - Batch coordinate transformation utility
```bash
python scripts/coord_convert.py 10.68 41.27 --from icrs --to galactic
python scripts/coord_convert.py --file coords.txt --from icrs --to galactic
```

**`cosmo_calc.py`** - Cosmological calculator
```bash
python scripts/cosmo_calc.py 0.5 1.0 1.5
python scripts/cosmo_calc.py --range 0 3 0.5 --cosmology Planck18
```

### references/

**`module_overview.md`** - Comprehensive reference of all astropy subpackages, classes, and methods. Consult this for detailed API information, available functions, and module capabilities.

**`common_workflows.md`** - Complete working examples for common astronomical data analysis tasks. Contains full code examples for FITS operations, coordinate transformations, cosmology, modeling, and complete analysis pipelines.

## Best Practices

1. **Use context managers for FITS files**:
   ```python
   with fits.open('file.fits') as hdul:
       # Work with file
   ```

2. **Prefer astropy.table over raw FITS tables** for better unit/metadata support

3. **Use SkyCoord for coordinates** (high-level interface) rather than low-level frame classes

4. **Always attach units** to quantities when possible for dimensional safety

5. **Use ECSV format** for saving tables when you want to preserve units and metadata

6. **Vectorize coordinate operations** rather than looping for performance

7. **Use memmap=True** when opening large FITS files to save memory

8. **Install Bottleneck** package for faster statistics operations

9. **Pre-compute composite units** for repeated operations to improve performance

10. **Consult `references/module_overview.md`** for detailed module information and `references/common_workflows.md`** for complete working examples

## Common Patterns

### Pattern: FITS → Process → FITS
```python
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

# Read
with fits.open('input.fits') as hdul:
    data = hdul[0].data
    header = hdul[0].header

    # Process
    mean, median, std = sigma_clipped_stats(data, sigma=3)
    processed = (data - median) / std

    # Write
    fits.writeto('output.fits', processed, header, overwrite=True)
```

### Pattern: Catalog Matching
```python
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u

# Load catalogs
cat1 = Table.read('catalog1.fits')
cat2 = Table.read('catalog2.fits')

# Create coordinate objects
coords1 = SkyCoord(ra=cat1['RA'], dec=cat1['DEC'], unit=u.degree)
coords2 = SkyCoord(ra=cat2['RA'], dec=cat2['DEC'], unit=u.degree)

# Match
idx, sep2d, dist3d = coords1.match_to_catalog_sky(coords2)

# Filter by separation
max_sep = 1 * u.arcsec
matched_mask = sep2d < max_sep

# Create matched catalog
matched_cat1 = cat1[matched_mask]
matched_cat2 = cat2[idx[matched_mask]]
```

### Pattern: Time Series Analysis
```python
from astropy.time import Time
from astropy.timeseries import TimeSeries
import astropy.units as u

# Create time series
times = Time(['2023-01-01', '2023-01-02', '2023-01-03'])
flux = [1.2, 2.3, 1.8] * u.Jy

ts = TimeSeries(time=times)
ts['flux'] = flux

# Fold on period
from astropy.timeseries import aggregate_downsample
period = 1.5 * u.day
folded = ts.fold(period=period)
```

### Pattern: Image Display with WCS
```python
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import ImageNormalize, SqrtStretch, PercentileInterval
import matplotlib.pyplot as plt

hdu = fits.open('image.fits')[0]
wcs = WCS(hdu.header)
data = hdu.data

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection=wcs)

norm = ImageNormalize(data, interval=PercentileInterval(99),
                     stretch=SqrtStretch())
im = ax.imshow(data, norm=norm, origin='lower', cmap='viridis')

ax.set_xlabel('RA')
ax.set_ylabel('Dec')
ax.coords.grid(color='white', alpha=0.5, linestyle='solid')
plt.colorbar(im, ax=ax)
```

## Installation Note

Ensure astropy is installed in the Python environment:
```bash
pip install astropy
```

For additional performance and features:
```bash
pip install astropy[all]  # Includes optional dependencies
```

## Additional Resources

- Official documentation: https://docs.astropy.org
- Tutorials: https://learn.astropy.org
- API reference: Consult `references/module_overview.md` in this skill
- Working examples: Consult `references/common_workflows.md` in this skill
