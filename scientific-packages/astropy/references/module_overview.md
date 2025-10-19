# Astropy Module Overview

This document provides a comprehensive reference of all major astropy subpackages and their capabilities.

## Core Data Structures

### astropy.units
**Purpose**: Handle physical units and dimensional analysis in computations.

**Key Classes**:
- `Quantity` - Combines numerical values with units
- `Unit` - Represents physical units

**Common Operations**:
```python
import astropy.units as u
distance = 5 * u.meter
time = 2 * u.second
velocity = distance / time  # Returns Quantity in m/s
wavelength = 500 * u.nm
frequency = wavelength.to(u.Hz, equivalencies=u.spectral())
```

**Equivalencies**:
- `u.spectral()` - Convert wavelength ↔ frequency
- `u.doppler_optical()`, `u.doppler_radio()` - Velocity conversions
- `u.temperature()` - Temperature unit conversions
- `u.pixel_scale()` - Pixel to physical units

### astropy.constants
**Purpose**: Provide physical and astronomical constants.

**Common Constants**:
- `c` - Speed of light
- `G` - Gravitational constant
- `h` - Planck constant
- `M_sun`, `R_sun`, `L_sun` - Solar mass, radius, luminosity
- `M_earth`, `R_earth` - Earth mass, radius
- `pc`, `au` - Parsec, astronomical unit

### astropy.time
**Purpose**: Represent and manipulate times and dates with astronomical precision.

**Time Scales**:
- `UTC` - Coordinated Universal Time
- `TAI` - International Atomic Time
- `TT` - Terrestrial Time
- `TCB`, `TCG` - Barycentric/Geocentric Coordinate Time
- `TDB` - Barycentric Dynamical Time
- `UT1` - Universal Time

**Common Formats**:
- `iso`, `isot` - ISO 8601 strings
- `jd`, `mjd` - Julian/Modified Julian Date
- `unix`, `gps` - Unix/GPS timestamps
- `datetime` - Python datetime objects

**Example**:
```python
from astropy.time import Time
t = Time('2023-01-01T00:00:00', format='isot', scale='utc')
print(t.mjd)  # Modified Julian Date
print(t.jd)   # Julian Date
print(t.tt)   # Convert to TT scale
```

### astropy.table
**Purpose**: Work with tabular data optimized for astronomical applications.

**Key Features**:
- Native support for astropy Quantity, Time, and SkyCoord columns
- Multi-dimensional columns
- Metadata preservation (units, descriptions, formats)
- Advanced operations: joins, grouping, binning
- File I/O via unified interface

**Example**:
```python
from astropy.table import Table
import astropy.units as u

t = Table()
t['name'] = ['Star1', 'Star2', 'Star3']
t['ra'] = [10.5, 11.2, 12.3] * u.degree
t['dec'] = [41.2, 42.1, 43.5] * u.degree
t['flux'] = [1.2, 2.3, 0.8] * u.Jy
```

## Coordinates and World Coordinate Systems

### astropy.coordinates
**Purpose**: Represent and transform celestial coordinates.

**Primary Interface**: `SkyCoord` - High-level class for sky positions

**Coordinate Frames**:
- `ICRS` - International Celestial Reference System (default)
- `FK5`, `FK4` - Fifth/Fourth Fundamental Katalog
- `Galactic`, `Supergalactic` - Galactic coordinates
- `AltAz` - Horizontal (altitude-azimuth) coordinates
- `GCRS`, `CIRS`, `ITRS` - Earth-based systems
- `BarycentricMeanEcliptic`, `HeliocentricMeanEcliptic`, `GeocentricMeanEcliptic` - Ecliptic coordinates

**Common Operations**:
```python
from astropy.coordinates import SkyCoord
import astropy.units as u

# Create coordinate
c = SkyCoord(ra=10.625*u.degree, dec=41.2*u.degree, frame='icrs')

# Transform to galactic
c_gal = c.galactic

# Calculate separation
c2 = SkyCoord(ra=11*u.degree, dec=42*u.degree, frame='icrs')
sep = c.separation(c2)

# Match catalogs
idx, sep2d, dist3d = c.match_to_catalog_sky(catalog_coords)
```

### astropy.wcs
**Purpose**: Handle World Coordinate System transformations for astronomical images.

**Key Class**: `WCS` - Maps between pixel and world coordinates

**Common Use Cases**:
- Convert pixel coordinates to RA/Dec
- Convert RA/Dec to pixel coordinates
- Handle distortion corrections (SIP, lookup tables)

**Example**:
```python
from astropy.wcs import WCS
from astropy.io import fits

hdu = fits.open('image.fits')[0]
wcs = WCS(hdu.header)

# Pixel to world
ra, dec = wcs.pixel_to_world_values(100, 200)

# World to pixel
x, y = wcs.world_to_pixel_values(ra, dec)
```

## File I/O

### astropy.io.fits
**Purpose**: Read and write FITS (Flexible Image Transport System) files.

**Key Classes**:
- `HDUList` - Container for all HDUs in a file
- `PrimaryHDU` - Primary header data unit
- `ImageHDU` - Image extension
- `BinTableHDU` - Binary table extension
- `Header` - FITS header keywords

**Common Operations**:
```python
from astropy.io import fits

# Read FITS file
with fits.open('file.fits') as hdul:
    hdul.info()  # Display structure
    header = hdul[0].header
    data = hdul[0].data

# Write FITS file
fits.writeto('output.fits', data, header)

# Update header keyword
fits.setval('file.fits', 'OBJECT', value='M31')
```

### astropy.io.ascii
**Purpose**: Read and write ASCII tables in various formats.

**Supported Formats**:
- Basic, CSV, tab-delimited
- CDS/MRT (Machine Readable Tables)
- IPAC, Daophot, SExtractor
- LaTeX tables
- HTML tables

### astropy.io.votable
**Purpose**: Handle Virtual Observatory (VO) table format.

### astropy.io.misc
**Purpose**: Additional formats including HDF5, Parquet, and YAML.

## Scientific Calculations

### astropy.cosmology
**Purpose**: Perform cosmological calculations.

**Common Models**:
- `FlatLambdaCDM` - Flat universe with cosmological constant (most common)
- `LambdaCDM` - Universe with cosmological constant
- `Planck18`, `Planck15`, `Planck13` - Pre-defined Planck parameters
- `WMAP9`, `WMAP7`, `WMAP5` - Pre-defined WMAP parameters

**Common Methods**:
```python
from astropy.cosmology import FlatLambdaCDM, Planck18
import astropy.units as u

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
# Or use built-in: cosmo = Planck18

z = 1.5
print(cosmo.age(z))  # Age of universe at z
print(cosmo.luminosity_distance(z))  # Luminosity distance
print(cosmo.angular_diameter_distance(z))  # Angular diameter distance
print(cosmo.comoving_distance(z))  # Comoving distance
print(cosmo.H(z))  # Hubble parameter at z
```

### astropy.modeling
**Purpose**: Framework for model evaluation and fitting.

**Model Categories**:
- 1D models: Gaussian1D, Lorentz1D, Voigt1D, Polynomial1D
- 2D models: Gaussian2D, Disk2D, Moffat2D
- Physical models: BlackBody, Drude1D, NFW
- Polynomial models: Chebyshev, Legendre

**Common Fitters**:
- `LinearLSQFitter` - Linear least squares
- `LevMarLSQFitter` - Levenberg-Marquardt
- `SimplexLSQFitter` - Downhill simplex

**Example**:
```python
from astropy.modeling import models, fitting

# Create model
g = models.Gaussian1D(amplitude=10, mean=5, stddev=1)

# Fit to data
fitter = fitting.LevMarLSQFitter()
fitted_model = fitter(g, x_data, y_data)
```

### astropy.convolution
**Purpose**: Convolve and filter astronomical data.

**Common Kernels**:
- `Gaussian2DKernel` - 2D Gaussian smoothing
- `Box2DKernel` - 2D boxcar smoothing
- `Tophat2DKernel` - 2D tophat filter
- Custom kernels via arrays

### astropy.stats
**Purpose**: Statistical tools for astronomical data analysis.

**Key Functions**:
- `sigma_clip()` - Remove outliers via sigma clipping
- `sigma_clipped_stats()` - Compute mean, median, std with clipping
- `mad_std()` - Median Absolute Deviation
- `biweight_location()`, `biweight_scale()` - Robust statistics
- `circmean()`, `circstd()` - Circular statistics

**Example**:
```python
from astropy.stats import sigma_clip, sigma_clipped_stats

# Remove outliers
filtered_data = sigma_clip(data, sigma=3, maxiters=5)

# Get statistics
mean, median, std = sigma_clipped_stats(data, sigma=3)
```

## Data Processing

### astropy.nddata
**Purpose**: Handle N-dimensional datasets with metadata.

**Key Class**: `NDData` - Container for array data with units, uncertainty, mask, and WCS

### astropy.timeseries
**Purpose**: Work with time series data.

**Key Classes**:
- `TimeSeries` - Time-indexed data table
- `BinnedTimeSeries` - Time-binned data

**Common Operations**:
- Period finding (Lomb-Scargle)
- Folding time series
- Binning data

### astropy.visualization
**Purpose**: Display astronomical data effectively.

**Key Features**:
- Image normalization (LogStretch, PowerStretch, SqrtStretch, etc.)
- Interval scaling (MinMaxInterval, PercentileInterval, ZScaleInterval)
- WCSAxes for plotting with coordinate overlays
- RGB image creation with stretching
- Astronomical colormaps

**Example**:
```python
from astropy.visualization import ImageNormalize, SqrtStretch, PercentileInterval
import matplotlib.pyplot as plt

norm = ImageNormalize(data, interval=PercentileInterval(99),
                     stretch=SqrtStretch())
plt.imshow(data, norm=norm, origin='lower')
```

## Utilities

### astropy.samp
**Purpose**: Simple Application Messaging Protocol for inter-application communication.

**Use Case**: Connect Python scripts with other astronomical tools (e.g., DS9, TOPCAT).

## Module Import Patterns

**Standard imports**:
```python
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.io import fits
from astropy.table import Table
from astropy import constants as const
```

## Performance Tips

1. **Pre-compute composite units** for repeated operations
2. **Use `<<` operator** for fast unit assignments: `array << u.m` instead of `array * u.m`
3. **Vectorize operations** rather than looping over coordinates/times
4. **Use memmap=True** when opening large FITS files
5. **Install Bottleneck** for faster stats operations
