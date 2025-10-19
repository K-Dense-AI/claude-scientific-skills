# Common Astropy Workflows

This document describes frequently used workflows when working with astronomical data using astropy.

## 1. Working with FITS Files

### Basic FITS Reading
```python
from astropy.io import fits
import numpy as np

# Open and examine structure
with fits.open('observation.fits') as hdul:
    hdul.info()

    # Access primary HDU
    primary_hdr = hdul[0].header
    primary_data = hdul[0].data

    # Access extension
    ext_data = hdul[1].data
    ext_hdr = hdul[1].header

    # Read specific header keywords
    object_name = primary_hdr['OBJECT']
    exposure = primary_hdr['EXPTIME']
```

### Writing FITS Files
```python
# Create new FITS file
from astropy.io import fits
import numpy as np

# Create data
data = np.random.random((100, 100))

# Create primary HDU
hdu = fits.PrimaryHDU(data)
hdu.header['OBJECT'] = 'M31'
hdu.header['EXPTIME'] = 300.0

# Write to file
hdu.writeto('output.fits', overwrite=True)

# Multi-extension FITS
hdul = fits.HDUList([
    fits.PrimaryHDU(data1),
    fits.ImageHDU(data2, name='SCI'),
    fits.ImageHDU(data3, name='ERR')
])
hdul.writeto('multi_ext.fits', overwrite=True)
```

### FITS Table Operations
```python
from astropy.io import fits

# Read binary table
with fits.open('catalog.fits') as hdul:
    table_data = hdul[1].data

    # Access columns
    ra = table_data['RA']
    dec = table_data['DEC']
    mag = table_data['MAG']

    # Filter data
    bright = table_data[table_data['MAG'] < 15]

# Write binary table
from astropy.table import Table
import astropy.units as u

t = Table([ra, dec, mag], names=['RA', 'DEC', 'MAG'])
t['RA'].unit = u.degree
t['DEC'].unit = u.degree
t.write('output_catalog.fits', format='fits', overwrite=True)
```

## 2. Coordinate Transformations

### Basic Coordinate Creation and Transformation
```python
from astropy.coordinates import SkyCoord
import astropy.units as u

# Create from RA/Dec
c = SkyCoord(ra=10.68458*u.degree, dec=41.26917*u.degree, frame='icrs')

# Alternative creation methods
c = SkyCoord('00:42:44.3 +41:16:09', unit=(u.hourangle, u.deg))
c = SkyCoord('00h42m44.3s +41d16m09s')

# Transform to different frames
c_gal = c.galactic
c_fk5 = c.fk5
print(f"Galactic: l={c_gal.l.deg}, b={c_gal.b.deg}")
```

### Coordinate Arrays and Separations
```python
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

# Create array of coordinates
ra_array = np.array([10.1, 10.2, 10.3]) * u.degree
dec_array = np.array([40.1, 40.2, 40.3]) * u.degree
coords = SkyCoord(ra=ra_array, dec=dec_array, frame='icrs')

# Calculate separations
c1 = SkyCoord(ra=10*u.degree, dec=40*u.degree)
c2 = SkyCoord(ra=11*u.degree, dec=41*u.degree)
sep = c1.separation(c2)
print(f"Separation: {sep.to(u.arcmin)}")

# Position angle
pa = c1.position_angle(c2)
```

### Catalog Matching
```python
from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy.units as u

# Two catalogs of coordinates
catalog1 = SkyCoord(ra=[10, 11, 12]*u.degree, dec=[40, 41, 42]*u.degree)
catalog2 = SkyCoord(ra=[10.01, 11.02, 13]*u.degree, dec=[40.01, 41.01, 43]*u.degree)

# Find nearest neighbors
idx, sep2d, dist3d = catalog1.match_to_catalog_sky(catalog2)

# Filter by separation threshold
max_sep = 1 * u.arcsec
matched = sep2d < max_sep
matching_indices = idx[matched]
```

### Horizontal Coordinates (Alt/Az)
```python
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u

# Observer location
location = EarthLocation(lat=40*u.deg, lon=-70*u.deg, height=300*u.m)

# Observation time
obstime = Time('2023-01-01 03:00:00')

# Target coordinate
target = SkyCoord(ra=10*u.degree, dec=40*u.degree, frame='icrs')

# Transform to Alt/Az
altaz_frame = AltAz(obstime=obstime, location=location)
target_altaz = target.transform_to(altaz_frame)

print(f"Altitude: {target_altaz.alt.deg}")
print(f"Azimuth: {target_altaz.az.deg}")
```

## 3. Units and Quantities

### Basic Unit Operations
```python
import astropy.units as u

# Create quantities
distance = 5.2 * u.parsec
time = 10 * u.year
velocity = 300 * u.km / u.s

# Unit conversion
distance_ly = distance.to(u.lightyear)
velocity_mps = velocity.to(u.m / u.s)

# Arithmetic with units
wavelength = 500 * u.nm
frequency = wavelength.to(u.Hz, equivalencies=u.spectral())

# Compose/decompose units
composite = (1 * u.kg * u.m**2 / u.s**2)
print(composite.decompose())  # Base SI units
print(composite.compose())     # Known compound units (Joule)
```

### Working with Arrays
```python
import numpy as np
import astropy.units as u

# Quantity arrays
wavelengths = np.array([400, 500, 600]) * u.nm
frequencies = wavelengths.to(u.THz, equivalencies=u.spectral())

# Mathematical operations preserve units
fluxes = np.array([1.2, 2.3, 1.8]) * u.Jy
luminosities = 4 * np.pi * (10*u.pc)**2 * fluxes
```

### Custom Units and Equivalencies
```python
import astropy.units as u

# Define custom unit
beam = u.def_unit('beam', 1.5e-10 * u.steradian)

# Register for session
u.add_enabled_units([beam])

# Use in calculations
flux_per_beam = 1.5 * u.Jy / beam

# Doppler equivalencies
rest_wavelength = 656.3 * u.nm  # H-alpha
observed = 656.5 * u.nm
velocity = observed.to(u.km/u.s,
                       equivalencies=u.doppler_optical(rest_wavelength))
```

## 4. Time Handling

### Time Creation and Conversion
```python
from astropy.time import Time
import astropy.units as u

# Create time objects
t1 = Time('2023-01-01T00:00:00', format='isot', scale='utc')
t2 = Time(2459945.5, format='jd', scale='utc')
t3 = Time(['2023-01-01', '2023-06-01'], format='iso')

# Convert formats
print(t1.jd)    # Julian Date
print(t1.mjd)   # Modified Julian Date
print(t1.unix)  # Unix timestamp
print(t1.iso)   # ISO format

# Convert time scales
print(t1.tai)   # Convert to TAI
print(t1.tt)    # Convert to TT
print(t1.tdb)   # Convert to TDB
```

### Time Arithmetic
```python
from astropy.time import Time, TimeDelta
import astropy.units as u

t1 = Time('2023-01-01T00:00:00')
dt = TimeDelta(1*u.day)

# Add time delta
t2 = t1 + dt

# Difference between times
diff = t2 - t1
print(diff.to(u.hour))

# Array of times
times = t1 + np.arange(10) * u.day
```

### Sidereal Time and Astronomical Calculations
```python
from astropy.time import Time
from astropy.coordinates import EarthLocation
import astropy.units as u

location = EarthLocation(lat=40*u.deg, lon=-70*u.deg)
t = Time('2023-01-01T00:00:00')

# Local sidereal time
lst = t.sidereal_time('apparent', longitude=location.lon)

# Light travel time correction
from astropy.coordinates import SkyCoord
target = SkyCoord(ra=10*u.deg, dec=40*u.deg)
ltt_bary = t.light_travel_time(target, location=location)
t_bary = t + ltt_bary
```

## 5. Tables and Data Management

### Creating and Manipulating Tables
```python
from astropy.table import Table, Column
import astropy.units as u
import numpy as np

# Create table
t = Table()
t['name'] = ['Star1', 'Star2', 'Star3']
t['ra'] = [10.5, 11.2, 12.3] * u.degree
t['dec'] = [41.2, 42.1, 43.5] * u.degree
t['flux'] = [1.2, 2.3, 0.8] * u.Jy

# Add column metadata
t['flux'].description = 'Flux at 1.4 GHz'
t['flux'].format = '.2f'

# Add new column
t['flux_mJy'] = t['flux'].to(u.mJy)

# Filter rows
bright = t[t['flux'] > 1.0 * u.Jy]

# Sort
t.sort('flux')
```

### Table I/O
```python
from astropy.table import Table

# Read various formats
t = Table.read('data.fits')
t = Table.read('data.csv', format='ascii.csv')
t = Table.read('data.ecsv', format='ascii.ecsv')  # Preserves units
t = Table.read('data.votable', format='votable')

# Write various formats
t.write('output.fits', overwrite=True)
t.write('output.csv', format='ascii.csv', overwrite=True)
t.write('output.ecsv', format='ascii.ecsv', overwrite=True)
t.write('output.votable', format='votable', overwrite=True)
```

### Advanced Table Operations
```python
from astropy.table import Table, join, vstack, hstack

# Join tables
t1 = Table([[1, 2], ['a', 'b']], names=['id', 'val1'])
t2 = Table([[1, 2], ['c', 'd']], names=['id', 'val2'])
joined = join(t1, t2, keys='id')

# Stack tables vertically
combined = vstack([t1, t2])

# Stack horizontally
combined = hstack([t1, t2])

# Grouping
t.group_by('category').groups.aggregate(np.mean)
```

### Tables with Astronomical Objects
```python
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

# Table with SkyCoord column
coords = SkyCoord(ra=[10, 11, 12]*u.deg, dec=[40, 41, 42]*u.deg)
times = Time(['2023-01-01', '2023-01-02', '2023-01-03'])

t = Table([coords, times], names=['coords', 'obstime'])

# Access individual coordinates
print(t['coords'][0].ra)
print(t['coords'][0].dec)
```

## 6. Cosmological Calculations

### Distance Calculations
```python
from astropy.cosmology import Planck18, FlatLambdaCDM
import astropy.units as u
import numpy as np

# Use built-in cosmology
cosmo = Planck18

# Redshifts
z = np.linspace(0, 5, 50)

# Calculate distances
comoving_dist = cosmo.comoving_distance(z)
angular_diam_dist = cosmo.angular_diameter_distance(z)
luminosity_dist = cosmo.luminosity_distance(z)

# Age of universe
age_at_z = cosmo.age(z)
lookback_time = cosmo.lookback_time(z)

# Hubble parameter
H_z = cosmo.H(z)
```

### Converting Observables
```python
from astropy.cosmology import Planck18
import astropy.units as u

cosmo = Planck18
z = 1.5

# Angular diameter distance
d_A = cosmo.angular_diameter_distance(z)

# Convert angular size to physical size
angular_size = 1 * u.arcsec
physical_size = (angular_size.to(u.radian) * d_A).to(u.kpc)

# Convert flux to luminosity
flux = 1e-17 * u.erg / u.s / u.cm**2
d_L = cosmo.luminosity_distance(z)
luminosity = flux * 4 * np.pi * d_L**2

# Find redshift for given distance
from astropy.cosmology import z_at_value
z_result = z_at_value(cosmo.luminosity_distance, 1000*u.Mpc)
```

### Custom Cosmology
```python
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

# Define custom cosmology
my_cosmo = FlatLambdaCDM(H0=70 * u.km/u.s/u.Mpc,
                         Om0=0.3,
                         Tcmb0=2.725 * u.K)

# Use it for calculations
print(my_cosmo.age(0))
print(my_cosmo.luminosity_distance(1.5))
```

## 7. Model Fitting

### Fitting 1D Models
```python
from astropy.modeling import models, fitting
import numpy as np
import matplotlib.pyplot as plt

# Generate data with noise
x = np.linspace(0, 10, 100)
true_model = models.Gaussian1D(amplitude=10, mean=5, stddev=1)
y = true_model(x) + np.random.normal(0, 0.5, x.shape)

# Create and fit model
g_init = models.Gaussian1D(amplitude=8, mean=4.5, stddev=0.8)
fitter = fitting.LevMarLSQFitter()
g_fit = fitter(g_init, x, y)

# Plot results
plt.plot(x, y, 'o', label='Data')
plt.plot(x, g_fit(x), label='Fit')
plt.legend()

# Get fitted parameters
print(f"Amplitude: {g_fit.amplitude.value}")
print(f"Mean: {g_fit.mean.value}")
print(f"Stddev: {g_fit.stddev.value}")
```

### Fitting with Constraints
```python
from astropy.modeling import models, fitting

# Set parameter bounds
g = models.Gaussian1D(amplitude=10, mean=5, stddev=1)
g.amplitude.bounds = (0, None)  # Positive only
g.mean.bounds = (4, 6)          # Constrain center
g.stddev.fixed = True           # Fix width

# Tie parameters (for multi-component models)
g1 = models.Gaussian1D(amplitude=10, mean=5, stddev=1, name='g1')
g2 = models.Gaussian1D(amplitude=5, mean=6, stddev=1, name='g2')
g2.stddev.tied = lambda model: model.g1.stddev

# Compound model
model = g1 + g2
```

### 2D Image Fitting
```python
from astropy.modeling import models, fitting
import numpy as np

# Create 2D data
y, x = np.mgrid[0:100, 0:100]
z = models.Gaussian2D(amplitude=100, x_mean=50, y_mean=50,
                     x_stddev=5, y_stddev=5)(x, y)
z += np.random.normal(0, 5, z.shape)

# Fit 2D Gaussian
g_init = models.Gaussian2D(amplitude=90, x_mean=48, y_mean=48,
                          x_stddev=4, y_stddev=4)
fitter = fitting.LevMarLSQFitter()
g_fit = fitter(g_init, x, y, z)

# Get parameters
print(f"Center: ({g_fit.x_mean.value}, {g_fit.y_mean.value})")
print(f"Width: ({g_fit.x_stddev.value}, {g_fit.y_stddev.value})")
```

## 8. Image Processing and Visualization

### Image Display with Proper Scaling
```python
from astropy.io import fits
from astropy.visualization import ImageNormalize, SqrtStretch, PercentileInterval
import matplotlib.pyplot as plt

# Read FITS image
data = fits.getdata('image.fits')

# Apply normalization
norm = ImageNormalize(data,
                     interval=PercentileInterval(99),
                     stretch=SqrtStretch())

# Display
plt.imshow(data, norm=norm, origin='lower', cmap='gray')
plt.colorbar()
```

### WCS Plotting
```python
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import ImageNormalize, LogStretch, PercentileInterval
import matplotlib.pyplot as plt

# Read FITS with WCS
hdu = fits.open('image.fits')[0]
wcs = WCS(hdu.header)
data = hdu.data

# Create figure with WCS projection
fig = plt.figure()
ax = fig.add_subplot(111, projection=wcs)

# Plot with coordinate grid
norm = ImageNormalize(data, interval=PercentileInterval(99.5),
                     stretch=LogStretch())
im = ax.imshow(data, norm=norm, origin='lower', cmap='viridis')

# Add coordinate labels
ax.set_xlabel('RA')
ax.set_ylabel('Dec')
ax.coords.grid(color='white', alpha=0.5)
plt.colorbar(im)
```

### Sigma Clipping and Statistics
```python
from astropy.stats import sigma_clip, sigma_clipped_stats
import numpy as np

# Data with outliers
data = np.random.normal(100, 15, 1000)
data[0:50] = np.random.normal(200, 10, 50)  # Add outliers

# Sigma clipping
clipped = sigma_clip(data, sigma=3, maxiters=5)

# Get statistics on clipped data
mean, median, std = sigma_clipped_stats(data, sigma=3)

print(f"Mean: {mean:.2f}")
print(f"Median: {median:.2f}")
print(f"Std: {std:.2f}")
print(f"Clipped {clipped.mask.sum()} values")
```

## 9. Complete Analysis Example

### Photometry Pipeline
```python
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ImageNormalize, LogStretch
import astropy.units as u
import numpy as np

# Read FITS file
hdu = fits.open('observation.fits')[0]
data = hdu.data
header = hdu.header
wcs = WCS(header)

# Calculate background statistics
mean, median, std = sigma_clipped_stats(data, sigma=3.0)
print(f"Background: {median:.2f} +/- {std:.2f}")

# Subtract background
data_sub = data - median

# Known source coordinates
source_coord = SkyCoord(ra='10:42:30', dec='+41:16:09', unit=(u.hourangle, u.deg))

# Convert to pixel coordinates
x_pix, y_pix = wcs.world_to_pixel(source_coord)

# Simple aperture photometry
aperture_radius = 10  # pixels
y, x = np.ogrid[:data.shape[0], :data.shape[1]]
mask = (x - x_pix)**2 + (y - y_pix)**2 <= aperture_radius**2

aperture_sum = np.sum(data_sub[mask])
npix = np.sum(mask)

print(f"Source position: ({x_pix:.1f}, {y_pix:.1f})")
print(f"Aperture sum: {aperture_sum:.2f}")
print(f"S/N: {aperture_sum / (std * np.sqrt(npix)):.2f}")
```

This workflow document provides practical examples for common astronomical data analysis tasks using astropy.
