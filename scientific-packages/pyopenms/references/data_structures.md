# pyOpenMS Data Structures Reference

This document provides comprehensive coverage of core data structures in pyOpenMS for representing mass spectrometry data.

## Core Hierarchy

```
MSExperiment                    # Top-level: Complete LC-MS/MS run
├── MSSpectrum[]               # Collection of mass spectra
│   ├── Peak1D[]              # Individual m/z, intensity pairs
│   └── SpectrumSettings      # Metadata (RT, MS level, precursor)
└── MSChromatogram[]          # Collection of chromatograms
    ├── ChromatogramPeak[]    # RT, intensity pairs
    └── ChromatogramSettings  # Metadata
```

## MSSpectrum

Represents a single mass spectrum (1-dimensional peak data).

### Creation and Basic Properties

```python
import pyopenms as oms

# Create empty spectrum
spectrum = oms.MSSpectrum()

# Set metadata
spectrum.setRT(123.45)           # Retention time in seconds
spectrum.setMSLevel(1)           # MS level (1 for MS1, 2 for MS2, etc.)
spectrum.setNativeID("scan=1234") # Native ID from file

# Additional metadata
spectrum.setDriftTime(15.2)      # Ion mobility drift time
spectrum.setName("MyScan")       # Optional name
```

### Peak Data Management

**Setting Peaks (Method 1 - Lists):**
```python
mz_values = [100.5, 200.3, 300.7, 400.2, 500.1]
intensity_values = [1000, 5000, 3000, 2000, 1500]

spectrum.set_peaks((mz_values, intensity_values))
```

**Setting Peaks (Method 2 - NumPy arrays):**
```python
import numpy as np

mz_array = np.array([100.5, 200.3, 300.7, 400.2, 500.1])
intensity_array = np.array([1000, 5000, 3000, 2000, 1500])

spectrum.set_peaks((mz_array, intensity_array))
```

**Retrieving Peaks:**
```python
# Get as numpy arrays (efficient)
mz_array, intensity_array = spectrum.get_peaks()

# Check number of peaks
n_peaks = spectrum.size()

# Get individual peak (slower)
for i in range(spectrum.size()):
    peak = spectrum[i]
    mz = peak.getMZ()
    intensity = peak.getIntensity()
```

### Precursor Information (for MS2/MSn spectra)

```python
# Create precursor
precursor = oms.Precursor()
precursor.setMZ(456.789)           # Precursor m/z
precursor.setCharge(2)             # Precursor charge
precursor.setIntensity(50000)      # Precursor intensity
precursor.setIsolationWindowLowerOffset(1.5)  # Lower isolation window
precursor.setIsolationWindowUpperOffset(1.5)  # Upper isolation window

# Set activation method
activation = oms.Activation()
activation.setActivationEnergy(35.0)  # Collision energy
activation.setMethod(oms.Activation.ActivationMethod.CID)
precursor.setActivation(activation)

# Assign to spectrum
spectrum.setPrecursors([precursor])

# Retrieve precursor information
precursors = spectrum.getPrecursors()
if len(precursors) > 0:
    prec = precursors[0]
    print(f"Precursor m/z: {prec.getMZ()}")
    print(f"Precursor charge: {prec.getCharge()}")
```

### Spectrum Metadata Access

```python
# Check if spectrum is sorted by m/z
is_sorted = spectrum.isSorted()

# Sort spectrum by m/z
spectrum.sortByPosition()

# Sort by intensity
spectrum.sortByIntensity()

# Clear all peaks
spectrum.clear(False)  # False = keep metadata, True = clear everything

# Get retention time
rt = spectrum.getRT()

# Get MS level
ms_level = spectrum.getMSLevel()
```

### Spectrum Types and Modes

```python
# Set spectrum type
spectrum.setType(oms.SpectrumSettings.SpectrumType.CENTROID)  # or PROFILE

# Get spectrum type
spec_type = spectrum.getType()
if spec_type == oms.SpectrumSettings.SpectrumType.CENTROID:
    print("Centroid spectrum")
elif spec_type == oms.SpectrumSettings.SpectrumType.PROFILE:
    print("Profile spectrum")
```

### Data Processing Annotations

```python
# Add processing information
processing = oms.DataProcessing()
processing.setMetaValue("smoothing", "gaussian")
spectrum.setDataProcessing([processing])
```

## MSExperiment

Represents a complete LC-MS/MS experiment containing multiple spectra and chromatograms.

### Creation and Population

```python
# Create empty experiment
exp = oms.MSExperiment()

# Add spectra
spectrum1 = oms.MSSpectrum()
spectrum1.setRT(100.0)
spectrum1.set_peaks(([100, 200], [1000, 2000]))

spectrum2 = oms.MSSpectrum()
spectrum2.setRT(200.0)
spectrum2.set_peaks(([100, 200], [1500, 2500]))

exp.addSpectrum(spectrum1)
exp.addSpectrum(spectrum2)

# Add chromatograms
chrom = oms.MSChromatogram()
chrom.set_peaks(([10.5, 11.0, 11.5], [1000, 5000, 3000]))
exp.addChromatogram(chrom)
```

### Accessing Spectra and Chromatograms

```python
# Get number of spectra and chromatograms
n_spectra = exp.getNrSpectra()
n_chroms = exp.getNrChromatograms()

# Access by index
first_spectrum = exp.getSpectrum(0)
last_spectrum = exp.getSpectrum(exp.getNrSpectra() - 1)

# Iterate over all spectra
for spectrum in exp:
    rt = spectrum.getRT()
    ms_level = spectrum.getMSLevel()
    n_peaks = spectrum.size()
    print(f"RT: {rt:.2f}s, MS{ms_level}, Peaks: {n_peaks}")

# Get all spectra as list
spectra = exp.getSpectra()

# Access chromatograms
chrom = exp.getChromatogram(0)
```

### Filtering Operations

```python
# Filter by MS level
exp.filterMSLevel(1)  # Keep only MS1 spectra
exp.filterMSLevel(2)  # Keep only MS2 spectra

# Filter by retention time range
exp.filterRT(100.0, 500.0)  # Keep RT between 100-500 seconds

# Filter by m/z range (all spectra)
exp.filterMZ(300.0, 1500.0)  # Keep m/z between 300-1500

# Filter by scan number
exp.filterScanNumber(100, 200)  # Keep scans 100-200
```

### Metadata and Properties

```python
# Set experiment metadata
exp.setMetaValue("operator", "John Doe")
exp.setMetaValue("instrument", "Q Exactive HF")

# Get metadata
operator = exp.getMetaValue("operator")

# Get RT range
rt_range = exp.getMinRT(), exp.getMaxRT()

# Get m/z range
mz_range = exp.getMinMZ(), exp.getMaxMZ()

# Clear all data
exp.clear(False)  # False = keep metadata
```

### Sorting and Organization

```python
# Sort spectra by retention time
exp.sortSpectra()

# Update ranges (call after modifications)
exp.updateRanges()

# Check if experiment is empty
is_empty = exp.empty()

# Reset (clear everything)
exp.reset()
```

## MSChromatogram

Represents an extracted or reconstructed chromatogram (retention time vs. intensity).

### Creation and Basic Usage

```python
# Create chromatogram
chrom = oms.MSChromatogram()

# Set peaks (RT, intensity pairs)
rt_values = [10.0, 10.5, 11.0, 11.5, 12.0]
intensity_values = [1000, 5000, 8000, 6000, 2000]
chrom.set_peaks((rt_values, intensity_values))

# Get peaks
rt_array, int_array = chrom.get_peaks()

# Get size
n_points = chrom.size()
```

### Chromatogram Types

```python
# Set chromatogram type
chrom.setChromatogramType(oms.ChromatogramSettings.ChromatogramType.SELECTED_ION_CURRENT_CHROMATOGRAM)

# Other types:
# - TOTAL_ION_CURRENT_CHROMATOGRAM
# - BASEPEAK_CHROMATOGRAM
# - SELECTED_ION_CURRENT_CHROMATOGRAM
# - SELECTED_REACTION_MONITORING_CHROMATOGRAM
```

### Metadata

```python
# Set native ID
chrom.setNativeID("TIC")

# Set name
chrom.setName("Total Ion Current")

# Access
native_id = chrom.getNativeID()
name = chrom.getName()
```

### Precursor and Product Information (for SRM/MRM)

```python
# For targeted experiments
precursor = oms.Precursor()
precursor.setMZ(456.7)
chrom.setPrecursor(precursor)

product = oms.Product()
product.setMZ(789.4)
chrom.setProduct(product)
```

## Peak1D and ChromatogramPeak

Individual peak data points.

### Peak1D (for mass spectra)

```python
# Create individual peak
peak = oms.Peak1D()
peak.setMZ(456.789)
peak.setIntensity(10000)

# Access
mz = peak.getMZ()
intensity = peak.getIntensity()

# Set position and intensity
peak.setPosition([456.789])
peak.setIntensity(10000)
```

### ChromatogramPeak (for chromatograms)

```python
# Create chromatogram peak
chrom_peak = oms.ChromatogramPeak()
chrom_peak.setRT(125.5)
chrom_peak.setIntensity(5000)

# Access
rt = chrom_peak.getRT()
intensity = chrom_peak.getIntensity()
```

## FeatureMap and Feature

For quantification results.

### Feature

Represents a detected LC-MS feature (peptide or metabolite signal).

```python
# Create feature
feature = oms.Feature()

# Set properties
feature.setMZ(456.789)
feature.setRT(123.45)
feature.setIntensity(1000000)
feature.setCharge(2)
feature.setWidth(15.0)  # RT width in seconds

# Set quality score
feature.setOverallQuality(0.95)

# Access
mz = feature.getMZ()
rt = feature.getRT()
intensity = feature.getIntensity()
charge = feature.getCharge()
```

### FeatureMap

Collection of features.

```python
# Create feature map
feature_map = oms.FeatureMap()

# Add features
feature1 = oms.Feature()
feature1.setMZ(456.789)
feature1.setRT(123.45)
feature1.setIntensity(1000000)

feature_map.push_back(feature1)

# Get size
n_features = feature_map.size()

# Iterate
for feature in feature_map:
    print(f"m/z: {feature.getMZ():.4f}, RT: {feature.getRT():.2f}")

# Access by index
first_feature = feature_map[0]

# Clear
feature_map.clear()
```

## PeptideIdentification and ProteinIdentification

For identification results.

### PeptideIdentification

```python
# Create peptide identification
pep_id = oms.PeptideIdentification()
pep_id.setRT(123.45)
pep_id.setMZ(456.789)

# Create peptide hit
hit = oms.PeptideHit()
hit.setSequence(oms.AASequence.fromString("PEPTIDE"))
hit.setCharge(2)
hit.setScore(25.5)
hit.setRank(1)

# Add to identification
pep_id.setHits([hit])
pep_id.setHigherScoreBetter(True)
pep_id.setScoreType("XCorr")

# Access
hits = pep_id.getHits()
for hit in hits:
    seq = hit.getSequence().toString()
    score = hit.getScore()
    print(f"Sequence: {seq}, Score: {score}")
```

### ProteinIdentification

```python
# Create protein identification
prot_id = oms.ProteinIdentification()

# Create protein hit
prot_hit = oms.ProteinHit()
prot_hit.setAccession("P12345")
prot_hit.setSequence("MKTAYIAKQRQISFVK...")
prot_hit.setScore(100.5)

# Add to identification
prot_id.setHits([prot_hit])
prot_id.setScoreType("Mascot Score")
prot_id.setHigherScoreBetter(True)

# Search parameters
search_params = oms.ProteinIdentification.SearchParameters()
search_params.db = "uniprot_human.fasta"
search_params.enzyme = "Trypsin"
prot_id.setSearchParameters(search_params)
```

## ConsensusMap and ConsensusFeature

For linking features across multiple samples.

### ConsensusFeature

```python
# Create consensus feature
cons_feature = oms.ConsensusFeature()
cons_feature.setMZ(456.789)
cons_feature.setRT(123.45)
cons_feature.setIntensity(5000000)  # Combined intensity

# Access linked features
for handle in cons_feature.getFeatureList():
    map_index = handle.getMapIndex()
    feature_index = handle.getIndex()
    intensity = handle.getIntensity()
```

### ConsensusMap

```python
# Create consensus map
consensus_map = oms.ConsensusMap()

# Add consensus features
consensus_map.push_back(cons_feature)

# Iterate
for cons_feat in consensus_map:
    mz = cons_feat.getMZ()
    rt = cons_feat.getRT()
    n_features = cons_feat.size()  # Number of linked features
```

## Best Practices

1. **Use numpy arrays** for peak data when possible - much faster than individual peak access
2. **Sort spectra** by position (m/z) before searching or filtering
3. **Update ranges** after modifying MSExperiment: `exp.updateRanges()`
4. **Check MS level** before processing - different algorithms for MS1 vs MS2
5. **Validate precursor info** for MS2 spectra - ensure charge and m/z are set
6. **Use appropriate containers** - MSExperiment for raw data, FeatureMap for quantification
7. **Clear metadata carefully** - use `clear(False)` to preserve metadata when clearing peaks

## Common Patterns

### Create MS2 Spectrum with Precursor

```python
spectrum = oms.MSSpectrum()
spectrum.setRT(205.2)
spectrum.setMSLevel(2)
spectrum.set_peaks(([100, 200, 300], [1000, 5000, 3000]))

precursor = oms.Precursor()
precursor.setMZ(450.5)
precursor.setCharge(2)
spectrum.setPrecursors([precursor])
```

### Extract MS1 Spectra from Experiment

```python
ms1_exp = oms.MSExperiment()
for spectrum in exp:
    if spectrum.getMSLevel() == 1:
        ms1_exp.addSpectrum(spectrum)
```

### Calculate Total Ion Current (TIC)

```python
tic_values = []
rt_values = []
for spectrum in exp:
    if spectrum.getMSLevel() == 1:
        mz, intensity = spectrum.get_peaks()
        tic = np.sum(intensity)
        tic_values.append(tic)
        rt_values.append(spectrum.getRT())
```

### Find Spectrum Closest to RT

```python
target_rt = 125.5
closest_spectrum = None
min_diff = float('inf')

for spectrum in exp:
    diff = abs(spectrum.getRT() - target_rt)
    if diff < min_diff:
        min_diff = diff
        closest_spectrum = spectrum
```
