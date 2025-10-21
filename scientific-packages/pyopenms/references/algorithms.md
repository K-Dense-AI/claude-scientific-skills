# pyOpenMS Algorithms Reference

This document provides comprehensive coverage of algorithms available in pyOpenMS for signal processing, feature detection, and quantification.

## Algorithm Usage Pattern

Most pyOpenMS algorithms follow a consistent pattern:

```python
import pyopenms as oms

# 1. Instantiate algorithm
algorithm = oms.AlgorithmName()

# 2. Get parameters
params = algorithm.getParameters()

# 3. Modify parameters
params.setValue("parameter_name", value)

# 4. Set parameters back
algorithm.setParameters(params)

# 5. Apply to data
algorithm.filterExperiment(exp)  # or .process(), .run(), etc.
```

## Signal Processing Algorithms

### Smoothing Filters

#### GaussFilter - Gaussian Smoothing

Applies Gaussian smoothing to reduce noise.

```python
gauss = oms.GaussFilter()

# Configure parameters
params = gauss.getParameters()
params.setValue("gaussian_width", 0.2)  # Gaussian width (larger = more smoothing)
params.setValue("ppm_tolerance", 10.0)  # PPM tolerance for spacing
params.setValue("use_ppm_tolerance", "true")
gauss.setParameters(params)

# Apply to experiment
gauss.filterExperiment(exp)

# Or apply to single spectrum
spectrum_smoothed = oms.MSSpectrum()
gauss.filter(spectrum, spectrum_smoothed)
```

**Key Parameters:**
- `gaussian_width`: Width of Gaussian kernel (default: 0.2 Da)
- `ppm_tolerance`: Tolerance in ppm for spacing
- `use_ppm_tolerance`: Whether to use ppm instead of absolute spacing

#### SavitzkyGolayFilter

Applies Savitzky-Golay smoothing (polynomial fitting).

```python
sg_filter = oms.SavitzkyGolayFilter()

params = sg_filter.getParameters()
params.setValue("frame_length", 11)  # Window size (must be odd)
params.setValue("polynomial_order", 3)  # Polynomial degree
sg_filter.setParameters(params)

sg_filter.filterExperiment(exp)
```

**Key Parameters:**
- `frame_length`: Size of smoothing window (must be odd)
- `polynomial_order`: Degree of polynomial (typically 2-4)

### Peak Filtering

#### NLargest - Keep Top N Peaks

Retains only the N most intense peaks per spectrum.

```python
n_largest = oms.NLargest()

params = n_largest.getParameters()
params.setValue("n", 100)  # Keep top 100 peaks
params.setValue("threshold", 0.0)  # Optional minimum intensity
n_largest.setParameters(params)

n_largest.filterExperiment(exp)
```

**Key Parameters:**
- `n`: Number of peaks to keep per spectrum
- `threshold`: Minimum absolute intensity threshold

#### ThresholdMower - Intensity Threshold Filtering

Removes peaks below a specified intensity threshold.

```python
threshold_filter = oms.ThresholdMower()

params = threshold_filter.getParameters()
params.setValue("threshold", 1000.0)  # Absolute intensity threshold
threshold_filter.setParameters(params)

threshold_filter.filterExperiment(exp)
```

**Key Parameters:**
- `threshold`: Absolute intensity cutoff

#### WindowMower - Window-Based Peak Selection

Divides m/z range into windows and keeps top N peaks per window.

```python
window_mower = oms.WindowMower()

params = window_mower.getParameters()
params.setValue("windowsize", 50.0)  # Window size in Da (or Thomson)
params.setValue("peakcount", 10)     # Peaks to keep per window
params.setValue("movetype", "jump")  # "jump" or "slide"
window_mower.setParameters(params)

window_mower.filterExperiment(exp)
```

**Key Parameters:**
- `windowsize`: Size of m/z window (Da)
- `peakcount`: Number of peaks to retain per window
- `movetype`: "jump" (non-overlapping) or "slide" (overlapping windows)

#### BernNorm - Bernoulli Normalization

Statistical normalization based on Bernoulli distribution.

```python
bern_norm = oms.BernNorm()

params = bern_norm.getParameters()
params.setValue("threshold", 0.7)  # Threshold for normalization
bern_norm.setParameters(params)

bern_norm.filterExperiment(exp)
```

### Spectrum Normalization

#### Normalizer

Normalizes spectrum intensities to unit total intensity or maximum intensity.

```python
normalizer = oms.Normalizer()

params = normalizer.getParameters()
params.setValue("method", "to_one")  # "to_one" or "to_TIC"
normalizer.setParameters(params)

normalizer.filterExperiment(exp)
```

**Methods:**
- `to_one`: Normalize max peak to 1.0
- `to_TIC`: Normalize to total ion current = 1.0

#### Scaler

Scales intensities by a constant factor.

```python
scaler = oms.Scaler()

params = scaler.getParameters()
params.setValue("scaling", 1000.0)  # Scaling factor
scaler.setParameters(params)

scaler.filterExperiment(exp)
```

## Centroiding and Peak Picking

### PeakPickerHiRes - High-Resolution Peak Picking

Converts profile spectra to centroid mode for high-resolution data.

```python
picker = oms.PeakPickerHiRes()

params = picker.getParameters()
params.setValue("signal_to_noise", 1.0)      # S/N threshold
params.setValue("spacing_difference", 1.5)    # Peak spacing factor
params.setValue("sn_win_len", 20.0)          # S/N window length
params.setValue("sn_bin_count", 30)          # Bins for S/N estimation
params.setValue("ms1_only", "false")         # Process only MS1
params.setValue("ms_levels", [1, 2])         # MS levels to process
picker.setParameters(params)

# Pick peaks
exp_centroided = oms.MSExperiment()
picker.pickExperiment(exp, exp_centroided)
```

**Key Parameters:**
- `signal_to_noise`: Minimum signal-to-noise ratio
- `spacing_difference`: Minimum spacing between peaks
- `ms_levels`: List of MS levels to process

### PeakPickerWavelet - Wavelet-Based Peak Picking

Uses continuous wavelet transform for peak detection.

```python
wavelet_picker = oms.PeakPickerWavelet()

params = wavelet_picker.getParameters()
params.setValue("signal_to_noise", 1.0)
params.setValue("peak_width", 0.15)  # Expected peak width (Da)
wavelet_picker.setParameters(params)

wavelet_picker.pickExperiment(exp, exp_centroided)
```

## Feature Detection

### FeatureFinder Algorithms

Feature finders detect 2D features (m/z and RT) in LC-MS data.

#### FeatureFinderMultiplex

For multiplex labeling experiments (SILAC, dimethyl labeling).

```python
ff = oms.FeatureFinderMultiplex()

params = ff.getParameters()
params.setValue("algorithm:labels", "[]")  # Empty for label-free
params.setValue("algorithm:charge", "2:4")  # Charge range
params.setValue("algorithm:rt_typical", 40.0)  # Expected feature RT width
params.setValue("algorithm:rt_min", 2.0)  # Minimum RT width
params.setValue("algorithm:mz_tolerance", 10.0)  # m/z tolerance (ppm)
params.setValue("algorithm:intensity_cutoff", 1000.0)  # Minimum intensity
ff.setParameters(params)

# Run feature detection
features = oms.FeatureMap()
ff.run(exp, features, oms.Param())

print(f"Found {features.size()} features")
```

**Key Parameters:**
- `algorithm:charge`: Charge state range to consider
- `algorithm:rt_typical`: Expected peak width in RT dimension
- `algorithm:mz_tolerance`: Mass tolerance in ppm
- `algorithm:intensity_cutoff`: Minimum intensity threshold

#### FeatureFinderCentroided

For centroided data, identifies isotope patterns and traces over RT.

```python
ff_centroided = oms.FeatureFinderCentroided()

params = ff_centroided.getParameters()
params.setValue("mass_trace:mz_tolerance", 10.0)  # ppm
params.setValue("mass_trace:min_spectra", 5)  # Min consecutive spectra
params.setValue("isotopic_pattern:charge_low", 1)
params.setValue("isotopic_pattern:charge_high", 4)
params.setValue("seed:min_score", 0.5)
ff_centroided.setParameters(params)

features = oms.FeatureMap()
seeds = oms.FeatureMap()  # Optional seed features
ff_centroided.run(exp, features, params, seeds)
```

#### FeatureFinderIdentification

Uses peptide identifications to guide feature detection.

```python
ff_id = oms.FeatureFinderIdentification()

params = ff_id.getParameters()
params.setValue("extract:mz_window", 10.0)  # ppm
params.setValue("extract:rt_window", 60.0)  # seconds
params.setValue("detect:peak_width", 30.0)  # Expected peak width
ff_id.setParameters(params)

# Requires peptide identifications
protein_ids = []
peptide_ids = []
features = oms.FeatureMap()

ff_id.run(exp, protein_ids, peptide_ids, features)
```

## Charge and Isotope Deconvolution

### Decharging and Charge State Deconvolution

#### FeatureDeconvolution

Resolves charge states and combines features.

```python
deconv = oms.FeatureDeconvolution()

params = deconv.getParameters()
params.setValue("charge_min", 1)
params.setValue("charge_max", 4)
params.setValue("q_value", 0.01)  # FDR threshold
deconv.setParameters(params)

features_deconv = oms.FeatureMap()
consensus_map = oms.ConsensusMap()
deconv.compute(features, features_deconv, consensus_map)
```

## Map Alignment

### MapAlignmentAlgorithm

Aligns retention times across multiple LC-MS runs.

#### MapAlignmentAlgorithmPoseClustering

Pose clustering-based RT alignment.

```python
aligner = oms.MapAlignmentAlgorithmPoseClustering()

params = aligner.getParameters()
params.setValue("max_num_peaks_considered", 1000)
params.setValue("pairfinder:distance_MZ:max_difference", 0.3)  # Da
params.setValue("pairfinder:distance_RT:max_difference", 60.0)  # seconds
aligner.setParameters(params)

# Align multiple feature maps
feature_maps = [features1, features2, features3]
transformations = []

# Create reference (e.g., use first map)
reference = oms.FeatureMap(feature_maps[0])

# Align others to reference
for fm in feature_maps[1:]:
    transformation = oms.TransformationDescription()
    aligner.align(fm, reference, transformation)
    transformations.append(transformation)

    # Apply transformation
    transformer = oms.MapAlignmentTransformer()
    transformer.transformRetentionTimes(fm, transformation)
```

## Feature Linking

### FeatureGroupingAlgorithm

Links features across samples to create consensus features.

#### FeatureGroupingAlgorithmQT

Quality threshold-based feature linking.

```python
grouper = oms.FeatureGroupingAlgorithmQT()

params = grouper.getParameters()
params.setValue("distance_RT:max_difference", 60.0)  # seconds
params.setValue("distance_MZ:max_difference", 10.0)  # ppm
params.setValue("distance_MZ:unit", "ppm")
grouper.setParameters(params)

# Create consensus map
consensus_map = oms.ConsensusMap()

# Group features from multiple samples
feature_maps = [features1, features2, features3]
grouper.group(feature_maps, consensus_map)

print(f"Created {consensus_map.size()} consensus features")
```

#### FeatureGroupingAlgorithmKD

KD-tree based linking (faster for large datasets).

```python
grouper_kd = oms.FeatureGroupingAlgorithmKD()

params = grouper_kd.getParameters()
params.setValue("mz_unit", "ppm")
params.setValue("mz_tolerance", 10.0)
params.setValue("rt_tolerance", 30.0)
grouper_kd.setParameters(params)

consensus_map = oms.ConsensusMap()
grouper_kd.group(feature_maps, consensus_map)
```

## Chromatographic Analysis

### ElutionPeakDetection

Detects elution peaks in chromatograms.

```python
epd = oms.ElutionPeakDetection()

params = epd.getParameters()
params.setValue("chrom_peak_snr", 3.0)  # Signal-to-noise threshold
params.setValue("chrom_fwhm", 5.0)  # Expected FWHM (seconds)
epd.setParameters(params)

# Apply to chromatograms
for chrom in exp.getChromatograms():
    peaks = epd.detectPeaks(chrom)
```

### MRMFeatureFinderScoring

Scoring and peak picking for targeted (MRM/SRM) experiments.

```python
mrm_finder = oms.MRMFeatureFinderScoring()

params = mrm_finder.getParameters()
params.setValue("TransitionGroupPicker:min_peak_width", 2.0)
params.setValue("TransitionGroupPicker:recalculate_peaks", "true")
params.setValue("TransitionGroupPicker:PeakPickerMRM:signal_to_noise", 1.0)
mrm_finder.setParameters(params)

# Requires chromatograms
features = oms.FeatureMap()
mrm_finder.pickExperiment(chrom_exp, features, targets, transformation, swath_maps)
```

## Quantification

### ProteinInference

Infers proteins from peptide identifications.

```python
protein_inference = oms.BasicProteinInferenceAlgorithm()

# Apply to identification results
protein_inference.run(peptide_ids, protein_ids)
```

### IsobaricQuantification

Quantification for isobaric labeling (TMT, iTRAQ).

```python
# For TMT/iTRAQ quantification
iso_quant = oms.IsobaricQuantification()

params = iso_quant.getParameters()
params.setValue("channel_116_description", "Sample1")
params.setValue("channel_117_description", "Sample2")
# ... configure all channels
iso_quant.setParameters(params)

# Run quantification
quant_method = oms.IsobaricQuantitationMethod.TMT_10PLEX
quant_info = oms.IsobaricQuantifierStatistics()
iso_quant.quantify(exp, quant_info)
```

## Data Processing

### BaselineFiltering

Removes baseline from spectra.

```python
baseline = oms.TopHatFilter()

params = baseline.getParameters()
params.setValue("struc_elem_length", 3.0)  # Structuring element size
params.setValue("struc_elem_unit", "Thomson")
baseline.setParameters(params)

baseline.filterExperiment(exp)
```

### SpectraMerger

Merges consecutive similar spectra.

```python
merger = oms.SpectraMerger()

params = merger.getParameters()
params.setValue("mz_binning_width", 0.05)  # Binning width (Da)
params.setValue("sort_blocks", "RT_ascending")
merger.setParameters(params)

merger.mergeSpectra(exp)
```

## Quality Control

### MzMLFileQuality

Analyzes mzML file quality.

```python
# Calculate basic QC metrics
def calculate_qc_metrics(exp):
    metrics = {
        'n_spectra': exp.getNrSpectra(),
        'n_ms1': sum(1 for s in exp if s.getMSLevel() == 1),
        'n_ms2': sum(1 for s in exp if s.getMSLevel() == 2),
        'rt_range': (exp.getMinRT(), exp.getMaxRT()),
        'mz_range': (exp.getMinMZ(), exp.getMaxMZ()),
    }

    # Calculate TIC
    tics = []
    for spectrum in exp:
        if spectrum.getMSLevel() == 1:
            mz, intensity = spectrum.get_peaks()
            tics.append(sum(intensity))

    metrics['median_tic'] = np.median(tics)
    metrics['mean_tic'] = np.mean(tics)

    return metrics
```

## FDR Control

### FalseDiscoveryRate

Estimates and controls false discovery rate.

```python
fdr = oms.FalseDiscoveryRate()

params = fdr.getParameters()
params.setValue("add_decoy_peptides", "false")
params.setValue("add_decoy_proteins", "false")
fdr.setParameters(params)

# Apply to identifications
fdr.apply(protein_ids, peptide_ids)

# Filter by FDR threshold
fdr_threshold = 0.01
filtered_peptides = [p for p in peptide_ids if p.getMetaValue("q-value") <= fdr_threshold]
```

## Algorithm Selection Guide

### When to Use Which Algorithm

**For Smoothing:**
- Use `GaussFilter` for general-purpose smoothing
- Use `SavitzkyGolayFilter` for preserving peak shapes

**For Peak Picking:**
- Use `PeakPickerHiRes` for high-resolution Orbitrap/FT-ICR data
- Use `PeakPickerWavelet` for lower-resolution TOF data

**For Feature Detection:**
- Use `FeatureFinderCentroided` for label-free proteomics (DDA)
- Use `FeatureFinderMultiplex` for SILAC/dimethyl labeling
- Use `FeatureFinderIdentification` when you have ID information
- Use `MRMFeatureFinderScoring` for targeted (MRM/SRM) experiments

**For Feature Linking:**
- Use `FeatureGroupingAlgorithmQT` for small-medium datasets (<10 samples)
- Use `FeatureGroupingAlgorithmKD` for large datasets (>10 samples)

## Parameter Tuning Tips

1. **S/N Thresholds**: Start with 1-3 for clean data, increase for noisy data
2. **m/z Tolerance**: Use 5-10 ppm for high-resolution instruments, 0.5-1 Da for low-res
3. **RT Tolerance**: Typically 30-60 seconds depending on chromatographic stability
4. **Peak Width**: Measure from real data - varies by instrument and gradient length
5. **Charge States**: Set based on expected analytes (1-2 for metabolites, 2-4 for peptides)

## Common Algorithm Workflows

### Complete Proteomics Workflow

```python
# 1. Load data
exp = oms.MSExperiment()
oms.MzMLFile().load("raw.mzML", exp)

# 2. Smooth
gauss = oms.GaussFilter()
gauss.filterExperiment(exp)

# 3. Peak picking
picker = oms.PeakPickerHiRes()
exp_centroid = oms.MSExperiment()
picker.pickExperiment(exp, exp_centroid)

# 4. Feature detection
ff = oms.FeatureFinderCentroided()
features = oms.FeatureMap()
ff.run(exp_centroid, features, oms.Param(), oms.FeatureMap())

# 5. Save results
oms.FeatureXMLFile().store("features.featureXML", features)
```

### Multi-Sample Quantification

```python
# Load multiple samples
feature_maps = []
for filename in ["sample1.mzML", "sample2.mzML", "sample3.mzML"]:
    exp = oms.MSExperiment()
    oms.MzMLFile().load(filename, exp)

    # Process and detect features
    features = detect_features(exp)  # Your processing function
    feature_maps.append(features)

# Align retention times
align_feature_maps(feature_maps)  # Implement alignment

# Link features
grouper = oms.FeatureGroupingAlgorithmQT()
consensus_map = oms.ConsensusMap()
grouper.group(feature_maps, consensus_map)

# Export quantification matrix
export_quant_matrix(consensus_map)
```
