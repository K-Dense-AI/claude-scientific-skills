---
name: pyopenms
description: "Mass spectrometry toolkit (OpenMS Python). Process mzML/mzXML, peak picking, feature detection, peptide ID, proteomics/metabolomics workflows, for LC-MS/MS analysis."
---

# pyOpenMS

## Overview

pyOpenMS is an open-source Python library for mass spectrometry data analysis in proteomics and metabolomics. Process LC-MS/MS data, perform peptide identification, detect and quantify features, and integrate with common proteomics tools (Comet, Mascot, MSGF+, Percolator, MSstats) using Python bindings to the OpenMS C++ library.

## When to Use This Skill

This skill should be used when:
- Processing mass spectrometry data (mzML, mzXML files)
- Performing peak picking and feature detection in LC-MS data
- Conducting peptide and protein identification workflows
- Quantifying metabolites or proteins
- Integrating proteomics or metabolomics tools into Python pipelines
- Working with OpenMS tools and file formats

## Core Capabilities

### 1. File I/O and Data Import/Export

Handle diverse mass spectrometry file formats efficiently:

**Supported Formats:**
- **mzML/mzXML**: Primary raw MS data formats (profile or centroid)
- **FASTA**: Protein/peptide sequence databases
- **mzTab**: Standardized reporting format for identification and quantification
- **mzIdentML**: Peptide and protein identification data
- **TraML**: Transition lists for targeted experiments
- **pepXML/protXML**: Search engine results

**Reading mzML Files:**
```python
import pyopenms as oms

# Load MS data
exp = oms.MSExperiment()
oms.MzMLFile().load("input_data.mzML", exp)

# Access basic information
print(f"Number of spectra: {exp.getNrSpectra()}")
print(f"Number of chromatograms: {exp.getNrChromatograms()}")
```

**Writing mzML Files:**
```python
# Save processed data
oms.MzMLFile().store("output_data.mzML", exp)
```

**File Encoding:** pyOpenMS automatically handles Base64 encoding, zlib compression, and Numpress compression internally.

### 2. MS Data Structures and Manipulation

Work with core mass spectrometry data structures. See `references/data_structures.md` for comprehensive details.

**MSSpectrum** - Individual mass spectrum:
```python
# Create spectrum with metadata
spectrum = oms.MSSpectrum()
spectrum.setRT(205.2)  # Retention time in seconds
spectrum.setMSLevel(2)  # MS2 spectrum

# Set peak data (m/z, intensity arrays)
mz_array = [100.5, 200.3, 300.7, 400.2]
intensity_array = [1000, 5000, 3000, 2000]
spectrum.set_peaks((mz_array, intensity_array))

# Add precursor information for MS2
precursor = oms.Precursor()
precursor.setMZ(450.5)
precursor.setCharge(2)
spectrum.setPrecursors([precursor])
```

**MSExperiment** - Complete LC-MS/MS run:
```python
# Create experiment and add spectra
exp = oms.MSExperiment()
exp.addSpectrum(spectrum)

# Access spectra
first_spectrum = exp.getSpectrum(0)
for spec in exp:
    print(f"RT: {spec.getRT()}, MS Level: {spec.getMSLevel()}")
```

**MSChromatogram** - Extracted ion chromatogram:
```python
# Create chromatogram
chrom = oms.MSChromatogram()
chrom.set_peaks(([10.5, 11.2, 11.8], [1000, 5000, 3000]))  # RT, intensity
exp.addChromatogram(chrom)
```

**Efficient Peak Access:**
```python
# Get peaks as numpy arrays for fast processing
mz_array, intensity_array = spectrum.get_peaks()

# Modify and set back
intensity_array *= 2  # Double all intensities
spectrum.set_peaks((mz_array, intensity_array))
```

### 3. Chemistry and Peptide Handling

Perform chemical calculations for proteomics and metabolomics. See `references/chemistry.md` for detailed examples.

**Molecular Formulas and Mass Calculations:**
```python
# Create empirical formula
formula = oms.EmpiricalFormula("C6H12O6")  # Glucose
print(f"Monoisotopic mass: {formula.getMonoWeight()}")
print(f"Average mass: {formula.getAverageWeight()}")

# Formula arithmetic
water = oms.EmpiricalFormula("H2O")
dehydrated = formula - water

# Isotope-specific formulas
heavy_carbon = oms.EmpiricalFormula("(13)C6H12O6")
```

**Isotopic Distributions:**
```python
# Generate coarse isotope pattern (unit mass resolution)
coarse_gen = oms.CoarseIsotopePatternGenerator()
pattern = coarse_gen.run(formula)

# Generate fine structure (high resolution)
fine_gen = oms.FineIsotopePatternGenerator(0.01)  # 0.01 Da resolution
fine_pattern = fine_gen.run(formula)
```

**Amino Acids and Residues:**
```python
# Access residue information
res_db = oms.ResidueDB()
leucine = res_db.getResidue("Leucine")
print(f"L monoisotopic mass: {leucine.getMonoWeight()}")
print(f"L formula: {leucine.getFormula()}")
print(f"L pKa: {leucine.getPka()}")
```

**Peptide Sequences:**
```python
# Create peptide sequence
peptide = oms.AASequence.fromString("PEPTIDE")
print(f"Peptide mass: {peptide.getMonoWeight()}")
print(f"Formula: {peptide.getFormula()}")

# Add modifications
modified = oms.AASequence.fromString("PEPTIDEM(Oxidation)")
print(f"Modified mass: {modified.getMonoWeight()}")

# Theoretical fragmentation
ions = []
for i in range(1, peptide.size()):
    b_ion = peptide.getPrefix(i)
    y_ion = peptide.getSuffix(i)
    ions.append(('b', i, b_ion.getMonoWeight()))
    ions.append(('y', i, y_ion.getMonoWeight()))
```

**Protein Digestion:**
```python
# Enzymatic digestion
dig = oms.ProteaseDigestion()
dig.setEnzyme("Trypsin")
dig.setMissedCleavages(2)

protein_seq = oms.AASequence.fromString("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK")
peptides = []
dig.digest(protein_seq, peptides)

for pep in peptides:
    print(f"{pep.toString()}: {pep.getMonoWeight():.2f} Da")
```

**Modifications:**
```python
# Access modification database
mod_db = oms.ModificationsDB()
oxidation = mod_db.getModification("Oxidation")
print(f"Oxidation mass diff: {oxidation.getDiffMonoMass()}")
print(f"Residues: {oxidation.getResidues()}")
```

### 4. Signal Processing and Filtering

Apply algorithms to process and filter MS data. See `references/algorithms.md` for comprehensive coverage.

**Spectral Smoothing:**
```python
# Gaussian smoothing
gauss_filter = oms.GaussFilter()
params = gauss_filter.getParameters()
params.setValue("gaussian_width", 0.2)
gauss_filter.setParameters(params)
gauss_filter.filterExperiment(exp)

# Savitzky-Golay filter
sg_filter = oms.SavitzkyGolayFilter()
sg_filter.filterExperiment(exp)
```

**Peak Filtering:**
```python
# Keep only N largest peaks per spectrum
n_largest = oms.NLargest()
params = n_largest.getParameters()
params.setValue("n", 100)  # Keep top 100 peaks
n_largest.setParameters(params)
n_largest.filterExperiment(exp)

# Threshold filtering
threshold_filter = oms.ThresholdMower()
params = threshold_filter.getParameters()
params.setValue("threshold", 1000.0)  # Remove peaks below 1000 intensity
threshold_filter.setParameters(params)
threshold_filter.filterExperiment(exp)

# Window-based filtering
window_filter = oms.WindowMower()
params = window_filter.getParameters()
params.setValue("windowsize", 50.0)  # 50 m/z windows
params.setValue("peakcount", 10)     # Keep 10 highest per window
window_filter.setParameters(params)
window_filter.filterExperiment(exp)
```

**Spectrum Normalization:**
```python
normalizer = oms.Normalizer()
normalizer.filterExperiment(exp)
```

**MS Level Filtering:**
```python
# Keep only MS2 spectra
exp.filterMSLevel(2)

# Filter by retention time range
exp.filterRT(100.0, 500.0)  # Keep RT between 100-500 seconds

# Filter by m/z range
exp.filterMZ(400.0, 1500.0)  # Keep m/z between 400-1500
```

### 5. Feature Detection and Quantification

Detect and quantify features in LC-MS data:

**Peak Picking (Centroiding):**
```python
# Convert profile data to centroid
picker = oms.PeakPickerHiRes()
params = picker.getParameters()
params.setValue("signal_to_noise", 1.0)
picker.setParameters(params)

exp_centroided = oms.MSExperiment()
picker.pickExperiment(exp, exp_centroided)
```

**Feature Detection:**
```python
# Detect features across LC-MS runs
feature_finder = oms.FeatureFinderMultiplex()

features = oms.FeatureMap()
feature_finder.run(exp, features, params)

print(f"Found {features.size()} features")
for feature in features:
    print(f"m/z: {feature.getMZ():.4f}, RT: {feature.getRT():.2f}, "
          f"Intensity: {feature.getIntensity():.0f}")
```

**Feature Linking (Map Alignment):**
```python
# Link features across multiple samples
feature_grouper = oms.FeatureGroupingAlgorithmQT()
consensus_map = oms.ConsensusMap()

# Provide multiple feature maps from different samples
feature_maps = [features1, features2, features3]
feature_grouper.group(feature_maps, consensus_map)
```

### 6. Peptide Identification Workflows

Integrate with search engines and process identification results:

**Database Searching:**
```python
# Prepare parameters for search engine
params = oms.Param()
params.setValue("database", "uniprot_human.fasta")
params.setValue("precursor_mass_tolerance", 10.0)  # ppm
params.setValue("fragment_mass_tolerance", 0.5)     # Da
params.setValue("enzyme", "Trypsin")
params.setValue("missed_cleavages", 2)

# Variable modifications
params.setValue("variable_modifications", ["Oxidation (M)", "Phospho (STY)"])

# Fixed modifications
params.setValue("fixed_modifications", ["Carbamidomethyl (C)"])
```

**FDR Control:**
```python
# False discovery rate estimation
fdr = oms.FalseDiscoveryRate()
fdr_threshold = 0.01  # 1% FDR

# Apply to peptide identifications
protein_ids = []
peptide_ids = []
oms.IdXMLFile().load("search_results.idXML", protein_ids, peptide_ids)

fdr.apply(protein_ids, peptide_ids)
```

### 7. Metabolomics Workflows

Analyze small molecule data:

**Adduct Detection:**
```python
# Common metabolite adducts
adducts = ["[M+H]+", "[M+Na]+", "[M+K]+", "[M-H]-", "[M+Cl]-"]

# Feature annotation with adducts
for feature in features:
    mz = feature.getMZ()
    # Calculate neutral mass for each adduct hypothesis
    for adduct in adducts:
        # Annotation logic
        pass
```

**Isotope Pattern Matching:**
```python
# Compare experimental to theoretical isotope patterns
experimental_pattern = []  # Extract from feature
theoretical = coarse_gen.run(formula)

# Calculate similarity score
similarity = compare_isotope_patterns(experimental_pattern, theoretical)
```

### 8. Quality Control and Visualization

Monitor data quality and visualize results:

**Basic Statistics:**
```python
# Calculate TIC (Total Ion Current)
tic_values = []
rt_values = []
for spectrum in exp:
    if spectrum.getMSLevel() == 1:
        tic = sum(spectrum.get_peaks()[1])  # Sum intensities
        tic_values.append(tic)
        rt_values.append(spectrum.getRT())

# Base peak chromatogram
bpc_values = []
for spectrum in exp:
    if spectrum.getMSLevel() == 1:
        max_intensity = max(spectrum.get_peaks()[1]) if spectrum.size() > 0 else 0
        bpc_values.append(max_intensity)
```

**Plotting (with pyopenms.plotting or matplotlib):**
```python
import matplotlib.pyplot as plt

# Plot TIC
plt.figure(figsize=(10, 4))
plt.plot(rt_values, tic_values)
plt.xlabel('Retention Time (s)')
plt.ylabel('Total Ion Current')
plt.title('TIC')
plt.show()

# Plot single spectrum
spectrum = exp.getSpectrum(0)
mz, intensity = spectrum.get_peaks()
plt.stem(mz, intensity, basefmt=' ')
plt.xlabel('m/z')
plt.ylabel('Intensity')
plt.title(f'Spectrum at RT {spectrum.getRT():.2f}s')
plt.show()
```

## Common Workflows

### Complete LC-MS/MS Processing Pipeline

```python
import pyopenms as oms

# 1. Load data
exp = oms.MSExperiment()
oms.MzMLFile().load("raw_data.mzML", exp)

# 2. Filter and smooth
exp.filterMSLevel(1)  # Keep only MS1 for feature detection
gauss = oms.GaussFilter()
gauss.filterExperiment(exp)

# 3. Peak picking
picker = oms.PeakPickerHiRes()
exp_centroid = oms.MSExperiment()
picker.pickExperiment(exp, exp_centroid)

# 4. Feature detection
ff = oms.FeatureFinderMultiplex()
features = oms.FeatureMap()
ff.run(exp_centroid, features, oms.Param())

# 5. Export results
oms.FeatureXMLFile().store("features.featureXML", features)
print(f"Detected {features.size()} features")
```

### Theoretical Peptide Mass Calculation

```python
# Calculate masses for peptide with modifications
peptide = oms.AASequence.fromString("PEPTIDEK")
print(f"Unmodified [M+H]+: {peptide.getMonoWeight() + 1.007276:.4f}")

# With modification
modified = oms.AASequence.fromString("PEPTIDEM(Oxidation)K")
print(f"Oxidized [M+H]+: {modified.getMonoWeight() + 1.007276:.4f}")

# Calculate for different charge states
for z in [1, 2, 3]:
    mz = (peptide.getMonoWeight() + z * 1.007276) / z
    print(f"[M+{z}H]^{z}+: {mz:.4f}")
```

## Installation

Ensure pyOpenMS is installed before using this skill:

```bash
# Via conda (recommended)
conda install -c bioconda pyopenms

# Via pip
pip install pyopenms
```

## Integration with Other Tools

pyOpenMS integrates seamlessly with:

- **Search Engines**: Comet, Mascot, MSGF+, MSFragger, Sage, SpectraST
- **Post-processing**: Percolator, MSstats, Epiphany
- **Metabolomics**: SIRIUS, CSI:FingerID
- **Data Analysis**: Pandas, NumPy, SciPy for downstream analysis
- **Visualization**: Matplotlib, Seaborn for plotting

## Resources

### references/

Detailed documentation on core concepts:

- **data_structures.md** - Comprehensive guide to MSExperiment, MSSpectrum, MSChromatogram, and peak data handling
- **algorithms.md** - Complete reference for signal processing, filtering, feature detection, and quantification algorithms
- **chemistry.md** - In-depth coverage of chemistry calculations, peptide handling, modifications, and isotope distributions

Load these references when needing detailed information about specific pyOpenMS capabilities.

## Best Practices

1. **File Format**: Always use mzML for raw MS data (standardized, well-supported)
2. **Peak Access**: Use `get_peaks()` and `set_peaks()` with numpy arrays for efficient processing
3. **Parameters**: Always check and configure algorithm parameters via `getParameters()` and `setParameters()`
4. **Memory**: For large datasets, process spectra iteratively rather than loading entire experiments
5. **Validation**: Check data integrity (MS levels, RT ordering, precursor information) after loading
6. **Modifications**: Use standard modification names from UniMod database
7. **Units**: RT in seconds, m/z in Thomson (Da/charge), intensity in arbitrary units

## Common Patterns

**Algorithm Application Pattern:**
```python
# 1. Instantiate algorithm
algorithm = oms.SomeAlgorithm()

# 2. Get and configure parameters
params = algorithm.getParameters()
params.setValue("parameter_name", value)
algorithm.setParameters(params)

# 3. Apply to data
algorithm.filterExperiment(exp)  # or .process(), .run(), depending on algorithm
```

**File I/O Pattern:**
```python
# Read
data_container = oms.DataContainer()  # MSExperiment, FeatureMap, etc.
oms.FileHandler().load("input.format", data_container)

# Process
# ... manipulate data_container ...

# Write
oms.FileHandler().store("output.format", data_container)
```

## Getting Help

- **Documentation**: https://pyopenms.readthedocs.io/
- **API Reference**: Browse class documentation for detailed method signatures
- **OpenMS Website**: https://www.openms.org/
- **GitHub Issues**: https://github.com/OpenMS/OpenMS/issues
