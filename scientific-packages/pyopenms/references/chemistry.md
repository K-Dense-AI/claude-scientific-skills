# pyOpenMS Chemistry Reference

This document provides comprehensive coverage of chemistry-related functionality in pyOpenMS, including elements, isotopes, molecular formulas, amino acids, peptides, proteins, and modifications.

## Elements and Isotopes

### ElementDB - Element Database

Access atomic and isotopic data for all elements.

```python
import pyopenms as oms

# Get element database instance
element_db = oms.ElementDB()

# Get element by symbol
carbon = element_db.getElement("C")
nitrogen = element_db.getElement("N")
oxygen = element_db.getElement("O")

# Element properties
print(f"Carbon monoisotopic weight: {carbon.getMonoWeight()}")
print(f"Carbon average weight: {carbon.getAverageWeight()}")
print(f"Atomic number: {carbon.getAtomicNumber()}")
print(f"Symbol: {carbon.getSymbol()}")
print(f"Name: {carbon.getName()}")
```

### Isotope Information

```python
# Get isotope distribution for an element
isotopes = carbon.getIsotopeDistribution()

# Access specific isotope
c12 = element_db.getElement("C", 12)  # Carbon-12
c13 = element_db.getElement("C", 13)  # Carbon-13

print(f"C-12 abundance: {isotopes.getContainer()[0].getIntensity()}")
print(f"C-13 abundance: {isotopes.getContainer()[1].getIntensity()}")

# Isotope mass
print(f"C-12 mass: {c12.getMonoWeight()}")
print(f"C-13 mass: {c13.getMonoWeight()}")
```

### Constants

```python
# Physical constants
avogadro = oms.Constants.AVOGADRO
electron_mass = oms.Constants.ELECTRON_MASS_U
proton_mass = oms.Constants.PROTON_MASS_U

print(f"Avogadro's number: {avogadro}")
print(f"Electron mass: {electron_mass} u")
print(f"Proton mass: {proton_mass} u")
```

## Empirical Formulas

### EmpiricalFormula - Molecular Formulas

Represent and manipulate molecular formulas.

#### Creating Formulas

```python
# From string
glucose = oms.EmpiricalFormula("C6H12O6")
water = oms.EmpiricalFormula("H2O")
ammonia = oms.EmpiricalFormula("NH3")

# From element composition
formula = oms.EmpiricalFormula()
formula.setCharge(1)  # Set charge state
```

#### Formula Arithmetic

```python
# Addition
sucrose = oms.EmpiricalFormula("C12H22O11")
hydrolyzed = sucrose + water  # Hydrolysis adds water

# Subtraction
dehydrated = glucose - water  # Dehydration removes water

# Multiplication
three_waters = water * 3  # 3 H2O = H6O3

# Division
formula_half = sucrose / 2  # Half the formula
```

#### Mass Calculations

```python
# Monoisotopic mass
mono_mass = glucose.getMonoWeight()
print(f"Glucose monoisotopic mass: {mono_mass:.6f} Da")

# Average mass
avg_mass = glucose.getAverageWeight()
print(f"Glucose average mass: {avg_mass:.6f} Da")

# Mass difference
mass_diff = (glucose - water).getMonoWeight()
```

#### Elemental Composition

```python
# Get element counts
formula = oms.EmpiricalFormula("C6H12O6")

# Access individual elements
n_carbon = formula.getNumberOf(element_db.getElement("C"))
n_hydrogen = formula.getNumberOf(element_db.getElement("H"))
n_oxygen = formula.getNumberOf(element_db.getElement("O"))

print(f"C: {n_carbon}, H: {n_hydrogen}, O: {n_oxygen}")

# String representation
print(f"Formula: {formula.toString()}")
```

#### Isotope-Specific Formulas

```python
# Specify specific isotopes using parentheses
labeled_glucose = oms.EmpiricalFormula("(13)C6H12O6")  # All carbons are C-13
partially_labeled = oms.EmpiricalFormula("C5(13)CH12O6")  # One C-13

# Deuterium labeling
deuterated = oms.EmpiricalFormula("C6D12O6")  # D2O instead of H2O
```

#### Charge States

```python
# Set charge
formula = oms.EmpiricalFormula("C6H12O6")
formula.setCharge(1)  # Positive charge

# Get charge
charge = formula.getCharge()

# Calculate m/z for charged molecule
mz = formula.getMonoWeight() / abs(charge) if charge != 0 else formula.getMonoWeight()
```

### Isotope Pattern Generation

Generate theoretical isotope patterns for formulas.

#### CoarseIsotopePatternGenerator

For unit mass resolution (low-resolution instruments).

```python
# Create generator
coarse_gen = oms.CoarseIsotopePatternGenerator()

# Generate pattern
formula = oms.EmpiricalFormula("C6H12O6")
pattern = coarse_gen.run(formula)

# Access isotope peaks
iso_dist = pattern.getContainer()
for peak in iso_dist:
    mass = peak.getMZ()
    abundance = peak.getIntensity()
    print(f"m/z: {mass:.4f}, Abundance: {abundance:.4f}")
```

#### FineIsotopePatternGenerator

For high-resolution instruments (hyperfine structure).

```python
# Create generator with resolution
fine_gen = oms.FineIsotopePatternGenerator(0.01)  # 0.01 Da resolution

# Generate fine pattern
fine_pattern = fine_gen.run(formula)

# Access fine isotope structure
for peak in fine_pattern.getContainer():
    print(f"m/z: {peak.getMZ():.6f}, Abundance: {peak.getIntensity():.6f}")
```

#### Isotope Pattern Matching

```python
# Compare experimental to theoretical
def compare_isotope_patterns(experimental_mz, experimental_int, formula):
    # Generate theoretical
    coarse_gen = oms.CoarseIsotopePatternGenerator()
    theoretical = coarse_gen.run(formula)

    # Extract theoretical peaks
    theo_peaks = theoretical.getContainer()
    theo_mz = [p.getMZ() for p in theo_peaks]
    theo_int = [p.getIntensity() for p in theo_peaks]

    # Normalize both patterns
    exp_int_norm = [i / max(experimental_int) for i in experimental_int]
    theo_int_norm = [i / max(theo_int) for i in theo_int]

    # Calculate similarity (e.g., cosine similarity)
    # ... implement similarity calculation
    return similarity_score
```

## Amino Acids and Residues

### Residue - Amino Acid Representation

Access properties of amino acids.

```python
# Get residue database
res_db = oms.ResidueDB()

# Get specific residue
leucine = res_db.getResidue("Leucine")
# Or by one-letter code
leu = res_db.getResidue("L")

# Residue properties
print(f"Name: {leucine.getName()}")
print(f"Three-letter code: {leucine.getThreeLetterCode()}")
print(f"One-letter code: {leucine.getOneLetterCode()}")
print(f"Monoisotopic mass: {leucine.getMonoWeight():.6f}")
print(f"Average mass: {leucine.getAverageWeight():.6f}")

# Chemical formula
formula = leucine.getFormula()
print(f"Formula: {formula.toString()}")

# pKa values
print(f"pKa (N-term): {leucine.getPka()}")
print(f"pKa (C-term): {leucine.getPkb()}")
print(f"pKa (side chain): {leucine.getPkc()}")

# Side chain basicity/acidity
print(f"Basicity: {leucine.getBasicity()}")
print(f"Hydrophobicity: {leucine.getHydrophobicity()}")
```

### All Standard Amino Acids

```python
# Iterate over all residues
for residue_name in ["Alanine", "Cysteine", "Aspartic acid", "Glutamic acid",
                     "Phenylalanine", "Glycine", "Histidine", "Isoleucine",
                     "Lysine", "Leucine", "Methionine", "Asparagine",
                     "Proline", "Glutamine", "Arginine", "Serine",
                     "Threonine", "Valine", "Tryptophan", "Tyrosine"]:
    res = res_db.getResidue(residue_name)
    print(f"{res.getOneLetterCode()}: {res.getMonoWeight():.4f} Da")
```

### Internal Residues vs. Termini

```python
# Get internal residue mass (no terminal groups)
internal_mass = leucine.getInternalToFull()

# Get residue with N-terminal modification
n_terminal = res_db.getResidue("L[1]")  # With NH2

# Get residue with C-terminal modification
c_terminal = res_db.getResidue("L[2]")  # With COOH
```

## Peptide Sequences

### AASequence - Amino Acid Sequences

Represent and manipulate peptide sequences.

#### Creating Sequences

```python
# From string
peptide = oms.AASequence.fromString("PEPTIDE")
longer = oms.AASequence.fromString("MKTAYIAKQRQISFVK")

# Empty sequence
empty_seq = oms.AASequence()
```

#### Sequence Properties

```python
peptide = oms.AASequence.fromString("PEPTIDE")

# Length
length = peptide.size()
print(f"Length: {length} residues")

# Mass
mono_mass = peptide.getMonoWeight()
avg_mass = peptide.getAverageWeight()
print(f"Monoisotopic mass: {mono_mass:.6f} Da")
print(f"Average mass: {avg_mass:.6f} Da")

# Formula
formula = peptide.getFormula()
print(f"Formula: {formula.toString()}")

# String representation
seq_str = peptide.toString()
print(f"Sequence: {seq_str}")
```

#### Accessing Individual Residues

```python
peptide = oms.AASequence.fromString("PEPTIDE")

# Access by index
first_aa = peptide[0]  # Returns Residue object
print(f"First amino acid: {first_aa.getOneLetterCode()}")

# Iterate
for i in range(peptide.size()):
    residue = peptide[i]
    print(f"Position {i}: {residue.getOneLetterCode()}")
```

#### Modifications

Add post-translational modifications (PTMs) to sequences.

```python
# Modifications in sequence string
# Format: AA(ModificationName)
oxidized_met = oms.AASequence.fromString("PEPTIDEM(Oxidation)")
phospho = oms.AASequence.fromString("PEPTIDES(Phospho)T(Phospho)")

# Multiple modifications
multi_mod = oms.AASequence.fromString("M(Oxidation)PEPTIDEK(Acetyl)")

# N-terminal modifications
n_term_acetyl = oms.AASequence.fromString("(Acetyl)PEPTIDE")

# C-terminal modifications
c_term_amide = oms.AASequence.fromString("PEPTIDE(Amidated)")

# Check mass change
unmodified = oms.AASequence.fromString("PEPTIDE")
modified = oms.AASequence.fromString("PEPTIDEM(Oxidation)")
mass_diff = modified.getMonoWeight() - unmodified.getMonoWeight()
print(f"Mass shift from oxidation: {mass_diff:.6f} Da")
```

#### Sequence Manipulation

```python
# Prefix (N-terminal fragment)
prefix = peptide.getPrefix(3)  # First 3 residues
print(f"Prefix: {prefix.toString()}")

# Suffix (C-terminal fragment)
suffix = peptide.getSuffix(3)  # Last 3 residues
print(f"Suffix: {suffix.toString()}")

# Subsequence
subseq = peptide.getSubsequence(2, 4)  # Residues 2-4
print(f"Subsequence: {subseq.toString()}")
```

#### Theoretical Fragmentation

Generate theoretical fragment ions for MS/MS.

```python
peptide = oms.AASequence.fromString("PEPTIDE")

# b-ions (N-terminal fragments)
b_ions = []
for i in range(1, peptide.size()):
    b_fragment = peptide.getPrefix(i)
    b_mass = b_fragment.getMonoWeight()
    b_ions.append(('b', i, b_mass))
    print(f"b{i}: {b_mass:.4f}")

# y-ions (C-terminal fragments)
y_ions = []
for i in range(1, peptide.size()):
    y_fragment = peptide.getSuffix(i)
    y_mass = y_fragment.getMonoWeight()
    y_ions.append(('y', i, y_mass))
    print(f"y{i}: {y_mass:.4f}")

# a-ions (b - CO)
a_ions = []
CO_mass = 27.994915  # CO loss
for ion_type, position, mass in b_ions:
    a_mass = mass - CO_mass
    a_ions.append(('a', position, a_mass))

# c-ions (b + NH3)
NH3_mass = 17.026549  # NH3 gain
c_ions = []
for ion_type, position, mass in b_ions:
    c_mass = mass + NH3_mass
    c_ions.append(('c', position, c_mass))

# z-ions (y - NH3)
z_ions = []
for ion_type, position, mass in y_ions:
    z_mass = mass - NH3_mass
    z_ions.append(('z', position, z_mass))
```

#### Calculate m/z for Charge States

```python
peptide = oms.AASequence.fromString("PEPTIDE")
proton_mass = 1.007276

# [M+H]+
mz_1 = peptide.getMonoWeight() + proton_mass
print(f"[M+H]+: {mz_1:.4f}")

# [M+2H]2+
mz_2 = (peptide.getMonoWeight() + 2 * proton_mass) / 2
print(f"[M+2H]2+: {mz_2:.4f}")

# [M+3H]3+
mz_3 = (peptide.getMonoWeight() + 3 * proton_mass) / 3
print(f"[M+3H]3+: {mz_3:.4f}")

# General formula for any charge
def calculate_mz(sequence, charge):
    proton_mass = 1.007276
    return (sequence.getMonoWeight() + charge * proton_mass) / charge

for z in range(1, 5):
    print(f"[M+{z}H]{z}+: {calculate_mz(peptide, z):.4f}")
```

## Protein Digestion

### ProteaseDigestion - Enzymatic Cleavage

Simulate enzymatic protein digestion.

#### Basic Digestion

```python
# Create digestion object
dig = oms.ProteaseDigestion()

# Set enzyme
dig.setEnzyme("Trypsin")  # Cleaves after K, R

# Other common enzymes:
# - "Trypsin" (K, R)
# - "Lys-C" (K)
# - "Arg-C" (R)
# - "Asp-N" (D)
# - "Glu-C" (E, D)
# - "Chymotrypsin" (F, Y, W, L)

# Set missed cleavages
dig.setMissedCleavages(0)  # No missed cleavages
dig.setMissedCleavages(2)  # Allow up to 2 missed cleavages

# Perform digestion
protein = oms.AASequence.fromString("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK")
peptides = []
dig.digest(protein, peptides)

# Print results
for pep in peptides:
    print(f"{pep.toString()}: {pep.getMonoWeight():.2f} Da")
```

#### Advanced Digestion Options

```python
# Get enzyme specificity
specificity = dig.getSpecificity()
# oms.EnzymaticDigestion.SPEC_FULL (both termini)
# oms.EnzymaticDigestion.SPEC_SEMI (one terminus)
# oms.EnzymaticDigestion.SPEC_NONE (no specificity)

# Set specificity for semi-tryptic search
dig.setSpecificity(oms.EnzymaticDigestion.SPEC_SEMI)

# Get cleavage sites
cleavage_residues = dig.getEnzyme().getCutAfterResidues()
restriction_residues = dig.getEnzyme().getRestriction()
```

#### Filter Peptides by Properties

```python
# Filter by mass range
min_mass = 600.0
max_mass = 4000.0
filtered = [p for p in peptides if min_mass <= p.getMonoWeight() <= max_mass]

# Filter by length
min_length = 6
max_length = 30
length_filtered = [p for p in peptides if min_length <= p.size() <= max_length]

# Combine filters
valid_peptides = [p for p in peptides
                  if min_mass <= p.getMonoWeight() <= max_mass
                  and min_length <= p.size() <= max_length]
```

## Modifications

### ModificationsDB - Modification Database

Access and apply post-translational modifications.

#### Accessing Modifications

```python
# Get modifications database
mod_db = oms.ModificationsDB()

# Get specific modification
oxidation = mod_db.getModification("Oxidation")
phospho = mod_db.getModification("Phospho")
acetyl = mod_db.getModification("Acetyl")

# Modification properties
print(f"Name: {oxidation.getFullName()}")
print(f"Mass difference: {oxidation.getDiffMonoMass():.6f} Da")
print(f"Formula: {oxidation.getDiffFormula().toString()}")

# Affected residues
print(f"Residues: {oxidation.getResidues()}")  # e.g., ['M']

# Specificity (N-term, C-term, anywhere)
print(f"Term specificity: {oxidation.getTermSpecificity()}")
```

#### Common Modifications

```python
# Oxidation (M)
oxidation = mod_db.getModification("Oxidation")
print(f"Oxidation: +{oxidation.getDiffMonoMass():.4f} Da")

# Phosphorylation (S, T, Y)
phospho = mod_db.getModification("Phospho")
print(f"Phospho: +{phospho.getDiffMonoMass():.4f} Da")

# Carbamidomethylation (C) - common alkylation
carbamido = mod_db.getModification("Carbamidomethyl")
print(f"Carbamidomethyl: +{carbamido.getDiffMonoMass():.4f} Da")

# Acetylation (K, N-term)
acetyl = mod_db.getModification("Acetyl")
print(f"Acetyl: +{acetyl.getDiffMonoMass():.4f} Da")

# Deamidation (N, Q)
deamid = mod_db.getModification("Deamidated")
print(f"Deamidation: +{deamid.getDiffMonoMass():.4f} Da")
```

#### Searching Modifications

```python
# Search modifications by mass
mass_tolerance = 0.01  # Da
target_mass = 15.9949  # Oxidation

# Get all modifications
all_mods = []
mod_db.getAllSearchModifications(all_mods)

# Find matching modifications
matching = []
for mod_name in all_mods:
    mod = mod_db.getModification(mod_name)
    if abs(mod.getDiffMonoMass() - target_mass) < mass_tolerance:
        matching.append(mod)
        print(f"Match: {mod.getFullName()} ({mod.getDiffMonoMass():.4f} Da)")
```

#### Variable vs. Fixed Modifications

```python
# In search engines, specify:
# Fixed modifications: applied to all occurrences
fixed_mods = ["Carbamidomethyl (C)"]

# Variable modifications: optionally present
variable_mods = ["Oxidation (M)", "Phospho (S)", "Phospho (T)", "Phospho (Y)"]
```

## Ribonucleotides (RNA)

### Ribonucleotide - RNA Building Blocks

```python
# Get ribonucleotide database
ribo_db = oms.RibonucleotideDB()

# Get specific ribonucleotide
adenine = ribo_db.getRibonucleotide("A")
uracil = ribo_db.getRibonucleotide("U")
guanine = ribo_db.getRibonucleotide("G")
cytosine = ribo_db.getRibonucleotide("C")

# Properties
print(f"Adenine mono mass: {adenine.getMonoWeight()}")
print(f"Formula: {adenine.getFormula().toString()}")

# Modified ribonucleotides
modified_ribo = ribo_db.getRibonucleotide("m6A")  # N6-methyladenosine
```

## Practical Examples

### Calculate Peptide Mass with Modifications

```python
def calculate_peptide_mz(sequence_str, charge):
    """Calculate m/z for a peptide sequence string with modifications."""
    peptide = oms.AASequence.fromString(sequence_str)
    proton_mass = 1.007276
    mz = (peptide.getMonoWeight() + charge * proton_mass) / charge
    return mz

# Examples
print(calculate_peptide_mz("PEPTIDE", 2))  # Unmodified [M+2H]2+
print(calculate_peptide_mz("PEPTIDEM(Oxidation)", 2))  # With oxidation
print(calculate_peptide_mz("(Acetyl)PEPTIDEK(Acetyl)", 2))  # Acetylated
```

### Generate Complete Fragment Ion Series

```python
def generate_fragment_ions(sequence_str, charge_states=[1, 2]):
    """Generate comprehensive fragment ion list."""
    peptide = oms.AASequence.fromString(sequence_str)
    proton_mass = 1.007276
    fragments = []

    for i in range(1, peptide.size()):
        # b and y ions
        b_frag = peptide.getPrefix(i)
        y_frag = peptide.getSuffix(i)

        for z in charge_states:
            b_mz = (b_frag.getMonoWeight() + z * proton_mass) / z
            y_mz = (y_frag.getMonoWeight() + z * proton_mass) / z

            fragments.append({
                'type': 'b',
                'position': i,
                'charge': z,
                'mz': b_mz
            })
            fragments.append({
                'type': 'y',
                'position': i,
                'charge': z,
                'mz': y_mz
            })

    return fragments

# Usage
ions = generate_fragment_ions("PEPTIDE", charge_states=[1, 2])
for ion in ions:
    print(f"{ion['type']}{ion['position']}^{ion['charge']}+: {ion['mz']:.4f}")
```

### Digest Protein and Calculate Peptide Masses

```python
def digest_and_calculate(protein_seq_str, enzyme="Trypsin", missed_cleavages=2,
                         min_mass=600, max_mass=4000):
    """Digest protein and return valid peptides with masses."""
    dig = oms.ProteaseDigestion()
    dig.setEnzyme(enzyme)
    dig.setMissedCleavages(missed_cleavages)

    protein = oms.AASequence.fromString(protein_seq_str)
    peptides = []
    dig.digest(protein, peptides)

    results = []
    for pep in peptides:
        mass = pep.getMonoWeight()
        if min_mass <= mass <= max_mass:
            results.append({
                'sequence': pep.toString(),
                'mass': mass,
                'length': pep.size()
            })

    return results

# Usage
protein = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
peptides = digest_and_calculate(protein)
for pep in peptides:
    print(f"{pep['sequence']}: {pep['mass']:.2f} Da ({pep['length']} aa)")
```
