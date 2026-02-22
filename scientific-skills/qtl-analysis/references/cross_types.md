# Cross Types Guide

Complete guide to experimental cross types for QTL mapping.

## Overview

| Cross | Mating | Progeny | Advantages | Best For |
|-------|--------|---------|------------|----------|
| **F2** | F1 × F1 | Segregating | Max recombination, simple | First-pass mapping |
| **BC** | F1 × Parent | Segregating | Homogeneous background | Confirming QTL |
| **RIL** | F2 × F2 (sib) | Inbred | Many breakpoints, repeatable | Fine mapping |
| **DO** | Multiple rounds | Outbred | High diversity, precision | Complex traits |
| **MPP** | Multiple parents | Segregating | High power, diverse | Systems genetics |

---

## F2 Intercross

### Design
```
Generation 0:  P1 (AA)  ×  P2 (aa)
                ↓
Generation 1:      F1 (Aa)
                     ↓
Generation 2:  F1 (Aa) × F1 (Aa)
                     ↓
Generation 3:      F2 Progeny
                   Genotypes: 1/4 AA, 1/2 Aa, 1/4 aa
```

### Genotype Frequencies

| Genotype | Frequency | Code |
|----------|-----------|------|
| AA (P1 homozygous) | 25% | 1 or "S" |
| Aa (Heterozygous) | 50% | 2 or "H" |
| aa (P2 homozygous) | 25% | 3 or "B" |

### Expected Genotype Probabilities

Given flanking markers A and B, the probability of Q being genotype g:

```
P(Q=AA) = P(A=AA) × P(B=AA) × (1-r)² / P(observed)
P(Q=Aa) = [P(A=AA) × P(B=Aa) + P(A=Aa) × P(B=AA)] × r(1-r) / P(observed)
...
```

Where r is the recombination fraction.

### When to Use
- **Best for**: Initial QTL discovery
- **Power**: Maximum for detecting dominance effects
- **Sample size**: 200-500 individuals recommended
- **Limitations**: Cannot distinguish maternal vs paternal effects

### R/qtl2 Implementation
```r
# Simulate F2 cross
map <- sim_map(len = c(100, 80), n_mar = c(50, 40))
cross <- sim_cross(map, n_ind = 200, type = "f2")

# Calculate genotype probabilities
pr <- calc_genoprob(cross, error_prob = 0.002)

# Scan for QTLs
out <- scan1(pr, cross$pheno[, 1])
```

---

## Backcross (BC)

### Design
```
Generation 0:  P1 (AA)  ×  P2 (aa)
                ↓
Generation 1:      F1 (Aa)
                     ↓
Generation 2:   F1 (Aa) × P1 (AA)
                     ↓
Generation 3:      BC Progeny
                   Genotypes: 1/2 AA, 1/2 Aa
```

### Genotype Frequencies

| Genotype | Frequency | Code |
|----------|-----------|------|
| AA | 50% | 1 or "A" |
| Aa | 50% | 2 or "H" |

### When to Use
- **Best for**: Confirming QTL in controlled background
- **Power**: Good for additive effects
- **Advantage**: Only two genotypes, simpler analysis
- **Limitation**: Cannot estimate dominance effects

### R/qtl2 Implementation
```r
# Simulate backcross
cross <- sim_cross(map, n_ind = 200, type = "bc")

# Calculate probabilities
pr <- calc_genoprob(cross, error_prob = 0.002)

# Scan (simpler - only 2 genotypes)
out <- scan1(pr, cross$pheno[, 1])
```

---

## Recombinant Inbred Lines (RIL)

### Types

#### RI by Selfing (RIL-self)
```
F2 → Individual plants self → F3 → ... → F8 → RILs
Result: Nearly homozygous lines (AA or aa)
```

#### RI by Sibmating (RIL-sib)
```
F2 → Brother-sister mating → F3 → ... → RILs
Result: Nearly homozygous lines
```

### Design
```
Generation 0:  P1 (AA)  ×  P2 (aa)
                ↓
Generation 1:      F1 (Aa)
                     ↓
Generation 2:      F2 (various)
                     ↓
Generations 3+:  Repeated selfing or sibmating
                     ↓
Final:           RIL (nearly AA or aa)
```

### Genotype Frequencies

| Genotype | Frequency | Code |
|----------|-----------|------|
| AA | ~50% | 1 or "S" |
| aa | ~50% | 3 or "B" |
| Heterozygous | <5% (residual) | 2 or "H" |

### Advantages
1. **Reproducible**: Same genotype in all replicates
2. **Many breakpoints**: ~2× F2 recombination
3. **Replication**: Can measure multiple times
4. **Environment studies**: Same genotype in different environments

### When to Use
- **Best for**: Fine mapping, multiple trait measurements
- **Generations**: 6-8 generations of inbreeding
- **Sample size**: 100-200 lines (but can be replicated)

### R/qtl2 Implementation
```r
# Simulate RIL (by selfing)
cross <- sim_cross(map, n_ind = 150, type = "riself")

# Or by sibmating
cross <- sim_cross(map, n_ind = 150, type = "risib")
```

---

## Diversity Outbred (DO)

### Design
```
Generation 0:  8 Founder strains (CC, AA, BB, DD, EE, FF, GG, HH)
                ↓
Generation 1:     F1 (2-way crosses)
                ↓
Generation 2:     F2 (4-way crosses)
                ↓
Generation 3:     F3 (8-way cross) - "G0"
                ↓
Generation 4+:   Random outbreeding for 20+ generations
                ↓
Final:          DO Population (highly recombined, outbred)
```

### Founder Strains (Collaborative Cross)
The 8 founder strains of the CC/DO:
1. A/J
2. C57BL/6J
3. 129S1/SvImJ
4. NOD/ShiLtJ
5. NZO/HlLtJ
6. CAST/EiJ
7. PWK/PhJ
8. WSB/EiJ

### Genotype States

| State | Haplotype Pair | Probability |
|-------|----------------|-------------|
| AA | A/J homozygous | ~1/64 |
| CC | C57BL/6J homozygous | ~1/64 |
| ... | ... | ... |
| AC | A/J + C57BL/6J | ~2/64 |
| etc. | 36 total genotypes | various |

### Advantages
1. **High resolution**: ~20 generations of recombination
2. **Mapping precision**: <1 Mb typical
3. **Population diversity**: Captures genetic variation from 8 founders
4. **Community resource**: Standardized populations available

### When to Use
- **Best for**: High-precision mapping, fine-mapping QTLs
- **Power**: Very high for common variants
- **Analysis**: Requires haplotype reconstruction (HMM)
- **Software**: R/qtl2, DOQTL

### R/qtl2 Implementation
```r
# DO crosses require special handling
# Use calc_genoprob with map and error probability
pr <- calc_genoprob(cross, map, error_prob = 0.002)

# For DO, use allele probabilities instead of genotype
apr <- genoprob_to_alleleprob(pr)

# Kinship is essential for DO
kinship <- calc_kinship(pr, type = "loco")

# Scan with kinship
out <- scan1(pr, cross$pheno[, 1], kinship)
```

---

## Multi-Parent Populations (MPP)

### Types

#### Nested Association Mapping (NAM)
```
Multiple diverse founders × Common parent
    ↓
Multiple RIL families with shared parent
```

#### Multiparent Advanced Generation InterCross (MAGIC)
```
Multiple founders → Multiple generations of intercrossing → RILs
```

#### Collaborative Cross (CC)
```
8 founders → F1 → 4-way → 8-way → Inbreeding → CC strains
```

### Design Comparison

| Population | Founders | Structure | Recombination | Power |
|------------|----------|-----------|---------------|-------|
| NAM | 25+ | Families | Moderate | High |
| MAGIC | 4-16 | RILs | High | Very High |
| CC | 8 | Inbred strains | Very High | High |
| DO | 8 | Outbred | Very High | Very High |

### Analysis Considerations

1. **Population structure**: Must account for relatedness
2. **Kinship**: Use LOCO (leave-one-chromosome-out) kinship
3. **Multiple alleles**: Not just biallelic
4. **Haplotypes**: Infer founder haplotypes, not just genotypes

### When to Use
- **Best for**: Complex trait dissection, systems genetics
- **Trade-off**: More complex analysis vs. more power
- **Requirements**: Specialized software, larger sample sizes

---

## Haploid and Doubled Haploid (DH)

### Design
```
F1 hybrid (Aa) → Gametes (A or a) → Haploid → Chromosome doubling → DH lines
Result: 50% AA, 50% aa (permanent, homozygous)
```

### Common in
- **Yeast**: Natural haploidy
- **Plants**: Anther culture, wide hybridization
- **Fish**: Gynogenesis, androgenesis

### Advantages
- Instant homozygosity
- No dominance effects to estimate
- Simple analysis

### R/qtl2 Implementation
```r
cross <- sim_cross(map, n_ind = 200, type = "dh")
# Or for haploid
cross <- sim_cross(map, n_ind = 200, type = "haploid")
```

---

## Cross Type Selection Guide

### Decision Tree

```
Starting QTL mapping?
    ├── Yes → F2 (maximum recombination, simple)
    └── No → Confirming specific QTL?
                ├── Yes → BC (controlled background)
                └── No → Need replication?
                            ├── Yes → RIL (reproducible, stable)
                            └── No → High resolution needed?
                                        ├── Yes → DO/MPP (>20 gen recomb)
                                        └── No → Model organism?
                                                    ├── Yeast → Haploid
                                                    └── Plant → DH possible?
```

### Recommended Sample Sizes

| Cross | Minimum | Adequate | High Power |
|-------|---------|----------|------------|
| F2 | 100 | 200-300 | 500+ |
| BC | 100 | 200 | 300+ |
| RIL | 50 | 100 | 200+ |
| DO | 100 | 200 | 500+ |
| MPP | 200 | 500 | 1000+ |

### Effect Size Detection

| Cross | Small Effect (5%) | Medium (10%) | Large (20%) |
|-------|-------------------|--------------|-------------|
| F2 | 500 | 200 | 50 |
| BC | 500 | 200 | 50 |
| RIL | 300 | 100 | 30 |
| DO | 200 | 75 | 25 |

---

## Software Compatibility

| Software | F2 | BC | RIL | DO | MPP |
|----------|----|-----|-----|-----|-----|
| R/qtl | ✓ | ✓ | ✓ | - | - |
| R/qtl2 | ✓ | ✓ | ✓ | ✓ | ✓ |
| MapManager | ✓ | ✓ | ✓ | - | - |
| QTL Cartographer | ✓ | ✓ | ✓ | - | - |
| DOQTL | - | - | - | ✓ | - |
| GEMMA | ✓ | ✓ | ✓ | ✓ | ✓ |
| tensorQTL | ✓ | ✓ | ✓ | ✓ | ✓ |

---

## Simulation Reference

### Recombination Fractions

```r
# Convert genetic (cM) to recombination fraction
# Haldane mapping function (no interference)
r <- 0.5 * (1 - exp(-d/50))

# Kosambi mapping function (some interference)
r <- 0.5 * tanh(2*d/100)

# Where d is distance in cM
```

### Expected Breakpoints

| Cross | Generations | Expected Breakpoints per Chromosome |
|-------|-------------|-------------------------------------|
| F2 | 1 | 1 |
| BC | 1 | 0.5 |
| RIL | 6 | ~3 |
| RIL | 8 | ~4 |
| DO | 20 | ~7-8 |
