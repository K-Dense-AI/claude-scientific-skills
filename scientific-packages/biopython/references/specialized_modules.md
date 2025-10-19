# BioPython Specialized Analysis Modules

This document covers BioPython's specialized modules for structural biology, phylogenetics, population genetics, and other advanced analyses.

## Structural Bioinformatics

### Bio.PDB - Protein Structure Analysis

Comprehensive tools for handling macromolecular crystal structures.

#### Structure Hierarchy

PDB structures are organized hierarchically:
- **Structure** → Models → Chains → Residues → Atoms

```python
from Bio.PDB import PDBParser

parser = PDBParser()
structure = parser.get_structure("protein", "1abc.pdb")

# Navigate hierarchy
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                print(atom.coord)  # xyz coordinates
```

#### Parsing Structure Files

**PDB format:**
```python
from Bio.PDB import PDBParser
parser = PDBParser(QUIET=True)
structure = parser.get_structure("myprotein", "structure.pdb")
```

**mmCIF format:**
```python
from Bio.PDB import MMCIFParser
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("myprotein", "structure.cif")
```

**Fast mmCIF parser:**
```python
from Bio.PDB import FastMMCIFParser
parser = FastMMCIFParser(QUIET=True)
structure = parser.get_structure("myprotein", "structure.cif")
```

**MMTF format:**
```python
from Bio.PDB import MMTFParser
parser = MMTFParser()
structure = parser.get_structure("structure.mmtf")
```

**Binary CIF:**
```python
from Bio.PDB.binary_cif import BinaryCIFParser
parser = BinaryCIFParser()
structure = parser.get_structure("structure.bcif")
```

#### Downloading Structures

```python
from Bio.PDB import PDBList
pdbl = PDBList()

# Download specific structure
pdbl.retrieve_pdb_file("1ABC", file_format="pdb", pdir="structures/")

# Download entire PDB (obsolete entries)
pdbl.download_obsolete_entries(pdir="obsolete/")

# Update local PDB mirror
pdbl.update_pdb()
```

#### Structure Selection and Filtering

```python
# Select specific chains
chain_A = structure[0]['A']

# Select specific residues
residue_10 = chain_A[10]

# Select specific atoms
ca_atom = residue_10['CA']

# Iterate over specific atom types
for atom in structure.get_atoms():
    if atom.name == 'CA':  # Alpha carbons only
        print(atom.coord)
```

**Structure selectors:**
```python
from Bio.PDB.Polypeptide import is_aa

# Filter by residue type
for residue in structure.get_residues():
    if is_aa(residue):
        print(f"Amino acid: {residue.resname}")
```

#### Secondary Structure Analysis

**DSSP integration:**
```python
from Bio.PDB import DSSP

# Requires DSSP program installed
model = structure[0]
dssp = DSSP(model, "structure.pdb")

# Access secondary structure
for key in dssp:
    secondary_structure = dssp[key][2]
    accessibility = dssp[key][3]
    print(f"Residue {key}: {secondary_structure}, accessible: {accessibility}")
```

DSSP codes:
- H: Alpha helix
- B: Beta bridge
- E: Extended strand (beta sheet)
- G: 3-10 helix
- I: Pi helix
- T: Turn
- S: Bend
- -: Coil

#### Solvent Accessibility

**Shrake-Rupley algorithm:**
```python
from Bio.PDB import ShrakeRupley

sr = ShrakeRupley()
sr.compute(structure, level="R")  # R=residue, A=atom, C=chain, M=model, S=structure

for residue in structure.get_residues():
    print(f"{residue.resname} {residue.id[1]}: {residue.sasa} Ų")
```

**NACCESS wrapper:**
```python
from Bio.PDB import NACCESS

# Requires NACCESS program
naccess = NACCESS("structure.pdb")
for residue_id, data in naccess.items():
    print(f"Residue {residue_id}: {data['all_atoms_abs']} Ų")
```

**Half-sphere exposure:**
```python
from Bio.PDB import HSExposure

# Requires DSSP
model = structure[0]
hse = HSExposure()
hse.calc_hs_exposure(model, "structure.pdb")

for chain in model:
    for residue in chain:
        if residue.has_id('EXP_HSE_A_U'):
            hse_up = residue.xtra['EXP_HSE_A_U']
            hse_down = residue.xtra['EXP_HSE_A_D']
```

#### Structural Alignment and Superimposition

**Standard superimposition:**
```python
from Bio.PDB import Superimposer

sup = Superimposer()
sup.set_atoms(ref_atoms, alt_atoms)  # Lists of atoms to align
sup.apply(structure2.get_atoms())  # Apply transformation

print(f"RMSD: {sup.rms}")
print(f"Rotation matrix: {sup.rotran[0]}")
print(f"Translation vector: {sup.rotran[1]}")
```

**QCP (Quaternion Characteristic Polynomial) method:**
```python
from Bio.PDB import QCPSuperimposer

qcp = QCPSuperimposer()
qcp.set(ref_coords, alt_coords)
qcp.run()
print(f"RMSD: {qcp.get_rms()}")
```

#### Geometric Calculations

**Distances and angles:**
```python
# Distance between atoms
from Bio.PDB import Vector
dist = atom1 - atom2  # Returns distance

# Angle between three atoms
from Bio.PDB import calc_angle
angle = calc_angle(atom1.coord, atom2.coord, atom3.coord)

# Dihedral angle
from Bio.PDB import calc_dihedral
dihedral = calc_dihedral(atom1.coord, atom2.coord, atom3.coord, atom4.coord)
```

**Vector operations:**
```python
from Bio.PDB.Vector import Vector

v1 = Vector(atom1.coord)
v2 = Vector(atom2.coord)

# Vector operations
v3 = v1 + v2
v4 = v1 - v2
dot_product = v1 * v2
cross_product = v1 ** v2
magnitude = v1.norm()
normalized = v1.normalized()
```

#### Internal Coordinates

Advanced residue geometry representation:
```python
from Bio.PDB import internal_coords

# Enable internal coordinates
structure.atom_to_internal_coordinates()

# Access phi, psi angles
for residue in structure.get_residues():
    if residue.internal_coord:
        print(f"Phi: {residue.internal_coord.get_angle('phi')}")
        print(f"Psi: {residue.internal_coord.get_angle('psi')}")
```

#### Writing Structures

```python
from Bio.PDB import PDBIO

io = PDBIO()
io.set_structure(structure)
io.save("output.pdb")

# Save specific selection
io.save("chain_A.pdb", select=ChainSelector("A"))
```

### Bio.SCOP - SCOP Database

Access to Structural Classification of Proteins database.

### Bio.KEGG - Pathway Analysis

Interface to KEGG (Kyoto Encyclopedia of Genes and Genomes) databases:

**Capabilities:**
- Access pathway maps
- Retrieve enzyme data
- Get compound information
- Query orthology relationships

## Phylogenetics

### Bio.Phylo - Phylogenetic Tree Analysis

Comprehensive phylogenetic tree manipulation and analysis.

#### Reading and Writing Trees

**Supported formats:**
- Newick: Simple, widely-used format
- NEXUS: Rich metadata format
- PhyloXML: XML-based with extensive annotations
- NeXML: Modern XML standard

```python
from Bio import Phylo

# Read tree
tree = Phylo.read("tree.nwk", "newick")

# Read multiple trees
trees = list(Phylo.parse("trees.nex", "nexus"))

# Write tree
Phylo.write(tree, "output.nwk", "newick")
```

#### Tree Visualization

**ASCII visualization:**
```python
Phylo.draw_ascii(tree)
```

**Matplotlib plotting:**
```python
import matplotlib.pyplot as plt
Phylo.draw(tree)
plt.show()

# With customization
fig, ax = plt.subplots(figsize=(10, 8))
Phylo.draw(tree, axes=ax, do_show=False)
ax.set_title("My Phylogenetic Tree")
plt.show()
```

#### Tree Navigation and Manipulation

**Find clades:**
```python
# Get all terminal nodes (leaves)
terminals = tree.get_terminals()

# Get all nonterminal nodes
nonterminals = tree.get_nonterminals()

# Find specific clade
target = tree.find_any(name="Species_A")

# Find all matching clades
matches = tree.find_clades(terminal=True)
```

**Tree properties:**
```python
# Count terminals
num_species = tree.count_terminals()

# Get total branch length
total_length = tree.total_branch_length()

# Check if tree is bifurcating
is_bifurcating = tree.is_bifurcating()

# Get maximum distance from root
max_dist = tree.distance(tree.root)
```

**Tree modification:**
```python
# Prune tree to specific taxa
keep_taxa = ["Species_A", "Species_B", "Species_C"]
tree.prune(keep_taxa)

# Collapse short branches
tree.collapse_all(lambda c: c.branch_length < 0.01)

# Ladderize (sort branches)
tree.ladderize()

# Root tree at midpoint
tree.root_at_midpoint()

# Root at specific clade
outgroup = tree.find_any(name="Outgroup_species")
tree.root_with_outgroup(outgroup)
```

**Calculate distances:**
```python
# Distance between two clades
dist = tree.distance(clade1, clade2)

# Distance from root
root_dist = tree.distance(tree.root, terminal_clade)
```

#### Tree Construction

**Distance-based methods:**
```python
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from Bio import AlignIO

# Load alignment
aln = AlignIO.read("alignment.fasta", "fasta")

# Calculate distance matrix
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(aln)

# Construct tree using UPGMA
constructor = DistanceTreeConstructor()
tree_upgma = constructor.upgma(dm)

# Or using Neighbor-Joining
tree_nj = constructor.nj(dm)
```

**Parsimony method:**
```python
from Bio.Phylo.TreeConstruction import ParsimonyScorer, NNITreeSearcher

scorer = ParsimonyScorer()
searcher = NNITreeSearcher(scorer)
tree = searcher.search(starting_tree, alignment)
```

**Distance calculators:**
- 'identity': Simple identity scoring
- 'blastn': BLAST nucleotide scoring
- 'blastp': BLAST protein scoring
- 'dnafull': EMBOSS DNA scoring matrix
- 'blosum62': BLOSUM62 protein matrix
- 'pam250': PAM250 protein matrix

#### Consensus Trees

```python
from Bio.Phylo.Consensus import majority_consensus, strict_consensus

# Strict consensus
consensus_strict = strict_consensus(trees)

# Majority rule consensus
consensus_majority = majority_consensus(trees, cutoff=0.5)

# Bootstrap consensus
from Bio.Phylo.Consensus import bootstrap_consensus
bootstrap_tree = bootstrap_consensus(trees, cutoff=0.7)
```

#### External Tool Wrappers

**PhyML:**
```python
from Bio.Phylo.Applications import PhymlCommandline

cmd = PhymlCommandline(input="alignment.phy", datatype="nt", model="HKY85", alpha="e", bootstrap=100)
stdout, stderr = cmd()
tree = Phylo.read("alignment.phy_phyml_tree.txt", "newick")
```

**RAxML:**
```python
from Bio.Phylo.Applications import RaxmlCommandline

cmd = RaxmlCommandline(
    sequences="alignment.phy",
    model="GTRGAMMA",
    name="mytree",
    parsimony_seed=12345
)
stdout, stderr = cmd()
```

**FastTree:**
```python
from Bio.Phylo.Applications import FastTreeCommandline

cmd = FastTreeCommandline(input="alignment.fasta", out="tree.nwk", gtr=True, gamma=True)
stdout, stderr = cmd()
```

### Bio.Phylo.PAML - Evolutionary Analysis

Interface to PAML (Phylogenetic Analysis by Maximum Likelihood):

**CODEML - Codon-based analysis:**
```python
from Bio.Phylo.PAML import codeml

cml = codeml.Codeml()
cml.alignment = "alignment.phy"
cml.tree = "tree.nwk"
cml.out_file = "results.out"
cml.working_dir = "./paml_wd"

# Set parameters
cml.set_options(
    seqtype=1,      # Codon sequences
    model=0,        # One omega ratio
    NSsites=[0, 1, 2],  # Test different models
    CodonFreq=2,    # F3x4 codon frequencies
)

results = cml.run()
```

**BaseML - Nucleotide-based analysis:**
```python
from Bio.Phylo.PAML import baseml

bml = baseml.Baseml()
bml.alignment = "alignment.phy"
bml.tree = "tree.nwk"
results = bml.run()
```

**YN00 - Yang-Nielsen method:**
```python
from Bio.Phylo.PAML import yn00

yn = yn00.Yn00()
yn.alignment = "alignment.phy"
results = yn.run()
```

## Population Genetics

### Bio.PopGen - Population Genetics Analysis

Tools for population-level genetic analysis.

**Capabilities:**
- Allele frequency calculations
- Hardy-Weinberg equilibrium testing
- Linkage disequilibrium analysis
- F-statistics (FST, FIS, FIT)
- Tajima's D
- Population structure analysis

## Clustering and Machine Learning

### Bio.Cluster - Clustering Algorithms

Statistical clustering for gene expression and other biological data:

**Hierarchical clustering:**
```python
from Bio.Cluster import treecluster

tree = treecluster(data, method='a', dist='e')
# method: 'a'=average, 's'=single, 'm'=maximum, 'c'=centroid
# dist: 'e'=Euclidean, 'c'=correlation, 'a'=absolute correlation
```

**k-means clustering:**
```python
from Bio.Cluster import kcluster

clusterid, error, nfound = kcluster(data, nclusters=5, npass=100)
```

**Self-Organizing Maps (SOM):**
```python
from Bio.Cluster import somcluster

clusterid, celldata = somcluster(data, nx=3, ny=3)
```

**Principal Component Analysis:**
```python
from Bio.Cluster import pca

columnmean, coordinates, components, eigenvalues = pca(data)
```

## Visualization

### Bio.Graphics - Genomic Visualization

Tools for creating publication-quality biological graphics.

**GenomeDiagram - Circular and linear genome maps:**
```python
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO

record = SeqIO.read("genome.gb", "genbank")

gd_diagram = GenomeDiagram.Diagram("Genome Map")
gd_track = gd_diagram.new_track(1, greytrack=True)
gd_feature_set = gd_track.new_set()

# Add features
for feature in record.features:
    if feature.type == "gene":
        gd_feature_set.add_feature(feature, color="blue", label=True)

gd_diagram.draw(format="linear", pagesize='A4', fragments=1)
gd_diagram.write("genome_map.pdf", "PDF")
```

**Chromosomes - Chromosome visualization:**
```python
from Bio.Graphics.BasicChromosome import Chromosome

chr = Chromosome("Chromosome 1")
chr.add("gene1", 1000, 2000, color="red")
chr.add("gene2", 3000, 4500, color="blue")
```

## Phenotype Analysis

### Bio.phenotype - Phenotypic Microarray Analysis

Tools for analyzing phenotypic microarray data (e.g., Biolog plates):

**Capabilities:**
- Parse PM plate data
- Growth curve analysis
- Compare phenotypic profiles
- Calculate similarity metrics
