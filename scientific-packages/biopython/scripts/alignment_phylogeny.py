#!/usr/bin/env python3
"""
Sequence alignment and phylogenetic analysis using BioPython.

This script demonstrates:
- Pairwise sequence alignment
- Multiple sequence alignment I/O
- Distance matrix calculation
- Phylogenetic tree construction
- Tree manipulation and visualization
"""

from Bio import Align, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import ParsimonyScorer, NNITreeSearcher
from Bio.Seq import Seq
import matplotlib.pyplot as plt


def pairwise_alignment_example():
    """Demonstrate pairwise sequence alignment."""

    print("Pairwise Sequence Alignment")
    print("=" * 60)

    # Create aligner
    aligner = Align.PairwiseAligner()

    # Set parameters
    aligner.mode = "global"  # or 'local' for local alignment
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    # Sequences to align
    seq1 = "ACGTACGTACGT"
    seq2 = "ACGTTACGTGT"

    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print()

    # Perform alignment
    alignments = aligner.align(seq1, seq2)

    # Show results
    print(f"Number of optimal alignments: {len(alignments)}")
    print(f"Best alignment score: {alignments.score:.1f}")
    print()

    # Display best alignment
    print("Best alignment:")
    print(alignments[0])
    print()


def local_alignment_example():
    """Demonstrate local alignment (Smith-Waterman)."""

    print("Local Sequence Alignment")
    print("=" * 60)

    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    seq1 = "AAAAACGTACGTACGTAAAAA"
    seq2 = "TTTTTTACGTACGTTTTTTT"

    print(f"Sequence 1: {seq1}")
    print(f"Sequence 2: {seq2}")
    print()

    alignments = aligner.align(seq1, seq2)

    print(f"Best local alignment score: {alignments.score:.1f}")
    print()
    print("Best local alignment:")
    print(alignments[0])
    print()


def read_and_analyze_alignment(alignment_file, format="fasta"):
    """Read and analyze a multiple sequence alignment."""

    print(f"Reading alignment from: {alignment_file}")
    print("-" * 60)

    # Read alignment
    alignment = AlignIO.read(alignment_file, format)

    print(f"Number of sequences: {len(alignment)}")
    print(f"Alignment length: {alignment.get_alignment_length()}")
    print()

    # Display alignment
    print("Alignment preview:")
    for record in alignment[:5]:  # Show first 5 sequences
        print(f"{record.id[:15]:15s} {record.seq[:50]}...")

    print()

    # Calculate some statistics
    analyze_alignment_statistics(alignment)

    return alignment


def analyze_alignment_statistics(alignment):
    """Calculate statistics for an alignment."""

    print("Alignment Statistics:")
    print("-" * 60)

    # Get alignment length
    length = alignment.get_alignment_length()

    # Count gaps
    total_gaps = sum(str(record.seq).count("-") for record in alignment)
    gap_percentage = (total_gaps / (length * len(alignment))) * 100

    print(f"Total positions: {length}")
    print(f"Number of sequences: {len(alignment)}")
    print(f"Total gaps: {total_gaps} ({gap_percentage:.1f}%)")
    print()

    # Calculate conservation at each position
    conserved_positions = 0
    for i in range(length):
        column = alignment[:, i]
        # Count most common residue
        if column.count(max(set(column), key=column.count)) == len(alignment):
            conserved_positions += 1

    conservation = (conserved_positions / length) * 100
    print(f"Fully conserved positions: {conserved_positions} ({conservation:.1f}%)")
    print()


def calculate_distance_matrix(alignment):
    """Calculate distance matrix from alignment."""

    print("Calculating Distance Matrix")
    print("-" * 60)

    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)

    print("Distance matrix:")
    print(dm)
    print()

    return dm


def build_upgma_tree(alignment):
    """Build phylogenetic tree using UPGMA."""

    print("Building UPGMA Tree")
    print("=" * 60)

    # Calculate distance matrix
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)

    # Construct tree
    constructor = DistanceTreeConstructor(calculator)
    tree = constructor.upgma(dm)

    print("UPGMA tree constructed")
    print(f"Number of terminals: {tree.count_terminals()}")
    print()

    return tree


def build_nj_tree(alignment):
    """Build phylogenetic tree using Neighbor-Joining."""

    print("Building Neighbor-Joining Tree")
    print("=" * 60)

    # Calculate distance matrix
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)

    # Construct tree
    constructor = DistanceTreeConstructor(calculator)
    tree = constructor.nj(dm)

    print("Neighbor-Joining tree constructed")
    print(f"Number of terminals: {tree.count_terminals()}")
    print()

    return tree


def visualize_tree(tree, title="Phylogenetic Tree"):
    """Visualize phylogenetic tree."""

    print("Visualizing tree...")
    print()

    # ASCII visualization
    print("ASCII tree:")
    Phylo.draw_ascii(tree)
    print()

    # Matplotlib visualization
    fig, ax = plt.subplots(figsize=(10, 8))
    Phylo.draw(tree, axes=ax, do_show=False)
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig("tree_visualization.png", dpi=300, bbox_inches="tight")
    print("Tree saved to tree_visualization.png")
    print()


def manipulate_tree(tree):
    """Demonstrate tree manipulation operations."""

    print("Tree Manipulation")
    print("=" * 60)

    # Get terminals
    terminals = tree.get_terminals()
    print(f"Terminal nodes: {[t.name for t in terminals]}")
    print()

    # Get nonterminals
    nonterminals = tree.get_nonterminals()
    print(f"Number of internal nodes: {len(nonterminals)}")
    print()

    # Calculate total branch length
    total_length = tree.total_branch_length()
    print(f"Total branch length: {total_length:.4f}")
    print()

    # Find specific clade
    if len(terminals) > 0:
        target_name = terminals[0].name
        found = tree.find_any(name=target_name)
        print(f"Found clade: {found.name}")
        print()

    # Ladderize tree (sort branches)
    tree.ladderize()
    print("Tree ladderized (branches sorted)")
    print()

    # Root at midpoint
    tree.root_at_midpoint()
    print("Tree rooted at midpoint")
    print()

    return tree


def read_and_analyze_tree(tree_file, format="newick"):
    """Read and analyze a phylogenetic tree."""

    print(f"Reading tree from: {tree_file}")
    print("-" * 60)

    tree = Phylo.read(tree_file, format)

    print(f"Tree format: {format}")
    print(f"Number of terminals: {tree.count_terminals()}")
    print(f"Is bifurcating: {tree.is_bifurcating()}")
    print(f"Total branch length: {tree.total_branch_length():.4f}")
    print()

    # Show tree structure
    print("Tree structure:")
    Phylo.draw_ascii(tree)
    print()

    return tree


def compare_trees(tree1, tree2):
    """Compare two phylogenetic trees."""

    print("Comparing Trees")
    print("=" * 60)

    # Get terminal names
    terminals1 = {t.name for t in tree1.get_terminals()}
    terminals2 = {t.name for t in tree2.get_terminals()}

    print(f"Tree 1 terminals: {len(terminals1)}")
    print(f"Tree 2 terminals: {len(terminals2)}")
    print(f"Shared terminals: {len(terminals1 & terminals2)}")
    print(f"Unique to tree 1: {len(terminals1 - terminals2)}")
    print(f"Unique to tree 2: {len(terminals2 - terminals1)}")
    print()


def create_example_alignment():
    """Create an example alignment for demonstration."""

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment

    sequences = [
        SeqRecord(Seq("ACTGCTAGCTAGCTAG"), id="seq1"),
        SeqRecord(Seq("ACTGCTAGCT-GCTAG"), id="seq2"),
        SeqRecord(Seq("ACTGCTAGCTAGCTGG"), id="seq3"),
        SeqRecord(Seq("ACTGCT-GCTAGCTAG"), id="seq4"),
    ]

    alignment = MultipleSeqAlignment(sequences)

    # Save alignment
    AlignIO.write(alignment, "example_alignment.fasta", "fasta")
    print("Created example alignment: example_alignment.fasta")
    print()

    return alignment


def example_workflow():
    """Demonstrate complete alignment and phylogeny workflow."""

    print("=" * 60)
    print("BioPython Alignment & Phylogeny Workflow")
    print("=" * 60)
    print()

    # Pairwise alignment examples
    pairwise_alignment_example()
    print()
    local_alignment_example()
    print()

    # Create example data
    alignment = create_example_alignment()

    # Analyze alignment
    analyze_alignment_statistics(alignment)

    # Calculate distance matrix
    dm = calculate_distance_matrix(alignment)

    # Build trees
    upgma_tree = build_upgma_tree(alignment)
    nj_tree = build_nj_tree(alignment)

    # Manipulate tree
    manipulate_tree(upgma_tree)

    # Visualize
    visualize_tree(upgma_tree, "UPGMA Tree")

    print("Workflow completed!")
    print()


if __name__ == "__main__":
    example_workflow()

    print("Note: For real analyses, use actual alignment files.")
    print("Supported alignment formats: clustal, phylip, stockholm, nexus, fasta")
    print("Supported tree formats: newick, nexus, phyloxml, nexml")
