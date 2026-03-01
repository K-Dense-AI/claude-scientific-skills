"""Static mapping of scientific skills to their primary domain agent."""

from __future__ import annotations

# Each skill appears in exactly one domain.
# Domain agents can still cross-search skills via the registry's skill_search() without domain filter.

DOMAIN_SKILLS: dict[str, list[str]] = {
    "bioinformatics": [
        # Genomics & sequencing
        "alphafold-database", "biopython", "ensembl-database", "gene-database",
        "gget", "pysam", "tiledbvcf", "ena-database", "gtars",
        # Single-cell & transcriptomics
        "scanpy", "anndata", "scvi-tools", "cellxgene-census", "pydeseq2", "arboreto",
        # Protein & structural biology
        "esm", "pdb-database", "uniprot-database", "string-database", "diffdock",
        # Pathway & metabolic modeling
        "cobrapy", "kegg-database", "reactome-database", "brenda-database",
        # Clinical genomics
        "clinvar-database", "cosmic-database", "gwas-database", "clinpgx-database",
        # Epigenomics & gene regulation
        "geniml", "deeptools",
        # Lab automation & integration
        "opentrons-integration", "pylabrobot", "benchling-integration",
        "protocolsio-integration", "labarchive-integration", "latchbio-integration",
        "dnanexus-integration", "omero-integration", "adaptyv",
        # Phylogenetics & biodiversity
        "etetoolkit", "scikit-bio",
        # Biomedical imaging
        "flowio", "histolab", "pathml", "pydicom", "imaging-data-commons",
        # Neuroscience
        "neuropixels-analysis", "neurokit2",
        # Molecular biology tools
        "primer-design",
        # Bio data management
        "lamindb", "geo-database", "bioservices", "zarr-python",
    ],
    "chemistry": [
        # Cheminformatics core
        "rdkit", "datamol", "molfeat",
        # Drug discovery & ML
        "deepchem", "pytdc", "torchdrug", "rowan",
        # Chemical databases
        "pubchem-database", "chembl-database", "drugbank-database",
        "zinc-database", "opentargets-database",
        # Medicinal chemistry
        "medchem",
        # Mass spectrometry & metabolomics
        "matchms", "pyopenms", "metabolomics-workbench-database", "hmdb-database",
        # Materials science
        "pymatgen",
        # Patents
        "uspto-database",
    ],
    "clinical": [
        # Clinical reports & documentation
        "clinical-reports", "clinical-decision-support", "treatment-plans",
        # Clinical trials & regulatory
        "clinicaltrials-database", "fda-database",
        # Digital health AI
        "pyhealth",
        # Quality & regulatory
        "iso-13485-certification",
    ],
    "computation": [
        # Machine learning
        "scikit-learn", "scikit-survival", "shap", "umap-learn", "hypogenic",
        # Deep learning
        "pytorch-lightning", "transformers", "torch_geometric",
        # Reinforcement learning
        "stable-baselines3", "pufferlib",
        # Statistics
        "statsmodels", "statistical-analysis", "pymc",
        # Optimization
        "pymoo",
        # Quantum computing
        "qiskit", "cirq", "pennylane", "qutip",
        # Simulation
        "simpy", "fluidsim",
        # Data processing
        "polars", "dask", "vaex",
        # Symbolic math
        "sympy",
        # Cloud compute
        "modal",
        # Graph analysis
        "networkx",
        # Astronomy & geospatial
        "astropy", "geopandas",
        # Time series
        "timesfm-forecasting", "aeon",
        # Resources & EDA
        "get-available-resources", "exploratory-data-analysis",
        # MATLAB
        "matlab",
        # Data commons
        "datacommons-client",
    ],
    "literature": [
        # Paper search & databases
        "pubmed-database", "openalex-database", "biorxiv-database",
        "bgpt-paper-search", "perplexity-search", "research-lookup",
        # Literature review
        "literature-review",
        # Scientific writing
        "scientific-writing", "manuscript-writer",
        # Research assistant (P1)
        "research-assistant", "research-commons",
        # Discussion (result interpretation)
        "research-discussion",
        # Citation management
        "citation-management", "pyzotero",
        # Review & evaluation
        "peer-review", "scholar-evaluation", "scientific-critical-thinking",
        # Grants & proposals
        "research-grants", "venue-templates",
        # Ideation
        "hypothesis-generation", "scientific-brainstorming",
        # Document conversion
        "markitdown", "paper-2-web",
    ],
    "visualization": [
        # Plotting libraries
        "matplotlib", "plotly", "seaborn",
        # Publication figures
        "scientific-visualization",
        # Diagrams & schematics
        "scientific-schematics", "markdown-mermaid-writing",
        # Image generation
        "generate-image", "infographics",
        # Presentations
        "scientific-slides", "journal-presentation-maker",
        # Posters
        "latex-posters", "pptx-posters",
    ],
    "workflow": [
        # Code implementation & quality (P2)
        "code-implementer", "code-excellence", "solid-principles",
        # Code validation (P3)
        "code-validator",
        # Improvement analysis (P5)
        "improvement-analyzer",
        # Git workflow
        "git-workflow-manager",
        # Reference surveyor (P1.5)
        "reference-surveyor",
        # Experiment management
        "experiment-hub",
        # Research automation
        "denario",
        # Monitoring & project management
        "monitor", "team-dashboard", "weekly-briefing", "open-notebook",
        "asana-extended-api",
        # Finance & market research
        "alpha-vantage", "fred-economic-data", "edgartools",
        "hedgefundmonitor", "usfiscaldata", "market-research-reports",
        # Misc
        "offer-k-dense-web",
    ],
}


def get_domain_for_skill(skill_name: str) -> str | None:
    """Look up which domain a skill belongs to."""
    for domain, skills in DOMAIN_SKILLS.items():
        if skill_name in skills:
            return domain
    return None


def get_all_skill_names() -> list[str]:
    """Return flat list of all mapped skill names."""
    return [s for skills in DOMAIN_SKILLS.values() for s in skills]
