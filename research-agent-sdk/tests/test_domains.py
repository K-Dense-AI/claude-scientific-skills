"""Test domain-to-skill mapping integrity."""

from __future__ import annotations

from research_agent.domains import DOMAIN_SKILLS, get_domain_for_skill, get_all_skill_names


def test_seven_domains_exist():
    """Exactly 7 domains defined."""
    assert len(DOMAIN_SKILLS) == 7
    expected = {"bioinformatics", "chemistry", "clinical", "computation",
                "literature", "visualization", "workflow"}
    assert set(DOMAIN_SKILLS.keys()) == expected


def test_no_duplicate_skills():
    """Each skill appears in exactly one domain."""
    seen: dict[str, str] = {}
    for domain, skills in DOMAIN_SKILLS.items():
        for skill in skills:
            assert skill not in seen, (
                f"Skill '{skill}' appears in both '{seen[skill]}' and '{domain}'"
            )
            seen[skill] = domain


def test_total_skill_count():
    """Should have 100+ skills mapped."""
    all_skills = get_all_skill_names()
    assert len(all_skills) >= 100, f"Only {len(all_skills)} skills mapped"


def test_get_domain_for_skill():
    """get_domain_for_skill returns correct domain."""
    assert get_domain_for_skill("scanpy") == "bioinformatics"
    assert get_domain_for_skill("rdkit") == "chemistry"
    assert get_domain_for_skill("clinical-reports") == "clinical"
    assert get_domain_for_skill("scikit-learn") == "computation"
    assert get_domain_for_skill("pubmed-database") == "literature"
    assert get_domain_for_skill("matplotlib") == "visualization"
    assert get_domain_for_skill("code-implementer") == "workflow"


def test_get_domain_for_unknown_skill():
    """Unknown skill returns None."""
    assert get_domain_for_skill("nonexistent_skill") is None


def test_key_skills_mapped():
    """Key scientific skills are in the mapping."""
    all_skills = get_all_skill_names()
    key_skills = [
        "biopython", "esm", "ensembl-database", "alphafold-database",
        "rdkit", "pubchem-database", "deepchem",
        "scikit-learn", "pytorch-lightning", "statsmodels",
        "pubmed-database", "scientific-writing",
        "matplotlib", "plotly",
        "code-implementer", "git-workflow-manager",
    ]
    for skill in key_skills:
        assert skill in all_skills, f"Key skill '{skill}' not in mapping"
