"""Test the SkillRegistry — parsing, indexing, search, load."""

from __future__ import annotations

from pathlib import Path

import pytest

from research_agent.domains import DOMAIN_SKILLS
from research_agent.skills.registry import SkillRegistry

SKILLS_DIR = Path(__file__).resolve().parents[2] / "scientific-skills"


@pytest.fixture(scope="module")
def registry():
    """Build registry once for all tests in this module."""
    return SkillRegistry(skills_dir=SKILLS_DIR, domain_map=DOMAIN_SKILLS)


class TestIndexBuilding:
    def test_index_has_skills(self, registry):
        """Should index a significant number of skills."""
        assert registry.size >= 100, f"Only {registry.size} skills indexed (expected 100+)"

    def test_all_domains_represented(self, registry):
        """Every domain should have at least one skill."""
        for domain in ["bioinformatics", "chemistry", "clinical", "computation",
                       "literature", "visualization", "workflow"]:
            skills = registry.list_skills(domain)
            assert len(skills) > 0, f"No skills in domain: {domain}"


class TestSearch:
    def test_exact_name_match(self, registry):
        """'scanpy' should return scanpy as first result."""
        results = registry.search("scanpy")
        assert len(results) > 0
        assert results[0].name == "scanpy"

    def test_keyword_match(self, registry):
        """'molecular docking' should find diffdock."""
        results = registry.search("molecular docking")
        names = [r.name for r in results]
        assert any("diffdock" in n or "dock" in n.lower() for n in names)

    def test_domain_filter(self, registry):
        """Search within a specific domain only returns that domain."""
        results = registry.search("analysis", domain="chemistry")
        assert all(r.domain == "chemistry" for r in results)

    def test_no_results(self, registry):
        """Nonsense query returns empty list."""
        results = registry.search("xyzzy_nonexistent_12345")
        assert len(results) == 0

    def test_limit(self, registry):
        """Limit parameter caps results."""
        results = registry.search("data", limit=3)
        assert len(results) <= 3


class TestLoad:
    def test_load_existing(self, registry):
        """Load a known skill."""
        content = registry.load("rdkit")
        assert content is not None
        assert len(content) > 100

    def test_load_nonexistent(self, registry):
        """Nonexistent skill returns None."""
        assert registry.load("nonexistent_skill_xyz") is None


class TestListSkills:
    def test_list_all(self, registry):
        """List all skills returns a non-empty list."""
        all_skills = registry.list_skills()
        assert len(all_skills) >= 100

    def test_list_by_domain(self, registry):
        """List by domain returns only that domain."""
        bio_skills = registry.list_skills("bioinformatics")
        assert len(bio_skills) > 10
        assert all(s.domain == "bioinformatics" for s in bio_skills)

    def test_get_domains(self, registry):
        """Should return 7 domains."""
        domains = registry.get_domains()
        assert len(domains) == 7
        assert "bioinformatics" in domains
        assert "literature" in domains
