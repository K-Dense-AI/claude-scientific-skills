"""Skill registry — index, search, and load scientific skills from SKILL.md files."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import yaml


@dataclass(frozen=True)
class SkillEntry:
    """Parsed skill metadata from SKILL.md YAML frontmatter."""

    name: str
    description: str
    domain: str
    path: Path
    license: str = ""
    metadata: dict = field(default_factory=dict)


class SkillRegistry:
    """Index and search scientific skills.

    Builds an in-memory index from SKILL.md files at construction time.
    Search is keyword-based against name + description fields.
    """

    def __init__(self, skills_dir: Path, domain_map: dict[str, list[str]]) -> None:
        self._skills_dir = skills_dir
        self._domain_map = domain_map
        self._reverse_map: dict[str, str] = {}
        self._entries: dict[str, SkillEntry] = {}
        self._build_reverse_map()
        self._build_index()

    def _build_reverse_map(self) -> None:
        """Invert domain_map to get skill_name -> domain lookup."""
        for domain, skills in self._domain_map.items():
            for skill_name in skills:
                self._reverse_map[skill_name] = domain

    def _build_index(self) -> None:
        """Scan skills_dir for SKILL.md files and parse frontmatter."""
        if not self._skills_dir.exists():
            return
        for skill_md in self._skills_dir.rglob("SKILL.md"):
            entry = self._parse_skill_md(skill_md)
            if entry:
                self._entries[entry.name] = entry

    def _parse_skill_md(self, path: Path) -> SkillEntry | None:
        """Parse YAML frontmatter from a SKILL.md file."""
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            return None

        frontmatter = self._extract_frontmatter(text)
        dir_name = path.parent.name

        if not frontmatter or "name" not in frontmatter:
            domain = self._reverse_map.get(dir_name)
            if domain:
                return SkillEntry(
                    name=dir_name,
                    description="",
                    domain=domain,
                    path=path,
                )
            return None

        name = frontmatter["name"]
        domain = self._reverse_map.get(name, self._reverse_map.get(dir_name, "uncategorized"))

        return SkillEntry(
            name=name,
            description=frontmatter.get("description", ""),
            domain=domain,
            path=path,
            license=frontmatter.get("license", ""),
            metadata=frontmatter.get("metadata", {}),
        )

    @staticmethod
    def _extract_frontmatter(text: str) -> dict | None:
        """Extract YAML frontmatter from text."""
        # Standard: file starts with ---
        match = re.match(r"^---\s*\n(.*?)\n---", text, re.DOTALL)
        if match:
            try:
                return yaml.safe_load(match.group(1))
            except yaml.YAMLError:
                return None

        # Preamble before ---
        match = re.search(r"---\s*\n(name:.*?)\n---", text, re.DOTALL)
        if match:
            try:
                return yaml.safe_load(match.group(1))
            except yaml.YAMLError:
                return None

        return None

    # --- Public API ---

    def search(self, query: str, domain: str | None = None, limit: int = 10) -> list[SkillEntry]:
        """Keyword search across skill names and descriptions."""
        terms = query.lower().split()
        candidates = list(self._entries.values())

        if domain:
            candidates = [e for e in candidates if e.domain == domain]

        scored: list[tuple[float, SkillEntry]] = []
        for entry in candidates:
            searchable = f"{entry.name} {entry.description}".lower()
            score = sum(1 for t in terms if t in searchable)
            if query.lower() == entry.name.lower():
                score += 10
            elif query.lower() in entry.name.lower():
                score += 3
            if score > 0:
                scored.append((score, entry))

        scored.sort(key=lambda x: (-x[0], x[1].name))
        return [entry for _, entry in scored[:limit]]

    def load(self, skill_name: str) -> str | None:
        """Load the full SKILL.md content for a given skill name."""
        entry = self._entries.get(skill_name)
        if not entry:
            return None
        try:
            return entry.path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            return None

    def list_skills(self, domain: str | None = None) -> list[SkillEntry]:
        """List all skills, optionally filtered by domain."""
        entries = list(self._entries.values())
        if domain:
            entries = [e for e in entries if e.domain == domain]
        entries.sort(key=lambda e: e.name)
        return entries

    def get_domains(self) -> list[str]:
        """Return list of all domain names."""
        return sorted(self._domain_map.keys())

    @property
    def size(self) -> int:
        """Total number of indexed skills."""
        return len(self._entries)
