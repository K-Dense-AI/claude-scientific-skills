#!/usr/bin/env python3
"""
Econ AI Skills Sync Script

This script synchronizes economics AI skills from the awesome-econ-ai-stuff
repository (https://github.com/meleantonio/awesome-econ-ai-stuff) to the
local skills directory.

Usage:
    python sync_econ_skills.py [--config CONFIG_PATH] [--dry-run] [--force]
"""

import argparse
import hashlib
import json
import logging
import os
import shutil
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import yaml

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
logger = logging.getLogger(__name__)

# Constants
GITHUB_API_BASE = "https://api.github.com"
GITHUB_RAW_BASE = "https://raw.githubusercontent.com"
DEFAULT_CONFIG_PATH = Path(__file__).parent.parent / "config.yml"


@dataclass
class SyncResult:
    """Result of a skill sync operation."""

    skill_name: str
    category: str
    action: str  # 'added', 'updated', 'skipped', 'failed'
    message: str
    old_hash: Optional[str] = None
    new_hash: Optional[str] = None


class EconSkillsSyncer:
    """Synchronizes economics AI skills from awesome-econ-ai-stuff repository."""

    def __init__(self, config_path: Path, dry_run: bool = False, force: bool = False):
        self.config_path = config_path
        self.dry_run = dry_run
        self.force = force
        self.config = self._load_config()
        self.results: list[SyncResult] = []
        self.github_token = os.environ.get("GITHUB_TOKEN")

    def _load_config(self) -> dict:
        """Load configuration from YAML file."""
        if not self.config_path.exists():
            raise FileNotFoundError(f"Config file not found: {self.config_path}")

        with open(self.config_path) as f:
            return yaml.safe_load(f)

    def _make_request(self, url: str) -> dict | list | bytes:
        """Make HTTP request to GitHub API or raw content."""
        headers = {"Accept": "application/vnd.github.v3+json"}
        if self.github_token:
            headers["Authorization"] = f"token {self.github_token}"

        request = Request(url, headers=headers)

        try:
            with urlopen(request, timeout=30) as response:
                content = response.read()
                content_type = response.headers.get("Content-Type", "")

                if "application/json" in content_type:
                    return json.loads(content.decode("utf-8"))
                return content
        except HTTPError as e:
            if e.code == 404:
                logger.warning(f"Resource not found: {url}")
                return {}
            raise
        except URLError as e:
            logger.error(f"Network error accessing {url}: {e}")
            raise

    def _get_file_hash(self, content: bytes) -> str:
        """Calculate SHA-256 hash of file content."""
        return hashlib.sha256(content).hexdigest()

    def _get_local_skill_hash(self, skill_path: Path) -> Optional[str]:
        """Get hash of locally installed skill."""
        skill_md = skill_path / "SKILL.md"
        if skill_md.exists():
            with open(skill_md, "rb") as f:
                return self._get_file_hash(f.read())
        return None

    def _list_remote_categories(self) -> list[str]:
        """List all skill categories in the remote repository."""
        repo = self.config["source"]["repository"]
        branch = self.config["source"]["branch"]
        skills_path = self.config["source"]["skills_path"]

        url = f"{GITHUB_API_BASE}/repos/{repo}/contents/{skills_path}?ref={branch}"
        contents = self._make_request(url)

        if not contents:
            return []

        categories = []
        for item in contents:
            if item.get("type") == "dir" and not item["name"].startswith("."):
                # Skip template files
                if item["name"] != "SKILL_TEMPLATE.md":
                    categories.append(item["name"])

        return categories

    def _list_skills_in_category(self, category: str) -> list[dict]:
        """List all skills in a given category."""
        repo = self.config["source"]["repository"]
        branch = self.config["source"]["branch"]
        skills_path = self.config["source"]["skills_path"]

        url = f"{GITHUB_API_BASE}/repos/{repo}/contents/{skills_path}/{category}?ref={branch}"
        contents = self._make_request(url)

        if not contents:
            return []

        skills = []
        for item in contents:
            if item.get("type") == "dir":
                skills.append({"name": item["name"], "path": item["path"], "sha": item.get("sha")})

        return skills

    def _download_skill(self, category: str, skill_name: str) -> dict[str, bytes]:
        """Download all files for a skill."""
        repo = self.config["source"]["repository"]
        branch = self.config["source"]["branch"]
        skills_path = self.config["source"]["skills_path"]

        skill_path = f"{skills_path}/{category}/{skill_name}"
        url = f"{GITHUB_API_BASE}/repos/{repo}/contents/{skill_path}?ref={branch}"
        contents = self._make_request(url)

        if not contents:
            return {}

        files = {}
        for item in contents:
            if item.get("type") == "file":
                file_url = f"{GITHUB_RAW_BASE}/{repo}/{branch}/{item['path']}"
                file_content = self._make_request(file_url)
                if isinstance(file_content, bytes):
                    files[item["name"]] = file_content

        return files

    def _install_skill(
        self, category: str, skill_name: str, files: dict[str, bytes]
    ) -> SyncResult:
        """Install a skill to the local directory."""
        target_base = Path(self.config["install"]["target_dir"])
        skill_dir = target_base / category / skill_name

        # Check if skill already exists
        old_hash = self._get_local_skill_hash(skill_dir)
        new_hash = self._get_file_hash(files.get("SKILL.md", b""))

        if old_hash and old_hash == new_hash and not self.force:
            return SyncResult(
                skill_name=skill_name,
                category=category,
                action="skipped",
                message="Skill is up to date",
                old_hash=old_hash,
                new_hash=new_hash,
            )

        if self.dry_run:
            action = "would_update" if old_hash else "would_add"
            return SyncResult(
                skill_name=skill_name,
                category=category,
                action=action,
                message=f"[DRY RUN] Would {'update' if old_hash else 'add'} skill",
                old_hash=old_hash,
                new_hash=new_hash,
            )

        # Backup existing skill if configured
        if old_hash and self.config["install"].get("backup_existing", True):
            backup_dir = target_base / ".backups" / f"{category}_{skill_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
            if skill_dir.exists():
                shutil.copytree(skill_dir, backup_dir)
                logger.info(f"Backed up existing skill to {backup_dir}")

        # Create skill directory
        skill_dir.mkdir(parents=True, exist_ok=True)

        # Write files
        for filename, content in files.items():
            file_path = skill_dir / filename
            with open(file_path, "wb") as f:
                f.write(content)

        action = "updated" if old_hash else "added"
        return SyncResult(
            skill_name=skill_name,
            category=category,
            action=action,
            message=f"Successfully {action} skill",
            old_hash=old_hash,
            new_hash=new_hash,
        )

    def _should_sync_category(self, category: str) -> bool:
        """Check if a category should be synced based on config."""
        categories_config = self.config.get("categories", {})
        return categories_config.get(category, True)

    def _should_sync_skill(self, skill_name: str) -> bool:
        """Check if a skill should be synced based on include/exclude lists."""
        exclude = self.config.get("exclude_skills", [])
        if skill_name in exclude:
            return False

        include = self.config.get("include_skills", [])
        if include and skill_name not in include:
            return False

        return True

    def sync(self) -> list[SyncResult]:
        """Perform the full sync operation."""
        logger.info("Starting Econ AI Skills sync...")
        logger.info(f"Source: {self.config['source']['repository']}")
        logger.info(f"Dry run: {self.dry_run}")

        # Get all categories
        categories = self._list_remote_categories()
        logger.info(f"Found {len(categories)} categories: {categories}")

        for category in categories:
            if not self._should_sync_category(category):
                logger.info(f"Skipping disabled category: {category}")
                continue

            logger.info(f"Processing category: {category}")
            skills = self._list_skills_in_category(category)

            for skill_info in skills:
                skill_name = skill_info["name"]

                if not self._should_sync_skill(skill_name):
                    logger.info(f"Skipping excluded skill: {skill_name}")
                    continue

                logger.info(f"Syncing skill: {category}/{skill_name}")

                try:
                    files = self._download_skill(category, skill_name)
                    if not files:
                        result = SyncResult(
                            skill_name=skill_name,
                            category=category,
                            action="failed",
                            message="No files found in skill directory",
                        )
                    else:
                        result = self._install_skill(category, skill_name, files)

                    self.results.append(result)
                    logger.info(f"  {result.action}: {result.message}")

                except Exception as e:
                    result = SyncResult(
                        skill_name=skill_name,
                        category=category,
                        action="failed",
                        message=str(e),
                    )
                    self.results.append(result)
                    logger.error(f"  Failed: {e}")

        return self.results

    def generate_report(self) -> str:
        """Generate a summary report of the sync operation."""
        added = [r for r in self.results if r.action == "added"]
        updated = [r for r in self.results if r.action == "updated"]
        skipped = [r for r in self.results if r.action == "skipped"]
        failed = [r for r in self.results if r.action == "failed"]

        report = []
        report.append("# Econ AI Skills Sync Report")
        report.append(f"\nGenerated: {datetime.now().isoformat()}")
        report.append(f"Source: {self.config['source']['repository']}")
        report.append(f"Dry run: {self.dry_run}")
        report.append("\n## Summary")
        report.append(f"- Added: {len(added)}")
        report.append(f"- Updated: {len(updated)}")
        report.append(f"- Skipped (up to date): {len(skipped)}")
        report.append(f"- Failed: {len(failed)}")

        if added:
            report.append("\n## Added Skills")
            for r in added:
                report.append(f"- {r.category}/{r.skill_name}")

        if updated:
            report.append("\n## Updated Skills")
            for r in updated:
                report.append(f"- {r.category}/{r.skill_name}")

        if failed:
            report.append("\n## Failed Skills")
            for r in failed:
                report.append(f"- {r.category}/{r.skill_name}: {r.message}")

        return "\n".join(report)

    def save_manifest(self, manifest_path: Path):
        """Save a manifest of installed skills."""
        manifest = {
            "last_sync": datetime.now().isoformat(),
            "source": self.config["source"],
            "skills": {},
        }

        for result in self.results:
            if result.action in ("added", "updated", "skipped"):
                category = result.category
                if category not in manifest["skills"]:
                    manifest["skills"][category] = {}
                manifest["skills"][category][result.skill_name] = {
                    "hash": result.new_hash,
                    "last_action": result.action,
                    "last_sync": datetime.now().isoformat(),
                }

        with open(manifest_path, "w") as f:
            json.dump(manifest, f, indent=2)


def main():
    parser = argparse.ArgumentParser(
        description="Sync economics AI skills from awesome-econ-ai-stuff repository"
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG_PATH,
        help="Path to configuration file",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without making changes",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force update even if skills are up to date",
    )
    parser.add_argument(
        "--output-report",
        type=Path,
        help="Path to save the sync report",
    )
    parser.add_argument(
        "--output-manifest",
        type=Path,
        help="Path to save the skills manifest",
    )

    args = parser.parse_args()

    try:
        syncer = EconSkillsSyncer(
            config_path=args.config,
            dry_run=args.dry_run,
            force=args.force,
        )

        results = syncer.sync()
        report = syncer.generate_report()

        print("\n" + report)

        if args.output_report:
            with open(args.output_report, "w") as f:
                f.write(report)
            logger.info(f"Report saved to {args.output_report}")

        if args.output_manifest:
            syncer.save_manifest(args.output_manifest)
            logger.info(f"Manifest saved to {args.output_manifest}")

        # Exit with error if any skills failed
        failed_count = len([r for r in results if r.action == "failed"])
        if failed_count > 0:
            sys.exit(1)

    except Exception as e:
        logger.error(f"Sync failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
