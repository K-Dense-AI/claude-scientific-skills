"""Configuration management for the systematic review pipeline."""

import os
from pathlib import Path
from typing import Any

import yaml
from pydantic import BaseModel, Field
from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    """Environment settings loaded from .env file."""

    ncbi_api_key: str = Field(default="", alias="NCBI_API_KEY")
    ncbi_email: str = Field(default="", alias="NCBI_EMAIL")
    epmc_email: str = Field(default="", alias="EPMC_EMAIL")
    s2_api_key: str = Field(default="", alias="S2_API_KEY")
    screening_db_path: str = Field(default="data/screening.db", alias="SCREENING_DB_PATH")
    extraction_db_path: str = Field(default="data/extraction.db", alias="EXTRACTION_DB_PATH")
    output_dir: str = Field(default="outputs", alias="OUTPUT_DIR")
    data_dir: str = Field(default="data", alias="DATA_DIR")

    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"
        extra = "ignore"


class ReviewConfig:
    """Loads and provides access to review configuration from YAML."""

    def __init__(self, config_path: str | Path = "config.yaml") -> None:
        """
        Initialize configuration loader.

        Args:
            config_path: Path to YAML configuration file
        """
        self.config_path = Path(config_path)
        if not self.config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_path}")

        with open(self.config_path, "r", encoding="utf-8") as f:
            self._config: dict[str, Any] = yaml.safe_load(f)

    def get(self, key: str, default: Any = None) -> Any:
        """
        Get configuration value by dot-notation key.

        Args:
            key: Configuration key (e.g., 'review.title')
            default: Default value if key not found

        Returns:
            Configuration value
        """
        keys = key.split(".")
        value = self._config
        for k in keys:
            if isinstance(value, dict):
                value = value.get(k, default)
            else:
                return default
        return value

    @property
    def review_title(self) -> str:
        """Get review title."""
        return self.get("review.title", "")

    @property
    def reviewers(self) -> list[dict[str, str]]:
        """Get list of reviewers."""
        return self.get("reviewers", [])

    @property
    def inclusion_criteria(self) -> list[str]:
        """Get inclusion criteria."""
        return self.get("inclusion_criteria", [])

    @property
    def exclusion_criteria(self) -> list[str]:
        """Get exclusion criteria."""
        return self.get("exclusion_criteria", [])

    @property
    def search_strategies(self) -> dict[str, Any]:
        """Get search strategies for all databases."""
        return self.get("search_strategies", {})

    @property
    def exclusion_reasons(self) -> list[dict[str, str]]:
        """Get exclusion reason codes and descriptions."""
        return self.get("screening.exclusion_reasons", [])

    @property
    def extraction_schema(self) -> dict[str, Any]:
        """Get data extraction schema."""
        return self.get("extraction", {})

    @property
    def meta_analysis_settings(self) -> dict[str, Any]:
        """Get meta-analysis settings."""
        return self.get("meta_analysis", {})

    def get_pubmed_query(self) -> str:
        """Get PubMed search query."""
        return self.get("search_strategies.pubmed.query", "").strip()

    def get_date_range(self, database: str = "pubmed") -> str:
        """Get date range for database search."""
        return self.get(f"search_strategies.{database}.date_range", "")


def get_settings() -> Settings:
    """Get application settings from environment."""
    return Settings()


def get_config(config_path: str | Path = "config.yaml") -> ReviewConfig:
    """
    Load review configuration.

    Args:
        config_path: Path to configuration file

    Returns:
        ReviewConfig instance
    """
    return ReviewConfig(config_path)
