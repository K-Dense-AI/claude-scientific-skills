"""Command-line interface for systematic review pipeline."""

import logging
from pathlib import Path

import click

from .config import get_config, get_settings
from .deduplication import Deduplicator, generate_dedup_report
from .exporters import export_to_csv, export_to_ris
from .extraction import (
    create_extraction_template,
    load_extraction_csv,
    save_extraction_json,
    validate_extraction_record,
)
from .importers import import_citations
from .meta import MetaAnalyzer, create_forest_plot, create_funnel_plot
from .meta.plots import create_sensitivity_plot
from .prisma import PRISMAFlowData, generate_prisma_diagram, generate_prisma_table
from .retrieval import PubMedRetriever, save_citations_jsonl, load_citations_jsonl
from .risk_of_bias import (
    create_nos_template,
    create_rob2_template,
    load_nos_assessments,
    load_rob2_assessments,
    generate_nos_summary_table,
    generate_rob2_summary_table,
)
from .screening.database import ScreeningDatabase

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


@click.group()
@click.option("--config", default="config.yaml", help="Path to configuration file")
@click.pass_context
def cli(ctx: click.Context, config: str) -> None:
    """Systematic Review Pipeline for Postoperative Delirium."""
    ctx.ensure_object(dict)
    ctx.obj["config_path"] = config
    ctx.obj["config"] = get_config(config)
    ctx.obj["settings"] = get_settings()


@cli.command()
@click.option("--max-results", default=10000, help="Maximum results to retrieve")
@click.option("--output-dir", default="data/raw", help="Output directory")
@click.pass_context
def retrieve(ctx: click.Context, max_results: int, output_dir: str) -> None:
    """Retrieve citations from PubMed."""
    config = ctx.obj["config"]
    settings = ctx.obj["settings"]

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Get PubMed query
    query = config.get_pubmed_query()
    logger.info(f"PubMed Query: {query}")

    # Retrieve
    retriever = PubMedRetriever(settings)
    citations = retriever.search_and_fetch(query, max_results)

    # Save
    jsonl_path = output_path / "pubmed_records.jsonl"
    save_citations_jsonl(citations, jsonl_path)

    # Export
    export_to_ris(citations, output_path / "pubmed_records.ris")
    export_to_csv(citations, output_path / "pubmed_records.csv")

    click.echo(f"✓ Retrieved {len(citations)} citations from PubMed")
    click.echo(f"✓ Saved to {output_path}")


@cli.command()
@click.argument("file_path", type=click.Path(exists=True))
@click.option("--format", type=click.Choice(["ris", "bibtex", "csv"]), help="File format")
@click.option("--output-dir", default="data/imported", help="Output directory")
@click.pass_context
def import_records(ctx: click.Context, file_path: str, format: str | None, output_dir: str) -> None:
    """Import citations from RIS, BibTeX, or CSV file."""
    file_path_obj = Path(file_path)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Import
    citations = import_citations(file_path_obj, format)

    # Save
    db_name = file_path_obj.stem
    jsonl_path = output_path / f"{db_name}_imported.jsonl"
    save_citations_jsonl(citations, jsonl_path)

    click.echo(f"✓ Imported {len(citations)} citations from {file_path}")
    click.echo(f"✓ Saved to {jsonl_path}")


@cli.command()
@click.option("--input-dir", default="data", help="Input directory containing JSONL files")
@click.option("--output-dir", default="data/processed", help="Output directory")
@click.pass_context
def dedup(ctx: click.Context, input_dir: str, output_dir: str) -> None:
    """Deduplicate citations."""
    config = ctx.obj["config"]

    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load all JSONL files
    all_citations = []
    for jsonl_file in input_path.rglob("*.jsonl"):
        citations = load_citations_jsonl(jsonl_file)
        all_citations.extend(citations)
        logger.info(f"Loaded {len(citations)} from {jsonl_file}")

    click.echo(f"Total citations before deduplication: {len(all_citations)}")

    # Deduplicate
    dedup_config = config.get("deduplication", {})
    deduplicator = Deduplicator(
        title_threshold=dedup_config.get("fuzzy_match", {}).get("title_threshold", 0.90),
        author_threshold=dedup_config.get("fuzzy_match", {}).get("author_threshold", 0.85),
    )

    unique_citations, matches = deduplicator.deduplicate(all_citations)

    # Save
    save_citations_jsonl(unique_citations, output_path / "deduplicated.jsonl")
    export_to_ris(unique_citations, output_path / "deduplicated.ris")
    export_to_csv(unique_citations, output_path / "deduplicated.csv")

    # Report
    report = generate_dedup_report(len(all_citations), len(unique_citations), matches)

    click.echo(f"✓ Unique citations: {len(unique_citations)}")
    click.echo(f"✓ Duplicates removed: {report['duplicates_removed']}")
    click.echo(f"✓ Duplicate rate: {report['duplicate_rate']:.1%}")
    click.echo(f"✓ Matches by type: {report['matches_by_type']}")


@cli.command()
@click.pass_context
def screen(ctx: click.Context) -> None:
    """Launch screening application."""
    config_path = ctx.obj["config_path"]

    click.echo("Launching Streamlit screening app...")
    click.echo("Open your browser to http://localhost:8501")

    from .screening import run_screening_app

    run_screening_app(config_path)


@cli.command()
@click.option("--input-file", required=True, help="Deduplicated citations JSONL file")
@click.pass_context
def load_screening(ctx: click.Context, input_file: str) -> None:
    """Load citations into screening database."""
    settings = ctx.obj["settings"]

    citations = load_citations_jsonl(Path(input_file))

    db = ScreeningDatabase(settings.screening_db_path)
    db.add_citations(citations)

    click.echo(f"✓ Loaded {len(citations)} citations into screening database")


@cli.command()
@click.option("--output-dir", default="outputs/prisma", help="Output directory")
@click.pass_context
def prisma(ctx: click.Context, output_dir: str) -> None:
    """Generate PRISMA flow diagram and table."""
    settings = ctx.obj["settings"]
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Get screening stats
    db = ScreeningDatabase(settings.screening_db_path)
    stats = db.get_screening_stats()

    # Build PRISMA flow data
    flow_data = PRISMAFlowData()
    # Populate from stats (simplified for now)
    flow_data.records_after_dedup = stats.get("total_citations", 0)

    # Generate diagram and table
    generate_prisma_diagram(flow_data, output_path / "prisma_flow.png")
    generate_prisma_table(flow_data, output_path / "prisma_counts.csv")

    click.echo(f"✓ PRISMA flow diagram saved to {output_path}")


@cli.command()
@click.option("--output-dir", default="outputs/extraction", help="Output directory")
@click.pass_context
def create_templates(ctx: click.Context, output_dir: str) -> None:
    """Create extraction and risk of bias templates."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    create_extraction_template(output_path / "extraction_template.csv")
    create_nos_template(output_path / "nos_template.csv")
    create_rob2_template(output_path / "rob2_template.csv")

    click.echo(f"✓ Templates created in {output_path}")


@cli.command()
@click.argument("extraction_file", type=click.Path(exists=True))
@click.argument("risk_factor")
@click.option("--output-dir", default="outputs/meta", help="Output directory")
@click.option("--effect-type", type=click.Choice(["OR", "RR", "HR"]), help="Filter by effect type")
@click.pass_context
def meta_analyze(
    ctx: click.Context, extraction_file: str, risk_factor: str, output_dir: str, effect_type: str | None
) -> None:
    """Perform meta-analysis for a risk factor."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load extraction data
    extraction_records = load_extraction_csv(Path(extraction_file))

    # Run meta-analysis
    analyzer = MetaAnalyzer()
    result = analyzer.analyze(extraction_records, risk_factor, effect_type)

    if result is None:
        click.echo(f"✗ Insufficient data for {risk_factor}")
        return

    # Create plots
    create_forest_plot(result, output_path / f"{risk_factor}_forest.png")

    if result.n_studies >= 10:
        create_funnel_plot(result, output_path / f"{risk_factor}_funnel.png")

    # Sensitivity analysis
    sensitivity = analyzer.leave_one_out_analysis(extraction_records, risk_factor)
    if sensitivity:
        create_sensitivity_plot(sensitivity, output_path / f"{risk_factor}_sensitivity.png")

    # Print results
    click.echo(f"\n✓ Meta-Analysis Results: {risk_factor}")
    click.echo(f"  Studies: {result.n_studies}")
    click.echo(f"  Pooled OR: {result.pooled_effect:.3f} (exp: {np.exp(result.pooled_effect):.3f})")
    click.echo(f"  95% CI: [{result.pooled_ci_lower:.3f}, {result.pooled_ci_upper:.3f}]")
    click.echo(f"  P-value: {result.pooled_p_value:.4f}")
    click.echo(f"  I²: {result.i_squared:.1f}%")
    click.echo(f"  τ²: {result.tau_squared:.4f}")
    if result.egger_p_value:
        click.echo(f"  Egger's test: p = {result.egger_p_value:.4f}")


@cli.command()
@click.option("--nos-file", type=click.Path(exists=True), help="NOS assessments CSV")
@click.option("--rob2-file", type=click.Path(exists=True), help="RoB2 assessments CSV")
@click.option("--output-dir", default="outputs/rob", help="Output directory")
def rob_summary(nos_file: str | None, rob2_file: str | None, output_dir: str) -> None:
    """Generate risk of bias summary tables."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    if nos_file:
        assessments = load_nos_assessments(Path(nos_file))
        generate_nos_summary_table(assessments, output_path / "nos_summary.csv")
        click.echo(f"✓ NOS summary saved")

    if rob2_file:
        assessments = load_rob2_assessments(Path(rob2_file))
        generate_rob2_summary_table(assessments, output_path / "rob2_summary.csv")
        click.echo(f"✓ RoB2 summary saved")


def main() -> None:
    """Main entry point."""
    import numpy as np  # Import for meta_analyze command

    cli(obj={})


if __name__ == "__main__":
    main()
