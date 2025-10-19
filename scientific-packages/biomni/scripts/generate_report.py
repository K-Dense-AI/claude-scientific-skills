#!/usr/bin/env python3
"""
Enhanced PDF Report Generation for Biomni

This script provides advanced PDF report generation with custom formatting,
styling, and metadata for Biomni analysis results.
"""

import argparse
import sys
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, Any


def generate_markdown_report(
    title: str,
    sections: list,
    metadata: Optional[Dict[str, Any]] = None,
    output_path: str = "report.md"
) -> str:
    """
    Generate a formatted markdown report.

    Args:
        title: Report title
        sections: List of dicts with 'heading' and 'content' keys
        metadata: Optional metadata dict (author, date, etc.)
        output_path: Path to save markdown file

    Returns:
        Path to generated markdown file
    """
    md_content = []

    # Title
    md_content.append(f"# {title}\n")

    # Metadata
    if metadata:
        md_content.append("---\n")
        for key, value in metadata.items():
            md_content.append(f"**{key}:** {value}  \n")
        md_content.append("---\n\n")

    # Sections
    for section in sections:
        heading = section.get('heading', 'Section')
        content = section.get('content', '')
        level = section.get('level', 2)  # Default to h2

        md_content.append(f"{'#' * level} {heading}\n\n")
        md_content.append(f"{content}\n\n")

    # Write to file
    output = Path(output_path)
    output.write_text('\n'.join(md_content))

    return str(output)


def convert_to_pdf_weasyprint(
    markdown_path: str,
    output_path: str,
    css_style: Optional[str] = None
) -> bool:
    """
    Convert markdown to PDF using WeasyPrint.

    Args:
        markdown_path: Path to markdown file
        output_path: Path for output PDF
        css_style: Optional CSS stylesheet path

    Returns:
        True if successful, False otherwise
    """
    try:
        import markdown
        from weasyprint import HTML, CSS

        # Read markdown
        with open(markdown_path, 'r') as f:
            md_content = f.read()

        # Convert to HTML
        html_content = markdown.markdown(
            md_content,
            extensions=['tables', 'fenced_code', 'codehilite']
        )

        # Wrap in HTML template
        html_template = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <meta charset="utf-8">
            <title>Biomni Report</title>
            <style>
                body {{
                    font-family: 'Helvetica', 'Arial', sans-serif;
                    line-height: 1.6;
                    color: #333;
                    max-width: 800px;
                    margin: 40px auto;
                    padding: 20px;
                }}
                h1 {{
                    color: #2c3e50;
                    border-bottom: 3px solid #3498db;
                    padding-bottom: 10px;
                }}
                h2 {{
                    color: #34495e;
                    margin-top: 30px;
                    border-bottom: 1px solid #bdc3c7;
                    padding-bottom: 5px;
                }}
                h3 {{
                    color: #7f8c8d;
                }}
                code {{
                    background-color: #f4f4f4;
                    padding: 2px 6px;
                    border-radius: 3px;
                    font-family: 'Courier New', monospace;
                }}
                pre {{
                    background-color: #f4f4f4;
                    padding: 15px;
                    border-radius: 5px;
                    overflow-x: auto;
                }}
                table {{
                    border-collapse: collapse;
                    width: 100%;
                    margin: 20px 0;
                }}
                th, td {{
                    border: 1px solid #ddd;
                    padding: 12px;
                    text-align: left;
                }}
                th {{
                    background-color: #3498db;
                    color: white;
                }}
                tr:nth-child(even) {{
                    background-color: #f9f9f9;
                }}
                .metadata {{
                    background-color: #ecf0f1;
                    padding: 15px;
                    border-radius: 5px;
                    margin: 20px 0;
                }}
            </style>
        </head>
        <body>
            {html_content}
        </body>
        </html>
        """

        # Generate PDF
        pdf = HTML(string=html_template)

        # Add custom CSS if provided
        stylesheets = []
        if css_style and Path(css_style).exists():
            stylesheets.append(CSS(filename=css_style))

        pdf.write_pdf(output_path, stylesheets=stylesheets)

        return True

    except ImportError:
        print("Error: WeasyPrint not installed. Install with: pip install weasyprint")
        return False
    except Exception as e:
        print(f"Error generating PDF: {e}")
        return False


def convert_to_pdf_pandoc(markdown_path: str, output_path: str) -> bool:
    """
    Convert markdown to PDF using Pandoc.

    Args:
        markdown_path: Path to markdown file
        output_path: Path for output PDF

    Returns:
        True if successful, False otherwise
    """
    try:
        import subprocess

        # Check if pandoc is installed
        result = subprocess.run(
            ['pandoc', '--version'],
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            print("Error: Pandoc not installed")
            return False

        # Convert with pandoc
        result = subprocess.run(
            [
                'pandoc',
                markdown_path,
                '-o', output_path,
                '--pdf-engine=pdflatex',
                '-V', 'geometry:margin=1in',
                '--toc'
            ],
            capture_output=True,
            text=True
        )

        if result.returncode != 0:
            print(f"Pandoc error: {result.stderr}")
            return False

        return True

    except FileNotFoundError:
        print("Error: Pandoc not found. Install from https://pandoc.org/")
        return False
    except Exception as e:
        print(f"Error: {e}")
        return False


def create_biomni_report(
    conversation_history: list,
    output_path: str = "biomni_report.pdf",
    method: str = "weasyprint"
) -> bool:
    """
    Create a formatted PDF report from Biomni conversation history.

    Args:
        conversation_history: List of conversation turns
        output_path: Output PDF path
        method: Conversion method ('weasyprint' or 'pandoc')

    Returns:
        True if successful
    """
    # Prepare report sections
    metadata = {
        'Date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'Tool': 'Biomni AI Agent',
        'Report Type': 'Analysis Summary'
    }

    sections = []

    # Executive Summary
    sections.append({
        'heading': 'Executive Summary',
        'level': 2,
        'content': 'This report contains the complete analysis workflow executed by the Biomni biomedical AI agent.'
    })

    # Conversation history
    for i, turn in enumerate(conversation_history, 1):
        sections.append({
            'heading': f'Task {i}: {turn.get("task", "Analysis")}',
            'level': 2,
            'content': f'**Input:**\n```\n{turn.get("input", "")}\n```\n\n**Output:**\n{turn.get("output", "")}'
        })

    # Generate markdown
    md_path = output_path.replace('.pdf', '.md')
    generate_markdown_report(
        title="Biomni Analysis Report",
        sections=sections,
        metadata=metadata,
        output_path=md_path
    )

    # Convert to PDF
    if method == 'weasyprint':
        success = convert_to_pdf_weasyprint(md_path, output_path)
    elif method == 'pandoc':
        success = convert_to_pdf_pandoc(md_path, output_path)
    else:
        print(f"Unknown method: {method}")
        return False

    if success:
        print(f"✓ Report generated: {output_path}")
        print(f"  Markdown: {md_path}")
    else:
        print("✗ Failed to generate PDF")
        print(f"  Markdown available: {md_path}")

    return success


def main():
    """CLI for report generation."""
    parser = argparse.ArgumentParser(
        description='Generate formatted PDF reports for Biomni analyses'
    )

    parser.add_argument(
        'input',
        type=str,
        help='Input markdown file or conversation history'
    )

    parser.add_argument(
        '-o', '--output',
        type=str,
        default='biomni_report.pdf',
        help='Output PDF path (default: biomni_report.pdf)'
    )

    parser.add_argument(
        '-m', '--method',
        type=str,
        choices=['weasyprint', 'pandoc'],
        default='weasyprint',
        help='Conversion method (default: weasyprint)'
    )

    parser.add_argument(
        '--css',
        type=str,
        help='Custom CSS stylesheet path'
    )

    args = parser.parse_args()

    # Check if input is markdown or conversation history
    input_path = Path(args.input)

    if not input_path.exists():
        print(f"Error: Input file not found: {args.input}")
        return 1

    # If input is markdown, convert directly
    if input_path.suffix == '.md':
        if args.method == 'weasyprint':
            success = convert_to_pdf_weasyprint(
                str(input_path),
                args.output,
                args.css
            )
        else:
            success = convert_to_pdf_pandoc(str(input_path), args.output)

        return 0 if success else 1

    # Otherwise, assume it's conversation history (JSON)
    try:
        import json
        with open(input_path) as f:
            history = json.load(f)

        success = create_biomni_report(
            history,
            args.output,
            args.method
        )

        return 0 if success else 1

    except json.JSONDecodeError:
        print("Error: Input file is not valid JSON or markdown")
        return 1


if __name__ == "__main__":
    sys.exit(main())
