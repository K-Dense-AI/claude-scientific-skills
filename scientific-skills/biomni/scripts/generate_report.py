#!/usr/bin/env python3
"""
Enhanced PDF report generation for biomni conversation histories.

This script provides additional customization options for biomni reports:
- Custom styling and branding
- Formatted code blocks
- Section organization
- Metadata inclusion
- Export format options (PDF, HTML, Markdown)

Usage:
    python generate_report.py --input conversation.json --output report.pdf
    python generate_report.py --agent-object agent --output report.pdf --format html
"""

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Any
from datetime import datetime


def format_conversation_history(
    messages: List[Dict[str, Any]],
    include_metadata: bool = True,
    include_code: bool = True,
    include_timestamps: bool = False
) -> str:
    """
    Format conversation history into structured markdown.

    Args:
        messages: List of conversation message dictionaries
        include_metadata: Include metadata section
        include_code: Include code blocks
        include_timestamps: Include message timestamps

    Returns:
        Formatted markdown string
    """
    sections = []

    # Header
    sections.append("# Biomni Analysis Report\n")

    # Metadata
    if include_metadata:
        sections.append("## Metadata\n")
        sections.append(f"- **Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        sections.append(f"- **Number of interactions**: {len(messages)}")
        sections.append("\n---\n")

    # Process messages
    sections.append("## Analysis\n")

    for i, msg in enumerate(messages, 1):
        role = msg.get('role', 'unknown')
        content = msg.get('content', '')

        if role == 'user':
            sections.append(f"### Task {i // 2 + 1}\n")
            sections.append(f"**Query:**\n```\n{content}\n```\n")

        elif role == 'assistant':
            sections.append(f"**Response:**\n")

            # Check if content contains code
            if include_code and ('```' in content or 'import ' in content):
                # Attempt to separate text and code
                parts = content.split('```')
                for j, part in enumerate(parts):
                    if j % 2 == 0:
                        # Text content
                        if part.strip():
                            sections.append(f"{part.strip()}\n")
                    else:
                        # Code content
                        # Check if language is specified
                        lines = part.split('\n', 1)
                        if len(lines) > 1 and lines[0].strip() in ['python', 'r', 'bash', 'sql']:
                            lang = lines[0].strip()
                            code = lines[1]
                        else:
                            lang = 'python'  # Default to python
                            code = part

                        sections.append(f"```{lang}\n{code}\n```\n")
            else:
                sections.append(f"{content}\n")

            sections.append("\n---\n")

    return '\n'.join(sections)


def markdown_to_html(markdown_content: str, title: str = "Biomni Report") -> str:
    """
    Convert markdown to styled HTML.

    Args:
        markdown_content: Markdown string
        title: HTML page title

    Returns:
        HTML string
    """
    # Simple markdown to HTML conversion
    # For production use, consider using a library like markdown or mistune

    html_template = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
            line-height: 1.6;
            max-width: 900px;
            margin: 0 auto;
            padding: 20px;
            color: #333;
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 30px;
            border-bottom: 2px solid #95a5a6;
            padding-bottom: 5px;
        }}
        h3 {{
            color: #555;
        }}
        code {{
            background-color: #f4f4f4;
            padding: 2px 6px;
            border-radius: 3px;
            font-family: 'Monaco', 'Menlo', 'Courier New', monospace;
        }}
        pre {{
            background-color: #f8f8f8;
            border: 1px solid #ddd;
            border-radius: 5px;
            padding: 15px;
            overflow-x: auto;
        }}
        pre code {{
            background-color: transparent;
            padding: 0;
        }}
        hr {{
            border: none;
            border-top: 1px solid #ddd;
            margin: 30px 0;
        }}
        .metadata {{
            background-color: #ecf0f1;
            padding: 15px;
            border-radius: 5px;
            margin-bottom: 20px;
        }}
        .task {{
            background-color: #e8f4f8;
            padding: 10px;
            border-left: 4px solid #3498db;
            margin: 20px 0;
        }}
        .footer {{
            margin-top: 50px;
            text-align: center;
            color: #7f8c8d;
            font-size: 0.9em;
        }}
    </style>
</head>
<body>
    <div class="content">
        {markdown_to_html_simple(markdown_content)}
    </div>
    <div class="footer">
        <p>Generated with Biomni | Stanford SNAP Lab</p>
        <p><a href="https://github.com/snap-stanford/biomni">github.com/snap-stanford/biomni</a></p>
    </div>
</body>
</html>
"""
    return html_template


def markdown_to_html_simple(md: str) -> str:
    """Simple markdown to HTML converter (basic implementation)."""
    lines = md.split('\n')
    html_lines = []
    in_code_block = False
    in_list = False

    for line in lines:
        # Code blocks
        if line.startswith('```'):
            if in_code_block:
                html_lines.append('</code></pre>')
                in_code_block = False
            else:
                lang = line[3:].strip()
                html_lines.append(f'<pre><code class="language-{lang}">')
                in_code_block = True
            continue

        if in_code_block:
            html_lines.append(line)
            continue

        # Headers
        if line.startswith('# '):
            html_lines.append(f'<h1>{line[2:]}</h1>')
        elif line.startswith('## '):
            html_lines.append(f'<h2>{line[3:]}</h2>')
        elif line.startswith('### '):
            html_lines.append(f'<h3>{line[4:]}</h3>')
        # Lists
        elif line.startswith('- '):
            if not in_list:
                html_lines.append('<ul>')
                in_list = True
            html_lines.append(f'<li>{line[2:]}</li>')
        else:
            if in_list:
                html_lines.append('</ul>')
                in_list = False

            # Horizontal rule
            if line.strip() == '---':
                html_lines.append('<hr>')
            # Bold
            elif '**' in line:
                line = line.replace('**', '<strong>', 1).replace('**', '</strong>', 1)
                html_lines.append(f'<p>{line}</p>')
            # Regular paragraph
            elif line.strip():
                html_lines.append(f'<p>{line}</p>')
            else:
                html_lines.append('<br>')

    if in_list:
        html_lines.append('</ul>')

    return '\n'.join(html_lines)


def generate_report(
    conversation_data: Dict[str, Any],
    output_path: Path,
    format: str = 'markdown',
    title: Optional[str] = None
):
    """
    Generate formatted report from conversation data.

    Args:
        conversation_data: Conversation history dictionary
        output_path: Output file path
        format: Output format ('markdown', 'html', or 'pdf')
        title: Report title
    """
    messages = conversation_data.get('messages', [])

    if not title:
        title = f"Biomni Analysis - {datetime.now().strftime('%Y-%m-%d')}"

    # Generate markdown
    markdown_content = format_conversation_history(messages)

    if format == 'markdown':
        output_path.write_text(markdown_content)
        print(f"✓ Markdown report saved to {output_path}")

    elif format == 'html':
        html_content = markdown_to_html(markdown_content, title)
        output_path.write_text(html_content)
        print(f"✓ HTML report saved to {output_path}")

    elif format == 'pdf':
        # For PDF generation, we'd typically use a library like weasyprint or reportlab
        # This is a placeholder implementation
        print("PDF generation requires additional dependencies (weasyprint or reportlab)")
        print("Falling back to HTML format...")

        html_path = output_path.with_suffix('.html')
        html_content = markdown_to_html(markdown_content, title)
        html_path.write_text(html_content)

        print(f"✓ HTML report saved to {html_path}")
        print("  To convert to PDF:")
        print(f"    1. Install weasyprint: pip install weasyprint")
        print(f"    2. Run: weasyprint {html_path} {output_path}")

    else:
        raise ValueError(f"Unsupported format: {format}")


def main():
    """Main entry point for CLI usage."""
    parser = argparse.ArgumentParser(
        description="Generate enhanced reports from biomni conversation histories"
    )

    parser.add_argument(
        '--input',
        type=Path,
        required=True,
        help='Input conversation history JSON file'
    )

    parser.add_argument(
        '--output',
        type=Path,
        required=True,
        help='Output report file path'
    )

    parser.add_argument(
        '--format',
        choices=['markdown', 'html', 'pdf'],
        default='markdown',
        help='Output format (default: markdown)'
    )

    parser.add_argument(
        '--title',
        type=str,
        help='Report title (optional)'
    )

    args = parser.parse_args()

    # Load conversation data
    try:
        with open(args.input, 'r') as f:
            conversation_data = json.load(f)
    except FileNotFoundError:
        print(f"❌ Input file not found: {args.input}")
        return 1
    except json.JSONDecodeError:
        print(f"❌ Invalid JSON in input file: {args.input}")
        return 1

    # Generate report
    try:
        generate_report(
            conversation_data,
            args.output,
            format=args.format,
            title=args.title
        )
        return 0
    except Exception as e:
        print(f"❌ Error generating report: {e}")
        return 1


if __name__ == '__main__':
    import sys
    sys.exit(main())
