#!/usr/bin/env python3
"""
Quick Document Generator - Helper for creating simple ReportLab documents

This script provides utility functions for quickly creating common document types
without writing boilerplate code.
"""

from reportlab.lib.pagesizes import letter, A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, PageBreak,
    Table, TableStyle, Image, KeepTogether
)
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_RIGHT, TA_JUSTIFY
from datetime import datetime


def create_simple_document(filename, title, author="", content_blocks=None, pagesize=letter):
    """
    Create a simple document with title and content blocks.

    Args:
        filename: Output PDF filename
        title: Document title
        author: Document author (optional)
        content_blocks: List of dicts with 'type' and 'content' keys
                       type can be: 'heading', 'paragraph', 'bullet', 'space'
        pagesize: Page size (default: letter)

    Example content_blocks:
        [
            {'type': 'heading', 'content': 'Introduction'},
            {'type': 'paragraph', 'content': 'This is a paragraph.'},
            {'type': 'bullet', 'content': 'Bullet point item'},
            {'type': 'space', 'height': 0.2},  # height in inches
        ]
    """
    if content_blocks is None:
        content_blocks = []

    # Create document
    doc = SimpleDocTemplate(
        filename,
        pagesize=pagesize,
        rightMargin=72,
        leftMargin=72,
        topMargin=72,
        bottomMargin=18,
        title=title,
        author=author
    )

    # Get styles
    styles = getSampleStyleSheet()
    story = []

    # Add title
    story.append(Paragraph(title, styles['Title']))
    story.append(Spacer(1, 0.3*inch))

    # Process content blocks
    for block in content_blocks:
        block_type = block.get('type', 'paragraph')
        content = block.get('content', '')

        if block_type == 'heading':
            story.append(Paragraph(content, styles['Heading1']))
            story.append(Spacer(1, 0.1*inch))

        elif block_type == 'heading2':
            story.append(Paragraph(content, styles['Heading2']))
            story.append(Spacer(1, 0.1*inch))

        elif block_type == 'paragraph':
            story.append(Paragraph(content, styles['BodyText']))
            story.append(Spacer(1, 0.1*inch))

        elif block_type == 'bullet':
            story.append(Paragraph(content, styles['Bullet']))

        elif block_type == 'space':
            height = block.get('height', 0.2)
            story.append(Spacer(1, height*inch))

        elif block_type == 'pagebreak':
            story.append(PageBreak())

    # Build PDF
    doc.build(story)
    return filename


def create_styled_table(data, col_widths=None, style_name='default'):
    """
    Create a styled table with common styling presets.

    Args:
        data: List of lists containing table data
        col_widths: List of column widths (None for auto)
        style_name: 'default', 'striped', 'minimal', 'report'

    Returns:
        Table object ready to add to story
    """
    table = Table(data, colWidths=col_widths)

    if style_name == 'striped':
        style = TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.darkblue),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 12),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.lightgrey]),
            ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
        ])

    elif style_name == 'minimal':
        style = TableStyle([
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('LINEABOVE', (0, 0), (-1, 0), 2, colors.black),
            ('LINEBELOW', (0, 0), (-1, 0), 1, colors.black),
            ('LINEBELOW', (0, -1), (-1, -1), 2, colors.black),
        ])

    elif style_name == 'report':
        style = TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.black),
            ('ALIGN', (0, 0), (-1, 0), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 11),
            ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
            ('GRID', (0, 0), (-1, -1), 1, colors.grey),
            ('LEFTPADDING', (0, 0), (-1, -1), 12),
            ('RIGHTPADDING', (0, 0), (-1, -1), 12),
        ])

    else:  # default
        style = TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 12),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.white),
            ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ])

    table.setStyle(style)
    return table


def add_header_footer(canvas, doc, header_text="", footer_text=""):
    """
    Callback function to add headers and footers to each page.

    Usage:
        from functools import partial
        callback = partial(add_header_footer, header_text="My Document", footer_text="Confidential")
        template = PageTemplate(id='normal', frames=[frame], onPage=callback)
    """
    canvas.saveState()

    # Header
    if header_text:
        canvas.setFont('Helvetica', 9)
        canvas.drawString(inch, doc.pagesize[1] - 0.5*inch, header_text)

    # Footer
    if footer_text:
        canvas.setFont('Helvetica', 9)
        canvas.drawString(inch, 0.5*inch, footer_text)

    # Page number
    canvas.drawRightString(doc.pagesize[0] - inch, 0.5*inch, f"Page {doc.page}")

    canvas.restoreState()


# Example usage
if __name__ == "__main__":
    # Example 1: Simple document
    content = [
        {'type': 'heading', 'content': 'Introduction'},
        {'type': 'paragraph', 'content': 'This is a sample paragraph with some text.'},
        {'type': 'space', 'height': 0.2},
        {'type': 'heading', 'content': 'Main Content'},
        {'type': 'paragraph', 'content': 'More content here with <b>bold</b> and <i>italic</i> text.'},
        {'type': 'bullet', 'content': 'First bullet point'},
        {'type': 'bullet', 'content': 'Second bullet point'},
    ]

    create_simple_document(
        "example_document.pdf",
        "Sample Document",
        author="John Doe",
        content_blocks=content
    )

    print("Created: example_document.pdf")

    # Example 2: Document with styled table
    doc = SimpleDocTemplate("table_example.pdf", pagesize=letter)
    story = []
    styles = getSampleStyleSheet()

    story.append(Paragraph("Sales Report", styles['Title']))
    story.append(Spacer(1, 0.3*inch))

    # Create table
    data = [
        ['Product', 'Q1', 'Q2', 'Q3', 'Q4'],
        ['Widget A', '100', '150', '130', '180'],
        ['Widget B', '80', '120', '110', '160'],
        ['Widget C', '90', '110', '100', '140'],
    ]

    table = create_styled_table(data, col_widths=[2*inch, 1*inch, 1*inch, 1*inch, 1*inch], style_name='striped')
    story.append(table)

    doc.build(story)
    print("Created: table_example.pdf")
