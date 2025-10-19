#!/usr/bin/env python3
"""
Report Template - Complete example of a professional multi-page report

This template demonstrates:
- Cover page
- Table of contents
- Multiple sections with headers
- Charts and graphs integration
- Tables with data
- Headers and footers
- Professional styling
"""

from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.platypus import (
    BaseDocTemplate, PageTemplate, Frame, Paragraph, Spacer,
    Table, TableStyle, PageBreak, KeepTogether, TableOfContents
)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_RIGHT, TA_CENTER, TA_JUSTIFY
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.barcharts import VerticalBarChart
from reportlab.graphics.charts.linecharts import HorizontalLineChart
from datetime import datetime


def header_footer(canvas, doc):
    """Draw header and footer on each page (except cover)"""
    canvas.saveState()

    # Skip header/footer on cover page (page 1)
    if doc.page > 1:
        # Header
        canvas.setFont('Helvetica', 9)
        canvas.setFillColor(colors.grey)
        canvas.drawString(inch, letter[1] - 0.5*inch, "Quarterly Business Report")
        canvas.line(inch, letter[1] - 0.55*inch, letter[0] - inch, letter[1] - 0.55*inch)

        # Footer
        canvas.drawString(inch, 0.5*inch, f"Generated: {datetime.now().strftime('%B %d, %Y')}")
        canvas.drawRightString(letter[0] - inch, 0.5*inch, f"Page {doc.page - 1}")

    canvas.restoreState()


def create_report(filename, report_data):
    """
    Create a comprehensive business report.

    Args:
        filename: Output PDF filename
        report_data: Dict containing report information
            {
                'title': 'Report Title',
                'subtitle': 'Report Subtitle',
                'author': 'Author Name',
                'date': 'Date',
                'sections': [
                    {
                        'title': 'Section Title',
                        'content': 'Section content...',
                        'subsections': [...],
                        'table': {...},
                        'chart': {...}
                    },
                    ...
                ]
            }
    """
    # Create document with custom page template
    doc = BaseDocTemplate(filename, pagesize=letter,
                         rightMargin=72, leftMargin=72,
                         topMargin=inch, bottomMargin=inch)

    # Define frame for content
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height - 0.5*inch, id='normal')

    # Create page template with header/footer
    template = PageTemplate(id='normal', frames=[frame], onPage=header_footer)
    doc.addPageTemplates([template])

    # Get styles
    styles = getSampleStyleSheet()

    # Custom styles
    title_style = ParagraphStyle(
        'ReportTitle',
        parent=styles['Title'],
        fontSize=28,
        textColor=colors.HexColor('#2C3E50'),
        spaceAfter=20,
        alignment=TA_CENTER,
    )

    subtitle_style = ParagraphStyle(
        'ReportSubtitle',
        parent=styles['Normal'],
        fontSize=14,
        textColor=colors.grey,
        alignment=TA_CENTER,
        spaceAfter=30,
    )

    heading1_style = ParagraphStyle(
        'CustomHeading1',
        parent=styles['Heading1'],
        fontSize=18,
        textColor=colors.HexColor('#2C3E50'),
        spaceAfter=12,
        spaceBefore=12,
    )

    heading2_style = ParagraphStyle(
        'CustomHeading2',
        parent=styles['Heading2'],
        fontSize=14,
        textColor=colors.HexColor('#34495E'),
        spaceAfter=10,
        spaceBefore=10,
    )

    body_style = ParagraphStyle(
        'ReportBody',
        parent=styles['BodyText'],
        fontSize=11,
        alignment=TA_JUSTIFY,
        spaceAfter=12,
        leading=14,
    )

    # Build story
    story = []

    # --- COVER PAGE ---
    story.append(Spacer(1, 2*inch))
    story.append(Paragraph(report_data['title'], title_style))
    story.append(Paragraph(report_data.get('subtitle', ''), subtitle_style))
    story.append(Spacer(1, inch))

    # Cover info table
    cover_info = [
        ['Prepared by:', report_data.get('author', '')],
        ['Date:', report_data.get('date', datetime.now().strftime('%B %d, %Y'))],
        ['Period:', report_data.get('period', 'Q4 2023')],
    ]

    cover_table = Table(cover_info, colWidths=[2*inch, 4*inch])
    cover_table.setStyle(TableStyle([
        ('ALIGN', (0, 0), (0, -1), 'RIGHT'),
        ('ALIGN', (1, 0), (1, -1), 'LEFT'),
        ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, -1), 11),
        ('TOPPADDING', (0, 0), (-1, -1), 6),
    ]))

    story.append(cover_table)
    story.append(PageBreak())

    # --- TABLE OF CONTENTS ---
    toc = TableOfContents()
    toc.levelStyles = [
        ParagraphStyle(name='TOCHeading1', fontSize=14, leftIndent=20, spaceBefore=10, spaceAfter=5),
        ParagraphStyle(name='TOCHeading2', fontSize=12, leftIndent=40, spaceBefore=3, spaceAfter=3),
    ]

    story.append(Paragraph("Table of Contents", heading1_style))
    story.append(toc)
    story.append(PageBreak())

    # --- SECTIONS ---
    for section in report_data.get('sections', []):
        # Section heading
        section_title = section['title']
        story.append(Paragraph(f'<a name="{section_title}"/>{section_title}', heading1_style))

        # Add to TOC
        toc.addEntry(0, section_title, doc.page)

        # Section content
        if 'content' in section:
            for para in section['content'].split('\n\n'):
                if para.strip():
                    story.append(Paragraph(para.strip(), body_style))

        story.append(Spacer(1, 0.2*inch))

        # Subsections
        for subsection in section.get('subsections', []):
            story.append(Paragraph(subsection['title'], heading2_style))

            if 'content' in subsection:
                story.append(Paragraph(subsection['content'], body_style))

            story.append(Spacer(1, 0.1*inch))

        # Add table if provided
        if 'table_data' in section:
            table = create_section_table(section['table_data'])
            story.append(table)
            story.append(Spacer(1, 0.2*inch))

        # Add chart if provided
        if 'chart_data' in section:
            chart = create_section_chart(section['chart_data'])
            story.append(chart)
            story.append(Spacer(1, 0.2*inch))

        story.append(Spacer(1, 0.3*inch))

    # Build PDF (twice for TOC to populate)
    doc.multiBuild(story)
    return filename


def create_section_table(table_data):
    """Create a styled table for report sections"""
    data = table_data['data']
    table = Table(data, colWidths=table_data.get('colWidths'))

    table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#34495E')),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 11),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('BACKGROUND', (0, 1), (-1, -1), colors.white),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.lightgrey]),
        ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
        ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
        ('FONTSIZE', (0, 1), (-1, -1), 10),
    ]))

    return table


def create_section_chart(chart_data):
    """Create a chart for report sections"""
    chart_type = chart_data.get('type', 'bar')
    drawing = Drawing(400, 200)

    if chart_type == 'bar':
        chart = VerticalBarChart()
        chart.x = 50
        chart.y = 30
        chart.width = 300
        chart.height = 150
        chart.data = chart_data['data']
        chart.categoryAxis.categoryNames = chart_data.get('categories', [])
        chart.valueAxis.valueMin = 0

        # Style bars
        for i in range(len(chart_data['data'])):
            chart.bars[i].fillColor = colors.HexColor(['#3498db', '#e74c3c', '#2ecc71'][i % 3])

    elif chart_type == 'line':
        chart = HorizontalLineChart()
        chart.x = 50
        chart.y = 30
        chart.width = 300
        chart.height = 150
        chart.data = chart_data['data']
        chart.categoryAxis.categoryNames = chart_data.get('categories', [])

        # Style lines
        for i in range(len(chart_data['data'])):
            chart.lines[i].strokeColor = colors.HexColor(['#3498db', '#e74c3c', '#2ecc71'][i % 3])
            chart.lines[i].strokeWidth = 2

    drawing.add(chart)
    return drawing


# Example usage
if __name__ == "__main__":
    report = {
        'title': 'Quarterly Business Report',
        'subtitle': 'Q4 2023 Performance Analysis',
        'author': 'Analytics Team',
        'date': 'January 15, 2024',
        'period': 'October - December 2023',
        'sections': [
            {
                'title': 'Executive Summary',
                'content': """
                This report provides a comprehensive analysis of our Q4 2023 performance.
                Overall, the quarter showed strong growth across all key metrics, with
                revenue increasing by 25% year-over-year and customer satisfaction
                scores reaching an all-time high of 4.8/5.0.

                Key highlights include the successful launch of three new products,
                expansion into two new markets, and the completion of our digital
                transformation initiative.
                """,
                'subsections': [
                    {
                        'title': 'Key Achievements',
                        'content': 'Successfully launched Product X with 10,000 units sold in first month.'
                    }
                ]
            },
            {
                'title': 'Financial Performance',
                'content': """
                The financial results for Q4 exceeded expectations across all categories.
                Revenue growth was driven primarily by strong product sales and increased
                market share in key regions.
                """,
                'table_data': {
                    'data': [
                        ['Metric', 'Q3 2023', 'Q4 2023', 'Change'],
                        ['Revenue', '$2.5M', '$3.1M', '+24%'],
                        ['Profit', '$500K', '$680K', '+36%'],
                        ['Expenses', '$2.0M', '$2.4M', '+20%'],
                    ],
                    'colWidths': [2*inch, 1.5*inch, 1.5*inch, 1*inch]
                },
                'chart_data': {
                    'type': 'bar',
                    'data': [[2.5, 3.1], [0.5, 0.68], [2.0, 2.4]],
                    'categories': ['Q3', 'Q4']
                }
            },
            {
                'title': 'Market Analysis',
                'content': """
                Market conditions remained favorable throughout the quarter, with
                strong consumer confidence and increasing demand for our products.
                """,
                'chart_data': {
                    'type': 'line',
                    'data': [[100, 120, 115, 140, 135, 150]],
                    'categories': ['Oct', 'Nov', 'Dec', 'Oct', 'Nov', 'Dec']
                }
            },
        ]
    }

    create_report("sample_report.pdf", report)
    print("Report created: sample_report.pdf")
