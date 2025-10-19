#!/usr/bin/env python3
"""
Invoice Template - Complete example of a professional invoice

This template demonstrates:
- Company header with logo placement
- Client information
- Invoice details table
- Calculations (subtotal, tax, total)
- Professional styling
- Terms and conditions footer
"""

from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image
)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_RIGHT, TA_CENTER
from datetime import datetime


def create_invoice(
    filename,
    invoice_number,
    invoice_date,
    due_date,
    company_info,
    client_info,
    items,
    tax_rate=0.0,
    notes="",
    terms="Payment due within 30 days.",
    logo_path=None
):
    """
    Create a professional invoice PDF.

    Args:
        filename: Output PDF filename
        invoice_number: Invoice number (e.g., "INV-2024-001")
        invoice_date: Date of invoice (datetime or string)
        due_date: Payment due date (datetime or string)
        company_info: Dict with company details
            {'name': 'Company Name', 'address': 'Address', 'phone': 'Phone', 'email': 'Email'}
        client_info: Dict with client details (same structure as company_info)
        items: List of dicts with item details
            [{'description': 'Item', 'quantity': 1, 'unit_price': 100.00}, ...]
        tax_rate: Tax rate as decimal (e.g., 0.08 for 8%)
        notes: Additional notes to client
        terms: Payment terms
        logo_path: Path to company logo image (optional)
    """
    # Create document
    doc = SimpleDocTemplate(filename, pagesize=letter,
                          rightMargin=0.5*inch, leftMargin=0.5*inch,
                          topMargin=0.5*inch, bottomMargin=0.5*inch)

    # Container for elements
    story = []
    styles = getSampleStyleSheet()

    # Create custom styles
    title_style = ParagraphStyle(
        'InvoiceTitle',
        parent=styles['Heading1'],
        fontSize=24,
        textColor=colors.HexColor('#2C3E50'),
        spaceAfter=12,
    )

    header_style = ParagraphStyle(
        'Header',
        parent=styles['Normal'],
        fontSize=10,
        textColor=colors.HexColor('#34495E'),
    )

    # --- HEADER SECTION ---
    header_data = []

    # Company info (left side)
    company_text = f"""
    <b><font size="14">{company_info['name']}</font></b><br/>
    {company_info.get('address', '')}<br/>
    Phone: {company_info.get('phone', '')}<br/>
    Email: {company_info.get('email', '')}
    """

    # Invoice title and number (right side)
    invoice_text = f"""
    <b><font size="16" color="#2C3E50">INVOICE</font></b><br/>
    <font size="10">Invoice #: {invoice_number}</font><br/>
    <font size="10">Date: {invoice_date}</font><br/>
    <font size="10">Due Date: {due_date}</font>
    """

    if logo_path:
        logo = Image(logo_path, width=1.5*inch, height=1*inch)
        header_data = [[logo, Paragraph(company_text, header_style), Paragraph(invoice_text, header_style)]]
        header_table = Table(header_data, colWidths=[1.5*inch, 3*inch, 2.5*inch])
    else:
        header_data = [[Paragraph(company_text, header_style), Paragraph(invoice_text, header_style)]]
        header_table = Table(header_data, colWidths=[4.5*inch, 2.5*inch])

    header_table.setStyle(TableStyle([
        ('VALIGN', (0, 0), (-1, -1), 'TOP'),
        ('ALIGN', (-1, 0), (-1, -1), 'RIGHT'),
    ]))

    story.append(header_table)
    story.append(Spacer(1, 0.3*inch))

    # --- CLIENT INFORMATION ---
    client_label = Paragraph("<b>Bill To:</b>", header_style)
    client_text = f"""
    <b>{client_info['name']}</b><br/>
    {client_info.get('address', '')}<br/>
    Phone: {client_info.get('phone', '')}<br/>
    Email: {client_info.get('email', '')}
    """
    client_para = Paragraph(client_text, header_style)

    client_table = Table([[client_label, client_para]], colWidths=[1*inch, 6*inch])
    story.append(client_table)
    story.append(Spacer(1, 0.3*inch))

    # --- ITEMS TABLE ---
    # Table header
    items_data = [['Description', 'Quantity', 'Unit Price', 'Amount']]

    # Calculate items
    subtotal = 0
    for item in items:
        desc = item['description']
        qty = item['quantity']
        price = item['unit_price']
        amount = qty * price
        subtotal += amount

        items_data.append([
            desc,
            str(qty),
            f"${price:,.2f}",
            f"${amount:,.2f}"
        ])

    # Create items table
    items_table = Table(items_data, colWidths=[3.5*inch, 1*inch, 1.5*inch, 1*inch])

    items_table.setStyle(TableStyle([
        # Header row
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#34495E')),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, 0), 'CENTER'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 11),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),

        # Data rows
        ('BACKGROUND', (0, 1), (-1, -1), colors.white),
        ('ALIGN', (1, 1), (-1, -1), 'RIGHT'),
        ('ALIGN', (0, 1), (0, -1), 'LEFT'),
        ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
        ('FONTSIZE', (0, 1), (-1, -1), 10),
        ('TOPPADDING', (0, 1), (-1, -1), 6),
        ('BOTTOMPADDING', (0, 1), (-1, -1), 6),
        ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
    ]))

    story.append(items_table)
    story.append(Spacer(1, 0.2*inch))

    # --- TOTALS SECTION ---
    tax_amount = subtotal * tax_rate
    total = subtotal + tax_amount

    totals_data = [
        ['Subtotal:', f"${subtotal:,.2f}"],
    ]

    if tax_rate > 0:
        totals_data.append([f'Tax ({tax_rate*100:.1f}%):', f"${tax_amount:,.2f}"])

    totals_data.append(['<b>Total:</b>', f"<b>${total:,.2f}</b>"])

    totals_table = Table(totals_data, colWidths=[5*inch, 2*inch])
    totals_table.setStyle(TableStyle([
        ('ALIGN', (0, 0), (-1, -1), 'RIGHT'),
        ('FONTNAME', (0, 0), (-1, -2), 'Helvetica'),
        ('FONTNAME', (0, -1), (-1, -1), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, -1), 11),
        ('TOPPADDING', (0, 0), (-1, -1), 6),
        ('LINEABOVE', (1, -1), (1, -1), 2, colors.HexColor('#34495E')),
    ]))

    story.append(totals_table)
    story.append(Spacer(1, 0.4*inch))

    # --- NOTES ---
    if notes:
        notes_style = ParagraphStyle('Notes', parent=styles['Normal'], fontSize=9)
        story.append(Paragraph(f"<b>Notes:</b><br/>{notes}", notes_style))
        story.append(Spacer(1, 0.2*inch))

    # --- TERMS ---
    terms_style = ParagraphStyle('Terms', parent=styles['Normal'],
                                fontSize=9, textColor=colors.grey)
    story.append(Paragraph(f"<b>Payment Terms:</b><br/>{terms}", terms_style))

    # Build PDF
    doc.build(story)
    return filename


# Example usage
if __name__ == "__main__":
    # Sample data
    company = {
        'name': 'Acme Corporation',
        'address': '123 Business St, Suite 100\nNew York, NY 10001',
        'phone': '(555) 123-4567',
        'email': 'info@acme.com'
    }

    client = {
        'name': 'John Doe',
        'address': '456 Client Ave\nLos Angeles, CA 90001',
        'phone': '(555) 987-6543',
        'email': 'john@example.com'
    }

    items = [
        {'description': 'Web Design Services', 'quantity': 1, 'unit_price': 2500.00},
        {'description': 'Content Writing (10 pages)', 'quantity': 10, 'unit_price': 50.00},
        {'description': 'SEO Optimization', 'quantity': 1, 'unit_price': 750.00},
        {'description': 'Hosting Setup', 'quantity': 1, 'unit_price': 200.00},
    ]

    create_invoice(
        filename="sample_invoice.pdf",
        invoice_number="INV-2024-001",
        invoice_date="January 15, 2024",
        due_date="February 15, 2024",
        company_info=company,
        client_info=client,
        items=items,
        tax_rate=0.08,
        notes="Thank you for your business! We appreciate your prompt payment.",
        terms="Payment due within 30 days. Late payments subject to 1.5% monthly fee.",
        logo_path=None  # Set to your logo path if available
    )

    print("Invoice created: sample_invoice.pdf")
