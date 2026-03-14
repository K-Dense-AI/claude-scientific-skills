---
name: browser-automation
description: Control browsers with Playwright for scientific web scraping, screenshots, and PDF generation. Use when you need to capture web content, scrape dynamic pages, or generate PDFs from HTML. For static API access use bioservices or gget instead.
license: MIT
metadata:
    skill-author: K-Dense Inc.
---

# Browser Automation: Playwright for Scientific Workflows

## Overview

Use Playwright to automate browser interactions for scientific data collection, documentation capture, and PDF generation. This skill covers headless browser control from Claude Code for tasks that require JavaScript rendering or visual capture.

## When to Use This Skill

- Scraping dynamic web pages that require JavaScript rendering
- Taking screenshots of scientific visualizations or dashboards
- Generating PDFs from HTML reports or papers
- Automating form submissions on scientific databases
- Capturing content from pages with lazy loading or infinite scroll
- Visual regression testing of scientific web apps

## When NOT to Use This Skill

- Static API queries → use `bioservices`, `gget`, or database-specific skills
- Simple HTTP requests → use `curl` or Python `requests`
- Data already available via MCP servers → use `biothings-mcp`, `gget-mcp`, etc.

## Prerequisites

```bash
# Install Playwright
npm install -D playwright
npx playwright install chromium

# Or with Python
uv pip install playwright
playwright install chromium
```

## Quick Start

### Screenshot Capture

```python
import asyncio
from playwright.async_api import async_playwright

async def screenshot(url: str, output: str = "screenshot.png"):
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page(viewport={"width": 1920, "height": 1080})
        await page.goto(url, wait_until="networkidle")
        await page.screenshot(path=output, full_page=True)
        await browser.close()
        print(f"Screenshot saved: {output}")

asyncio.run(screenshot("https://pubmed.ncbi.nlm.nih.gov/?term=BRCA1"))
```

### PDF Generation from HTML

```python
async def html_to_pdf(html_content: str, output: str = "report.pdf"):
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()
        await page.set_content(html_content)
        await page.pdf(
            path=output,
            format="A4",
            margin={"top": "1cm", "right": "1cm", "bottom": "1cm", "left": "1cm"},
            print_background=True
        )
        await browser.close()
        print(f"PDF saved: {output}")
```

### Scientific Web Scraping

```python
async def scrape_pubmed(query: str, max_results: int = 10):
    """Scrape PubMed search results with full JavaScript rendering."""
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()
        await page.goto(f"https://pubmed.ncbi.nlm.nih.gov/?term={query}")
        await page.wait_for_selector(".docsum-title")

        results = await page.evaluate("""
            () => {
                return Array.from(document.querySelectorAll('.docsum-content')).map(el => ({
                    title: el.querySelector('.docsum-title')?.textContent?.trim(),
                    authors: el.querySelector('.docsum-authors')?.textContent?.trim(),
                    journal: el.querySelector('.docsum-journal-citation')?.textContent?.trim(),
                    pmid: el.closest('article')?.getAttribute('data-pmid')
                }));
            }
        """)
        await browser.close()
        return results[:max_results]
```

## Core Patterns

### Pattern 1: Batch Screenshot Pipeline

```python
async def batch_screenshots(urls: list[str], output_dir: str = "./screenshots"):
    """Take screenshots of multiple URLs in parallel."""
    import os
    os.makedirs(output_dir, exist_ok=True)

    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)

        async def capture(url, index):
            page = await browser.new_page()
            await page.goto(url, wait_until="networkidle", timeout=30000)
            filename = f"{output_dir}/page_{index:03d}.png"
            await page.screenshot(path=filename, full_page=True)
            await page.close()
            return filename

        tasks = [capture(url, i) for i, url in enumerate(urls)]
        results = await asyncio.gather(*tasks, return_exceptions=True)
        await browser.close()
        return results
```

### Pattern 2: Dynamic Data Extraction

```python
async def extract_table_data(url: str, table_selector: str = "table"):
    """Extract table data from a dynamically rendered page."""
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()
        await page.goto(url, wait_until="networkidle")
        await page.wait_for_selector(table_selector)

        data = await page.evaluate(f"""
            () => {{
                const table = document.querySelector('{table_selector}');
                const headers = Array.from(table.querySelectorAll('th')).map(th => th.textContent.trim());
                const rows = Array.from(table.querySelectorAll('tbody tr')).map(tr =>
                    Array.from(tr.querySelectorAll('td')).map(td => td.textContent.trim())
                );
                return {{ headers, rows }};
            }}
        """)
        await browser.close()
        return data
```

### Pattern 3: Scientific Report PDF

```python
async def generate_report_pdf(title: str, sections: list[dict], output: str = "report.pdf"):
    """Generate a formatted PDF report from structured data."""
    html = f"""
    <html>
    <head>
        <style>
            body {{ font-family: 'Helvetica Neue', sans-serif; margin: 2cm; color: #333; }}
            h1 {{ color: #1a5276; border-bottom: 2px solid #1a5276; padding-bottom: 10px; }}
            h2 {{ color: #2980b9; margin-top: 30px; }}
            table {{ border-collapse: collapse; width: 100%; margin: 15px 0; }}
            th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            th {{ background-color: #f4f6f7; }}
            .timestamp {{ color: #888; font-size: 0.9em; }}
        </style>
    </head>
    <body>
        <h1>{title}</h1>
        <p class="timestamp">Generated: <script>document.write(new Date().toISOString())</script></p>
    """
    for section in sections:
        html += f"<h2>{section['title']}</h2><p>{section['content']}</p>"
    html += "</body></html>"

    await html_to_pdf(html, output)
```

## Integration with Claude Code Scripts

### Headless Screenshot from CLI
```bash
# Take a screenshot via Claude Code headless mode
claude -p "Use Playwright to screenshot https://example.com and save to ./screenshot.png" --output-format text

# Batch screenshots
claude -p "Use Playwright to screenshot these URLs and save to ./screenshots/: url1, url2, url3" --output-format text
```

### PDF Generation Pipeline
```bash
# Generate PDF from analysis results
claude -p "Analyze the data in results.csv, create an HTML report with charts, then use Playwright to convert it to PDF" --output-format text
```

## Troubleshooting

- **Timeout on page load**: Increase timeout or use `wait_until="domcontentloaded"` instead of `"networkidle"`
- **Missing elements**: Use `page.wait_for_selector()` before extraction
- **Headless rendering issues**: Try `headless=False` for debugging, then switch back
- **Large pages**: Use `page.set_viewport_size()` for consistent dimensions
- **CAPTCHA/bot detection**: This skill respects bot detection -- do not attempt to bypass
