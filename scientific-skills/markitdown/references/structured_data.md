# Structured Data Handling Reference

This document provides detailed information about converting structured data formats (CSV, JSON, XML) to Markdown.

## CSV Files

Convert CSV (Comma-Separated Values) files to Markdown tables.

### Basic CSV Conversion

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("data.csv")
print(result.text_content)
```

### CSV to Markdown Table

CSV files are automatically converted to Markdown table format:

**Input CSV (`data.csv`):**
```csv
Name,Age,City
Alice,30,New York
Bob,25,Los Angeles
Charlie,35,Chicago
```

**Output Markdown:**
```markdown
| Name    | Age | City        |
|---------|-----|-------------|
| Alice   | 30  | New York    |
| Bob     | 25  | Los Angeles |
| Charlie | 35  | Chicago     |
```

### CSV Conversion Features

**What's preserved:**
- All column headers
- All data rows
- Cell values (text and numbers)
- Column structure

**Formatting:**
- Headers are bolded (Markdown table format)
- Columns are aligned
- Empty cells are preserved
- Special characters are escaped

### Large CSV Files

For large CSV files:

```python
from markitdown import MarkItDown

md = MarkItDown()

# Convert large CSV
result = md.convert("large_dataset.csv")

# Save to file instead of printing
with open("output.md", "w") as f:
    f.write(result.text_content)
```

**Performance considerations:**
- Very large files may take time to process
- Consider previewing first few rows for testing
- Memory usage scales with file size
- Very wide tables may not display well in all Markdown viewers

### CSV with Special Characters

CSV files containing special characters are handled automatically:

```python
from markitdown import MarkItDown

md = MarkItDown()

# Handles UTF-8, special characters, quotes, etc.
result = md.convert("international_data.csv")
```

### CSV Delimiters

Standard CSV delimiters are supported:
- Comma (`,`) - standard
- Semicolon (`;`) - common in European formats
- Tab (`\t`) - TSV files

### Command-Line CSV Conversion

```bash
# Basic conversion
markitdown data.csv -o data.md

# Multiple CSV files
for file in *.csv; do
    markitdown "$file" -o "${file%.csv}.md"
done
```

## JSON Files

Convert JSON data to readable Markdown format.

### Basic JSON Conversion

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("data.json")
print(result.text_content)
```

### JSON Formatting

JSON is converted to a readable, structured Markdown format:

**Input JSON (`config.json`):**
```json
{
  "name": "MyApp",
  "version": "1.0.0",
  "dependencies": {
    "library1": "^2.0.0",
    "library2": "^3.1.0"
  },
  "features": ["auth", "api", "database"]
}
```

**Output Markdown:**
```markdown
## Configuration

**name:** MyApp
**version:** 1.0.0

### dependencies
- **library1:** ^2.0.0
- **library2:** ^3.1.0

### features
- auth
- api
- database
```

### JSON Array Handling

JSON arrays are converted to lists or tables:

**Array of objects:**
```json
[
  {"id": 1, "name": "Alice", "active": true},
  {"id": 2, "name": "Bob", "active": false}
]
```

**Converted to table:**
```markdown
| id | name  | active |
|----|-------|--------|
| 1  | Alice | true   |
| 2  | Bob   | false  |
```

### Nested JSON Structures

Nested JSON is converted with appropriate indentation and hierarchy:

```python
from markitdown import MarkItDown

md = MarkItDown()

# Handles deeply nested structures
result = md.convert("complex_config.json")
print(result.text_content)
```

### JSON Lines (JSONL)

For JSON Lines format (one JSON object per line):

```python
from markitdown import MarkItDown
import json

md = MarkItDown()

# Read JSONL file
with open("data.jsonl", "r") as f:
    for line in f:
        obj = json.loads(line)

        # Convert to JSON temporarily
        with open("temp.json", "w") as temp:
            json.dump(obj, temp)

        result = md.convert("temp.json")
        print(result.text_content)
        print("\n---\n")
```

### Large JSON Files

For large JSON files:

```python
from markitdown import MarkItDown

md = MarkItDown()

# Convert large JSON
result = md.convert("large_data.json")

# Save to file
with open("output.md", "w") as f:
    f.write(result.text_content)
```

## XML Files

Convert XML documents to structured Markdown.

### Basic XML Conversion

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("data.xml")
print(result.text_content)
```

### XML Structure Preservation

XML is converted to Markdown maintaining hierarchical structure:

**Input XML (`book.xml`):**
```xml
<?xml version="1.0"?>
<book>
  <title>Example Book</title>
  <author>John Doe</author>
  <chapters>
    <chapter id="1">
      <title>Introduction</title>
      <content>Chapter 1 content...</content>
    </chapter>
    <chapter id="2">
      <title>Background</title>
      <content>Chapter 2 content...</content>
    </chapter>
  </chapters>
</book>
```

**Output Markdown:**
```markdown
# book

## title
Example Book

## author
John Doe

## chapters

### chapter (id: 1)
#### title
Introduction

#### content
Chapter 1 content...

### chapter (id: 2)
#### title
Background

#### content
Chapter 2 content...
```

### XML Attributes

XML attributes are preserved in the conversion:

```python
from markitdown import MarkItDown

md = MarkItDown()
result = md.convert("data.xml")
# Attributes shown as (attr: value) in headings
```

### XML Namespaces

XML namespaces are handled:

```python
from markitdown import MarkItDown

md = MarkItDown()

# Handles xmlns and namespaced elements
result = md.convert("namespaced.xml")
```

### XML Use Cases

**Configuration files:**
- Convert XML configs to readable format
- Document system configurations
- Compare configuration files

**Data interchange:**
- Convert XML APIs responses
- Process XML data feeds
- Transform between formats

**Document processing:**
- Convert DocBook to Markdown
- Process SVG descriptions
- Extract structured data

## Structured Data Workflows

### CSV Data Analysis Pipeline

```python
from markitdown import MarkItDown
import pandas as pd

md = MarkItDown()

# Read CSV for analysis
df = pd.read_csv("data.csv")

# Do analysis
summary = df.describe()

# Convert both to Markdown
original = md.convert("data.csv")

# Save summary as CSV then convert
summary.to_csv("summary.csv")
summary_md = md.convert("summary.csv")

print("## Original Data\n")
print(original.text_content)
print("\n## Statistical Summary\n")
print(summary_md.text_content)
```

### JSON API Documentation

```python
from markitdown import MarkItDown
import requests
import json

md = MarkItDown()

# Fetch JSON from API
response = requests.get("https://api.example.com/data")
data = response.json()

# Save as JSON
with open("api_response.json", "w") as f:
    json.dump(data, f, indent=2)

# Convert to Markdown
result = md.convert("api_response.json")

# Create documentation
doc = f"""# API Response Documentation

## Endpoint
GET https://api.example.com/data

## Response
{result.text_content}
"""

with open("api_docs.md", "w") as f:
    f.write(doc)
```

### XML to Markdown Documentation

```python
from markitdown import MarkItDown

md = MarkItDown()

# Convert XML documentation
xml_files = ["config.xml", "schema.xml", "data.xml"]

for xml_file in xml_files:
    result = md.convert(xml_file)

    output_name = xml_file.replace('.xml', '.md')
    with open(f"docs/{output_name}", "w") as f:
        f.write(result.text_content)
```

### Multi-Format Data Processing

```python
from markitdown import MarkItDown
import os

md = MarkItDown()

def convert_structured_data(directory):
    """Convert all structured data files in directory."""
    extensions = {'.csv', '.json', '.xml'}

    for filename in os.listdir(directory):
        ext = os.path.splitext(filename)[1]

        if ext in extensions:
            input_path = os.path.join(directory, filename)
            result = md.convert(input_path)

            # Save Markdown
            output_name = filename.replace(ext, '.md')
            output_path = os.path.join("markdown", output_name)

            with open(output_path, 'w') as f:
                f.write(result.text_content)

            print(f"Converted: {filename} → {output_name}")

# Process all structured data
convert_structured_data("data")
```

### CSV to JSON to Markdown

```python
import pandas as pd
from markitdown import MarkItDown
import json

md = MarkItDown()

# Read CSV
df = pd.read_csv("data.csv")

# Convert to JSON
json_data = df.to_dict(orient='records')
with open("temp.json", "w") as f:
    json.dump(json_data, f, indent=2)

# Convert JSON to Markdown
result = md.convert("temp.json")
print(result.text_content)
```

### Database Export to Markdown

```python
from markitdown import MarkItDown
import sqlite3
import csv

md = MarkItDown()

# Export database query to CSV
conn = sqlite3.connect("database.db")
cursor = conn.execute("SELECT * FROM users")

with open("users.csv", "w", newline='') as f:
    writer = csv.writer(f)
    writer.writerow([description[0] for description in cursor.description])
    writer.writerows(cursor.fetchall())

# Convert to Markdown
result = md.convert("users.csv")
print(result.text_content)
```

## Error Handling

### CSV Errors

```python
from markitdown import MarkItDown

md = MarkItDown()

try:
    result = md.convert("data.csv")
    print(result.text_content)
except FileNotFoundError:
    print("CSV file not found")
except Exception as e:
    print(f"CSV conversion error: {e}")
    # Common issues: encoding problems, malformed CSV, delimiter issues
```

### JSON Errors

```python
from markitdown import MarkItDown

md = MarkItDown()

try:
    result = md.convert("data.json")
    print(result.text_content)
except Exception as e:
    print(f"JSON conversion error: {e}")
    # Common issues: invalid JSON syntax, encoding issues
```

### XML Errors

```python
from markitdown import MarkItDown

md = MarkItDown()

try:
    result = md.convert("data.xml")
    print(result.text_content)
except Exception as e:
    print(f"XML conversion error: {e}")
    # Common issues: malformed XML, encoding problems, namespace issues
```

## Best Practices

### CSV Processing
- Check delimiter before conversion
- Verify encoding (UTF-8 recommended)
- Handle large files with streaming if needed
- Preview output for very wide tables

### JSON Processing
- Validate JSON before conversion
- Consider pretty-printing complex structures
- Handle circular references appropriately
- Be aware of large array performance

### XML Processing
- Validate XML structure first
- Handle namespaces consistently
- Consider XPath for selective extraction
- Be mindful of very deep nesting

### Data Quality
- Clean data before conversion when possible
- Handle missing values appropriately
- Verify special character handling
- Test with representative samples

### Performance
- Process large files in batches
- Use streaming for very large datasets
- Monitor memory usage
- Cache converted results when appropriate
