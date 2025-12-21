# Pie Charts

Pie charts visualize proportional data and percentages.

## Basic Syntax

```mermaid
pie
    title Browser Usage
    "Chrome" : 60
    "Firefox" : 20
    "Safari" : 15
    "Edge" : 5
```

## Simple Data

```mermaid
pie
    "Category A" : 30
    "Category B" : 25
    "Category C" : 45
```

## Common Patterns

### Budget Distribution
```mermaid
pie
    title Annual Budget Allocation
    "Development" : 35
    "Marketing" : 25
    "Operations" : 20
    "Sales" : 15
    "Administration" : 5
```

### Market Share
```mermaid
pie
    title Smartphone Market Share 2024
    "Apple" : 28
    "Samsung" : 22
    "Xiaomi" : 13
    "Oppo" : 10
    "Vivo" : 9
    "Others" : 18
```

### Time Allocation
```mermaid
pie
    title Developer Time Distribution
    "Coding" : 40
    "Meetings" : 20
    "Code Review" : 15
    "Testing" : 10
    "Documentation" : 10
    "Other" : 5
```

### Survey Results
```mermaid
pie
    title Customer Satisfaction Survey
    "Very Satisfied" : 45
    "Satisfied" : 35
    "Neutral" : 12
    "Unsatisfied" : 6
    "Very Unsatisfied" : 2
```

### Technology Stack Usage
```mermaid
pie
    title Backend Technologies Used
    "Node.js" : 35
    "Python" : 30
    "Java" : 18
    "Go" : 10
    "Ruby" : 5
    "PHP" : 2
```

### Traffic Sources
```mermaid
pie
    title Website Traffic Sources
    "Organic Search" : 42
    "Direct" : 28
    "Social Media" : 15
    "Referral" : 10
    "Email" : 5
```

### Product Sales
```mermaid
pie
    title Q1 Product Sales Distribution
    "Product A" : 350
    "Product B" : 280
    "Product C" : 190
    "Product D" : 120
    "Product E" : 60
```

### Employee Distribution
```mermaid
pie
    title Company Departments
    "Engineering" : 120
    "Sales" : 80
    "Marketing" : 45
    "HR" : 25
    "Finance" : 20
    "Operations" : 35
```

### Cloud Provider Usage
```mermaid
pie
    title Infrastructure Distribution
    "AWS" : 55
    "Azure" : 25
    "Google Cloud" : 15
    "On-Premise" : 5
```

### Issue Types
```mermaid
pie
    title Support Ticket Categories
    "Technical Issues" : 40
    "Billing Questions" : 25
    "Feature Requests" : 20
    "Account Issues" : 10
    "Other" : 5
```

## Tips

- Values are relative (will be converted to percentages)
- Use quotes around labels with spaces
- Keep number of slices reasonable (5-8 max for readability)
- Order slices by size (largest to smallest) for clarity
- Group small categories into "Others" if needed
- Use meaningful, concise labels
- Values can be any positive number
- Percentages are calculated automatically
