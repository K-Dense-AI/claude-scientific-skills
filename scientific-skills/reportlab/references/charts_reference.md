# Charts and Graphics Reference

Comprehensive guide to creating charts and data visualizations in ReportLab.

## Graphics Architecture

ReportLab's graphics system provides platform-independent drawing:

- **Drawings** - Container for shapes and charts
- **Shapes** - Primitives (rectangles, circles, lines, polygons, paths)
- **Renderers** - Convert to PDF, PostScript, SVG, or bitmaps (PNG, GIF, JPG)
- **Coordinate System** - Y-axis points upward (like PDF, unlike web graphics)

## Quick Start

```python
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.barcharts import VerticalBarChart
from reportlab.graphics import renderPDF

# Create drawing (canvas for chart)
drawing = Drawing(400, 200)

# Create chart
chart = VerticalBarChart()
chart.x = 50
chart.y = 50
chart.width = 300
chart.height = 125
chart.data = [[100, 150, 130, 180]]
chart.categoryAxis.categoryNames = ['Q1', 'Q2', 'Q3', 'Q4']

# Add chart to drawing
drawing.add(chart)

# Render to PDF
renderPDF.drawToFile(drawing, 'chart.pdf', 'Chart Title')

# Or add as flowable to Platypus document
story.append(drawing)
```

## Available Chart Types

### Bar Charts

```python
from reportlab.graphics.charts.barcharts import (
    VerticalBarChart,
    HorizontalBarChart,
)

# Vertical bar chart
chart = VerticalBarChart()
chart.x = 50
chart.y = 50
chart.width = 300
chart.height = 150

# Single series
chart.data = [[100, 150, 130, 180, 140]]

# Multiple series (grouped bars)
chart.data = [
    [100, 150, 130, 180],  # Series 1
    [80, 120, 110, 160],   # Series 2
]

# Categories
chart.categoryAxis.categoryNames = ['Q1', 'Q2', 'Q3', 'Q4']

# Colors for each series
chart.bars[0].fillColor = colors.blue
chart.bars[1].fillColor = colors.red

# Bar spacing
chart.barWidth = 10
chart.groupSpacing = 10
chart.barSpacing = 2
```

### Stacked Bar Charts

```python
from reportlab.graphics.charts.barcharts import VerticalBarChart

chart = VerticalBarChart()
# ... set position and size ...

chart.data = [
    [100, 150, 130, 180],  # Bottom layer
    [50, 70, 60, 90],      # Top layer
]
chart.categoryAxis.categoryNames = ['Q1', 'Q2', 'Q3', 'Q4']

# Enable stacking
chart.barLabelFormat = 'values'
chart.valueAxis.visible = 1
```

### Horizontal Bar Charts

```python
from reportlab.graphics.charts.barcharts import HorizontalBarChart

chart = HorizontalBarChart()
chart.x = 50
chart.y = 50
chart.width = 300
chart.height = 150

chart.data = [[100, 150, 130, 180]]
chart.categoryAxis.categoryNames = ['Product A', 'Product B', 'Product C', 'Product D']

# Horizontal charts use valueAxis horizontally
chart.valueAxis.valueMin = 0
chart.valueAxis.valueMax = 200
```

### Line Charts

```python
from reportlab.graphics.charts.linecharts import HorizontalLineChart

chart = HorizontalLineChart()
chart.x = 50
chart.y = 50
chart.width = 300
chart.height = 150

# Multiple lines
chart.data = [
    [100, 150, 130, 180, 140],  # Line 1
    [80, 120, 110, 160, 130],   # Line 2
]

chart.categoryAxis.categoryNames = ['Jan', 'Feb', 'Mar', 'Apr', 'May']

# Line styling
chart.lines[0].strokeColor = colors.blue
chart.lines[0].strokeWidth = 2
chart.lines[1].strokeColor = colors.red
chart.lines[1].strokeWidth = 2

# Show/hide points
chart.lines[0].symbol = None  # No symbols
# Or use symbols from makeMarker()
```

### Line Plots (X-Y Plots)

```python
from reportlab.graphics.charts.lineplots import LinePlot

chart = LinePlot()
chart.x = 50
chart.y = 50
chart.width = 300
chart.height = 150

# Data as (x, y) tuples
chart.data = [
    [(0, 0), (1, 1), (2, 4), (3, 9), (4, 16)],  # y = x^2
    [(0, 0), (1, 2), (2, 4), (3, 6), (4, 8)],   # y = 2x
]

# Both axes are value axes (not category)
chart.xValueAxis.valueMin = 0
chart.xValueAxis.valueMax = 5
chart.yValueAxis.valueMin = 0
chart.yValueAxis.valueMax = 20

# Line styling
chart.lines[0].strokeColor = colors.blue
chart.lines[1].strokeColor = colors.red
```

### Pie Charts

```python
from reportlab.graphics.charts.piecharts import Pie

chart = Pie()
chart.x = 100
chart.y = 50
chart.width = 200
chart.height = 200

chart.data = [25, 35, 20, 20]
chart.labels = ['Q1', 'Q2', 'Q3', 'Q4']

# Slice colors
chart.slices[0].fillColor = colors.blue
chart.slices[1].fillColor = colors.red
chart.slices[2].fillColor = colors.green
chart.slices[3].fillColor = colors.yellow

# Pop out a slice
chart.slices[1].popout = 10

# Label positioning
chart.slices.strokeColor = colors.white
chart.slices.strokeWidth = 2
```

### Pie Chart with Side Labels

```python
from reportlab.graphics.charts.piecharts import Pie

chart = Pie()
# ... set position, data, labels ...

# Side label mode (labels in columns beside pie)
chart.sideLabels = 1
chart.sideLabelsOffset = 0.1  # Distance from pie

# Simple labels (not fancy layout)
chart.simpleLabels = 1
```

### Area Charts

```python
from reportlab.graphics.charts.areacharts import HorizontalAreaChart

chart = HorizontalAreaChart()
chart.x = 50
chart.y = 50
chart.width = 300
chart.height = 150

# Areas stack on top of each other
chart.data = [
    [100, 150, 130, 180],  # Bottom area
    [50, 70, 60, 90],      # Top area
]

chart.categoryAxis.categoryNames = ['Q1', 'Q2', 'Q3', 'Q4']

# Area colors
chart.strands[0].fillColor = colors.lightblue
chart.strands[1].fillColor = colors.pink
```

### Scatter Charts

```python
from reportlab.graphics.charts.lineplots import ScatterPlot

chart = ScatterPlot()
chart.x = 50
chart.y = 50
chart.width = 300
chart.height = 150

# Data points
chart.data = [
    [(1, 2), (2, 3), (3, 5), (4, 4), (5, 6)],  # Series 1
    [(1, 1), (2, 2), (3, 3), (4, 3), (5, 4)],  # Series 2
]

# Hide lines, show points only
chart.lines[0].strokeColor = None
chart.lines[1].strokeColor = None

# Marker symbols
from reportlab.graphics.widgets.markers import makeMarker
chart.lines[0].symbol = makeMarker('Circle')
chart.lines[1].symbol = makeMarker('Square')
```

## Axes Configuration

### Category Axis (XCategoryAxis)

For categorical data (labels, not numbers):

```python
# Access via chart
axis = chart.categoryAxis

# Labels
axis.categoryNames = ['Jan', 'Feb', 'Mar', 'Apr']

# Label angle (for long labels)
axis.labels.angle = 45
axis.labels.dx = 0
axis.labels.dy = -5

# Label formatting
axis.labels.fontSize = 10
axis.labels.fontName = 'Helvetica'

# Visibility
axis.visible = 1
```

### Value Axis (YValueAxis)

For numeric data:

```python
# Access via chart
axis = chart.valueAxis

# Range
axis.valueMin = 0
axis.valueMax = 200
axis.valueStep = 50  # Tick interval

# Or auto-configure
axis.valueSteps = [0, 50, 100, 150, 200]  # Explicit steps

# Label formatting
axis.labels.fontSize = 10
axis.labelTextFormat = '%d%%'  # Add percentage sign

# Grid lines
axis.strokeWidth = 1
axis.strokeColor = colors.black
```

## Styling and Customization

### Colors

```python
from reportlab.lib import colors

# Named colors
colors.blue, colors.red, colors.green, colors.yellow

# RGB
colors.Color(0.5, 0.5, 0.5)  # Grey

# With alpha
colors.Color(1, 0, 0, alpha=0.5)  # Semi-transparent red

# Hex colors
colors.HexColor('#FF5733')
```

### Line Styling

```python
# For line charts
chart.lines[0].strokeColor = colors.blue
chart.lines[0].strokeWidth = 2
chart.lines[0].strokeDashArray = [2, 2]  # Dashed line
```

### Bar Labels

```python
# Show values on bars
chart.barLabels.nudge = 5  # Offset from bar top
chart.barLabels.fontSize = 8
chart.barLabelFormat = '%d'  # Number format

# For negative values
chart.barLabels.dy = -5  # Position below bar
```

## Legends

Charts can have associated legends:

```python
from reportlab.graphics.charts.legends import Legend

# Create legend
legend = Legend()
legend.x = 350
legend.y = 150
legend.columnMaximum = 10

# Link to chart (share colors)
legend.colorNamePairs = [
    (chart.bars[0].fillColor, 'Series 1'),
    (chart.bars[1].fillColor, 'Series 2'),
]

# Add to drawing
drawing.add(legend)
```

## Drawing Shapes

### Basic Shapes

```python
from reportlab.graphics.shapes import (
    Drawing, Rect, Circle, Ellipse, Line, Polygon, String
)
from reportlab.lib import colors

drawing = Drawing(400, 200)

# Rectangle
rect = Rect(50, 50, 100, 50)
rect.fillColor = colors.blue
rect.strokeColor = colors.black
rect.strokeWidth = 1
drawing.add(rect)

# Circle
circle = Circle(200, 100, 30)
circle.fillColor = colors.red
drawing.add(circle)

# Line
line = Line(50, 150, 350, 150)
line.strokeColor = colors.black
line.strokeWidth = 2
drawing.add(line)

# Text
text = String(50, 175, "Label Text")
text.fontSize = 12
text.fontName = 'Helvetica'
drawing.add(text)
```

### Paths (Complex Shapes)

```python
from reportlab.graphics.shapes import Path

path = Path()
path.moveTo(50, 50)
path.lineTo(100, 100)
path.curveTo(120, 120, 140, 100, 150, 50)
path.closePath()

path.fillColor = colors.lightblue
path.strokeColor = colors.blue
path.strokeWidth = 2

drawing.add(path)
```

## Rendering Options

### Render to PDF

```python
from reportlab.graphics import renderPDF

# Direct to file
renderPDF.drawToFile(drawing, 'output.pdf', 'Chart Title')

# As flowable in Platypus
story.append(drawing)
```

### Render to Image

```python
from reportlab.graphics import renderPM

# PNG
renderPM.drawToFile(drawing, 'chart.png', fmt='PNG')

# GIF
renderPM.drawToFile(drawing, 'chart.gif', fmt='GIF')

# JPG
renderPM.drawToFile(drawing, 'chart.jpg', fmt='JPG')

# With specific DPI
renderPM.drawToFile(drawing, 'chart.png', fmt='PNG', dpi=150)
```

### Render to SVG

```python
from reportlab.graphics import renderSVG

renderSVG.drawToFile(drawing, 'chart.svg')
```

## Advanced Customization

### Inspect Properties

```python
# List all properties
print(chart.getProperties())

# Dump properties (for debugging)
chart.dumpProperties()

# Set multiple properties
chart.setProperties({
    'width': 400,
    'height': 200,
    'data': [[100, 150, 130]],
})
```

### Custom Colors for Series

```python
# Define color scheme
from reportlab.lib.colors import PCMYKColor

colors_list = [
    PCMYKColor(100, 67, 0, 23),   # Blue
    PCMYKColor(0, 100, 100, 0),   # Red
    PCMYKColor(66, 13, 0, 22),    # Green
]

# Apply to chart
for i, color in enumerate(colors_list):
    chart.bars[i].fillColor = color
```

## Complete Examples

### Sales Report Bar Chart

```python
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.barcharts import VerticalBarChart
from reportlab.graphics.charts.legends import Legend
from reportlab.lib import colors

drawing = Drawing(400, 250)

# Create chart
chart = VerticalBarChart()
chart.x = 50
chart.y = 50
chart.width = 300
chart.height = 150

# Data
chart.data = [
    [120, 150, 180, 200],  # 2023
    [100, 130, 160, 190],  # 2022
]
chart.categoryAxis.categoryNames = ['Q1', 'Q2', 'Q3', 'Q4']

# Styling
chart.bars[0].fillColor = colors.HexColor('#3498db')
chart.bars[1].fillColor = colors.HexColor('#e74c3c')
chart.valueAxis.valueMin = 0
chart.valueAxis.valueMax = 250
chart.categoryAxis.labels.fontSize = 10
chart.valueAxis.labels.fontSize = 10

# Add legend
legend = Legend()
legend.x = 325
legend.y = 200
legend.columnMaximum = 2
legend.colorNamePairs = [
    (chart.bars[0].fillColor, '2023'),
    (chart.bars[1].fillColor, '2022'),
]

drawing.add(chart)
drawing.add(legend)

# Add to story or save
story.append(drawing)
```

### Multi-Line Trend Chart

```python
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.linecharts import HorizontalLineChart
from reportlab.lib import colors

drawing = Drawing(400, 250)

chart = HorizontalLineChart()
chart.x = 50
chart.y = 50
chart.width = 320
chart.height = 170

# Data
chart.data = [
    [10, 15, 12, 18, 20, 25],  # Product A
    [8, 10, 14, 16, 18, 22],   # Product B
    [12, 11, 13, 15, 17, 19],  # Product C
]

chart.categoryAxis.categoryNames = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun']

# Line styling
chart.lines[0].strokeColor = colors.blue
chart.lines[0].strokeWidth = 2
chart.lines[1].strokeColor = colors.red
chart.lines[1].strokeWidth = 2
chart.lines[2].strokeColor = colors.green
chart.lines[2].strokeWidth = 2

# Axes
chart.valueAxis.valueMin = 0
chart.valueAxis.valueMax = 30
chart.categoryAxis.labels.angle = 0
chart.categoryAxis.labels.fontSize = 9
chart.valueAxis.labels.fontSize = 9

drawing.add(chart)
story.append(drawing)
```

## Best Practices

1. **Set explicit dimensions** for Drawing to ensure consistent sizing
2. **Position charts** with enough margin (x, y at least 30-50 from edge)
3. **Use consistent color schemes** throughout document
4. **Set valueMin and valueMax** explicitly for consistent scales
5. **Test with realistic data** to ensure labels fit and don't overlap
6. **Add legends** for multi-series charts
7. **Angle category labels** if they're long (45° works well)
8. **Keep it simple** - fewer data series are easier to read
9. **Use appropriate chart types** - bars for comparisons, lines for trends, pies for proportions
10. **Consider colorblind-friendly palettes** - avoid red/green combinations
