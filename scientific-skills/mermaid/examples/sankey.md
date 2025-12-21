# Sankey Diagrams

Sankey diagrams visualize flow and distribution of quantities between nodes.

## Basic Syntax

```mermaid
sankey-beta

Source,Target,Value
A,B,10
A,C,20
B,D,15
C,D,25
```

## Common Patterns

### Website Traffic Flow
```mermaid
sankey-beta

Source,Target,Value
Homepage,Products,5000
Homepage,Blog,3000
Homepage,About,1000
Homepage,Contact,500
Products,Product A,2000
Products,Product B,1500
Products,Product C,1000
Products,Cart,500
Blog,Products,1000
Blog,Subscribe,500
Cart,Checkout,300
Cart,Exit,200
Checkout,Purchase,250
Checkout,Exit,50
```

### Energy Distribution
```mermaid
sankey-beta

Source,Target,Value
Solar,Battery,400
Solar,Grid,600
Wind,Battery,300
Wind,Grid,700
Battery,Residential,500
Battery,Commercial,200
Grid,Residential,800
Grid,Commercial,500
Grid,Industrial,1000
```

### Budget Allocation
```mermaid
sankey-beta

Source,Target,Value
Total Budget,Engineering,500000
Total Budget,Marketing,300000
Total Budget,Sales,250000
Total Budget,Operations,200000
Engineering,Development,300000
Engineering,QA,120000
Engineering,DevOps,80000
Marketing,Digital,180000
Marketing,Events,70000
Marketing,Content,50000
Sales,Team Salaries,150000
Sales,Tools,60000
Sales,Travel,40000
```

### Customer Journey Conversion
```mermaid
sankey-beta

Source,Target,Value
Visitors,Homepage,10000
Visitors,Landing Page,5000
Homepage,Signup,1500
Homepage,Exit,8500
Landing Page,Signup,2000
Landing Page,Exit,3000
Signup,Trial,3000
Signup,Exit,500
Trial,Paid,1200
Trial,Churned,1800
```

### Supply Chain Flow
```mermaid
sankey-beta

Source,Target,Value
Supplier A,Warehouse 1,800
Supplier B,Warehouse 1,600
Supplier C,Warehouse 2,900
Warehouse 1,Distribution Center,1400
Warehouse 2,Distribution Center,900
Distribution Center,Store A,700
Distribution Center,Store B,600
Distribution Center,Store C,500
Distribution Center,Online Orders,500
```

### Revenue Streams
```mermaid
sankey-beta

Source,Target,Value
Total Revenue,Subscriptions,6000000
Total Revenue,Services,3000000
Total Revenue,Products,1000000
Subscriptions,Basic Plan,2000000
Subscriptions,Pro Plan,3000000
Subscriptions,Enterprise,1000000
Services,Consulting,1800000
Services,Training,800000
Services,Support,400000
Products,Hardware,600000
Products,Merchandise,400000
```

### Manufacturing Process
```mermaid
sankey-beta

Source,Target,Value
Raw Materials,Processing Plant A,5000
Raw Materials,Processing Plant B,3000
Processing Plant A,Assembly Line 1,3000
Processing Plant A,Assembly Line 2,2000
Processing Plant B,Assembly Line 1,1500
Processing Plant B,Assembly Line 2,1500
Assembly Line 1,Quality Check,4500
Assembly Line 2,Quality Check,3500
Quality Check,Packaging,7500
Quality Check,Defects,500
Packaging,Shipping,7500
```

### User Acquisition Channels
```mermaid
sankey-beta

Source,Target,Value
Organic Search,Landing Page,15000
Paid Ads,Landing Page,8000
Social Media,Landing Page,6000
Referral,Landing Page,4000
Landing Page,Signup,12000
Landing Page,Bounce,21000
Signup,Active User,8000
Signup,Inactive,4000
Active User,Paying Customer,3000
Active User,Free Tier,5000
```

### Cost Breakdown
```mermaid
sankey-beta

Source,Target,Value
Operating Costs,Personnel,450000
Operating Costs,Infrastructure,250000
Operating Costs,Marketing,150000
Operating Costs,Other,100000
Personnel,Salaries,350000
Personnel,Benefits,100000
Infrastructure,Cloud Services,150000
Infrastructure,Office,80000
Infrastructure,Equipment,20000
Marketing,Advertising,100000
Marketing,Events,30000
Marketing,Content,20000
```

### Data Flow Architecture
```mermaid
sankey-beta

Source,Target,Value
Data Sources,API Gateway,50000
Data Sources,Direct DB,20000
API Gateway,Service A,25000
API Gateway,Service B,15000
API Gateway,Service C,10000
Direct DB,Service A,10000
Direct DB,Service C,10000
Service A,Cache,20000
Service A,Database,15000
Service B,Database,15000
Service C,Database,20000
```

### Educational Path
```mermaid
sankey-beta

Source,Target,Value
Enrolled Students,Computer Science,500
Enrolled Students,Engineering,400
Enrolled Students,Business,300
Computer Science,Software Dev,300
Computer Science,Data Science,150
Computer Science,Other,50
Engineering,Mechanical,200
Engineering,Electrical,150
Engineering,Other,50
Business,Finance,150
Business,Marketing,100
Business,Other,50
```

### Food Waste Analysis
```mermaid
sankey-beta

Source,Target,Value
Total Food,Consumed,7000
Total Food,Wasted,3000
Wasted,Preparation Waste,800
Wasted,Spoilage,1200
Wasted,Plate Waste,1000
Preparation Waste,Compost,400
Preparation Waste,Landfill,400
Spoilage,Compost,600
Spoilage,Landfill,600
Plate Waste,Compost,300
Plate Waste,Landfill,700
```

## Data Format

The CSV-like format requires three columns:
- **Source**: Starting node
- **Target**: Ending node  
- **Value**: Flow amount (positive number)

## Tips

- Values represent flow magnitude/quantity
- Node names are case-sensitive
- Flows are directional (Source → Target)
- Larger values = thicker flows
- Use meaningful node names
- Group related flows for clarity
- Keep node count reasonable (< 20 for readability)
- Values should be proportional and meaningful
- Perfect for: conversions, distributions, processes, journeys
- Show where quantities split, merge, or transform
