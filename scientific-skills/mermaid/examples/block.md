# Block Diagrams

Block diagrams visualize system components and their connections.

**Note:** Block diagram support in Mermaid is new and experimental.

## Basic Syntax

```mermaid
block-beta
  columns 3
  A["Component A"] B["Component B"] C["Component C"]
  D["Component D"] E["Component E"] F["Component F"]
  
  A --> B
  B --> C
  D --> E
  E --> F
```

## Common Patterns

### Web Application Architecture
```mermaid
block-beta
  columns 3
  
  Client["Web Browser"]
  space
  space
  
  space
  LB["Load Balancer"]
  space
  
  Server1["App Server 1"]
  Server2["App Server 2"]
  Server3["App Server 3"]
  
  space
  Cache["Redis Cache"]
  space
  
  space
  DB["Database"]
  space
  
  Client --> LB
  LB --> Server1
  LB --> Server2
  LB --> Server3
  Server1 --> Cache
  Server2 --> Cache
  Server3 --> Cache
  Server1 --> DB
  Server2 --> DB
  Server3 --> DB
```

### Data Processing Pipeline
```mermaid
block-beta
  columns 4
  
  Source["Data Source"]
  Ingest["Ingestion Layer"]
  Process["Processing"]
  Store["Storage"]
  
  API["REST API"]
  Queue["Message Queue"]
  Transform["ETL"]
  Warehouse["Data Warehouse"]
  
  Stream["Stream"]
  Buffer["Buffer"]
  Enrich["Enrichment"]
  Lake["Data Lake"]
  
  Source --> Ingest
  Ingest --> Process
  Process --> Store
  API --> Queue
  Queue --> Transform
  Transform --> Warehouse
  Stream --> Buffer
  Buffer --> Enrich
  Enrich --> Lake
```

### Microservices Architecture
```mermaid
block-beta
  columns 3
  
  Gateway["API Gateway"]
  space
  space
  
  Auth["Auth Service"]
  User["User Service"]
  Order["Order Service"]
  
  Payment["Payment Service"]
  Inventory["Inventory Service"]
  Notification["Notification Service"]
  
  space
  Database["Shared Database"]
  space
  
  Gateway --> Auth
  Gateway --> User
  Gateway --> Order
  Auth --> Database
  User --> Database
  Order --> Payment
  Order --> Inventory
  Inventory --> Notification
```

### Network Topology
```mermaid
block-beta
  columns 4
  
  Internet["Internet"]
  space
  space
  space
  
  space
  Firewall["Firewall"]
  space
  space
  
  space
  Router["Router"]
  space
  space
  
  Switch1["Switch 1"]
  Switch2["Switch 2"]
  Switch3["Switch 3"]
  space
  
  PC1["PC 1"]
  PC2["PC 2"]
  Server1["Server 1"]
  Server2["Server 2"]
  
  Internet --> Firewall
  Firewall --> Router
  Router --> Switch1
  Router --> Switch2
  Router --> Switch3
  Switch1 --> PC1
  Switch1 --> PC2
  Switch2 --> Server1
  Switch3 --> Server2
```

### CI/CD Pipeline
```mermaid
block-beta
  columns 5
  
  Code["Source Code"]
  Build["Build"]
  Test["Test"]
  Deploy["Deploy"]
  Monitor["Monitor"]
  
  Git["Git Push"]
  Compile["Compile"]
  Unit["Unit Tests"]
  Staging["Staging"]
  Logs["Logs"]
  
  space
  Package["Package"]
  Integration["Integration Tests"]
  Production["Production"]
  Metrics["Metrics"]
  
  Code --> Build
  Build --> Test
  Test --> Deploy
  Deploy --> Monitor
  Git --> Compile
  Compile --> Package
  Package --> Unit
  Unit --> Integration
  Integration --> Staging
  Staging --> Production
  Production --> Logs
  Production --> Metrics
```

### IoT System
```mermaid
block-beta
  columns 3
  
  Sensors["IoT Sensors"]
  space
  space
  
  Gateway["IoT Gateway"]
  space
  space
  
  space
  Cloud["Cloud Platform"]
  space
  
  Analytics["Analytics"]
  Storage["Storage"]
  API["API Layer"]
  
  Dashboard["Dashboard"]
  Mobile["Mobile App"]
  Alerts["Alert System"]
  
  Sensors --> Gateway
  Gateway --> Cloud
  Cloud --> Analytics
  Cloud --> Storage
  Cloud --> API
  API --> Dashboard
  API --> Mobile
  Analytics --> Alerts
```

### E-Commerce System
```mermaid
block-beta
  columns 4
  
  Customer["Customer"]
  space
  space
  space
  
  space
  Web["Web App"]
  space
  space
  
  Catalog["Product Catalog"]
  Cart["Shopping Cart"]
  Checkout["Checkout"]
  space
  
  Inventory["Inventory"]
  Payment["Payment Gateway"]
  Shipping["Shipping Service"]
  Email["Email Service"]
  
  Customer --> Web
  Web --> Catalog
  Web --> Cart
  Cart --> Checkout
  Checkout --> Inventory
  Checkout --> Payment
  Checkout --> Shipping
  Checkout --> Email
```

### Video Streaming Platform
```mermaid
block-beta
  columns 3
  
  Upload["Video Upload"]
  space
  space
  
  space
  Transcode["Transcoding Service"]
  space
  
  Storage["Object Storage"]
  CDN["CDN"]
  space
  
  space
  Player["Video Player"]
  space
  
  space
  Analytics["Analytics"]
  space
  
  Upload --> Transcode
  Transcode --> Storage
  Storage --> CDN
  CDN --> Player
  Player --> Analytics
```

### Authentication System
```mermaid
block-beta
  columns 3
  
  User["User"]
  space
  space
  
  space
  Login["Login Page"]
  space
  
  Auth["Auth Service"]
  Session["Session Manager"]
  space
  
  Database["User DB"]
  Cache["Session Cache"]
  MFA["MFA Service"]
  
  space
  App["Application"]
  space
  
  User --> Login
  Login --> Auth
  Auth --> Database
  Auth --> Session
  Auth --> MFA
  Session --> Cache
  Session --> App
```

### Machine Learning Pipeline
```mermaid
block-beta
  columns 4
  
  Data["Raw Data"]
  Clean["Data Cleaning"]
  Features["Feature Engineering"]
  Split["Train/Test Split"]
  
  space
  Train["Model Training"]
  Validate["Validation"]
  space
  
  space
  Model["Trained Model"]
  space
  space
  
  space
  Deploy["Model Deployment"]
  space
  space
  
  space
  Predict["Prediction Service"]
  Monitor["Monitoring"]
  space
  
  Data --> Clean
  Clean --> Features
  Features --> Split
  Split --> Train
  Train --> Validate
  Validate --> Model
  Model --> Deploy
  Deploy --> Predict
  Predict --> Monitor
```

## Tips

- Use `columns` to set grid layout
- Use `space` for empty cells
- Blocks are defined with labels in quotes
- Arrows show connections
- Organize logically top to bottom or left to right
- Group related components
- Show data flow with arrows
- Keep layout clean and readable
- Use consistent block sizing
- Label all components clearly
