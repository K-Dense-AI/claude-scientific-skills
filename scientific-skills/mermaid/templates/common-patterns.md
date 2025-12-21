# Common Mermaid Diagram Patterns

This file contains reusable templates and patterns for frequently used diagram scenarios.

## Software Development Workflows

### Feature Development Flow
```mermaid
flowchart TD
    Start([Feature Request]) --> Analysis[Analyze Requirements]
    Analysis --> Design[Create Design]
    Design --> Review{Design Review}
    Review -->|Approved| Dev[Development]
    Review -->|Changes Needed| Design
    Dev --> CodeReview[Code Review]
    CodeReview --> Tests{Tests Pass?}
    Tests -->|No| Dev
    Tests -->|Yes| Deploy[Deploy to Staging]
    Deploy --> QA[QA Testing]
    QA --> Prod{Ready for Prod?}
    Prod -->|No| Dev
    Prod -->|Yes| Release[Production Release]
    Release --> End([Done])
```

### Bug Fix Workflow
```mermaid
flowchart LR
    Report[Bug Reported] --> Triage{Severity}
    Triage -->|Critical| Immediate[Immediate Fix]
    Triage -->|High| Sprint[Add to Sprint]
    Triage -->|Low| Backlog[Add to Backlog]
    Immediate --> Fix[Develop Fix]
    Sprint --> Fix
    Fix --> Test[Test Fix]
    Test --> Deploy[Deploy]
    Deploy --> Verify[Verify Resolution]
    Verify --> Close[Close Ticket]
```

## Authentication Patterns

### OAuth 2.0 Flow
```mermaid
sequenceDiagram
    actor User
    participant App
    participant AuthServer
    participant ResourceServer
    
    User->>+App: Click Login
    App->>+AuthServer: Authorization Request
    AuthServer->>User: Login Page
    User->>AuthServer: Credentials
    AuthServer->>-App: Authorization Code
    App->>+AuthServer: Exchange Code for Token
    AuthServer->>-App: Access Token
    App->>+ResourceServer: API Request + Token
    ResourceServer->>ResourceServer: Validate Token
    ResourceServer->>-App: Protected Resource
    App->>-User: Display Data
```

### JWT Authentication
```mermaid
sequenceDiagram
    participant Client
    participant API
    participant Auth
    participant DB
    
    Client->>+API: POST /login (credentials)
    API->>+Auth: Validate credentials
    Auth->>+DB: Query user
    DB->>-Auth: User data
    
    alt Valid credentials
        Auth->>Auth: Generate JWT
        Auth->>-API: JWT token
        API->>-Client: 200 OK + token
    else Invalid
        Auth->>-API: Invalid credentials
        API->>-Client: 401 Unauthorized
    end
    
    Note over Client,API: Subsequent requests
    Client->>+API: GET /protected (+ JWT)
    API->>API: Verify JWT
    API->>-Client: Protected data
```

## Database Patterns

### E-Commerce Database Schema
```mermaid
erDiagram
    CUSTOMER ||--o{ ORDER : places
    CUSTOMER {
        int id PK
        string email UK
        string name
        datetime created_at
    }
    
    ORDER ||--|{ ORDER_ITEM : contains
    ORDER {
        int id PK
        int customer_id FK
        datetime order_date
        decimal total
        string status
    }
    
    PRODUCT ||--o{ ORDER_ITEM : "ordered in"
    PRODUCT {
        int id PK
        string sku UK
        string name
        decimal price
        int stock
    }
    
    ORDER_ITEM {
        int id PK
        int order_id FK
        int product_id FK
        int quantity
        decimal unit_price
    }
```

### User Management Schema
```mermaid
erDiagram
    USER ||--o{ SESSION : has
    USER ||--o{ ROLE_ASSIGNMENT : has
    ROLE ||--o{ ROLE_ASSIGNMENT : assigned
    ROLE ||--o{ PERMISSION : has
    
    USER {
        uuid id PK
        string email UK
        string password_hash
        boolean is_active
        datetime last_login
    }
    
    SESSION {
        uuid id PK
        uuid user_id FK
        string token
        datetime expires_at
    }
    
    ROLE {
        int id PK
        string name UK
        string description
    }
    
    ROLE_ASSIGNMENT {
        uuid user_id FK
        int role_id FK
        datetime assigned_at
    }
    
    PERMISSION {
        int id PK
        int role_id FK
        string resource
        string action
    }
```

## Architecture Patterns

### Microservices Architecture
```mermaid
flowchart TB
    subgraph Client Layer
        Web[Web App]
        Mobile[Mobile App]
    end
    
    subgraph API Layer
        Gateway[API Gateway]
    end
    
    subgraph Services
        Auth[Auth Service]
        User[User Service]
        Product[Product Service]
        Order[Order Service]
        Payment[Payment Service]
    end
    
    subgraph Data Layer
        AuthDB[(Auth DB)]
        UserDB[(User DB)]
        ProductDB[(Product DB)]
        OrderDB[(Order DB)]
    end
    
    subgraph Infrastructure
        Cache[(Redis)]
        Queue[Message Queue]
    end
    
    Web --> Gateway
    Mobile --> Gateway
    Gateway --> Auth
    Gateway --> User
    Gateway --> Product
    Gateway --> Order
    
    Auth --> AuthDB
    User --> UserDB
    Product --> ProductDB
    Order --> OrderDB
    
    Order --> Payment
    Order --> Queue
    Product --> Cache
```

### Layered Architecture
```mermaid
flowchart TD
    subgraph Presentation Layer
        UI[User Interface]
        API[REST API]
    end
    
    subgraph Business Layer
        BL[Business Logic]
        Validation[Validation]
        Rules[Business Rules]
    end
    
    subgraph Data Access Layer
        Repo[Repositories]
        ORM[ORM/Data Mappers]
    end
    
    subgraph Database
        DB[(Database)]
    end
    
    UI --> BL
    API --> BL
    BL --> Validation
    BL --> Rules
    BL --> Repo
    Repo --> ORM
    ORM --> DB
```

## State Management

### Order Lifecycle
```mermaid
stateDiagram-v2
    [*] --> Created
    Created --> Pending : Submit Order
    Pending --> Processing : Payment Confirmed
    Pending --> Cancelled : Cancel Order
    
    Processing --> Shipped : Items Shipped
    Processing --> Cancelled : Out of Stock
    
    Shipped --> Delivered : Delivery Confirmed
    Shipped --> Returned : Return Initiated
    
    Delivered --> Returned : Return Request
    Delivered --> Completed : Return Window Closed
    
    Returned --> Refunded : Refund Processed
    Cancelled --> Refunded : Refund Processed
    
    Completed --> [*]
    Refunded --> [*]
```

### User Account States
```mermaid
stateDiagram-v2
    [*] --> Registered
    Registered --> Active : Email Verified
    Registered --> Pending : Awaiting Verification
    
    Pending --> Active : Verify Email
    Pending --> Expired : Verification Timeout
    
    Active --> Suspended : Policy Violation
    Active --> Locked : Too Many Login Attempts
    Active --> Inactive : No Activity
    
    Suspended --> Active : Appeal Approved
    Locked --> Active : Password Reset
    Inactive --> Active : User Login
    
    Active --> Deleted : User Request
    Suspended --> Deleted : Admin Action
    Expired --> Deleted : Cleanup Job
    
    Deleted --> [*]
```

## Class Diagram Patterns

### Repository Pattern
```mermaid
classDiagram
    class IRepository~T~ {
        <<interface>>
        +findById(id) T
        +findAll() List~T~
        +save(entity T) void
        +delete(id) void
    }
    
    class BaseRepository~T~ {
        <<abstract>>
        #db Database
        +findById(id) T
        +findAll() List~T~
        +save(entity T) void
        +delete(id) void
    }
    
    class UserRepository {
        +findByEmail(email) User
        +findByRole(role) List~User~
    }
    
    class ProductRepository {
        +findBySku(sku) Product
        +findInStock() List~Product~
    }
    
    class User {
        +String id
        +String email
        +String name
    }
    
    class Product {
        +String id
        +String sku
        +decimal price
    }
    
    IRepository~T~ <|.. BaseRepository~T~
    BaseRepository~T~ <|-- UserRepository
    BaseRepository~T~ <|-- ProductRepository
    UserRepository --> User
    ProductRepository --> Product
```

### Service Layer Pattern
```mermaid
classDiagram
    class Controller {
        -service Service
        +handleRequest(req) Response
    }
    
    class Service {
        -repository Repository
        -validator Validator
        +createEntity(data) Entity
        +updateEntity(id, data) Entity
        +deleteEntity(id) void
    }
    
    class Repository {
        -db Database
        +save(entity) void
        +findById(id) Entity
        +delete(id) void
    }
    
    class Validator {
        +validate(data) ValidationResult
    }
    
    class Entity {
        +String id
        +Map attributes
    }
    
    Controller --> Service
    Service --> Repository
    Service --> Validator
    Repository --> Entity
```

## CI/CD Pipeline

### Standard Pipeline
```mermaid
flowchart LR
    Commit[Git Commit] --> Build[Build]
    Build --> UnitTest[Unit Tests]
    UnitTest --> Lint[Linting]
    Lint --> IntTest[Integration Tests]
    IntTest --> Security[Security Scan]
    Security --> Package[Package]
    Package --> StageDeploy[Deploy to Staging]
    StageDeploy --> E2E[E2E Tests]
    E2E --> Approve{Manual Approval}
    Approve -->|Yes| ProdDeploy[Deploy to Production]
    Approve -->|No| End([End])
    ProdDeploy --> Monitor[Monitor]
    Monitor --> End
```

## API Patterns

### REST API Standard Flow
```mermaid
sequenceDiagram
    participant Client
    participant Gateway
    participant Service
    participant Cache
    participant Database
    
    Client->>+Gateway: HTTP Request
    Gateway->>Gateway: Authenticate
    Gateway->>Gateway: Rate Limit Check
    
    Gateway->>+Service: Forward Request
    Service->>+Cache: Check Cache
    
    alt Cache Hit
        Cache->>-Service: Cached Data
    else Cache Miss
        Service->>+Database: Query
        Database->>-Service: Data
        Service->>Cache: Update Cache
    end
    
    Service->>-Gateway: Response
    Gateway->>-Client: HTTP Response
```

## Tips for Using Templates

1. **Customize**: Adjust templates to match your specific needs
2. **Simplify**: Remove unnecessary elements for clarity
3. **Extend**: Add domain-specific elements as needed
4. **Validate**: Always validate with mmdc after customization
5. **Document**: Add comments for complex logic
6. **Consistent**: Use consistent naming conventions
7. **Focused**: Keep diagrams focused on one aspect
8. **Updated**: Maintain templates as patterns evolve
