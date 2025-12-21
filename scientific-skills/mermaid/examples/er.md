# Entity Relationship Diagrams

ER diagrams model database schemas and relationships between entities.

## Basic Syntax

```mermaid
erDiagram
    CUSTOMER {
        string id PK
        string name
        string email UK
    }
    
    ORDER {
        string id PK
        string customer_id FK
        date order_date
        decimal total
    }
    
    CUSTOMER ||--o{ ORDER : places
```

## Relationship Types

```mermaid
erDiagram
    A ||--|| B : "one to one"
    C ||--o{ D : "one to many"
    E }o--o{ F : "many to many"
    G }|--|{ H : "one or more to one or more"
```

## Cardinality Symbols

- `||` - Exactly one
- `o|` - Zero or one
- `}o` - Zero or more
- `}|` - One or more

## Attribute Types

```mermaid
erDiagram
    USER {
        uuid id PK "Primary Key"
        string username UK "Unique constraint"
        string email UK
        string password_hash
        datetime created_at
        datetime updated_at
        boolean is_active "Default true"
    }
```

## Common Key Constraints

- `PK` - Primary Key
- `FK` - Foreign Key
- `UK` - Unique Key

## Common Patterns

### E-Commerce Database
```mermaid
erDiagram
    CUSTOMER {
        int id PK
        string email UK
        string name
        string phone
        datetime created_at
    }
    
    ORDER {
        int id PK
        int customer_id FK
        datetime order_date
        string status
        decimal total_amount
    }
    
    ORDER_ITEM {
        int id PK
        int order_id FK
        int product_id FK
        int quantity
        decimal unit_price
    }
    
    PRODUCT {
        int id PK
        string sku UK
        string name
        text description
        decimal price
        int stock_quantity
    }
    
    CATEGORY {
        int id PK
        string name
        int parent_id FK
    }
    
    CUSTOMER ||--o{ ORDER : places
    ORDER ||--|{ ORDER_ITEM : contains
    PRODUCT ||--o{ ORDER_ITEM : "ordered in"
    CATEGORY ||--o{ PRODUCT : contains
    CATEGORY ||--o{ CATEGORY : "has subcategory"
```

### Blog Platform
```mermaid
erDiagram
    USER {
        uuid id PK
        string username UK
        string email UK
        string password_hash
        datetime created_at
    }
    
    POST {
        uuid id PK
        uuid author_id FK
        string title
        text content
        string slug UK
        datetime published_at
        boolean is_published
    }
    
    COMMENT {
        uuid id PK
        uuid post_id FK
        uuid user_id FK
        text content
        datetime created_at
    }
    
    TAG {
        uuid id PK
        string name UK
    }
    
    POST_TAG {
        uuid post_id FK
        uuid tag_id FK
    }
    
    USER ||--o{ POST : writes
    USER ||--o{ COMMENT : writes
    POST ||--o{ COMMENT : has
    POST ||--o{ POST_TAG : tagged
    TAG ||--o{ POST_TAG : tags
```

### Social Network
```mermaid
erDiagram
    USER {
        bigint id PK
        string username UK
        string email UK
        string bio
        string avatar_url
        datetime joined_at
    }
    
    POST {
        bigint id PK
        bigint user_id FK
        text content
        datetime created_at
        int likes_count
    }
    
    FRIENDSHIP {
        bigint user_id FK
        bigint friend_id FK
        datetime created_at
        string status
    }
    
    LIKE {
        bigint id PK
        bigint user_id FK
        bigint post_id FK
        datetime created_at
    }
    
    MESSAGE {
        bigint id PK
        bigint sender_id FK
        bigint receiver_id FK
        text content
        datetime sent_at
        boolean is_read
    }
    
    USER ||--o{ POST : creates
    USER ||--o{ FRIENDSHIP : "has friend"
    USER ||--o{ LIKE : gives
    POST ||--o{ LIKE : receives
    USER ||--o{ MESSAGE : sends
    USER ||--o{ MESSAGE : receives
```

### SaaS Application
```mermaid
erDiagram
    ORGANIZATION {
        uuid id PK
        string name
        string subdomain UK
        datetime created_at
    }
    
    USER {
        uuid id PK
        string email UK
        string name
        datetime created_at
    }
    
    ORGANIZATION_MEMBER {
        uuid id PK
        uuid organization_id FK
        uuid user_id FK
        string role
        datetime joined_at
    }
    
    PROJECT {
        uuid id PK
        uuid organization_id FK
        string name
        text description
        datetime created_at
    }
    
    TASK {
        uuid id PK
        uuid project_id FK
        uuid assigned_to FK
        string title
        text description
        string status
        datetime due_date
    }
    
    ORGANIZATION ||--|{ ORGANIZATION_MEMBER : has
    USER ||--o{ ORGANIZATION_MEMBER : "member of"
    ORGANIZATION ||--o{ PROJECT : owns
    PROJECT ||--o{ TASK : contains
    USER ||--o{ TASK : "assigned to"
```

### Inventory Management
```mermaid
erDiagram
    WAREHOUSE {
        int id PK
        string name
        string location
        string manager
    }
    
    PRODUCT {
        int id PK
        string sku UK
        string name
        decimal price
    }
    
    INVENTORY {
        int id PK
        int warehouse_id FK
        int product_id FK
        int quantity
        datetime last_updated
    }
    
    SUPPLIER {
        int id PK
        string name
        string contact_email
        string phone
    }
    
    PURCHASE_ORDER {
        int id PK
        int supplier_id FK
        int warehouse_id FK
        datetime order_date
        datetime expected_delivery
        string status
    }
    
    PURCHASE_ORDER_ITEM {
        int id PK
        int purchase_order_id FK
        int product_id FK
        int quantity
        decimal unit_cost
    }
    
    WAREHOUSE ||--|{ INVENTORY : stores
    PRODUCT ||--o{ INVENTORY : "stocked in"
    SUPPLIER ||--o{ PURCHASE_ORDER : supplies
    WAREHOUSE ||--o{ PURCHASE_ORDER : receives
    PURCHASE_ORDER ||--|{ PURCHASE_ORDER_ITEM : contains
    PRODUCT ||--o{ PURCHASE_ORDER_ITEM : "ordered as"
```

## Advanced Features

### Self-Referencing Relationship
```mermaid
erDiagram
    EMPLOYEE {
        int id PK
        string name
        int manager_id FK
    }
    
    EMPLOYEE ||--o{ EMPLOYEE : manages
```

### Multiple Relationships
```mermaid
erDiagram
    USER {
        int id PK
        string name
    }
    
    DOCUMENT {
        int id PK
        int created_by FK
        int updated_by FK
    }
    
    USER ||--o{ DOCUMENT : creates
    USER ||--o{ DOCUMENT : updates
```

## Tips

- Use meaningful entity names (UPPERCASE by convention)
- Include all foreign keys with FK annotation
- Mark unique constraints with UK
- Add comments for complex attributes
- Show cardinality accurately
- Keep entity count manageable (split large diagrams)
- Use junction tables for many-to-many relationships
- Consider normalization when designing
