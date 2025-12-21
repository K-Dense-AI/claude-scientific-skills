# Sequence Diagrams

Sequence diagrams show interactions between participants over time.

## Basic Syntax

```mermaid
sequenceDiagram
    participant A as Alice
    participant B as Bob
    A->>B: Hello Bob!
    B->>A: Hi Alice!
```

## Participants

```mermaid
sequenceDiagram
    participant Client
    participant API
    participant DB as Database
    actor User
```

## Arrow Types

```mermaid
sequenceDiagram
    A->>B: Solid arrow (message)
    A-->>B: Dotted arrow (response)
    A-)B: Open arrow
    A--)B: Dotted open arrow
    A-xB: Cross (failure)
    A--xB: Dotted cross
```

## Activation Boxes

```mermaid
sequenceDiagram
    Client->>+API: Request
    API->>+Database: Query
    Database-->>-API: Result
    API-->>-Client: Response
```

## Notes

```mermaid
sequenceDiagram
    participant A
    participant B
    Note left of A: Note on left
    Note right of B: Note on right
    Note over A: Note over one
    Note over A,B: Note spanning both
```

## Loops

```mermaid
sequenceDiagram
    Client->>API: Start
    loop Every 5 seconds
        API->>Server: Poll
        Server-->>API: Status
    end
    API-->>Client: Complete
```

## Alternatives (if/else)

```mermaid
sequenceDiagram
    User->>API: Login Request
    alt Valid Credentials
        API-->>User: Success Token
    else Invalid Credentials
        API-->>User: Error Message
    end
```

## Optional Blocks

```mermaid
sequenceDiagram
    User->>System: Action
    opt Cache Available
        System->>Cache: Get Data
        Cache-->>System: Cached Data
    end
    System-->>User: Response
```

## Parallel Processing

```mermaid
sequenceDiagram
    Client->>Server: Request
    par Fetch User Data
        Server->>DB1: Query Users
    and Fetch Settings
        Server->>DB2: Query Settings
    end
    DB1-->>Server: User Data
    DB2-->>Server: Settings
    Server-->>Client: Combined Response
```

## Common Patterns

### Authentication Flow
```mermaid
sequenceDiagram
    actor User
    participant App
    participant Auth
    participant DB
    
    User->>+App: Login(username, password)
    App->>+Auth: Authenticate
    Auth->>+DB: Query User
    DB-->>-Auth: User Record
    
    alt Valid Password
        Auth->>Auth: Generate Token
        Auth-->>App: Token
        App-->>User: Success
    else Invalid Password
        Auth-->>App: Error
        App-->>-User: Invalid Credentials
    end
```

### API Request Flow
```mermaid
sequenceDiagram
    participant Client
    participant API
    participant Cache
    participant DB
    
    Client->>+API: GET /users/123
    API->>+Cache: Check Cache
    
    alt Cache Hit
        Cache-->>API: User Data
    else Cache Miss
        Cache-->>API: Not Found
        API->>+DB: SELECT * FROM users
        DB-->>-API: User Record
        API->>Cache: Store in Cache
    end
    
    API-->>-Client: User Response
```

### Microservice Communication
```mermaid
sequenceDiagram
    participant Client
    participant Gateway
    participant UserService
    participant OrderService
    participant PaymentService
    
    Client->>+Gateway: Place Order
    Gateway->>+UserService: Validate User
    UserService-->>-Gateway: User Valid
    
    Gateway->>+OrderService: Create Order
    OrderService-->>Gateway: Order ID
    
    par Process Payment
        Gateway->>+PaymentService: Charge Card
        PaymentService-->>-Gateway: Payment Confirmed
    and Update Inventory
        OrderService->>OrderService: Reserve Items
    end
    
    Gateway-->>-Client: Order Complete
```

## Background Highlighting

```mermaid
sequenceDiagram
    participant A
    participant B
    
    rect rgb(200, 220, 250)
        note right of A: Important Section
        A->>B: Message 1
        B->>A: Response 1
    end
    
    rect rgba(255, 200, 200, 0.5)
        note right of A: Error Handling
        A->>B: Message 2
        B-xA: Error
    end
```

## Autonumbering

```mermaid
sequenceDiagram
    autonumber
    Alice->>Bob: Message 1
    Bob->>Alice: Message 2
    Alice->>Bob: Message 3
```

## Tips

- Use `participant` to control order of appearance
- Use `actor` for human participants
- Activate (+) and deactivate (-) to show when components are processing
- Use `alt/else` for conditional flows
- Use `par/and` for concurrent operations
- Add `autonumber` for step numbering
