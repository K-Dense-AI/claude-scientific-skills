# Flowchart Diagrams

Flowcharts visualize processes, algorithms, and decision flows.

## Syntax

```mermaid
flowchart TD
    Start --> Process
    Process --> Decision{Question?}
    Decision -->|Yes| ActionA
    Decision -->|No| ActionB
    ActionA --> End
    ActionB --> End
```

## Orientations

- `flowchart TD` or `flowchart TB` - Top to Bottom
- `flowchart BT` - Bottom to Top
- `flowchart LR` - Left to Right
- `flowchart RL` - Right to Left

## Node Shapes

```mermaid
flowchart LR
    A[Rectangle]
    B(Rounded)
    C([Stadium])
    D[[Subroutine]]
    E[(Database)]
    F((Circle))
    G>Asymmetric]
    H{Diamond}
    I{{Hexagon}}
    J[/Parallelogram/]
    K[\Reverse Parallelogram\]
    L[/Trapezoid\]
    M[\Reverse Trapezoid/]
```

## Arrow Types

```mermaid
flowchart LR
    A --> B     %% Solid arrow
    C -.-> D    %% Dotted arrow
    E ==> F     %% Thick arrow
    G --- H     %% Solid line (no arrow)
    I -.- J     %% Dotted line
    K === L     %% Thick line
```

## Arrow Labels

```mermaid
flowchart LR
    A -->|Label text| B
    C -.->|Optional| D
    E ==>|Required| F
```

## Subgraphs

```mermaid
flowchart TB
    subgraph Frontend
        A[React App]
        B[Vue App]
    end
    subgraph Backend
        C[API Server]
        D[Database]
    end
    A --> C
    B --> C
    C --> D
```

## Common Patterns

### Simple Process Flow
```mermaid
flowchart TD
    Start([Start]) --> Input[Get User Input]
    Input --> Validate{Valid?}
    Validate -->|Yes| Process[Process Data]
    Validate -->|No| Error[Show Error]
    Error --> Input
    Process --> Save[(Save to DB)]
    Save --> End([End])
```

### Decision Tree
```mermaid
flowchart TD
    Start --> Q1{Is User Logged In?}
    Q1 -->|Yes| Q2{Has Permission?}
    Q1 -->|No| Login[Redirect to Login]
    Q2 -->|Yes| Access[Grant Access]
    Q2 -->|No| Deny[Access Denied]
```

### Error Handling Flow
```mermaid
flowchart LR
    Try[Try Operation] --> Check{Success?}
    Check -->|Yes| Continue[Continue]
    Check -->|No| Log[Log Error]
    Log --> Retry{Retry?}
    Retry -->|Yes| Try
    Retry -->|No| Fail[Fail Gracefully]
```

## Styling

```mermaid
flowchart LR
    A[Normal]
    B[Styled]
    style B fill:#f9f,stroke:#333,stroke-width:4px
```

### Class Definitions
```mermaid
flowchart TD
    A[Start]:::important --> B[Process]:::normal
    B --> C[End]:::important
    
    classDef important fill:#f96,stroke:#333,stroke-width:2px
    classDef normal fill:#bbf,stroke:#333,stroke-width:1px
```

## Tips

- Use meaningful IDs (not just A, B, C)
- Quote labels with special characters: `Node["Text with spaces"]`
- Maximum clarity: one decision per diamond
- Group related nodes with subgraphs
- Keep flows unidirectional when possible
