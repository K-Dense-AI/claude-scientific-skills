# State Diagrams

State diagrams model state machines and lifecycle transitions.

## Basic Syntax

```mermaid
stateDiagram-v2
    [*] --> Idle
    Idle --> Processing
    Processing --> Idle
    Processing --> [*]
```

## States with Descriptions

```mermaid
stateDiagram-v2
    [*] --> New
    New --> Active : Activate
    Active --> Suspended : Suspend
    Suspended --> Active : Resume
    Active --> Closed : Close
    Closed --> [*]
```

## Composite States

```mermaid
stateDiagram-v2
    [*] --> Active
    
    state Active {
        [*] --> Running
        Running --> Paused
        Paused --> Running
    }
    
    Active --> Stopped
    Stopped --> [*]
```

## Concurrent States

```mermaid
stateDiagram-v2
    [*] --> Working
    
    state Working {
        [*] --> TaskA
        --
        [*] --> TaskB
    }
    
    Working --> Done
    Done --> [*]
```

## Choice (Conditional)

```mermaid
stateDiagram-v2
    [*] --> Check
    Check --> state1 : if x > 0
    Check --> state2 : if x <= 0
    state1 --> [*]
    state2 --> [*]
```

## Fork and Join

```mermaid
stateDiagram-v2
    [*] --> Start
    Start --> fork_state
    
    fork_state --> Task1
    fork_state --> Task2
    
    Task1 --> join_state
    Task2 --> join_state
    
    join_state --> End
    End --> [*]
```

## Notes

```mermaid
stateDiagram-v2
    [*] --> Active
    Active --> Inactive
    
    note right of Active
        This is an important state
        where processing occurs
    end note
    
    note left of Inactive : System is idle
```

## Common Patterns

### Order Lifecycle
```mermaid
stateDiagram-v2
    [*] --> Created
    Created --> Pending : Submit
    Pending --> Processing : Confirm
    Processing --> Shipped : Ship
    Shipped --> Delivered : Deliver
    Delivered --> [*]
    
    Pending --> Cancelled : Cancel
    Processing --> Cancelled : Cancel
    Cancelled --> [*]
```

### User Authentication States
```mermaid
stateDiagram-v2
    [*] --> LoggedOut
    LoggedOut --> LoggingIn : Login Attempt
    
    LoggingIn --> LoggedIn : Success
    LoggingIn --> LoggedOut : Failure
    
    LoggedIn --> RefreshingToken : Token Expired
    RefreshingToken --> LoggedIn : Success
    RefreshingToken --> LoggedOut : Failure
    
    LoggedIn --> LoggedOut : Logout
    LoggedOut --> [*]
```

### Connection States
```mermaid
stateDiagram-v2
    [*] --> Disconnected
    Disconnected --> Connecting : Connect
    
    Connecting --> Connected : Success
    Connecting --> Failed : Timeout/Error
    
    Failed --> Reconnecting : Retry
    Reconnecting --> Connected : Success
    Reconnecting --> Failed : Error
    
    Connected --> Disconnecting : Close
    Disconnecting --> Disconnected
    
    Disconnected --> [*]
```

### Document Workflow
```mermaid
stateDiagram-v2
    [*] --> Draft
    Draft --> Review : Submit
    
    state Review {
        [*] --> Pending
        Pending --> InReview : Assign Reviewer
        InReview --> Approved : Approve
        InReview --> Rejected : Reject
    }
    
    Review --> Draft : if Rejected
    Review --> Published : if Approved
    Published --> Archived
    Archived --> [*]
```

### Media Player States
```mermaid
stateDiagram-v2
    [*] --> Stopped
    Stopped --> Playing : Play
    Playing --> Paused : Pause
    Paused --> Playing : Resume
    Playing --> Stopped : Stop
    Paused --> Stopped : Stop
    
    state Playing {
        [*] --> Loading
        Loading --> Buffering : if buffer low
        Loading --> Streaming : if ready
        Buffering --> Streaming : buffer filled
        Streaming --> Buffering : if buffer low
    }
```

### Payment Processing
```mermaid
stateDiagram-v2
    [*] --> Initiated
    Initiated --> Authorizing : Authorize
    
    Authorizing --> Authorized : Success
    Authorizing --> Failed : Declined
    
    Authorized --> Capturing : Capture
    Capturing --> Completed : Success
    Capturing --> Failed : Error
    
    Completed --> Refunding : Refund Request
    Refunding --> Refunded : Success
    Refunding --> Failed : Error
    
    Failed --> [*]
    Refunded --> [*]
```

## State Actions

```mermaid
stateDiagram-v2
    [*] --> Active
    Active : entry / initialize()
    Active : do / processData()
    Active : exit / cleanup()
    Active --> Inactive
    Inactive --> [*]
```

## Direction

```mermaid
stateDiagram-v2
    direction LR
    [*] --> A
    A --> B
    B --> C
    C --> [*]
```

## Styling

```mermaid
stateDiagram-v2
    [*] --> Important
    Important --> Normal
    Normal --> [*]
    
    style Important fill:#f96,stroke:#333,stroke-width:4px
```

## Tips

- `[*]` represents start/end states
- Use composite states for sub-states
- Use `--` for concurrent regions
- Transition labels show triggering events/conditions
- Keep state names concise and descriptive
- Model only relevant states (avoid over-complication)
- Use choice nodes for conditional transitions
