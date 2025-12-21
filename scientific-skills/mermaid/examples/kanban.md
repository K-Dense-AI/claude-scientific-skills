# Kanban Diagrams

Kanban boards visualize workflow and task status.

**Note:** Kanban diagram support in Mermaid is new and experimental.

## Basic Syntax

```mermaid
kanban
  Todo
    Task 1
    Task 2
  In Progress
    Task 3
  Done
    Task 4
    Task 5
```

## Common Patterns

### Software Development Board
```mermaid
kanban
  Backlog
    User authentication
    Payment integration
    Mobile app design
    API documentation
  To Do
    Dashboard redesign
    Bug fixes batch 1
  In Progress
    User profile page
    Database optimization
  Code Review
    Shopping cart feature
  Testing
    Email notifications
  Done
    Login page
    Homepage redesign
    Search functionality
```

### Content Production
```mermaid
kanban
  Ideas
    Blog post: AI trends
    Video: Product tutorial
    Podcast: Industry interview
  Writing
    Blog post: DevOps guide
    Case study: Customer X
  Review
    Blog post: Cloud migration
  Editing
    Video: Getting started
  Scheduled
    Social media campaign
    Newsletter issue 45
  Published
    Blog post: Web performance
    Video: Feature overview
```

### Marketing Campaign
```mermaid
kanban
  Planning
    Q2 campaign strategy
    Budget allocation
    Target audience research
  Design
    Email templates
    Landing page mockups
  Development
    Campaign landing page
    Email automation
  Review
    Ad copy review
    Creative approval
  Launched
    Email campaign wave 1
    Social media ads
  Completed
    Q1 campaign wrap-up
    Webinar series
```

### Support Ticket Board
```mermaid
kanban
  New
    Login issues - User A
    Feature request - User B
    Billing question - User C
  Assigned
    App crashes on iOS
    Email not sending
  In Progress
    Performance degradation
    Data export bug
  Waiting on Customer
    Clarification needed
  Resolved
    Password reset issue
    API timeout fixed
    UI alignment bug
  Closed
    How-to question
    Feature explanation
```

### Product Roadmap
```mermaid
kanban
  Ideas
    AI-powered search
    Dark mode
    Collaborative editing
    Mobile offline mode
  Planned
    Advanced analytics
    SSO integration
  In Development
    API v2
    New dashboard
  Testing
    File upload feature
    Notification system
  Releasing
    Performance improvements
  Released
    User profiles
    Comments feature
    Export functionality
```

### Sprint Board
```mermaid
kanban
  Sprint Backlog
    Story: User onboarding
    Story: Payment flow
    Story: Admin panel
  To Do
    Task: Design login screen
    Task: Create API endpoints
  In Progress
    Task: Implement auth
    Task: Write unit tests
  Code Review
    Task: User dashboard
    Task: Settings page
  QA
    Task: Checkout flow
  Done
    Task: Homepage layout
    Task: Database schema
```

### Sales Pipeline
```mermaid
kanban
  Leads
    Company A
    Company B
    Company C
  Qualified
    Company D
    Company E
  Demo Scheduled
    Company F
    Company G
  Proposal Sent
    Company H
  Negotiation
    Company I
  Closed Won
    Company J
    Company K
  Closed Lost
    Company L
```

### Bug Tracking
```mermaid
kanban
  Reported
    Critical: Data loss
    High: Login broken
    Medium: UI glitch
    Low: Typo in message
  Triaged
    High: API timeout
    Medium: Slow query
  Assigned
    Critical: Security issue
    High: Payment failure
  In Progress
    High: Memory leak
  Testing
    Medium: Form validation
  Verified
    High: Cache bug
  Closed
    Medium: CSS issue
    Low: Console warning
```

### Event Planning
```mermaid
kanban
  Ideas
    Annual conference
    Workshop series
    Networking mixer
  Planning
    Tech meetup Q2
    Webinar on AI
  Vendor Selection
    Catering quotes
    Venue booking
  Confirmed
    Speaker confirmations
    Sponsor agreements
  Marketing
    Event page live
    Email invitations
  Completed
    Q1 Workshop
    Product launch event
```

### HR Recruitment
```mermaid
kanban
  Sourcing
    Senior Developer
    Product Designer
    Marketing Manager
  Phone Screen
    Candidate A - Dev
    Candidate B - Design
  Interview
    Candidate C - Dev
    Candidate D - PM
  Offer
    Candidate E - Design
  Onboarding
    New Developer
  Hired
    Frontend Engineer
    UX Designer
```

## Tips

- Columns represent workflow stages
- Cards represent work items/tasks
- Move cards left to right as work progresses
- Keep columns focused on clear stages
- Limit work in progress (WIP limits)
- Use for visual workflow management
- Update regularly to reflect current status
- Common columns: To Do, In Progress, Done
- Customize columns for your workflow
- Include all relevant stages
