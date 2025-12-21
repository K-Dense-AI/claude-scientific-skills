# User Journey Diagrams

User journey maps track user experience across tasks and touchpoints.

## Basic Syntax

```mermaid
journey
    title My Working Day
    section Go to work
      Make tea: 5: Me
      Go upstairs: 3: Me
      Do work: 1: Me, Cat
    section Go home
      Go downstairs: 5: Me
      Sit down: 5: Me
```

## Satisfaction Scores

- `1` - Very Dissatisfied
- `2` - Dissatisfied
- `3` - Neutral
- `4` - Satisfied
- `5` - Very Satisfied

## Common Patterns

### E-Commerce Purchase Journey
```mermaid
journey
    title Online Shopping Experience
    section Discovery
      See social media ad: 4: Customer
      Click through to website: 4: Customer
      Browse products: 5: Customer
    section Evaluation
      Read product details: 4: Customer
      Check reviews: 5: Customer
      Compare prices: 3: Customer
    section Purchase
      Add to cart: 5: Customer
      Enter shipping info: 3: Customer
      Payment process: 2: Customer
      Confirm order: 4: Customer
    section Post-Purchase
      Receive confirmation email: 5: Customer, System
      Track shipment: 4: Customer
      Receive product: 5: Customer
      Leave review: 4: Customer
```

### SaaS Onboarding
```mermaid
journey
    title New User Onboarding
    section Signup
      Land on homepage: 4: User
      Click signup: 5: User
      Fill registration form: 3: User
      Verify email: 2: User
    section Setup
      Welcome tutorial: 4: User, System
      Configure preferences: 3: User
      Connect integrations: 2: User
    section First Use
      Create first project: 5: User
      Invite team members: 4: User
      Complete sample task: 5: User
    section Activation
      Daily usage: 4: User
      Explore features: 5: User
      Contact support: 3: User, Support
      Upgrade to paid: 5: User
```

### Customer Support Journey
```mermaid
journey
    title Support Ticket Resolution
    section Problem Discovery
      Encounter issue: 1: Customer
      Search help docs: 2: Customer
      Can't find solution: 1: Customer
    section Contact Support
      Find contact form: 3: Customer
      Submit ticket: 3: Customer
      Receive auto-reply: 3: Customer, System
    section Resolution
      Agent reviews ticket: 4: Agent
      Agent responds: 4: Agent, Customer
      Back and forth: 3: Customer, Agent
      Problem solved: 5: Customer, Agent
    section Follow-up
      Satisfaction survey: 4: Customer
      Close ticket: 5: Agent
```

### Job Application Process
```mermaid
journey
    title Candidate Application Journey
    section Application
      Find job posting: 4: Candidate
      Read description: 5: Candidate
      Prepare resume: 3: Candidate
      Submit application: 4: Candidate
    section Initial Review
      Application reviewed: 3: Recruiter
      Initial screening call: 4: Candidate, Recruiter
      Schedule interview: 5: Candidate, Recruiter
    section Interview
      Technical interview: 3: Candidate, Engineer
      Culture fit interview: 4: Candidate, Manager
      Final round: 3: Candidate, Team
    section Decision
      Wait for decision: 2: Candidate
      Receive offer: 5: Candidate, Recruiter
      Negotiate terms: 4: Candidate, Recruiter
      Accept offer: 5: Candidate
```

### Restaurant Dining Experience
```mermaid
journey
    title Restaurant Visit
    section Arrival
      Make reservation: 5: Customer
      Arrive at restaurant: 4: Customer
      Wait for table: 2: Customer
      Get seated: 5: Customer, Host
    section Ordering
      Review menu: 4: Customer
      Ask questions: 5: Customer, Waiter
      Place order: 5: Customer, Waiter
      Receive drinks: 4: Customer, Waiter
    section Dining
      Food arrives: 5: Customer, Waiter
      Enjoy meal: 5: Customer
      Request check: 4: Customer, Waiter
    section Payment
      Review bill: 4: Customer
      Pay: 4: Customer
      Leave tip: 5: Customer
      Exit: 5: Customer
```

### Mobile App User Flow
```mermaid
journey
    title First Time App User
    section Download
      See app in store: 4: User
      Read reviews: 5: User
      Download app: 5: User
      Install: 3: User, System
    section Onboarding
      Open app: 5: User
      See splash screen: 4: User
      Skip tutorial: 2: User
      Grant permissions: 2: User
    section First Session
      Explore interface: 3: User
      Try main feature: 4: User
      Get confused: 2: User
      Find help: 3: User
    section Retention
      Return next day: 4: User
      Complete task: 5: User
      Enable notifications: 3: User
      Recommend to friend: 5: User
```

### Banking Customer Journey
```mermaid
journey
    title Opening Bank Account
    section Research
      Compare banks: 3: Customer
      Check requirements: 4: Customer
      Read reviews: 4: Customer
    section Application
      Start online application: 4: Customer
      Upload documents: 2: Customer
      Identity verification: 2: Customer, System
      Submit application: 3: Customer
    section Approval
      Application review: 4: Bank
      Request additional info: 2: Customer, Bank
      Approval notification: 5: Customer, System
    section Setup
      Choose account type: 5: Customer
      Set up online banking: 3: Customer
      Order debit card: 4: Customer
      Make initial deposit: 5: Customer
    section Usage
      First transaction: 5: Customer
      Set up alerts: 4: Customer
      Contact support: 3: Customer, Support
```

### Healthcare Appointment
```mermaid
journey
    title Doctor's Appointment
    section Booking
      Feel symptoms: 1: Patient
      Search for doctor: 3: Patient
      Call to book: 2: Patient, Receptionist
      Receive confirmation: 4: Patient, System
    section Preparation
      Receive reminder: 4: Patient, System
      Fill pre-visit forms: 2: Patient
      Travel to clinic: 3: Patient
    section Visit
      Check in: 3: Patient, Receptionist
      Wait in lobby: 2: Patient
      Consultation: 4: Patient, Doctor
      Receive prescription: 5: Patient, Doctor
    section Follow-up
      Pay bill: 2: Patient, Billing
      Pick up medication: 4: Patient
      Schedule follow-up: 4: Patient, Receptionist
```

### Event Registration
```mermaid
journey
    title Conference Registration
    section Awareness
      See announcement: 4: Attendee
      Visit website: 5: Attendee
      Review agenda: 5: Attendee
    section Registration
      Select ticket type: 4: Attendee
      Fill registration form: 3: Attendee
      Payment: 3: Attendee
      Confirmation email: 5: Attendee, System
    section Preparation
      Add to calendar: 5: Attendee
      Book travel: 2: Attendee
      Review speakers: 5: Attendee
    section Event Day
      Check in: 4: Attendee, Staff
      Attend sessions: 5: Attendee
      Network: 5: Attendee, Other
      Collect swag: 5: Attendee
    section Post-Event
      Receive recording links: 5: Attendee, System
      Complete survey: 4: Attendee
      Connect on LinkedIn: 5: Attendee
```

## Multiple Actors

```mermaid
journey
    title Product Development
    section Planning
      Define requirements: 4: PM, Team
      Create designs: 5: Designer
      Review designs: 4: PM, Developer
    section Development
      Write code: 4: Developer
      Code review: 5: Developer, Lead
      QA testing: 3: QA, Developer
    section Launch
      Deploy: 4: DevOps
      Monitor: 4: DevOps, Developer
      Support users: 3: Support, Developer
```

## Tips

- Scores range from 1 (worst) to 5 (best)
- Include multiple actors when relevant
- Break journey into logical sections
- Focus on key touchpoints
- Identify pain points (low scores)
- Highlight moments of delight (high scores)
- Use realistic scenarios
- Consider the entire end-to-end experience
