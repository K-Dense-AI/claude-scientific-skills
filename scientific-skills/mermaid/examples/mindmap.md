# Mindmaps

Mindmaps organize ideas hierarchically around a central concept.

## Basic Syntax

```mermaid
mindmap
  root((Central Idea))
    Topic 1
      Subtopic 1.1
      Subtopic 1.2
    Topic 2
      Subtopic 2.1
      Subtopic 2.2
```

## Node Shapes

```mermaid
mindmap
  root((Root - Circle))
    [Square]
    (Rounded)
    ((Circle))
    ))Cloud((
    {{Hexagon}}
```

## Common Patterns

### Project Planning
```mermaid
mindmap
  root((New Product Launch))
    Research
      Market Analysis
        Competitors
        Target Audience
        Trends
      User Needs
        Surveys
        Interviews
        Analytics
    Development
      Backend
        Database Design
        API Development
        Security
      Frontend
        UI Design
        User Experience
        Mobile App
      Testing
        Unit Tests
        Integration Tests
        User Testing
    Marketing
      Branding
        Logo
        Color Scheme
        Messaging
      Channels
        Social Media
        Email Campaign
        Content Marketing
      Launch Event
    Operations
      Support
        Documentation
        Training
        Help Desk
      Infrastructure
        Hosting
        Monitoring
        Backups
```

### Learning Path
```mermaid
mindmap
  root((Web Development))
    Frontend
      HTML
        Semantic Tags
        Forms
        Accessibility
      CSS
        Flexbox
        Grid
        Animations
        Preprocessors
          Sass
          Less
      JavaScript
        ES6+
        DOM Manipulation
        Async/Await
      Frameworks
        React
        Vue
        Angular
    Backend
      Languages
        Node.js
        Python
        Java
        Go
      Databases
        SQL
          PostgreSQL
          MySQL
        NoSQL
          MongoDB
          Redis
      APIs
        REST
        GraphQL
        WebSockets
    DevOps
      Version Control
        Git
        GitHub
      CI/CD
        Jenkins
        GitHub Actions
      Cloud
        AWS
        Azure
        Docker
```

### Business Strategy
```mermaid
mindmap
  root((Business Growth))
    Revenue Streams
      Products
        SaaS Platform
        Mobile App
        Enterprise Solution
      Services
        Consulting
        Training
        Support
    Customer Acquisition
      Marketing
        Content
        SEO
        Paid Ads
      Sales
        Outbound
        Partnerships
        Referrals
    Operations
      Team
        Hiring
        Training
        Culture
      Processes
        Workflows
        Automation
        Quality Control
    Technology
      Infrastructure
        Scalability
        Security
        Monitoring
      Innovation
        R&D
        New Features
        AI/ML
```

### Personal Goals
```mermaid
mindmap
  root((2024 Goals))
    Career
      Skills
        Programming
          Python
          TypeScript
        Soft Skills
          Leadership
          Communication
      Achievements
        Promotion
        Certification
        Side Project
    Health
      Fitness
        Gym 3x/week
        Running
        Yoga
      Nutrition
        Meal Prep
        Water Intake
        Reduce Sugar
    Finance
      Income
        Salary Increase
        Freelancing
        Investments
      Savings
        Emergency Fund
        Retirement
        Travel Fund
    Personal
      Learning
        Books
        Courses
        Languages
      Hobbies
        Photography
        Music
        Cooking
```

### Software Architecture
```mermaid
mindmap
  root((E-Commerce Platform))
    Frontend
      Web App
        React
        Redux
        TypeScript
      Mobile
        React Native
        iOS
        Android
    Backend
      API Gateway
        Authentication
        Rate Limiting
        Routing
      Microservices
        User Service
        Product Service
        Order Service
        Payment Service
      Databases
        PostgreSQL
        MongoDB
        Redis Cache
    Infrastructure
      Cloud Provider
        AWS
        Load Balancers
        Auto Scaling
      Monitoring
        Logs
        Metrics
        Alerts
      Security
        SSL/TLS
        Firewall
        Encryption
```

### Content Strategy
```mermaid
mindmap
  root((Content Marketing))
    Platforms
      Blog
        SEO Articles
        Tutorials
        Case Studies
      Social Media
        Twitter
        LinkedIn
        YouTube
      Email
        Newsletter
        Drip Campaigns
        Announcements
    Content Types
      Educational
        How-to Guides
        Webinars
        Courses
      Promotional
        Product Updates
        Success Stories
        Comparisons
      Engagement
        Polls
        Q&A
        Community
    Metrics
      Traffic
        Page Views
        Unique Visitors
        Bounce Rate
      Engagement
        Time on Page
        Comments
        Shares
      Conversion
        Signups
        Downloads
        Sales
```

### Event Planning
```mermaid
mindmap
  root((Annual Conference))
    Venue
      Location
        City Selection
        Accessibility
        Hotels
      Facilities
        Main Hall
        Breakout Rooms
        Catering
    Program
      Speakers
        Keynotes
        Workshops
        Panels
      Schedule
        Day 1
        Day 2
        Networking
    Logistics
      Registration
        Tickets
        Check-in
        Badges
      Technology
        AV Equipment
        WiFi
        Event App
      Catering
        Breakfast
        Lunch
        Coffee Breaks
    Marketing
      Promotion
        Email Campaign
        Social Media
        Partners
      Sponsorship
        Tiers
        Benefits
        Activation
```

### Problem Solving
```mermaid
mindmap
  root((Performance Issue))
    Symptoms
      Slow Load Time
        Homepage: 8s
        API: 3s
        Database: 2s
      High CPU Usage
        Peak: 95%
        Average: 75%
    Analysis
      Database
        Slow Queries
        Missing Indexes
        N+1 Problems
      Application
        Memory Leaks
        Inefficient Code
        Large Payloads
      Infrastructure
        Server Resources
        Network Latency
        CDN Issues
    Solutions
      Quick Wins
        Add Indexes
        Enable Caching
        Optimize Images
      Long Term
        Refactor Code
        Scale Infrastructure
        Migrate Database
```

## Tips

- Use indentation to show hierarchy (2 or 4 spaces)
- Central node uses `root((text))` syntax
- Keep node text concise
- Use shapes to categorize types of information
- Limit depth to 3-4 levels for readability
- Organize clockwise or by importance
- Group related concepts together
- Use mindmaps for brainstorming and planning
