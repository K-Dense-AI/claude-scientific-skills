# Gantt Charts

Gantt charts visualize project schedules and timelines.

## Basic Syntax

```mermaid
gantt
    title Project Schedule
    dateFormat YYYY-MM-DD
    
    section Planning
    Requirements    :a1, 2024-01-01, 7d
    Design         :a2, after a1, 5d
    
    section Development
    Backend        :b1, 2024-01-15, 14d
    Frontend       :b2, after b1, 10d
    
    section Testing
    QA Testing     :c1, after b2, 7d
```

## Date Formats

```mermaid
gantt
    dateFormat YYYY-MM-DD
    title Different Date Formats
    
    section Example
    Task 1 :2024-01-01, 5d
    Task 2 :2024-01-06, 3d
```

Common formats:
- `YYYY-MM-DD` - 2024-01-15
- `DD-MM-YYYY` - 15-01-2024
- `YYYY-MM-DD HH:mm` - 2024-01-15 14:30

## Task Dependencies

```mermaid
gantt
    dateFormat YYYY-MM-DD
    
    section Phase 1
    Task A :a1, 2024-01-01, 5d
    Task B :a2, after a1, 3d
    
    section Phase 2
    Task C :b1, after a2, 4d
    Task D :b2, after a1, 6d
```

## Task States

```mermaid
gantt
    dateFormat YYYY-MM-DD
    
    section Status Examples
    Completed Task  :done, t1, 2024-01-01, 3d
    Active Task     :active, t2, after t1, 2d
    Critical Task   :crit, t3, after t2, 3d
    Future Task     :t4, after t3, 2d
    
    section Combinations
    Critical Done   :crit, done, t5, 2024-01-01, 2d
    Critical Active :crit, active, t6, after t5, 3d
```

## Milestones

```mermaid
gantt
    dateFormat YYYY-MM-DD
    
    section Development
    Coding :2024-01-01, 10d
    
    section Milestones
    Alpha Release :milestone, 2024-01-11, 0d
    Beta Release  :milestone, 2024-01-25, 0d
    Final Release :milestone, crit, 2024-02-15, 0d
```

## Common Patterns

### Software Development Sprint
```mermaid
gantt
    title Sprint 1 - User Authentication
    dateFormat YYYY-MM-DD
    
    section Planning
    Sprint Planning     :done, plan1, 2024-01-02, 1d
    
    section Backend
    Database Schema     :done, db1, after plan1, 2d
    API Endpoints       :active, api1, after db1, 3d
    Authentication      :auth1, after api1, 2d
    
    section Frontend
    Login UI            :ui1, 2024-01-05, 3d
    Registration UI     :ui2, after ui1, 2d
    Integration         :int1, after auth1, 2d
    
    section Testing
    Unit Tests          :test1, after auth1, 2d
    Integration Tests   :test2, after int1, 2d
    
    section Review
    Code Review         :review1, after test2, 1d
    Sprint Demo         :milestone, demo1, after review1, 0d
```

### Product Launch
```mermaid
gantt
    title Product Launch Plan
    dateFormat YYYY-MM-DD
    
    section Research
    Market Analysis     :done, r1, 2024-01-01, 14d
    Competitor Study    :done, r2, 2024-01-01, 14d
    
    section Development
    MVP Development     :done, d1, after r1, 30d
    Beta Features       :active, d2, after d1, 21d
    Bug Fixes           :crit, d3, after d2, 7d
    
    section Marketing
    Brand Strategy      :m1, 2024-01-15, 20d
    Content Creation    :m2, after m1, 15d
    Ad Campaign         :crit, m3, after m2, 10d
    
    section Launch
    Soft Launch         :milestone, l1, after d3, 0d
    Public Launch       :milestone, crit, l2, after m3, 0d
```

### Construction Project
```mermaid
gantt
    title House Construction Timeline
    dateFormat YYYY-MM-DD
    
    section Planning
    Permits & Approvals :crit, p1, 2024-01-01, 30d
    
    section Foundation
    Site Preparation    :f1, after p1, 7d
    Foundation Pour     :crit, f2, after f1, 14d
    
    section Structure
    Framing            :s1, after f2, 21d
    Roofing            :s2, after s1, 10d
    Windows & Doors    :s3, after s2, 7d
    
    section Systems
    Electrical         :sys1, after s3, 14d
    Plumbing           :sys2, after s3, 14d
    HVAC               :sys3, after s3, 14d
    
    section Finishing
    Drywall            :fin1, after sys1, 10d
    Painting           :fin2, after fin1, 7d
    Flooring           :fin3, after fin2, 10d
    
    section Final
    Inspection         :milestone, crit, i1, after fin3, 0d
    Handover           :milestone, h1, after i1, 0d
```

### Marketing Campaign
```mermaid
gantt
    title Q1 Marketing Campaign
    dateFormat YYYY-MM-DD
    excludes weekends
    
    section Strategy
    Campaign Planning   :done, s1, 2024-01-02, 5d
    Budget Approval     :done, s2, after s1, 2d
    
    section Content
    Blog Posts          :active, c1, after s2, 15d
    Video Production    :c2, after s2, 20d
    Social Media        :c3, after s2, 30d
    
    section Ads
    Ad Design           :a1, after s2, 7d
    Ad Placement        :a2, after a1, 3d
    Campaign Launch     :milestone, crit, launch, after a2, 0d
    
    section Analysis
    A/B Testing         :t1, after launch, 14d
    Analytics Review    :t2, after t1, 5d
    Optimization        :t3, after t2, 10d
```

### Event Planning
```mermaid
gantt
    title Conference Organization
    dateFormat YYYY-MM-DD
    
    section Preparation
    Venue Booking       :done, v1, 2024-01-01, 7d
    Speaker Invites     :done, v2, 2024-01-01, 14d
    Sponsor Outreach    :active, v3, 2024-01-01, 30d
    
    section Marketing
    Website Launch      :m1, after v1, 10d
    Ticket Sales Start  :milestone, m2, after m1, 0d
    Social Campaign     :m3, after m2, 45d
    
    section Logistics
    Catering Contract   :l1, 2024-02-01, 7d
    AV Equipment        :l2, 2024-02-01, 7d
    Signage & Materials :l3, 2024-02-15, 10d
    
    section Final Week
    Venue Setup         :crit, f1, 2024-03-01, 2d
    Tech Rehearsal      :crit, f2, after f1, 1d
    Event Day           :milestone, crit, event, after f2, 0d
    
    section Post-Event
    Follow-up Emails    :p1, after event, 3d
    Feedback Analysis   :p2, after p1, 5d
```

## Excluding Days

```mermaid
gantt
    title Project with Holidays
    dateFormat YYYY-MM-DD
    excludes weekends 2024-01-15
    
    section Tasks
    Task 1 :2024-01-01, 10d
    Task 2 :after Task 1, 5d
```

## Axis Format

```mermaid
gantt
    title Custom Time Format
    dateFormat YYYY-MM-DD
    axisFormat %d/%m
    
    section Phase 1
    Development :2024-01-01, 30d
```

Format options:
- `%Y` - 4-digit year (2024)
- `%y` - 2-digit year (24)
- `%m` - Month number (01-12)
- `%b` - Month abbreviation (Jan)
- `%B` - Full month (January)
- `%d` - Day of month (01-31)
- `%a` - Day abbreviation (Mon)
- `%A` - Full day (Monday)
- `%H:%M` - Time (14:30)

## Today Marker

```mermaid
gantt
    title Project Progress
    dateFormat YYYY-MM-DD
    todayMarker stroke-width:3px,stroke:#f00,opacity:0.5
    
    section Tasks
    Completed :done, 2024-01-01, 5d
    In Progress :active, 2024-01-06, 5d
    Upcoming :2024-01-11, 5d
```

## Tips

- Use `excludes weekends` for business day planning
- Mark critical path with `crit` tag
- Use milestones (0d duration) for key dates
- Combine tags: `crit, done` or `crit, active`
- Use `after taskId` for dependencies
- Keep section names concise
- Use meaningful task IDs for dependencies
- Consider time zones for global teams
