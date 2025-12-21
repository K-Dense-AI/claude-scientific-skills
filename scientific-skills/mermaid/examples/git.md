# Git Graphs

Git graphs visualize branching strategies and commit histories.

## Basic Syntax

```mermaid
gitGraph
    commit
    commit
    branch develop
    checkout develop
    commit
    commit
    checkout main
    merge develop
```

## Commit with Messages

```mermaid
gitGraph
    commit id: "Initial commit"
    commit id: "Add feature"
    branch feature
    checkout feature
    commit id: "Work in progress"
    commit id: "Complete feature"
    checkout main
    merge feature
    commit id: "Release v1.0"
```

## Tags

```mermaid
gitGraph
    commit
    commit tag: "v1.0"
    branch develop
    checkout develop
    commit
    checkout main
    merge develop
    commit tag: "v1.1"
```

## Common Patterns

### Git Flow
```mermaid
gitGraph
    commit id: "Initial commit"
    branch develop
    checkout develop
    commit id: "Setup project"
    
    branch feature/user-auth
    checkout feature/user-auth
    commit id: "Add login"
    commit id: "Add registration"
    checkout develop
    merge feature/user-auth
    
    branch feature/dashboard
    checkout feature/dashboard
    commit id: "Create dashboard"
    commit id: "Add widgets"
    
    checkout develop
    commit id: "Fix bug"
    merge feature/dashboard
    
    checkout main
    merge develop tag: "v1.0.0"
    
    checkout develop
    branch hotfix/security
    checkout hotfix/security
    commit id: "Fix security issue"
    checkout main
    merge hotfix/security tag: "v1.0.1"
    checkout develop
    merge hotfix/security
```

### Feature Branch Workflow
```mermaid
gitGraph
    commit id: "Initial"
    commit id: "Setup"
    
    branch feature-a
    checkout feature-a
    commit id: "Feature A: Start"
    commit id: "Feature A: WIP"
    
    checkout main
    branch feature-b
    checkout feature-b
    commit id: "Feature B: Start"
    
    checkout feature-a
    commit id: "Feature A: Done"
    
    checkout main
    merge feature-a
    
    checkout feature-b
    commit id: "Feature B: Done"
    checkout main
    merge feature-b
```

### Release Process
```mermaid
gitGraph
    commit id: "0.1.0" tag: "v0.1.0"
    branch develop
    checkout develop
    commit id: "Add features"
    commit id: "More features"
    
    branch release/1.0
    checkout release/1.0
    commit id: "Prepare release"
    commit id: "Update docs"
    commit id: "Bug fixes"
    
    checkout main
    merge release/1.0 tag: "v1.0.0"
    
    checkout develop
    merge release/1.0
    commit id: "Continue development"
```

### Hotfix Workflow
```mermaid
gitGraph
    commit id: "1.0.0" tag: "v1.0.0"
    branch develop
    checkout develop
    commit id: "New features"
    commit id: "More work"
    
    checkout main
    branch hotfix/critical-bug
    checkout hotfix/critical-bug
    commit id: "Fix critical bug"
    commit id: "Add tests"
    
    checkout main
    merge hotfix/critical-bug tag: "v1.0.1"
    
    checkout develop
    merge hotfix/critical-bug
    commit id: "Continue dev"
```

### Multiple Feature Development
```mermaid
gitGraph
    commit id: "Start"
    branch develop
    checkout develop
    commit id: "Base setup"
    
    branch feature/auth
    checkout feature/auth
    commit id: "Auth: OAuth"
    commit id: "Auth: JWT"
    
    checkout develop
    branch feature/api
    checkout feature/api
    commit id: "API: Endpoints"
    
    checkout develop
    branch feature/ui
    checkout feature/ui
    commit id: "UI: Components"
    
    checkout feature/auth
    commit id: "Auth: Complete"
    checkout develop
    merge feature/auth
    
    checkout feature/api
    commit id: "API: Complete"
    checkout develop
    merge feature/api
    
    checkout feature/ui
    commit id: "UI: Complete"
    checkout develop
    merge feature/ui
    
    checkout main
    merge develop tag: "v2.0.0"
```

### Trunk-Based Development
```mermaid
gitGraph
    commit id: "Init"
    commit id: "Feature 1"
    
    branch short-lived-1
    checkout short-lived-1
    commit id: "Quick fix"
    checkout main
    merge short-lived-1
    
    commit id: "Feature 2"
    
    branch short-lived-2
    checkout short-lived-2
    commit id: "Small feature"
    checkout main
    merge short-lived-2
    
    commit id: "Release" tag: "v1.0"
```

### Parallel Release Branches
```mermaid
gitGraph
    commit id: "1.0.0" tag: "v1.0.0"
    branch release-1.x
    checkout release-1.x
    commit id: "1.1 features"
    
    checkout main
    commit id: "2.0 breaking changes"
    branch release-2.x
    checkout release-2.x
    commit id: "2.1 features"
    
    checkout release-1.x
    commit id: "1.1.0" tag: "v1.1.0"
    
    checkout release-2.x
    commit id: "2.1.0" tag: "v2.1.0"
    
    checkout main
    merge release-2.x
```

### Revert and Fix
```mermaid
gitGraph
    commit id: "Feature A"
    commit id: "Feature B"
    commit id: "Bad commit"
    commit id: "Revert bad commit"
    
    branch fix
    checkout fix
    commit id: "Proper fix"
    commit id: "Add tests"
    
    checkout main
    merge fix
```

## Branch Ordering

```mermaid
gitGraph
    commit
    branch feature-1 order: 3
    branch feature-2 order: 2
    branch feature-3 order: 1
    
    checkout feature-1
    commit
    checkout feature-2
    commit
    checkout feature-3
    commit
```

## Cherry-Pick Pattern

```mermaid
gitGraph
    commit id: "A"
    branch develop
    checkout develop
    commit id: "B"
    commit id: "C - Important fix"
    commit id: "D"
    
    checkout main
    commit id: "Cherry-pick C"
    commit tag: "hotfix"
```

## Tips

- Use meaningful commit messages
- Tag releases with version numbers
- Show merge direction clearly
- Keep graph readable (limit branches shown)
- Use consistent branch naming
- Order branches logically with `order:` parameter
- Illustrate specific workflows, not entire history
- Consider chronological flow (top to bottom)
