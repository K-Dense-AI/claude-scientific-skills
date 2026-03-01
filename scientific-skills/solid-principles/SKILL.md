---
name: solid-principles
description: Automatically apply SOLID design principles when writing or reviewing code. Activated for code design, refactoring, architecture discussions, design patterns, or class structure tasks.
---

# SOLID Principles Skill

## Purpose
Automatically apply SOLID principles to all code writing, review, and refactoring tasks to ensure maintainable, extensible, and robust software design.

## Instructions

### 1. Always Check for SOLID Violations

When writing or reviewing code, systematically check for these violations:

| Principle | Violation Signs | Quick Fix |
|-----------|----------------|-----------|
| **S**ingle Responsibility | • Class > 100 lines<br>• Multiple reasons to change<br>• Mixed concerns (data + UI + logic) | Extract Class pattern |
| **O**pen/Closed | • Long if/elif chains<br>• Switch statements on type<br>• Modifying existing code for new features | Strategy/Factory pattern |
| **L**iskov Substitution | • Empty method overrides<br>• Throwing NotImplementedError<br>• Subclass breaks parent contract | Redesign inheritance hierarchy |
| **I**nterface Segregation | • Fat interfaces with unused methods<br>• Clients forced to implement irrelevant methods | Split into focused interfaces |
| **D**ependency Inversion | • Direct instantiation in classes<br>• Concrete class dependencies<br>• No abstraction layer | Dependency Injection + Protocols |

### 2. Refactoring Patterns

When violations are detected, apply these patterns:

**Single Responsibility Violations:**
```python
# Before: God class
class UserManager:
    def save_to_db(self): ...
    def send_email(self): ...
    def validate(self): ...

# After: Separated concerns
class User: ...  # Data model
class UserRepository: ...  # Persistence
class EmailService: ...  # Notifications
class UserValidator: ...  # Validation
```

**Open/Closed Violations:**
```python
# Before: Modification required for new types
def process_payment(method):
    if method == "credit": ...
    elif method == "paypal": ...
    # Adding new method requires modifying this

# After: Extension via strategy pattern
class PaymentProcessor(Protocol):
    def process(self): ...

class CreditCardProcessor: ...
class PayPalProcessor: ...
```

**Dependency Inversion Violations:**
```python
# Before: Direct dependency
class OrderService:
    def __init__(self):
        self.db = MySQLDatabase()  # Concrete dependency

# After: Dependency injection with protocol
class Database(Protocol):
    def save(self, data): ...

class OrderService:
    def __init__(self, db: Database):
        self.db = db  # Abstract dependency
```

### 3. Project Structure

Organize code in layers that respect dependency inversion:

```
project/
├── domain/          # Business logic (no dependencies)
│   ├── entities/
│   └── protocols/
├── application/     # Use cases (depends on domain)
│   └── services/
└── infrastructure/  # Implementation details (depends on domain)
    ├── repositories/
    └── external/
```

### 4. Analysis Workflow

When asked to review or write code:

1. **Identify violations** using `scripts/solid_analyzer.py`
2. **Suggest refactoring** with specific patterns
3. **Generate templates** using `scripts/solid_templates.py`
4. **Validate** that refactored code follows all 5 principles

### 5. Code Generation Rules

When writing new code:

- Start with protocols/interfaces (abstractions first)
- Keep classes under 100 lines
- Use composition over inheritance
- Inject dependencies via `__init__`
- Create factory functions for complex object creation
- Avoid concrete dependencies in business logic

## Examples of When to Activate

- "Design a payment processing system"
- "Refactor this class"
- "Review this code for design issues"
- "How should I structure this feature?"
- "This class is getting too large"
- "Add a new notification type"

## Tools Available

- `scripts/solid_analyzer.py` - Analyze code for violations
- `scripts/solid_templates.py` - Generate refactoring templates
- `reference.md` - Detailed examples and patterns

## Expected Outcomes

After applying this skill, code should:
- Have clear single responsibilities
- Be extensible without modification
- Use proper abstraction layers
- Have focused, cohesive interfaces
- Depend on abstractions, not concretions
