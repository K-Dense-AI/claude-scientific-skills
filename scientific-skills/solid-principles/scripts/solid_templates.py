#!/usr/bin/env python3
"""
SOLID Templates Generator
Generates code templates following SOLID principles.
"""

import sys
from typing import Dict, List


class TemplateGenerator:
    """Generate SOLID-compliant code templates."""

    @staticmethod
    def protocol_template(name: str, methods: List[Dict[str, str]]) -> str:
        """Generate a Protocol (interface) template."""
        method_defs = []
        for method in methods:
            params = method.get('params', '')
            return_type = method.get('return_type', 'None')
            method_defs.append(f'    def {method["name"]}({params}) -> {return_type}: ...')

        methods_str = '\n'.join(method_defs)

        return f'''from typing import Protocol


class {name}(Protocol):
    """Protocol defining the interface for {name.lower().replace("protocol", "")}."""

{methods_str}
'''

    @staticmethod
    def dependency_injection_template(class_name: str, dependencies: List[str]) -> str:
        """Generate a class with dependency injection."""
        dep_params = ', '.join([f'{dep.lower()}: {dep}' for dep in dependencies])
        dep_assignments = '\n        '.join([f'self.{dep.lower()} = {dep.lower()}' for dep in dependencies])

        return f'''class {class_name}:
    """Business logic class with injected dependencies."""

    def __init__(self, {dep_params}):
        """{class_name} constructor with dependency injection.

        Args:
            {chr(10).join([f'{dep.lower()}: Implementation of {dep} protocol' for dep in dependencies])}
        """
        {dep_assignments}

    def execute(self):
        """Execute the main business logic."""
        # Use injected dependencies
        # Example: self.repository.save(data)
        raise NotImplementedError("Implement business logic here")
'''

    @staticmethod
    def strategy_pattern_template(context_name: str, strategy_name: str) -> str:
        """Generate Strategy pattern template."""
        return f'''from typing import Protocol


class {strategy_name}(Protocol):
    """Strategy interface for different implementations."""

    def execute(self, data: any) -> any:
        """Execute the strategy."""
        ...


class Concrete{strategy_name}A:
    """First concrete strategy implementation."""

    def execute(self, data: any) -> any:
        """Execute strategy A."""
        # Implementation A
        pass


class Concrete{strategy_name}B:
    """Second concrete strategy implementation."""

    def execute(self, data: any) -> any:
        """Execute strategy B."""
        # Implementation B
        pass


class {context_name}:
    """Context that uses a strategy."""

    def __init__(self, strategy: {strategy_name}):
        """Initialize with a strategy.

        Args:
            strategy: The strategy to use
        """
        self._strategy = strategy

    def set_strategy(self, strategy: {strategy_name}):
        """Change the strategy at runtime."""
        self._strategy = strategy

    def do_something(self, data: any) -> any:
        """Execute using the current strategy."""
        return self._strategy.execute(data)
'''

    @staticmethod
    def repository_pattern_template(entity_name: str) -> str:
        """Generate Repository pattern template."""
        return f'''from typing import Protocol, List, Optional


class {entity_name}:
    """Domain entity."""

    def __init__(self, id: str, **attributes):
        self.id = id
        self.__dict__.update(attributes)


class {entity_name}Repository(Protocol):
    """Repository protocol for {entity_name}."""

    def get_by_id(self, id: str) -> Optional[{entity_name}]:
        """Get entity by ID."""
        ...

    def get_all(self) -> List[{entity_name}]:
        """Get all entities."""
        ...

    def save(self, entity: {entity_name}) -> None:
        """Save entity."""
        ...

    def delete(self, id: str) -> None:
        """Delete entity by ID."""
        ...


class InMemory{entity_name}Repository:
    """In-memory implementation of {entity_name}Repository."""

    def __init__(self):
        self._storage: Dict[str, {entity_name}] = {{}}

    def get_by_id(self, id: str) -> Optional[{entity_name}]:
        return self._storage.get(id)

    def get_all(self) -> List[{entity_name}]:
        return list(self._storage.values())

    def save(self, entity: {entity_name}) -> None:
        self._storage[entity.id] = entity

    def delete(self, id: str) -> None:
        self._storage.pop(id, None)


class Database{entity_name}Repository:
    """Database implementation of {entity_name}Repository."""

    def __init__(self, connection):
        self._connection = connection

    def get_by_id(self, id: str) -> Optional[{entity_name}]:
        # Database query implementation
        raise NotImplementedError()

    def get_all(self) -> List[{entity_name}]:
        # Database query implementation
        raise NotImplementedError()

    def save(self, entity: {entity_name}) -> None:
        # Database save implementation
        raise NotImplementedError()

    def delete(self, id: str) -> None:
        # Database delete implementation
        raise NotImplementedError()
'''

    @staticmethod
    def factory_pattern_template(product_name: str) -> str:
        """Generate Factory pattern template."""
        return f'''from typing import Protocol
from enum import Enum


class {product_name}(Protocol):
    """Product interface."""

    def operation(self) -> str:
        """Perform the operation."""
        ...


class Concrete{product_name}A:
    """Concrete product A."""

    def operation(self) -> str:
        return "Result from Product A"


class Concrete{product_name}B:
    """Concrete product B."""

    def operation(self) -> str:
        return "Result from Product B"


class {product_name}Type(Enum):
    """Types of products."""
    TYPE_A = "a"
    TYPE_B = "b"


class {product_name}Factory:
    """Factory for creating {product_name} instances."""

    @staticmethod
    def create(product_type: {product_name}Type) -> {product_name}:
        """Create a product based on type.

        Args:
            product_type: The type of product to create

        Returns:
            A concrete product instance
        """
        if product_type == {product_name}Type.TYPE_A:
            return Concrete{product_name}A()
        elif product_type == {product_name}Type.TYPE_B:
            return Concrete{product_name}B()
        else:
            raise ValueError(f"Unknown product type: {{product_type}}")


# Usage example:
# factory = {product_name}Factory()
# product = factory.create({product_name}Type.TYPE_A)
# result = product.operation()
'''

    @staticmethod
    def layered_architecture_template() -> str:
        """Generate layered architecture project structure."""
        return '''# Layered Architecture Structure

```
project/
├── domain/                    # Business logic layer (no dependencies)
│   ├── __init__.py
│   ├── entities/             # Business entities
│   │   ├── __init__.py
│   │   └── user.py
│   └── protocols/            # Interfaces/contracts
│       ├── __init__.py
│       └── repository.py
│
├── application/               # Use cases layer (depends on domain)
│   ├── __init__.py
│   └── services/
│       ├── __init__.py
│       └── user_service.py
│
└── infrastructure/            # Implementation layer (depends on domain)
    ├── __init__.py
    ├── repositories/         # Data access implementations
    │   ├── __init__.py
    │   └── user_repository.py
    └── external/             # External service integrations
        ├── __init__.py
        └── email_service.py
```

## domain/entities/user.py
```python
from dataclasses import dataclass


@dataclass
class User:
    """User entity - pure business object."""
    id: str
    name: str
    email: str
```

## domain/protocols/repository.py
```python
from typing import Protocol, Optional, List
from domain.entities.user import User


class UserRepository(Protocol):
    """Repository protocol - defines the contract."""

    def get_by_id(self, user_id: str) -> Optional[User]: ...
    def get_all(self) -> List[User]: ...
    def save(self, user: User) -> None: ...
    def delete(self, user_id: str) -> None: ...
```

## application/services/user_service.py
```python
from domain.entities.user import User
from domain.protocols.repository import UserRepository


class UserService:
    """Application service - orchestrates business logic."""

    def __init__(self, repository: UserRepository):
        """Inject repository dependency."""
        self._repository = repository

    def register_user(self, name: str, email: str) -> User:
        """Business logic for user registration."""
        # Validation
        if not email or '@' not in email:
            raise ValueError("Invalid email")

        # Create entity
        user = User(id=self._generate_id(), name=name, email=email)

        # Persist via repository
        self._repository.save(user)

        return user

    def _generate_id(self) -> str:
        import uuid
        return str(uuid.uuid4())
```

## infrastructure/repositories/user_repository.py
```python
from typing import Optional, List, Dict
from domain.entities.user import User


class InMemoryUserRepository:
    """Concrete implementation of UserRepository."""

    def __init__(self):
        self._storage: Dict[str, User] = {}

    def get_by_id(self, user_id: str) -> Optional[User]:
        return self._storage.get(user_id)

    def get_all(self) -> List[User]:
        return list(self._storage.values())

    def save(self, user: User) -> None:
        self._storage[user.id] = user

    def delete(self, user_id: str) -> None:
        self._storage.pop(user_id, None)
```

## Dependency Flow

```
infrastructure → domain ← application
     ↓              ↑
     └──────────────┘
     (implements)  (uses)
```

**Key principles:**
- Domain layer has NO dependencies (pure business logic)
- Application layer depends ONLY on domain abstractions
- Infrastructure layer implements domain protocols
- Dependency injection connects the layers at runtime
'''


def main():
    """Main CLI for template generation."""
    generator = TemplateGenerator()

    templates = {
        '1': ('Protocol', lambda: generator.protocol_template(
            input("Protocol name: "),
            [{'name': input("Method name: "), 'params': input("Parameters: "), 'return_type': input("Return type: ")}]
        )),
        '2': ('Dependency Injection', lambda: generator.dependency_injection_template(
            input("Class name: "),
            input("Dependencies (comma-separated): ").split(',')
        )),
        '3': ('Strategy Pattern', lambda: generator.strategy_pattern_template(
            input("Context class name: "),
            input("Strategy name: ")
        )),
        '4': ('Repository Pattern', lambda: generator.repository_pattern_template(
            input("Entity name: ")
        )),
        '5': ('Factory Pattern', lambda: generator.factory_pattern_template(
            input("Product name: ")
        )),
        '6': ('Layered Architecture', lambda: generator.layered_architecture_template()),
    }

    print("\n=== SOLID Template Generator ===\n")
    print("Choose a template:")
    for key, (name, _) in templates.items():
        print(f"  {key}. {name}")
    print()

    choice = input("Enter choice (1-6): ").strip()

    if choice in templates:
        name, generator_func = templates[choice]
        print(f"\n--- {name} Template ---\n")
        template = generator_func()
        print(template)
    else:
        print("Invalid choice")
        sys.exit(1)


if __name__ == '__main__':
    main()
