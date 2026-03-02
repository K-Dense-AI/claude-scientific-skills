# SOLID Principles - Comprehensive Reference

## Table of Contents
1. [Single Responsibility Principle (SRP)](#single-responsibility-principle)
2. [Open/Closed Principle (OCP)](#openclosed-principle)
3. [Liskov Substitution Principle (LSP)](#liskov-substitution-principle)
4. [Interface Segregation Principle (ISP)](#interface-segregation-principle)
5. [Dependency Inversion Principle (DIP)](#dependency-inversion-principle)
6. [Real-World Examples](#real-world-examples)

---

## Single Responsibility Principle

**Definition:** A class should have one, and only one, reason to change.

### ❌ Bad Example

```python
class UserManager:
    """God class that does everything - violates SRP."""

    def __init__(self):
        self.db_connection = self._connect_to_db()

    def register_user(self, name, email, password):
        # Validation
        if not self._validate_email(email):
            raise ValueError("Invalid email")

        # Hash password
        hashed = self._hash_password(password)

        # Save to database
        self.db_connection.execute(
            "INSERT INTO users VALUES (?, ?, ?)",
            (name, email, hashed)
        )

        # Send welcome email
        self._send_email(email, "Welcome!", "Thanks for joining!")

        # Log activity
        self._log("User registered: " + email)

    def _validate_email(self, email): ...
    def _hash_password(self, password): ...
    def _send_email(self, to, subject, body): ...
    def _log(self, message): ...
    def _connect_to_db(self): ...
```

**Problems:**
- Changes to email validation require changing UserManager
- Changes to database schema require changing UserManager
- Changes to email service require changing UserManager
- Changes to logging require changing UserManager

### ✅ Good Example

```python
from dataclasses import dataclass


@dataclass
class User:
    """Entity - just data."""
    name: str
    email: str
    password_hash: str


class EmailValidator:
    """Single responsibility: email validation."""

    @staticmethod
    def validate(email: str) -> bool:
        return '@' in email and '.' in email


class PasswordHasher:
    """Single responsibility: password hashing."""

    @staticmethod
    def hash(password: str) -> str:
        import hashlib
        return hashlib.sha256(password.encode()).hexdigest()


class UserRepository:
    """Single responsibility: user persistence."""

    def __init__(self, db_connection):
        self._db = db_connection

    def save(self, user: User) -> None:
        self._db.execute(
            "INSERT INTO users VALUES (?, ?, ?)",
            (user.name, user.email, user.password_hash)
        )


class EmailService:
    """Single responsibility: sending emails."""

    def send_welcome_email(self, email: str) -> None:
        self._send(email, "Welcome!", "Thanks for joining!")

    def _send(self, to, subject, body):
        # Email sending logic
        pass


class ActivityLogger:
    """Single responsibility: logging."""

    def log_user_registration(self, email: str) -> None:
        self.log(f"User registered: {email}")

    def log(self, message: str) -> None:
        print(f"[LOG] {message}")


class UserRegistrationService:
    """Orchestrates user registration using focused components."""

    def __init__(
        self,
        repository: UserRepository,
        email_service: EmailService,
        logger: ActivityLogger
    ):
        self.repository = repository
        self.email_service = email_service
        self.logger = logger
        self.validator = EmailValidator()
        self.hasher = PasswordHasher()

    def register(self, name: str, email: str, password: str) -> User:
        # Validate
        if not self.validator.validate(email):
            raise ValueError("Invalid email")

        # Create user
        user = User(
            name=name,
            email=email,
            password_hash=self.hasher.hash(password)
        )

        # Save
        self.repository.save(user)

        # Send email
        self.email_service.send_welcome_email(email)

        # Log
        self.logger.log_user_registration(email)

        return user
```

**Benefits:**
- Each class has ONE reason to change
- Easy to test in isolation
- Easy to reuse components
- Clear responsibilities

---

## Open/Closed Principle

**Definition:** Software entities should be open for extension but closed for modification.

### ❌ Bad Example

```python
class PaymentProcessor:
    """Closed for extension - must modify to add new payment methods."""

    def process_payment(self, amount: float, method: str):
        if method == "credit_card":
            # Process credit card
            print(f"Processing ${amount} via credit card")
        elif method == "paypal":
            # Process PayPal
            print(f"Processing ${amount} via PayPal")
        elif method == "bitcoin":  # New method requires modifying this class!
            # Process Bitcoin
            print(f"Processing ${amount} via Bitcoin")
        else:
            raise ValueError(f"Unknown payment method: {method}")
```

**Problem:** Every new payment method requires modifying this class.

### ✅ Good Example

```python
from typing import Protocol


class PaymentMethod(Protocol):
    """Protocol defining payment method interface."""

    def process(self, amount: float) -> None:
        """Process a payment."""
        ...


class CreditCardPayment:
    """Concrete payment method."""

    def process(self, amount: float) -> None:
        print(f"Processing ${amount} via credit card")


class PayPalPayment:
    """Concrete payment method."""

    def process(self, amount: float) -> None:
        print(f"Processing ${amount} via PayPal")


class BitcoinPayment:
    """New payment method - no modification to existing code!"""

    def process(self, amount: float) -> None:
        print(f"Processing ${amount} via Bitcoin")


class PaymentProcessor:
    """Open for extension via new PaymentMethod implementations."""

    def __init__(self, payment_method: PaymentMethod):
        self._payment_method = payment_method

    def process_payment(self, amount: float) -> None:
        self._payment_method.process(amount)


# Usage:
processor = PaymentProcessor(CreditCardPayment())
processor.process_payment(100.0)

# Add new method without changing PaymentProcessor
processor = PaymentProcessor(BitcoinPayment())
processor.process_payment(50.0)
```

**Benefits:**
- Add new payment methods without modifying existing code
- Existing code is stable and tested
- Easy to add features via extension

---

## Liskov Substitution Principle

**Definition:** Subclasses should be substitutable for their base classes without altering program correctness.

### ❌ Bad Example

```python
class Bird:
    """Base class."""

    def fly(self):
        print("Flying...")


class Sparrow(Bird):
    """Sparrow can fly - OK."""
    pass


class Penguin(Bird):
    """Penguin cannot fly - violates LSP!"""

    def fly(self):
        raise NotImplementedError("Penguins can't fly!")


# Code that expects any Bird to fly:
def make_bird_fly(bird: Bird):
    bird.fly()  # Crashes if bird is a Penguin!


# This breaks:
penguin = Penguin()
make_bird_fly(penguin)  # NotImplementedError!
```

**Problem:** Penguin breaks the contract that all Birds can fly.

### ✅ Good Example

```python
from typing import Protocol


class Bird:
    """Base for all birds."""
    pass


class FlyingBird(Protocol):
    """Protocol for birds that can fly."""

    def fly(self) -> None:
        """Fly."""
        ...


class Sparrow(Bird):
    """Sparrow is a bird that can fly."""

    def fly(self) -> None:
        print("Sparrow flying...")


class Penguin(Bird):
    """Penguin is a bird but doesn't implement FlyingBird."""

    def swim(self) -> None:
        print("Penguin swimming...")


# Code that needs flying behavior:
def make_bird_fly(bird: FlyingBird):
    bird.fly()  # Only accepts birds that can actually fly


# This works:
sparrow = Sparrow()
make_bird_fly(sparrow)

# This won't compile (type error):
# penguin = Penguin()
# make_bird_fly(penguin)  # Type error - Penguin doesn't have fly()
```

**Benefits:**
- Subclasses don't break parent contracts
- Type system catches violations
- No surprising runtime errors

---

## Interface Segregation Principle

**Definition:** Clients should not be forced to depend on interfaces they don't use.

### ❌ Bad Example

```python
from typing import Protocol


class Worker(Protocol):
    """Fat interface - forces all workers to implement everything."""

    def work(self) -> None: ...
    def eat(self) -> None: ...
    def sleep(self) -> None: ...


class HumanWorker:
    """Human needs all three."""

    def work(self): print("Working...")
    def eat(self): print("Eating...")
    def sleep(self): print("Sleeping...")


class RobotWorker:
    """Robot doesn't eat or sleep - forced to implement useless methods!"""

    def work(self): print("Working...")

    def eat(self):
        # Robot doesn't eat!
        pass  # Violates ISP - forced empty implementation

    def sleep(self):
        # Robot doesn't sleep!
        pass  # Violates ISP - forced empty implementation
```

**Problem:** RobotWorker is forced to implement methods it doesn't need.

### ✅ Good Example

```python
from typing import Protocol


class Workable(Protocol):
    """Focused interface for work capability."""

    def work(self) -> None: ...


class Eatable(Protocol):
    """Focused interface for eating capability."""

    def eat(self) -> None: ...


class Sleepable(Protocol):
    """Focused interface for sleeping capability."""

    def sleep(self) -> None: ...


class HumanWorker:
    """Implements all three interfaces."""

    def work(self): print("Working...")
    def eat(self): print("Eating...")
    def sleep(self): print("Sleeping...")


class RobotWorker:
    """Only implements Workable - clean!"""

    def work(self): print("Working...")


# Code that only needs work:
def manage_work(worker: Workable):
    worker.work()  # Only requires work() method


# Both work:
manage_work(HumanWorker())
manage_work(RobotWorker())
```

**Benefits:**
- Clients only depend on methods they use
- No forced empty implementations
- More flexible composition

---

## Dependency Inversion Principle

**Definition:** High-level modules should not depend on low-level modules. Both should depend on abstractions.

### ❌ Bad Example

```python
class MySQLDatabase:
    """Concrete low-level database implementation."""

    def save(self, data):
        print(f"Saving to MySQL: {data}")


class UserService:
    """High-level business logic depends on concrete MySQL!"""

    def __init__(self):
        # Direct dependency on concrete implementation
        self.db = MySQLDatabase()  # Violates DIP!

    def register_user(self, name):
        user_data = {"name": name}
        self.db.save(user_data)  # Tightly coupled to MySQL


# Problems:
# - Cannot test UserService without MySQL
# - Cannot switch to PostgreSQL without changing UserService
# - Cannot use in-memory DB for tests
```

**Problem:** High-level UserService directly depends on low-level MySQLDatabase.

### ✅ Good Example

```python
from typing import Protocol


class Database(Protocol):
    """Abstraction that both high and low-level modules depend on."""

    def save(self, data: dict) -> None:
        """Save data."""
        ...


class MySQLDatabase:
    """Low-level module implements the abstraction."""

    def save(self, data: dict) -> None:
        print(f"Saving to MySQL: {data}")


class PostgreSQLDatabase:
    """Another low-level implementation."""

    def save(self, data: dict) -> None:
        print(f"Saving to PostgreSQL: {data}")


class InMemoryDatabase:
    """In-memory implementation for testing."""

    def __init__(self):
        self.storage = []

    def save(self, data: dict) -> None:
        self.storage.append(data)


class UserService:
    """High-level module depends on abstraction."""

    def __init__(self, db: Database):
        # Depends on abstraction, not concrete implementation
        self.db = db

    def register_user(self, name: str):
        user_data = {"name": name}
        self.db.save(user_data)


# Usage - inject any implementation:
service = UserService(MySQLDatabase())
service.register_user("Alice")

# Easy to switch:
service = UserService(PostgreSQLDatabase())
service.register_user("Bob")

# Easy to test:
service = UserService(InMemoryDatabase())
service.register_user("Test User")
```

**Benefits:**
- Easy to test (inject mock/fake)
- Easy to switch implementations
- High-level logic independent of low-level details
- Both depend on stable abstractions

---

## Real-World Examples

### Example 1: E-Commerce Order System

```python
from typing import Protocol, List
from dataclasses import dataclass
from decimal import Decimal


# Domain Layer (no dependencies)
@dataclass
class OrderItem:
    product_id: str
    quantity: int
    price: Decimal


@dataclass
class Order:
    order_id: str
    customer_id: str
    items: List[OrderItem]

    @property
    def total(self) -> Decimal:
        return sum(item.price * item.quantity for item in self.items)


# Protocols (abstractions)
class OrderRepository(Protocol):
    def save(self, order: Order) -> None: ...
    def get_by_id(self, order_id: str) -> Order: ...


class PaymentGateway(Protocol):
    def charge(self, amount: Decimal, customer_id: str) -> bool: ...


class NotificationService(Protocol):
    def send_order_confirmation(self, order: Order) -> None: ...


# Application Layer (business logic)
class OrderService:
    """Orchestrates order processing - depends only on abstractions."""

    def __init__(
        self,
        repository: OrderRepository,
        payment: PaymentGateway,
        notification: NotificationService
    ):
        self._repository = repository
        self._payment = payment
        self._notification = notification

    def place_order(self, order: Order) -> bool:
        """Place an order (use case)."""
        # Charge payment
        if not self._payment.charge(order.total, order.customer_id):
            return False

        # Save order
        self._repository.save(order)

        # Send confirmation
        self._notification.send_order_confirmation(order)

        return True


# Infrastructure Layer (implementations)
class SQLOrderRepository:
    """Concrete repository implementation."""

    def __init__(self, db_connection):
        self._db = db_connection

    def save(self, order: Order) -> None:
        # Save to SQL database
        pass

    def get_by_id(self, order_id: str) -> Order:
        # Query from SQL database
        pass


class StripePaymentGateway:
    """Concrete payment implementation."""

    def charge(self, amount: Decimal, customer_id: str) -> bool:
        # Charge via Stripe API
        return True


class EmailNotificationService:
    """Concrete notification implementation."""

    def send_order_confirmation(self, order: Order) -> None:
        # Send email
        print(f"Sending confirmation for order {order.order_id}")


# Dependency Injection (wiring)
def create_order_service(db_connection) -> OrderService:
    """Factory function to wire dependencies."""
    repository = SQLOrderRepository(db_connection)
    payment = StripePaymentGateway()
    notification = EmailNotificationService()

    return OrderService(repository, payment, notification)
```

### Example 2: Data Processing Pipeline

```python
from typing import Protocol, Any, List


# Strategy pattern for Open/Closed
class DataProcessor(Protocol):
    """Strategy for processing data."""

    def process(self, data: Any) -> Any: ...


class CSVProcessor:
    """Process CSV data."""

    def process(self, data: str) -> List[dict]:
        # Parse CSV
        return [{"row": 1}, {"row": 2}]


class JSONProcessor:
    """Process JSON data."""

    def process(self, data: str) -> dict:
        import json
        return json.loads(data)


class XMLProcessor:
    """Process XML data."""

    def process(self, data: str) -> dict:
        # Parse XML
        return {"xml": "data"}


class DataPipeline:
    """Open for extension via new processors."""

    def __init__(self, processor: DataProcessor):
        self._processor = processor

    def run(self, raw_data: str) -> Any:
        # Pre-processing
        cleaned = self._clean(raw_data)

        # Process (delegated to strategy)
        processed = self._processor.process(cleaned)

        # Post-processing
        return self._validate(processed)

    def _clean(self, data: str) -> str:
        return data.strip()

    def _validate(self, data: Any) -> Any:
        # Validation logic
        return data


# Usage - easy to extend with new formats
pipeline = DataPipeline(CSVProcessor())
result = pipeline.run("data,data,data")

pipeline = DataPipeline(JSONProcessor())
result = pipeline.run('{"key": "value"}')

# Add XML support without changing DataPipeline
pipeline = DataPipeline(XMLProcessor())
result = pipeline.run('<xml>data</xml>')
```

---

## Summary

| Principle | Violation Sign | Fix Pattern |
|-----------|---------------|-------------|
| **S** - Single Responsibility | 100+ line class, multiple concerns | Extract Class |
| **O** - Open/Closed | if/elif chains, modifying for new features | Strategy/Factory |
| **L** - Liskov Substitution | Empty overrides, NotImplementedError | Redesign hierarchy |
| **I** - Interface Segregation | Fat interfaces, unused methods | Split interfaces |
| **D** - Dependency Inversion | Direct instantiation, concrete deps | Dependency Injection + Protocols |

**Golden Rule:** Always code to abstractions (Protocols), inject dependencies, and keep classes focused.
