# Class Diagrams

Class diagrams show object-oriented structure with classes, attributes, methods, and relationships.

## Basic Syntax

```mermaid
classDiagram
    class Animal {
        +String name
        +int age
        +makeSound()
    }
```

## Visibility Modifiers

```mermaid
classDiagram
    class Example {
        +public attribute
        -private attribute
        #protected attribute
        ~package attribute
        +publicMethod()
        -privateMethod()
    }
```

## Data Types

```mermaid
classDiagram
    class User {
        +String username
        +String email
        +int age
        +boolean isActive
        +Date createdAt
        +login(String password) boolean
        +updateProfile(User data) void
    }
```

## Relationships

```mermaid
classDiagram
    ClassA --|> ClassB : Inheritance
    ClassC --* ClassD : Composition
    ClassE --o ClassF : Aggregation
    ClassG --> ClassH : Association
    ClassI -- ClassJ : Link
    ClassK ..> ClassL : Dependency
    ClassM ..|> ClassN : Realization
```

### Relationship Cardinality

```mermaid
classDiagram
    Customer "1" --> "*" Order : places
    Order "*" --> "1..*" Product : contains
    User "1" --o "0..1" Profile : has
```

## Abstract Classes & Interfaces

```mermaid
classDiagram
    class Animal {
        <<abstract>>
        +String name
        +makeSound()* void
    }
    
    class Flyable {
        <<interface>>
        +fly() void
    }
    
    class Bird {
        +int wingspan
    }
    
    Bird --|> Animal
    Bird ..|> Flyable
```

## Enumerations

```mermaid
classDiagram
    class Status {
        <<enumeration>>
        PENDING
        ACTIVE
        COMPLETED
        CANCELLED
    }
    
    class Order {
        +Status status
    }
    
    Order --> Status
```

## Generic Types

```mermaid
classDiagram
    class List~T~ {
        +add(T item)
        +get(int index) T
    }
    
    class StringList {
        +String items
    }
    
    StringList --|> List~String~
```

## Common Patterns

### Simple Inheritance
```mermaid
classDiagram
    class Vehicle {
        +String brand
        +int year
        +start()
        +stop()
    }
    
    class Car {
        +int doors
        +openTrunk()
    }
    
    class Motorcycle {
        +boolean hasSidecar
    }
    
    Car --|> Vehicle
    Motorcycle --|> Vehicle
```

### Repository Pattern
```mermaid
classDiagram
    class IRepository~T~ {
        <<interface>>
        +find(id) T
        +findAll() List~T~
        +save(T entity) void
        +delete(id) void
    }
    
    class UserRepository {
        -Database db
        +findByEmail(email) User
    }
    
    class User {
        +String id
        +String email
        +String name
    }
    
    UserRepository ..|> IRepository~User~
    UserRepository --> User
```

### MVC Pattern
```mermaid
classDiagram
    class Model {
        +Object data
        +getData()
        +setData(Object data)
        +notifyObservers()
    }
    
    class View {
        +render()
        +update()
    }
    
    class Controller {
        -Model model
        -View view
        +handleInput()
        +updateModel()
    }
    
    Controller --> Model
    Controller --> View
    View ..> Model : observes
```

### Service Layer Architecture
```mermaid
classDiagram
    class Controller {
        -UserService service
        +register(data) Response
        +login(credentials) Response
    }
    
    class UserService {
        -UserRepository repo
        -EmailService email
        +createUser(data) User
        +authenticate(credentials) Token
    }
    
    class UserRepository {
        -Database db
        +save(user) void
        +findByEmail(email) User
    }
    
    class EmailService {
        +sendWelcomeEmail(user) void
    }
    
    class User {
        +String id
        +String email
        +String password
    }
    
    Controller --> UserService
    UserService --> UserRepository
    UserService --> EmailService
    UserRepository --> User
```

### Factory Pattern
```mermaid
classDiagram
    class ShapeFactory {
        <<factory>>
        +createShape(type) Shape
    }
    
    class Shape {
        <<interface>>
        +draw() void
        +area() double
    }
    
    class Circle {
        +int radius
        +draw() void
        +area() double
    }
    
    class Rectangle {
        +int width
        +int height
        +draw() void
        +area() double
    }
    
    ShapeFactory ..> Shape : creates
    Circle ..|> Shape
    Rectangle ..|> Shape
```

## Annotations

```mermaid
classDiagram
    class User {
        <<Entity>>
        +String id
    }
    
    class UserDTO {
        <<DataTransferObject>>
        +String username
    }
```

## Notes

```mermaid
classDiagram
    class User {
        +String email
    }
    note for User "This is the main user entity.\nIt handles authentication."
```

## Styling

```mermaid
classDiagram
    class ImportantClass {
        +criticalMethod()
    }
    
    style ImportantClass fill:#f96,stroke:#333,stroke-width:4px
```

## Tips

- Use meaningful class names (PascalCase)
- Show only relevant attributes/methods
- Use interfaces for contracts
- Mark abstract classes and methods
- Show cardinality on associations
- Group related classes visually
- Use composition over inheritance when modeling "has-a" relationships
