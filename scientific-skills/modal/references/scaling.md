# Scaling Out on Modal

## Automatic Autoscaling

Every Modal Function corresponds to an autoscaling pool of containers. Modal's autoscaler:
- Spins up containers when no capacity available
- Spins down containers when resources idle
- Scales to zero by default when no inputs to process

Autoscaling decisions are made quickly and frequently.

## Parallel Execution with `.map()`

Run function repeatedly with different inputs in parallel:

```python
@app.function()
def evaluate_model(x):
    return x ** 2

@app.local_entrypoint()
def main():
    inputs = list(range(100))
    # Runs 100 inputs in parallel across containers
    for result in evaluate_model.map(inputs):
        print(result)
```

### Multiple Arguments with `.starmap()`

For functions with multiple arguments:

```python
@app.function()
def add(a, b):
    return a + b

@app.local_entrypoint()
def main():
    results = list(add.starmap([(1, 2), (3, 4)]))
    # [3, 7]
```

### Exception Handling

```python
@app.function()
def may_fail(a):
    if a == 2:
        raise Exception("error")
    return a ** 2

@app.local_entrypoint()
def main():
    results = list(may_fail.map(
        range(3),
        return_exceptions=True,
        wrap_returned_exceptions=False
    ))
    # [0, 1, Exception('error')]
```

## Autoscaling Configuration

Configure autoscaler behavior with parameters:

```python
@app.function(
    max_containers=100,      # Upper limit on containers
    min_containers=2,        # Keep warm even when inactive
    buffer_containers=5,     # Maintain buffer while active
    scaledown_window=60,     # Max idle time before scaling down (seconds)
)
def my_function():
    ...
```

Parameters:
- **max_containers**: Upper limit on total containers
- **min_containers**: Minimum kept warm even when inactive
- **buffer_containers**: Buffer size while function active (additional inputs won't need to queue)
- **scaledown_window**: Maximum idle duration before scale down (seconds)

Trade-offs:
- Larger warm pool/buffer → Higher cost, lower latency
- Longer scaledown window → Less churn for infrequent requests

## Dynamic Autoscaler Updates

Update autoscaler settings without redeployment:

```python
f = modal.Function.from_name("my-app", "f")
f.update_autoscaler(max_containers=100)
```

Settings revert to decorator configuration on next deploy, or are overridden by further updates:

```python
f.update_autoscaler(min_containers=2, max_containers=10)
f.update_autoscaler(min_containers=4)  # max_containers=10 still in effect
```

### Time-Based Scaling

Adjust warm pool based on time of day:

```python
@app.function()
def inference_server():
    ...

@app.function(schedule=modal.Cron("0 6 * * *", timezone="America/New_York"))
def increase_warm_pool():
    inference_server.update_autoscaler(min_containers=4)

@app.function(schedule=modal.Cron("0 22 * * *", timezone="America/New_York"))
def decrease_warm_pool():
    inference_server.update_autoscaler(min_containers=0)
```

### For Classes

Update autoscaler for specific parameter instances:

```python
MyClass = modal.Cls.from_name("my-app", "MyClass")
obj = MyClass(model_version="3.5")
obj.update_autoscaler(buffer_containers=2)  # type: ignore
```

## Input Concurrency

Process multiple inputs per container with `@modal.concurrent`:

```python
@app.function()
@modal.concurrent(max_inputs=100)
def my_function(input: str):
    # Container can handle up to 100 concurrent inputs
    ...
```

Ideal for I/O-bound workloads:
- Database queries
- External API requests
- Remote Modal Function calls

### Concurrency Mechanisms

**Synchronous Functions**: Separate threads (must be thread-safe)

```python
@app.function()
@modal.concurrent(max_inputs=10)
def sync_function():
    time.sleep(1)  # Must be thread-safe
```

**Async Functions**: Separate asyncio tasks (must not block event loop)

```python
@app.function()
@modal.concurrent(max_inputs=10)
async def async_function():
    await asyncio.sleep(1)  # Must not block event loop
```

### Target vs Max Inputs

```python
@app.function()
@modal.concurrent(
    max_inputs=120,    # Hard limit
    target_inputs=100  # Autoscaler target
)
def my_function(input: str):
    # Allow 20% burst above target
    ...
```

Autoscaler aims for `target_inputs`, but containers can burst to `max_inputs` during scale-up.

## Scaling Limits

Modal enforces limits per function:
- 2,000 pending inputs (not yet assigned to containers)
- 25,000 total inputs (running + pending)

For `.spawn()` async jobs: up to 1 million pending inputs.

Exceeding limits returns `Resource Exhausted` error - retry later.

Each `.map()` invocation: max 1,000 concurrent inputs.

## Async Usage

Use async APIs for arbitrary parallel execution patterns:

```python
@app.function()
async def async_task(x):
    await asyncio.sleep(1)
    return x * 2

@app.local_entrypoint()
async def main():
    tasks = [async_task.remote.aio(i) for i in range(100)]
    results = await asyncio.gather(*tasks)
```

## Common Gotchas

**Incorrect**: Using Python's builtin map (runs sequentially)
```python
# DON'T DO THIS
results = map(evaluate_model, inputs)
```

**Incorrect**: Calling function first
```python
# DON'T DO THIS
results = evaluate_model(inputs).map()
```

**Correct**: Call .map() on Modal function object
```python
# DO THIS
results = evaluate_model.map(inputs)
```
