# Modal Functions

## Basic Function Definition

Decorate Python functions with `@app.function()`:

```python
import modal

app = modal.App(name="my-app")

@app.function()
def my_function():
    print("Hello from Modal!")
    return "result"
```

## Calling Functions

### Remote Execution

Call `.remote()` to run on Modal:

```python
@app.local_entrypoint()
def main():
    result = my_function.remote()
    print(result)
```

### Local Execution

Call `.local()` to run locally (useful for testing):

```python
result = my_function.local()
```

## Function Parameters

Functions accept standard Python arguments:

```python
@app.function()
def process(x: int, y: str):
    return f"{y}: {x * 2}"

@app.local_entrypoint()
def main():
    result = process.remote(42, "answer")
```

## Deployment

### Ephemeral Apps

Run temporarily:
```bash
modal run script.py
```

### Deployed Apps

Deploy persistently:
```bash
modal deploy script.py
```

Access deployed functions from other code:

```python
f = modal.Function.from_name("my-app", "my_function")
result = f.remote(args)
```

## Entrypoints

### Local Entrypoint

Code that runs on local machine:

```python
@app.local_entrypoint()
def main():
    result = my_function.remote()
    print(result)
```

### Remote Entrypoint

Use `@app.function()` without local_entrypoint - runs entirely on Modal:

```python
@app.function()
def train_model():
    # All code runs in Modal
    ...
```

Invoke with:
```bash
modal run script.py::app.train_model
```

## Argument Parsing

Entrypoints with primitive type arguments get automatic CLI parsing:

```python
@app.local_entrypoint()
def main(foo: int, bar: str):
    some_function.remote(foo, bar)
```

Run with:
```bash
modal run script.py --foo 1 --bar "hello"
```

For custom parsing, accept variable-length arguments:

```python
import argparse

@app.function()
def train(*arglist):
    parser = argparse.ArgumentParser()
    parser.add_argument("--foo", type=int)
    args = parser.parse_args(args=arglist)
```

## Function Configuration

Common parameters:

```python
@app.function(
    image=my_image,           # Custom environment
    gpu="A100",               # GPU type
    cpu=2.0,                  # CPU cores
    memory=4096,              # Memory in MB
    timeout=3600,             # Timeout in seconds
    retries=3,                # Number of retries
    secrets=[my_secret],      # Environment secrets
    volumes={"/data": vol},   # Persistent storage
)
def my_function():
    ...
```

## Parallel Execution

### Map

Run function on multiple inputs in parallel:

```python
@app.function()
def evaluate_model(x):
    return x ** 2

@app.local_entrypoint()
def main():
    inputs = list(range(100))
    for result in evaluate_model.map(inputs):
        print(result)
```

### Starmap

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
results = my_func.map(
    range(3),
    return_exceptions=True,
    wrap_returned_exceptions=False
)
# [0, 1, Exception('error')]
```

## Async Functions

Define async functions:

```python
@app.function()
async def async_function(x: int):
    await asyncio.sleep(1)
    return x * 2

@app.local_entrypoint()
async def main():
    result = await async_function.remote.aio(42)
```

## Generator Functions

Return iterators for streaming results:

```python
@app.function()
def generate_data():
    for i in range(10):
        yield i

@app.local_entrypoint()
def main():
    for value in generate_data.remote_gen():
        print(value)
```

## Spawning Functions

Submit functions for background execution:

```python
@app.function()
def process_job(data):
    # Long-running job
    return result

@app.local_entrypoint()
def main():
    # Spawn without waiting
    call = process_job.spawn(data)

    # Get result later
    result = call.get(timeout=60)
```

## Programmatic Execution

Run apps programmatically:

```python
def main():
    with modal.enable_output():
        with app.run():
            result = some_function.remote()
```

## Specifying Entrypoint

With multiple functions, specify which to run:

```python
@app.function()
def f():
    print("Function f")

@app.function()
def g():
    print("Function g")
```

Run specific function:
```bash
modal run script.py::app.f
modal run script.py::app.g
```
