# Web Endpoints

## Quick Start

Create web endpoint with single decorator:

```python
image = modal.Image.debian_slim().pip_install("fastapi[standard]")

@app.function(image=image)
@modal.fastapi_endpoint()
def hello():
    return "Hello world!"
```

## Development and Deployment

### Development with `modal serve`

```bash
modal serve server.py
```

Creates ephemeral app with live-reloading. Changes to endpoints appear almost immediately.

### Deployment with `modal deploy`

```bash
modal deploy server.py
```

Creates persistent endpoint with stable URL.

## Simple Endpoints

### Query Parameters

```python
@app.function(image=image)
@modal.fastapi_endpoint()
def square(x: int):
    return {"square": x**2}
```

Call with:
```bash
curl "https://workspace--app-square.modal.run?x=42"
```

### POST Requests

```python
@app.function(image=image)
@modal.fastapi_endpoint(method="POST")
def square(item: dict):
    return {"square": item['x']**2}
```

Call with:
```bash
curl -X POST -H 'Content-Type: application/json' \
  --data '{"x": 42}' \
  https://workspace--app-square.modal.run
```

### Pydantic Models

```python
from pydantic import BaseModel

class Item(BaseModel):
    name: str
    qty: int = 42

@app.function()
@modal.fastapi_endpoint(method="POST")
def process(item: Item):
    return {"processed": item.name, "quantity": item.qty}
```

## ASGI Apps (FastAPI, Starlette, FastHTML)

Serve full ASGI applications:

```python
image = modal.Image.debian_slim().pip_install("fastapi[standard]")

@app.function(image=image)
@modal.concurrent(max_inputs=100)
@modal.asgi_app()
def fastapi_app():
    from fastapi import FastAPI

    web_app = FastAPI()

    @web_app.get("/")
    async def root():
        return {"message": "Hello"}

    @web_app.post("/echo")
    async def echo(request: Request):
        body = await request.json()
        return body

    return web_app
```

## WSGI Apps (Flask, Django)

Serve synchronous web frameworks:

```python
image = modal.Image.debian_slim().pip_install("flask")

@app.function(image=image)
@modal.concurrent(max_inputs=100)
@modal.wsgi_app()
def flask_app():
    from flask import Flask, request

    web_app = Flask(__name__)

    @web_app.post("/echo")
    def echo():
        return request.json

    return web_app
```

## Non-ASGI Web Servers

For frameworks with custom network binding:

```python
@app.function()
@modal.concurrent(max_inputs=100)
@modal.web_server(8000)
def my_server():
    import subprocess
    # Must bind to 0.0.0.0, not 127.0.0.1
    subprocess.Popen("python -m http.server -d / 8000", shell=True)
```

## Streaming Responses

Use FastAPI's `StreamingResponse`:

```python
import time

def event_generator():
    for i in range(10):
        yield f"data: event {i}\n\n".encode()
        time.sleep(0.5)

@app.function(image=modal.Image.debian_slim().pip_install("fastapi[standard]"))
@modal.fastapi_endpoint()
def stream():
    from fastapi.responses import StreamingResponse
    return StreamingResponse(
        event_generator(),
        media_type="text/event-stream"
    )
```

### Streaming from Modal Functions

```python
@app.function(gpu="any")
def process_gpu():
    for i in range(10):
        yield f"data: result {i}\n\n".encode()
        time.sleep(1)

@app.function(image=modal.Image.debian_slim().pip_install("fastapi[standard]"))
@modal.fastapi_endpoint()
def hook():
    from fastapi.responses import StreamingResponse
    return StreamingResponse(
        process_gpu.remote_gen(),
        media_type="text/event-stream"
    )
```

### With .map()

```python
@app.function()
def process_segment(i):
    return f"segment {i}\n"

@app.function(image=modal.Image.debian_slim().pip_install("fastapi[standard]"))
@modal.fastapi_endpoint()
def stream_parallel():
    from fastapi.responses import StreamingResponse
    return StreamingResponse(
        process_segment.map(range(10)),
        media_type="text/plain"
    )
```

## WebSockets

Supported with `@web_server`, `@asgi_app`, and `@wsgi_app`. Maintains single function call per connection. Use with `@modal.concurrent` for multiple simultaneous connections.

Full WebSocket protocol (RFC 6455) supported. Messages up to 2 MiB each.

## Authentication

### Proxy Auth Tokens

First-class authentication via Modal:

```python
@app.function()
@modal.fastapi_endpoint()
def protected():
    return "authenticated!"
```

Protect with tokens in settings, pass in headers:
- `Modal-Key`
- `Modal-Secret`

### Bearer Token Authentication

```python
from fastapi import Depends, HTTPException, status
from fastapi.security import HTTPBearer, HTTPAuthorizationCredentials

auth_scheme = HTTPBearer()

@app.function(secrets=[modal.Secret.from_name("auth-token")])
@modal.fastapi_endpoint()
async def protected(token: HTTPAuthorizationCredentials = Depends(auth_scheme)):
    import os
    if token.credentials != os.environ["AUTH_TOKEN"]:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid token"
        )
    return "success!"
```

### Client IP Address

```python
from fastapi import Request

@app.function()
@modal.fastapi_endpoint()
def get_ip(request: Request):
    return f"Your IP: {request.client.host}"
```

## Web Endpoint URLs

### Auto-Generated URLs

Format: `https://<workspace>--<app>-<function>.modal.run`

With environment suffix: `https://<workspace>-<suffix>--<app>-<function>.modal.run`

### Custom Labels

```python
@app.function()
@modal.fastapi_endpoint(label="api")
def handler():
    ...
# URL: https://workspace--api.modal.run
```

### Programmatic URL Retrieval

```python
@app.function()
@modal.fastapi_endpoint()
def my_endpoint():
    url = my_endpoint.get_web_url()
    return {"url": url}

# From deployed function
f = modal.Function.from_name("app-name", "my_endpoint")
url = f.get_web_url()
```

### Custom Domains

Available on Team and Enterprise plans:

```python
@app.function()
@modal.fastapi_endpoint(custom_domains=["api.example.com"])
def hello(message: str):
    return {"message": f"hello {message}"}
```

Multiple domains:
```python
@modal.fastapi_endpoint(custom_domains=["api.example.com", "api.example.net"])
```

Wildcard domains:
```python
@modal.fastapi_endpoint(custom_domains=["*.example.com"])
```

TLS certificates automatically generated and renewed.

## Performance

### Cold Starts

First request may experience cold start (few seconds). Modal keeps containers alive for subsequent requests.

### Scaling

- Autoscaling based on traffic
- Use `@modal.concurrent` for multiple requests per container
- Beyond concurrency limit, additional containers spin up
- Requests queue when at max containers

### Rate Limits

Default: 200 requests/second with 5-second burst multiplier
- Excess returns 429 status code
- Contact support to increase limits

### Size Limits

- Request body: up to 4 GiB
- Response body: unlimited
- WebSocket messages: up to 2 MiB
