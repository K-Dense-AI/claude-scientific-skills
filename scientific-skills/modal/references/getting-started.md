# Getting Started with Modal

## Sign Up

Sign up for free at https://modal.com and get $30/month of credits.

## Authentication

Set up authentication using the Modal CLI:

```bash
modal token new
```

This creates credentials in `~/.modal.toml`. Alternatively, set environment variables:
- `MODAL_TOKEN_ID`
- `MODAL_TOKEN_SECRET`

## Basic Concepts

### Modal is Serverless

Modal is a serverless platform - only pay for resources used and spin up containers on demand in seconds.

### Core Components

**App**: Represents an application running on Modal, grouping one or more Functions for atomic deployment.

**Function**: Acts as an independent unit that scales up and down independently. No containers run (and no charges) when there are no live inputs.

**Image**: The environment code runs in - a container snapshot with dependencies installed.

## First Modal App

Create a file `hello_modal.py`:

```python
import modal

app = modal.App(name="hello-modal")

@app.function()
def hello():
    print("Hello from Modal!")
    return "success"

@app.local_entrypoint()
def main():
    hello.remote()
```

Run with:
```bash
modal run hello_modal.py
```

## Running Apps

### Ephemeral Apps (Development)

Run temporarily with `modal run`:
```bash
modal run script.py
```

The app stops when the script exits. Use `--detach` to keep running after client exits.

### Deployed Apps (Production)

Deploy persistently with `modal deploy`:
```bash
modal deploy script.py
```

View deployed apps at https://modal.com/apps or with:
```bash
modal app list
```

Stop deployed apps:
```bash
modal app stop app-name
```

## Key Features

- **Fast prototyping**: Write Python, run on GPUs in seconds
- **Serverless APIs**: Create web endpoints with a decorator
- **Scheduled jobs**: Run cron jobs in the cloud
- **GPU inference**: Access T4, L4, A10, A100, H100, H200, B200 GPUs
- **Distributed volumes**: Persistent storage for ML models
- **Sandboxes**: Secure containers for untrusted code
