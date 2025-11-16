---
name: modal
description: Run Python code in the cloud with serverless containers, GPUs, and autoscaling. Use when deploying ML models, running batch processing jobs, scheduling compute-intensive tasks, or serving APIs that require GPU acceleration or dynamic scaling.
---

# Modal

## Overview

Modal is a serverless platform for running Python code in the cloud with minimal configuration. Execute functions on powerful GPUs, scale automatically to thousands of containers, and pay only for compute used.

Modal is particularly suited for AI/ML workloads, high-performance batch processing, scheduled jobs, GPU inference, and serverless APIs. Sign up for free at https://modal.com and receive $30/month in credits.

## When to Use This Skill

Use Modal for:
- Deploying and serving ML models (LLMs, image generation, embedding models)
- Running GPU-accelerated computation (training, inference, rendering)
- Batch processing large datasets in parallel
- Scheduling compute-intensive jobs (daily data processing, model training)
- Building serverless APIs that need automatic scaling
- Scientific computing requiring distributed compute or specialized hardware

## Authentication and Setup

Modal requires authentication via API token.

### Initial Setup

```bash
# Install Modal
uv uv pip install modal

# Authenticate (opens browser for login)
modal token new
```

This creates a token stored in `~/.modal.toml`. The token authenticates all Modal operations.

### Verify Setup

```python
import modal

app = modal.App("test-app")

@app.function()
def hello():
    print("Modal is working!")
```

Run with: `modal run script.py`

## Core Capabilities

Modal provides serverless Python execution through Functions that run in containers. Define compute requirements, dependencies, and scaling behavior declaratively.

### 1. Define Container Images

Specify dependencies and environment for functions using Modal Images.

```python
import modal

# Basic image with Python packages
image = (
    modal.Image.debian_slim(python_version="3.12")
    .uv_pip_install("torch", "transformers", "numpy")
)

app = modal.App("ml-app", image=image)
```

**Common patterns:**
- Install Python packages: `.uv_pip_install("pandas", "scikit-learn")`
- Install system packages: `.apt_install("ffmpeg", "git")`
- Use existing Docker images: `modal.Image.from_registry("nvidia/cuda:12.1.0-base")`
- Add local code: `.add_local_python_source("my_module")`

See `references/images.md` for comprehensive image building documentation.

### 2. Create Functions

Define functions that run in the cloud with the `@app.function()` decorator.

```python
@app.function()
def process_data(file_path: str):
    import pandas as pd
    df = pd.read_csv(file_path)
    return df.describe()
```

**Call functions:**
```python
# From local entrypoint
@app.local_entrypoint()
def main():
    result = process_data.remote("data.csv")
    print(result)
```

Run with: `modal run script.py`

See `references/functions.md` for function patterns, deployment, and parameter handling.

### 3. Request GPUs

Attach GPUs to functions for accelerated computation.

```python
@app.function(gpu="H100")
def train_model():
    import torch
    assert torch.cuda.is_available()
    # GPU-accelerated code here
```

**Available GPU types:**
- `T4`, `L4` - Cost-effective inference
- `A10`, `A100`, `A100-80GB` - Standard training/inference
- `L40S` - Excellent cost/performance balance (48GB)
- `H100`, `H200` - High-performance training
- `B200` - Flagship performance (most powerful)

**Request multiple GPUs:**
```python
@app.function(gpu="H100:8")  # 8x H100 GPUs
def train_large_model():
    pass
```

See `references/gpu.md` for GPU selection guidance, CUDA setup, and multi-GPU configuration.

### 4. Configure Resources

Request CPU cores, memory, and disk for functions.

```python
@app.function(
    cpu=8.0,           # 8 physical cores
    memory=32768,      # 32 GiB RAM
    ephemeral_disk=10240  # 10 GiB disk
)
def memory_intensive_task():
    pass
```

Default allocation: 0.125 CPU cores, 128 MiB memory. Billing based on reservation or actual usage, whichever is higher.

See `references/resources.md` for resource limits and billing details.

### 5. Scale Automatically

Modal autoscales functions from zero to thousands of containers based on demand.

**Process inputs in parallel:**
```python
@app.function()
def analyze_sample(sample_id: int):
    # Process single sample
    return result

@app.local_entrypoint()
def main():
    sample_ids = range(1000)
    # Automatically parallelized across containers
    results = list(analyze_sample.map(sample_ids))
```

**Configure autoscaling:**
```python
@app.function(
    max_containers=100,      # Upper limit
    min_containers=2,        # Keep warm
    buffer_containers=5      # Idle buffer for bursts
)
def inference():
    pass
```

See `references/scaling.md` for autoscaling configuration, concurrency, and scaling limits.

### 6. Store Data Persistently

Use Volumes for persistent storage across function invocations.

```python
volume = modal.Volume.from_name("my-data", create_if_missing=True)

@app.function(volumes={"/data": volume})
def save_results(data):
    with open("/data/results.txt", "w") as f:
        f.write(data)
    volume.commit()  # Persist changes
```

Volumes persist data between runs, store model weights, cache datasets, and share data between functions.

See `references/volumes.md` for volume management, commits, and caching patterns.

### 7. Manage Secrets

Store API keys and credentials securely using Modal Secrets.

```python
@app.function(secrets=[modal.Secret.from_name("huggingface")])
def download_model():
    import os
    token = os.environ["HF_TOKEN"]
    # Use token for authentication
```

**Create secrets in Modal dashboard or via CLI:**
```bash
modal secret create my-secret KEY=value API_TOKEN=xyz
```

See `references/secrets.md` for secret management and authentication patterns.

### 8. Deploy Web Endpoints

Serve HTTP endpoints, APIs, and webhooks with `@modal.web_endpoint()`.

```python
@app.function()
@modal.web_endpoint(method="POST")
def predict(data: dict):
    # Process request
    result = model.predict(data["input"])
    return {"prediction": result}
```

**Deploy with:**
```bash
modal deploy script.py
```

Modal provides HTTPS URL for the endpoint.

See `references/web-endpoints.md` for FastAPI integration, streaming, authentication, and WebSocket support.

### 9. Schedule Jobs

Run functions on a schedule with cron expressions.

```python
@app.function(schedule=modal.Cron("0 2 * * *"))  # Daily at 2 AM
def daily_backup():
    # Backup data
    pass

@app.function(schedule=modal.Period(hours=4))  # Every 4 hours
def refresh_cache():
    # Update cache
    pass
```

Scheduled functions run automatically without manual invocation.

See `references/scheduled-jobs.md` for cron syntax, timezone configuration, and monitoring.

## Common Workflows

### Deploy ML Model for Inference

```python
import modal

# Define dependencies
image = modal.Image.debian_slim().uv_pip_install("torch", "transformers")
app = modal.App("llm-inference", image=image)

# Download model at build time
@app.function()
def download_model():
    from transformers import AutoModel
    AutoModel.from_pretrained("bert-base-uncased")

# Serve model
@app.cls(gpu="L40S")
class Model:
    @modal.enter()
    def load_model(self):
        from transformers import pipeline
        self.pipe = pipeline("text-classification", device="cuda")

    @modal.method()
    def predict(self, text: str):
        return self.pipe(text)

@app.local_entrypoint()
def main():
    model = Model()
    result = model.predict.remote("Modal is great!")
    print(result)
```

### Batch Process Large Dataset

```python
@app.function(cpu=2.0, memory=4096)
def process_file(file_path: str):
    import pandas as pd
    df = pd.read_csv(file_path)
    # Process data
    return df.shape[0]

@app.local_entrypoint()
def main():
    files = ["file1.csv", "file2.csv", ...]  # 1000s of files
    # Automatically parallelized across containers
    for count in process_file.map(files):
        print(f"Processed {count} rows")
```

### Train Model on GPU

```python
@app.function(
    gpu="A100:2",      # 2x A100 GPUs
    timeout=3600       # 1 hour timeout
)
def train_model(config: dict):
    import torch
    # Multi-GPU training code
    model = create_model(config)
    train(model)
    return metrics
```

## Reference Documentation

Detailed documentation for specific features:

- **`references/getting-started.md`** - Authentication, setup, basic concepts
- **`references/images.md`** - Image building, dependencies, Dockerfiles
- **`references/functions.md`** - Function patterns, deployment, parameters
- **`references/gpu.md`** - GPU types, CUDA, multi-GPU configuration
- **`references/resources.md`** - CPU, memory, disk management
- **`references/scaling.md`** - Autoscaling, parallel execution, concurrency
- **`references/volumes.md`** - Persistent storage, data management
- **`references/secrets.md`** - Environment variables, authentication
- **`references/web-endpoints.md`** - APIs, webhooks, endpoints
- **`references/scheduled-jobs.md`** - Cron jobs, periodic tasks
- **`references/examples.md`** - Common patterns for scientific computing

## Best Practices

1. **Pin dependencies** in `.uv_pip_install()` for reproducible builds
2. **Use appropriate GPU types** - L40S for inference, H100/A100 for training
3. **Leverage caching** - Use Volumes for model weights and datasets
4. **Configure autoscaling** - Set `max_containers` and `min_containers` based on workload
5. **Import packages in function body** if not available locally
6. **Use `.map()` for parallel processing** instead of sequential loops
7. **Store secrets securely** - Never hardcode API keys
8. **Monitor costs** - Check Modal dashboard for usage and billing

## Troubleshooting

**"Module not found" errors:**
- Add packages to image with `.uv_pip_install("package-name")`
- Import packages inside function body if not available locally

**GPU not detected:**
- Verify GPU specification: `@app.function(gpu="A100")`
- Check CUDA availability: `torch.cuda.is_available()`

**Function timeout:**
- Increase timeout: `@app.function(timeout=3600)`
- Default timeout is 5 minutes

**Volume changes not persisting:**
- Call `volume.commit()` after writing files
- Verify volume mounted correctly in function decorator

For additional help, see Modal documentation at https://modal.com/docs or join Modal Slack community.
