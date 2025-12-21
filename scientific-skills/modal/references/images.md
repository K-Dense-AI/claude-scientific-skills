# Modal Images

## Overview

Modal Images define the environment code runs in - containers with dependencies installed. Images are built from method chains starting from a base image.

## Base Images

Start with a base image and chain methods:

```python
image = (
    modal.Image.debian_slim(python_version="3.13")
    .apt_install("git")
    .uv_pip_install("torch<3")
    .env({"HALT_AND_CATCH_FIRE": "0"})
    .run_commands("git clone https://github.com/modal-labs/agi")
)
```

Available base images:
- `Image.debian_slim()` - Debian Linux with Python
- `Image.micromamba()` - Base with Micromamba package manager
- `Image.from_registry()` - Pull from Docker Hub, ECR, etc.
- `Image.from_dockerfile()` - Build from existing Dockerfile

## Installing Python Packages

### With uv (Recommended)

Use `.uv_pip_install()` for fast package installation:

```python
image = (
    modal.Image.debian_slim()
    .uv_pip_install("pandas==2.2.0", "numpy")
)
```

### With pip

Fallback to standard pip if needed:

```python
image = (
    modal.Image.debian_slim(python_version="3.13")
    .pip_install("pandas==2.2.0", "numpy")
)
```

Pin dependencies tightly (e.g., `"torch==2.8.0"`) for reproducibility.

## Installing System Packages

Install Linux packages with apt:

```python
image = modal.Image.debian_slim().apt_install("git", "curl")
```

## Setting Environment Variables

Pass a dictionary to `.env()`:

```python
image = modal.Image.debian_slim().env({"PORT": "6443"})
```

## Running Shell Commands

Execute commands during image build:

```python
image = (
    modal.Image.debian_slim()
    .apt_install("git")
    .run_commands("git clone https://github.com/modal-labs/gpu-glossary")
)
```

## Running Python Functions at Build Time

Download model weights or perform setup:

```python
def download_models():
    import diffusers
    model_name = "segmind/small-sd"
    pipe = diffusers.StableDiffusionPipeline.from_pretrained(model_name)

hf_cache = modal.Volume.from_name("hf-cache")

image = (
    modal.Image.debian_slim()
    .pip_install("diffusers[torch]", "transformers")
    .run_function(
        download_models,
        secrets=[modal.Secret.from_name("huggingface-secret")],
        volumes={"/root/.cache/huggingface": hf_cache},
    )
)
```

## Adding Local Files

### Add Files or Directories

```python
image = modal.Image.debian_slim().add_local_dir(
    "/user/erikbern/.aws",
    remote_path="/root/.aws"
)
```

By default, files are added at container startup. Use `copy=True` to include in built image.

### Add Python Source

Add importable Python modules:

```python
image = modal.Image.debian_slim().add_local_python_source("local_module")

@app.function(image=image)
def f():
    import local_module
    local_module.do_stuff()
```

## Using Existing Container Images

### From Public Registry

```python
sklearn_image = modal.Image.from_registry("huanjason/scikit-learn")

@app.function(image=sklearn_image)
def fit_knn():
    from sklearn.neighbors import KNeighborsClassifier
    ...
```

Can pull from Docker Hub, Nvidia NGC, AWS ECR, GitHub ghcr.io.

### From Private Registry

Use Modal Secrets for authentication:

**Docker Hub**:
```python
secret = modal.Secret.from_name("my-docker-secret")
image = modal.Image.from_registry(
    "private-repo/image:tag",
    secret=secret
)
```

**AWS ECR**:
```python
aws_secret = modal.Secret.from_name("my-aws-secret")
image = modal.Image.from_aws_ecr(
    "000000000000.dkr.ecr.us-east-1.amazonaws.com/my-private-registry:latest",
    secret=aws_secret,
)
```

### From Dockerfile

```python
image = modal.Image.from_dockerfile("Dockerfile")

@app.function(image=image)
def fit():
    import sklearn
    ...
```

Can still extend with other image methods after importing.

## Using Micromamba

For coordinated installation of Python and system packages:

```python
numpyro_pymc_image = (
    modal.Image.micromamba()
    .micromamba_install("pymc==5.10.4", "numpyro==0.13.2", channels=["conda-forge"])
)
```

## GPU Support at Build Time

Run build steps on GPU instances:

```python
image = (
    modal.Image.debian_slim()
    .pip_install("bitsandbytes", gpu="H100")
)
```

## Image Caching

Images are cached per layer. Breaking cache on one layer causes cascading rebuilds for subsequent layers.

Define frequently-changing layers last to maximize cache reuse.

### Force Rebuild

```python
image = (
    modal.Image.debian_slim()
    .apt_install("git")
    .pip_install("slack-sdk", force_build=True)
)
```

Or set environment variable:
```bash
MODAL_FORCE_BUILD=1 modal run ...
```

## Handling Different Local/Remote Packages

Import packages only available remotely inside function bodies:

```python
@app.function(image=image)
def my_function():
    import pandas as pd  # Only imported remotely
    df = pd.DataFrame()
    ...
```

Or use the imports context manager:

```python
pandas_image = modal.Image.debian_slim().pip_install("pandas")

with pandas_image.imports():
    import pandas as pd

@app.function(image=pandas_image)
def my_function():
    df = pd.DataFrame()
```

## Fast Pull from Registry with eStargz

Improve pull performance with eStargz compression:

```bash
docker buildx build --tag "<registry>/<namespace>/<repo>:<version>" \
  --output type=registry,compression=estargz,force-compression=true,oci-mediatypes=true \
  .
```

Supported registries:
- AWS ECR
- Docker Hub
- Google Artifact Registry
