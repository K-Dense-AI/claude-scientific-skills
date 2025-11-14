# GPU Acceleration on Modal

## Quick Start

Run functions on GPUs with the `gpu` parameter:

```python
import modal

image = modal.Image.debian_slim().pip_install("torch")
app = modal.App(image=image)

@app.function(gpu="A100")
def run():
    import torch
    assert torch.cuda.is_available()
```

## Available GPU Types

Modal supports the following GPUs:

- `T4` - Entry-level GPU
- `L4` - Balanced performance and cost
- `A10` - Up to 4 GPUs, 96 GB total
- `A100` - 40GB or 80GB variants
- `A100-40GB` - Specific 40GB variant
- `A100-80GB` - Specific 80GB variant
- `L40S` - 48 GB, excellent for inference
- `H100` / `H100!` - Top-tier Hopper architecture
- `H200` - Improved Hopper with more memory
- `B200` - Latest Blackwell architecture

See https://modal.com/pricing for pricing.

## GPU Count

Request multiple GPUs per container with `:n` syntax:

```python
@app.function(gpu="H100:8")
def run_llama_405b():
    # 8 H100 GPUs available
    ...
```

Supported counts:
- B200, H200, H100, A100, L4, T4, L40S: up to 8 GPUs (up to 1,536 GB)
- A10: up to 4 GPUs (up to 96 GB)

Note: Requesting >2 GPUs may result in longer wait times.

## GPU Selection Guide

**For Inference (Recommended)**: Start with L40S
- Excellent cost/performance
- 48 GB memory
- Good for LLaMA, Stable Diffusion, etc.

**For Training**: Consider H100 or A100
- High compute throughput
- Large memory for batch processing

**For Memory-Bound Tasks**: H200 or A100-80GB
- More memory capacity
- Better for large models

## B200 GPUs

NVIDIA's flagship Blackwell chip:

```python
@app.function(gpu="B200:8")
def run_deepseek():
    # Most powerful option
    ...
```

## H200 and H100 GPUs

Hopper architecture GPUs with excellent software support:

```python
@app.function(gpu="H100")
def train():
    ...
```

### Automatic H200 Upgrades

Modal may upgrade `gpu="H100"` to H200 at no extra cost. H200 provides:
- 141 GB memory (vs 80 GB for H100)
- 4.8 TB/s bandwidth (vs 3.35 TB/s)

To avoid automatic upgrades (e.g., for benchmarking):
```python
@app.function(gpu="H100!")
def benchmark():
    ...
```

## A100 GPUs

Ampere architecture with 40GB or 80GB variants:

```python
# May be automatically upgraded to 80GB
@app.function(gpu="A100")
def qwen_7b():
    ...

# Specific variants
@app.function(gpu="A100-40GB")
def model_40gb():
    ...

@app.function(gpu="A100-80GB")
def llama_70b():
    ...
```

## GPU Fallbacks

Specify multiple GPU types with fallback:

```python
@app.function(gpu=["H100", "A100-40GB:2"])
def run_on_80gb():
    # Tries H100 first, falls back to 2x A100-40GB
    ...
```

Modal respects ordering and allocates most preferred available GPU.

## Multi-GPU Training

Modal supports multi-GPU training on a single node. Multi-node training is in closed beta.

### PyTorch Example

For frameworks that re-execute entrypoints, use subprocess or specific strategies:

```python
@app.function(gpu="A100:2")
def train():
    import subprocess
    import sys
    subprocess.run(
        ["python", "train.py"],
        stdout=sys.stdout,
        stderr=sys.stderr,
        check=True,
    )
```

For PyTorch Lightning, set strategy to `ddp_spawn` or `ddp_notebook`.

## Performance Considerations

**Memory-Bound vs Compute-Bound**:
- Running models with small batch sizes is memory-bound
- Newer GPUs have faster arithmetic than memory access
- Speedup from newer hardware may not justify cost for memory-bound workloads

**Optimization**:
- Use batching when possible
- Consider L40S before jumping to H100/B200
- Profile to identify bottlenecks
