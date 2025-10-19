# Distributed and Model Parallel Training

Comprehensive guide for distributed training strategies in PyTorch Lightning.

## Overview

PyTorch Lightning provides seamless distributed training across multiple GPUs, machines, and TPUs with minimal code changes. The framework automatically handles the complexity of distributed training while keeping code device-agnostic and readable.

## Training Strategies

### Data Parallel (DDP - DistributedDataParallel)

**Best for:** Most models (< 500M parameters) where the full model fits in GPU memory.

**How it works:** Each GPU holds a complete copy of the model and trains on a different batch subset. Gradients are synchronized across GPUs during backward pass.

```python
# Single-node, multi-GPU
trainer = Trainer(
    accelerator='gpu',
    devices=4,  # Use 4 GPUs
    strategy='ddp',
)

# Multi-node, multi-GPU
trainer = Trainer(
    accelerator='gpu',
    devices=4,  # GPUs per node
    num_nodes=2,  # Number of nodes
    strategy='ddp',
)
```

**Advantages:**
- Most widely used and tested
- Works with most PyTorch code
- Good scaling efficiency
- No code changes required in LightningModule

**When to use:** Default choice for most distributed training scenarios.

### FSDP (Fully Sharded Data Parallel)

**Best for:** Large models (500M+ parameters) that don't fit in single GPU memory.

**How it works:** Shards model parameters, gradients, and optimizer states across GPUs. Each GPU only stores a subset of the model.

```python
trainer = Trainer(
    accelerator='gpu',
    devices=4,
    strategy='fsdp',
)

# With configuration
from lightning.pytorch.strategies import FSDPStrategy

strategy = FSDPStrategy(
    sharding_strategy="FULL_SHARD",  # Full sharding
    cpu_offload=False,  # Offload to CPU
    mixed_precision=torch.float16,
)

trainer = Trainer(
    accelerator='gpu',
    devices=4,
    strategy=strategy,
)
```

**Sharding Strategies:**
- `FULL_SHARD` - Shard parameters, gradients, and optimizer states
- `SHARD_GRAD_OP` - Shard only gradients and optimizer states
- `NO_SHARD` - DDP-like (no sharding)
- `HYBRID_SHARD` - Shard within node, DDP across nodes

**Advanced FSDP Configuration:**
```python
from lightning.pytorch.strategies import FSDPStrategy

strategy = FSDPStrategy(
    sharding_strategy="FULL_SHARD",
    activation_checkpointing=True,  # Save memory
    cpu_offload=True,  # Offload parameters to CPU
    backward_prefetch="BACKWARD_PRE",  # Prefetch strategy
    forward_prefetch=True,
    limit_all_gathers=True,
)
```

**When to use:**
- Models > 500M parameters
- Limited GPU memory
- Native PyTorch solution preferred
- Migrating from standalone PyTorch FSDP

### DeepSpeed

**Best for:** Cutting-edge features, massive models, or existing DeepSpeed users.

**How it works:** Comprehensive optimization library with multiple stages of memory and compute optimization.

```python
# Basic DeepSpeed
trainer = Trainer(
    accelerator='gpu',
    devices=4,
    strategy='deepspeed',
    precision='16-mixed',
)

# With configuration
from lightning.pytorch.strategies import DeepSpeedStrategy

strategy = DeepSpeedStrategy(
    stage=2,  # ZeRO Stage (1, 2, or 3)
    offload_optimizer=True,
    offload_parameters=True,
)

trainer = Trainer(
    accelerator='gpu',
    devices=4,
    strategy=strategy,
)
```

**ZeRO Stages:**
- **Stage 1:** Shard optimizer states
- **Stage 2:** Shard optimizer states + gradients
- **Stage 3:** Shard optimizer states + gradients + parameters (like FSDP)

**With DeepSpeed Config File:**
```python
strategy = DeepSpeedStrategy(config="deepspeed_config.json")
```

Example `deepspeed_config.json`:
```json
{
  "zero_optimization": {
    "stage": 2,
    "offload_optimizer": {
      "device": "cpu",
      "pin_memory": true
    },
    "allgather_bucket_size": 2e8,
    "reduce_bucket_size": 2e8
  },
  "activation_checkpointing": {
    "partition_activations": true,
    "cpu_checkpointing": true
  },
  "fp16": {
    "enabled": true
  },
  "gradient_clipping": 1.0
}
```

**When to use:**
- Need specific DeepSpeed features
- Maximum memory efficiency required
- Already familiar with DeepSpeed
- Training extremely large models

### DDP Spawn

**Note:** Generally avoid using `ddp_spawn`. Use `ddp` instead.

```python
trainer = Trainer(strategy='ddp_spawn')  # Not recommended
```

**Issues with ddp_spawn:**
- Cannot return values from `.fit()`
- Pickling issues with unpicklable objects
- Slower than `ddp`
- More memory overhead

**When to use:** Only for debugging or if `ddp` doesn't work on your system.

## Multi-Node Training

### Basic Multi-Node Setup

```python
# On each node, run the same command
trainer = Trainer(
    accelerator='gpu',
    devices=4,  # GPUs per node
    num_nodes=8,  # Total number of nodes
    strategy='ddp',
)
```

### SLURM Cluster

Lightning automatically detects SLURM environment:

```python
trainer = Trainer(
    accelerator='gpu',
    devices=4,
    num_nodes=8,
    strategy='ddp',
)
```

**SLURM Submit Script:**
```bash
#!/bin/bash
#SBATCH --nodes=8
#SBATCH --gres=gpu:4
#SBATCH --ntasks-per-node=4
#SBATCH --job-name=lightning_training

python train.py
```

### Manual Cluster Setup

```python
from lightning.pytorch.strategies import DDPStrategy

strategy = DDPStrategy(
    cluster_environment='TorchElastic',  # or 'SLURM', 'LSF', 'Kubeflow'
)

trainer = Trainer(
    accelerator='gpu',
    devices=4,
    num_nodes=8,
    strategy=strategy,
)
```

## Memory Optimization Techniques

### Gradient Accumulation

Simulate larger batch sizes without increasing memory:

```python
trainer = Trainer(
    accumulate_grad_batches=4,  # Accumulate 4 batches before optimizer step
)

# Variable accumulation by epoch
trainer = Trainer(
    accumulate_grad_batches={
        0: 8,  # Epochs 0-4: accumulate 8 batches
        5: 4,  # Epochs 5+: accumulate 4 batches
    }
)
```

### Activation Checkpointing

Trade computation for memory by recomputing activations during backward pass:

```python
# FSDP
from torch.distributed.algorithms._checkpoint.checkpoint_wrapper import (
    checkpoint_wrapper,
    CheckpointImpl,
    apply_activation_checkpointing,
)

class MyModule(L.LightningModule):
    def configure_model(self):
        # Wrap specific layers for activation checkpointing
        self.model = MyTransformer()
        apply_activation_checkpointing(
            self.model,
            checkpoint_wrapper_fn=lambda m: checkpoint_wrapper(m, CheckpointImpl.NO_REENTRANT),
            check_fn=lambda m: isinstance(m, TransformerBlock),
        )
```

### Mixed Precision Training

Reduce memory usage and increase speed with mixed precision:

```python
# 16-bit mixed precision
trainer = Trainer(precision='16-mixed')

# BFloat16 mixed precision (more stable, requires newer GPUs)
trainer = Trainer(precision='bf16-mixed')
```

### CPU Offloading

Offload parameters or optimizer states to CPU:

```python
# FSDP with CPU offload
from lightning.pytorch.strategies import FSDPStrategy

strategy = FSDPStrategy(
    cpu_offload=True,  # Offload parameters to CPU
)

# DeepSpeed with CPU offload
from lightning.pytorch.strategies import DeepSpeedStrategy

strategy = DeepSpeedStrategy(
    stage=3,
    offload_optimizer=True,
    offload_parameters=True,
)
```

## Performance Optimization

### Synchronize Batch Normalization

Synchronize batch norm statistics across GPUs:

```python
trainer = Trainer(
    accelerator='gpu',
    devices=4,
    strategy='ddp',
    sync_batchnorm=True,  # Sync batch norm across GPUs
)
```

### Find Optimal Batch Size

```python
from lightning.pytorch.tuner import Tuner

trainer = Trainer()
tuner = Tuner(trainer)

# Auto-scale batch size
tuner.scale_batch_size(model, mode="power")  # or "binsearch"
```

### Gradient Clipping

Prevent gradient explosion in distributed training:

```python
trainer = Trainer(
    gradient_clip_val=1.0,
    gradient_clip_algorithm='norm',  # or 'value'
)
```

### Benchmark Mode

Enable cudnn.benchmark for consistent input sizes:

```python
trainer = Trainer(
    benchmark=True,  # Optimize for consistent input sizes
)
```

## Distributed Data Loading

### Automatic Distributed Sampling

Lightning automatically handles distributed sampling:

```python
# No changes needed - Lightning handles this automatically
def train_dataloader(self):
    return DataLoader(
        self.train_dataset,
        batch_size=32,
        shuffle=True,  # Lightning converts to DistributedSampler
    )
```

### Manual Control

```python
# Disable automatic distributed sampler
trainer = Trainer(
    use_distributed_sampler=False,
)

# Manual distributed sampler
from torch.utils.data.distributed import DistributedSampler

def train_dataloader(self):
    sampler = DistributedSampler(self.train_dataset)
    return DataLoader(
        self.train_dataset,
        batch_size=32,
        sampler=sampler,
    )
```

### Data Loading Best Practices

```python
def train_dataloader(self):
    return DataLoader(
        self.train_dataset,
        batch_size=32,
        num_workers=4,  # Use multiple workers
        pin_memory=True,  # Faster CPU-GPU transfer
        persistent_workers=True,  # Keep workers alive between epochs
    )
```

## Common Patterns

### Logging in Distributed Training

```python
def training_step(self, batch, batch_idx):
    loss = self.compute_loss(batch)

    # Automatically syncs across processes
    self.log('train_loss', loss, sync_dist=True)

    return loss
```

### Rank-Specific Operations

```python
def training_step(self, batch, batch_idx):
    # Run only on rank 0 (main process)
    if self.trainer.is_global_zero:
        print("This only prints once across all processes")

    # Get current rank
    rank = self.trainer.global_rank
    world_size = self.trainer.world_size

    return loss
```

### Barrier Synchronization

```python
def on_train_epoch_end(self):
    # Wait for all processes
    self.trainer.strategy.barrier()

    # Now all processes are synchronized
    if self.trainer.is_global_zero:
        # Save something only once
        self.save_artifacts()
```

## Troubleshooting

### Common Issues

**1. Out of Memory:**
- Reduce batch size
- Enable gradient accumulation
- Use FSDP or DeepSpeed
- Enable activation checkpointing
- Use mixed precision

**2. Slow Training:**
- Check data loading (use `num_workers > 0`)
- Enable `pin_memory=True` and `persistent_workers=True`
- Use `benchmark=True` for consistent input sizes
- Profile with `profiler='simple'`

**3. Hanging:**
- Ensure all processes execute same collectives
- Check for `if` statements that differ across ranks
- Use barrier synchronization when needed

**4. Inconsistent Results:**
- Set `deterministic=True`
- Use `seed_everything()`
- Ensure proper gradient synchronization

### Debugging Distributed Training

```python
# Test with single GPU first
trainer = Trainer(accelerator='gpu', devices=1)

# Then test with 2 GPUs
trainer = Trainer(accelerator='gpu', devices=2, strategy='ddp')

# Use fast_dev_run for quick testing
trainer = Trainer(
    accelerator='gpu',
    devices=2,
    strategy='ddp',
    fast_dev_run=10,  # Run 10 batches only
)
```

## Strategy Selection Guide

| Model Size | Available Memory | Recommended Strategy |
|-----------|------------------|---------------------|
| < 500M params | Fits in 1 GPU | Single GPU |
| < 500M params | Fits across GPUs | DDP |
| 500M - 3B params | Limited memory | FSDP or DeepSpeed Stage 2 |
| 3B+ params | Very limited memory | FSDP or DeepSpeed Stage 3 |
| Any size | Maximum efficiency | DeepSpeed with offloading |
| Multiple nodes | Any | DDP (< 500M) or FSDP/DeepSpeed (> 500M) |
