---
name: pytorch-lightning
description: PyTorch Lightning deep learning framework skill for organizing PyTorch code and automating training workflows. Use this skill for: creating LightningModules with training_step/validation_step hooks, implementing DataModules for data loading and preprocessing, configuring Trainer with accelerators/devices/strategies, setting up distributed training (DDP/FSDP/DeepSpeed), implementing callbacks (ModelCheckpoint/EarlyStopping), configuring loggers (TensorBoard/WandB/MLflow), converting PyTorch code to Lightning format, optimizing performance with mixed precision/gradient accumulation, debugging with fast_dev_run/overfit_batches, checkpointing and resuming training, hyperparameter tuning with Tuner, handling multi-GPU/multi-node training, memory optimization for large models, experiment tracking and reproducibility, custom training loops, validation/testing workflows, prediction pipelines, and production deployment. Includes templates, API references, distributed training guides, and best practices for efficient deep learning development.
---

# PyTorch Lightning

## Overview

PyTorch Lightning is a deep learning framework that organizes PyTorch code to decouple research from engineering. It automates training loop complexity (multi-GPU, mixed precision, checkpointing, logging) while maintaining full flexibility over model architecture and training logic.

**Core Philosophy:** Separate concerns
- **LightningModule** - Research code (model architecture, training logic)
- **Trainer** - Engineering automation (hardware, optimization, logging)
- **DataModule** - Data processing (downloading, loading, transforms)
- **Callbacks** - Non-essential functionality (checkpointing, early stopping)

## When to Use This Skill

Use this skill when:
- Building or training deep learning models with PyTorch
- Converting existing PyTorch code to Lightning structure
- Setting up distributed training across multiple GPUs or nodes
- Implementing custom training loops with validation and testing
- Organizing data processing pipelines
- Configuring experiment logging and model checkpointing
- Optimizing training performance and memory usage
- Working with large models requiring model parallelism

## Quick Start

### Basic Lightning Workflow

1. **Define a LightningModule** (organize your model)
2. **Create a DataModule or DataLoaders** (organize your data)
3. **Configure a Trainer** (automate training)
4. **Train** with `trainer.fit()`

### Minimal Example

```python
import lightning as L
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset

# 1. Define LightningModule
class SimpleModel(L.LightningModule):
    def __init__(self, input_dim, output_dim):
        super().__init__()
        self.save_hyperparameters()
        self.model = nn.Linear(input_dim, output_dim)

    def forward(self, x):
        return self.model(x)

    def training_step(self, batch, batch_idx):
        x, y = batch
        y_hat = self(x)
        loss = nn.functional.mse_loss(y_hat, y)
        self.log('train_loss', loss)
        return loss

    def configure_optimizers(self):
        return torch.optim.Adam(self.parameters(), lr=1e-3)

# 2. Prepare data
train_data = TensorDataset(torch.randn(1000, 10), torch.randn(1000, 1))
train_loader = DataLoader(train_data, batch_size=32)

# 3. Create Trainer
trainer = L.Trainer(max_epochs=10, accelerator='auto')

# 4. Train
model = SimpleModel(input_dim=10, output_dim=1)
trainer.fit(model, train_loader)
```

## Core Workflows

### 1. Creating a LightningModule

Structure model code by implementing essential hooks:

**Template:** Use `scripts/template_lightning_module.py` as a starting point.

```python
class MyLightningModule(L.LightningModule):
    def __init__(self, hyperparameters):
        super().__init__()
        self.save_hyperparameters()  # Save for checkpointing
        self.model = YourModel()

    def forward(self, x):
        """Inference forward pass."""
        return self.model(x)

    def training_step(self, batch, batch_idx):
        """Define training loop logic."""
        x, y = batch
        y_hat = self(x)
        loss = self.compute_loss(y_hat, y)
        self.log('train_loss', loss, on_step=True, on_epoch=True)
        return loss

    def validation_step(self, batch, batch_idx):
        """Define validation logic."""
        x, y = batch
        y_hat = self(x)
        loss = self.compute_loss(y_hat, y)
        self.log('val_loss', loss, on_epoch=True)
        return loss

    def configure_optimizers(self):
        """Return optimizer and optional scheduler."""
        optimizer = torch.optim.Adam(self.parameters(), lr=self.hparams.lr)
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer)
        return {
            "optimizer": optimizer,
            "lr_scheduler": {
                "scheduler": scheduler,
                "monitor": "val_loss",
            }
        }
```

**Key Points:**
- Use `self.save_hyperparameters()` to automatically save init args
- Use `self.log()` to track metrics across loggers
- Return loss from training_step for automatic optimization
- Keep model architecture separate from training logic

### 2. Creating a DataModule

Organize all data processing in a reusable module:

**Template:** Use `scripts/template_datamodule.py` as a starting point.

```python
class MyDataModule(L.LightningDataModule):
    def __init__(self, data_dir, batch_size=32):
        super().__init__()
        self.save_hyperparameters()

    def prepare_data(self):
        """Download data (called once, single process)."""
        # Download datasets, tokenize, etc.
        pass

    def setup(self, stage=None):
        """Create datasets (called on every process)."""
        if stage == 'fit' or stage is None:
            # Create train/val datasets
            self.train_dataset = ...
            self.val_dataset = ...

        if stage == 'test' or stage is None:
            # Create test dataset
            self.test_dataset = ...

    def train_dataloader(self):
        return DataLoader(self.train_dataset, batch_size=self.hparams.batch_size)

    def val_dataloader(self):
        return DataLoader(self.val_dataset, batch_size=self.hparams.batch_size)

    def test_dataloader(self):
        return DataLoader(self.test_dataset, batch_size=self.hparams.batch_size)
```

**Key Points:**
- `prepare_data()` for downloading (single process)
- `setup()` for creating datasets (every process)
- Use `stage` parameter to separate fit/test logic
- Makes data code reusable across projects

### 3. Configuring the Trainer

The Trainer automates training complexity:

**Helper:** Use `scripts/quick_trainer_setup.py` for preset configurations.

```python
from lightning.pytorch.callbacks import ModelCheckpoint, EarlyStopping

trainer = L.Trainer(
    # Training duration
    max_epochs=100,

    # Hardware
    accelerator='auto',  # 'cpu', 'gpu', 'tpu'
    devices=1,  # Number of devices or specific IDs

    # Optimization
    precision='16-mixed',  # Mixed precision training
    gradient_clip_val=1.0,
    accumulate_grad_batches=4,  # Gradient accumulation

    # Validation
    check_val_every_n_epoch=1,
    val_check_interval=1.0,  # Validate every epoch

    # Logging
    log_every_n_steps=50,
    logger=TensorBoardLogger('logs/'),

    # Callbacks
    callbacks=[
        ModelCheckpoint(monitor='val_loss', mode='min'),
        EarlyStopping(monitor='val_loss', patience=10),
    ],

    # Debugging
    fast_dev_run=False,  # Quick test with few batches
    enable_progress_bar=True,
)
```

**Common Presets:**

```python
from scripts.quick_trainer_setup import create_trainer

# Development preset (fast debugging)
trainer = create_trainer(preset='fast_dev', max_epochs=3)

# Production preset (full features)
trainer = create_trainer(preset='production', max_epochs=100)

# Distributed preset (multi-GPU)
trainer = create_trainer(preset='distributed', devices=4)
```

### 4. Training and Evaluation

```python
# Training
trainer.fit(model, datamodule=dm)
# Or with dataloaders
trainer.fit(model, train_loader, val_loader)

# Resume from checkpoint
trainer.fit(model, datamodule=dm, ckpt_path='checkpoint.ckpt')

# Testing
trainer.test(model, datamodule=dm)
# Or load best checkpoint
trainer.test(ckpt_path='best', datamodule=dm)

# Prediction
predictions = trainer.predict(model, predict_loader)

# Validation only
trainer.validate(model, datamodule=dm)
```

### 5. Distributed Training

Lightning handles distributed training automatically:

```python
# Single machine, multiple GPUs (Data Parallel)
trainer = L.Trainer(
    accelerator='gpu',
    devices=4,
    strategy='ddp',  # DistributedDataParallel
)

# Multiple machines, multiple GPUs
trainer = L.Trainer(
    accelerator='gpu',
    devices=4,  # GPUs per node
    num_nodes=8,  # Number of machines
    strategy='ddp',
)

# Large models (Model Parallel with FSDP)
trainer = L.Trainer(
    accelerator='gpu',
    devices=4,
    strategy='fsdp',  # Fully Sharded Data Parallel
)

# Large models (Model Parallel with DeepSpeed)
trainer = L.Trainer(
    accelerator='gpu',
    devices=4,
    strategy='deepspeed_stage_2',
    precision='16-mixed',
)
```

**For detailed distributed training guide, see:** `references/distributed_training.md`

**Strategy Selection:**
- Models < 500M params → Use `ddp`
- Models > 500M params → Use `fsdp` or `deepspeed`
- Maximum memory efficiency → Use DeepSpeed Stage 3 with offloading
- Native PyTorch → Use `fsdp`
- Cutting-edge features → Use `deepspeed`

### 6. Callbacks

Extend training with modular functionality:

```python
from lightning.pytorch.callbacks import (
    ModelCheckpoint,
    EarlyStopping,
    LearningRateMonitor,
    RichProgressBar,
)

callbacks = [
    # Save best models
    ModelCheckpoint(
        monitor='val_loss',
        mode='min',
        save_top_k=3,
        filename='{epoch}-{val_loss:.2f}',
    ),

    # Stop when no improvement
    EarlyStopping(
        monitor='val_loss',
        patience=10,
        mode='min',
    ),

    # Log learning rate
    LearningRateMonitor(logging_interval='epoch'),

    # Rich progress bar
    RichProgressBar(),
]

trainer = L.Trainer(callbacks=callbacks)
```

**Custom Callbacks:**

```python
from lightning.pytorch.callbacks import Callback

class MyCustomCallback(Callback):
    def on_train_epoch_end(self, trainer, pl_module):
        # Custom logic at end of each epoch
        print(f"Epoch {trainer.current_epoch} completed")

    def on_validation_end(self, trainer, pl_module):
        val_loss = trainer.callback_metrics.get('val_loss')
        # Custom validation logic
        pass
```

### 7. Logging

Track experiments with various loggers:

```python
from lightning.pytorch.loggers import (
    TensorBoardLogger,
    WandbLogger,
    CSVLogger,
    MLFlowLogger,
)

# Single logger
logger = TensorBoardLogger('logs/', name='my_experiment')

# Multiple loggers
loggers = [
    TensorBoardLogger('logs/'),
    WandbLogger(project='my_project'),
    CSVLogger('logs/'),
]

trainer = L.Trainer(logger=loggers)
```

**Logging in LightningModule:**

```python
def training_step(self, batch, batch_idx):
    loss = self.compute_loss(batch)

    # Log single metric
    self.log('train_loss', loss, on_step=True, on_epoch=True, prog_bar=True)

    # Log multiple metrics
    metrics = {'loss': loss, 'acc': acc, 'f1': f1}
    self.log_dict(metrics, on_step=True, on_epoch=True)

    return loss
```

## Converting Existing PyTorch Code

### Standard PyTorch → Lightning

**Before (PyTorch):**
```python
model = MyModel()
optimizer = torch.optim.Adam(model.parameters())

for epoch in range(num_epochs):
    for batch in train_loader:
        optimizer.zero_grad()
        x, y = batch
        y_hat = model(x)
        loss = F.cross_entropy(y_hat, y)
        loss.backward()
        optimizer.step()
```

**After (Lightning):**
```python
class MyLightningModel(L.LightningModule):
    def __init__(self):
        super().__init__()
        self.model = MyModel()

    def training_step(self, batch, batch_idx):
        x, y = batch
        y_hat = self.model(x)
        loss = F.cross_entropy(y_hat, y)
        self.log('train_loss', loss)
        return loss

    def configure_optimizers(self):
        return torch.optim.Adam(self.parameters())

trainer = L.Trainer(max_epochs=num_epochs)
trainer.fit(model, train_loader)
```

**Key Changes:**
1. Wrap model in LightningModule
2. Move training loop logic to `training_step()`
3. Move optimizer setup to `configure_optimizers()`
4. Replace manual loop with `trainer.fit()`
5. Lightning handles: `.zero_grad()`, `.backward()`, `.step()`, device placement

## Common Patterns

### Reproducibility

```python
from lightning.pytorch import seed_everything

# Set seed for reproducibility
seed_everything(42, workers=True)

trainer = L.Trainer(deterministic=True)
```

### Mixed Precision Training

```python
# 16-bit mixed precision
trainer = L.Trainer(precision='16-mixed')

# BFloat16 mixed precision (more stable)
trainer = L.Trainer(precision='bf16-mixed')
```

### Gradient Accumulation

```python
# Effective batch size = 4x actual batch size
trainer = L.Trainer(accumulate_grad_batches=4)
```

### Learning Rate Finding

```python
from lightning.pytorch.tuner import Tuner

trainer = L.Trainer()
tuner = Tuner(trainer)

# Find optimal learning rate
lr_finder = tuner.lr_find(model, train_dataloader)
model.hparams.learning_rate = lr_finder.suggestion()

# Find optimal batch size
tuner.scale_batch_size(model, mode="power")
```

### Checkpointing and Loading

```python
# Save checkpoint
trainer.fit(model, datamodule=dm)
# Checkpoint automatically saved to checkpoints/

# Load from checkpoint
model = MyLightningModule.load_from_checkpoint('path/to/checkpoint.ckpt')

# Resume training
trainer.fit(model, datamodule=dm, ckpt_path='checkpoint.ckpt')

# Test from checkpoint
trainer.test(ckpt_path='best', datamodule=dm)
```

### Debugging

```python
# Quick test with few batches
trainer = L.Trainer(fast_dev_run=10)

# Overfit on small data (debug model)
trainer = L.Trainer(overfit_batches=100)

# Limit batches for quick iteration
trainer = L.Trainer(
    limit_train_batches=100,
    limit_val_batches=50,
)

# Profile training
trainer = L.Trainer(profiler='simple')  # or 'advanced'
```

## Best Practices

### Code Organization

1. **Separate concerns:**
   - Model architecture in `__init__()`
   - Training logic in `training_step()`
   - Validation logic in `validation_step()`
   - Data processing in DataModule

2. **Use `save_hyperparameters()`:**
   ```python
   def __init__(self, lr, hidden_dim, dropout):
       super().__init__()
       self.save_hyperparameters()  # Automatically saves all args
   ```

3. **Device-agnostic code:**
   ```python
   # Avoid manual device placement
   # BAD: tensor.cuda()
   # GOOD: Lightning handles this automatically

   # Create tensors on model's device
   new_tensor = torch.zeros(10, device=self.device)
   ```

4. **Log comprehensively:**
   ```python
   self.log('metric', value, on_step=True, on_epoch=True, prog_bar=True)
   ```

### Performance Optimization

1. **Use DataLoader best practices:**
   ```python
   DataLoader(
       dataset,
       batch_size=32,
       num_workers=4,  # Multiple workers
       pin_memory=True,  # Faster GPU transfer
       persistent_workers=True,  # Keep workers alive
   )
   ```

2. **Enable benchmark mode for fixed input sizes:**
   ```python
   trainer = L.Trainer(benchmark=True)
   ```

3. **Use gradient clipping:**
   ```python
   trainer = L.Trainer(gradient_clip_val=1.0)
   ```

4. **Enable mixed precision:**
   ```python
   trainer = L.Trainer(precision='16-mixed')
   ```

### Distributed Training

1. **Sync metrics across devices:**
   ```python
   self.log('metric', value, sync_dist=True)
   ```

2. **Rank-specific operations:**
   ```python
   if self.trainer.is_global_zero:
       # Only run on main process
       self.save_artifacts()
   ```

3. **Use appropriate strategy:**
   - Small models → `ddp`
   - Large models → `fsdp` or `deepspeed`

## Resources

### Scripts

Executable templates for quick implementation:

- **`template_lightning_module.py`** - Complete LightningModule template with all hooks, logging, and optimization patterns
- **`template_datamodule.py`** - Complete DataModule template with data loading, splitting, and transformation patterns
- **`quick_trainer_setup.py`** - Helper functions to create Trainers with preset configurations (development, production, distributed)

### References

Comprehensive documentation for deep-dive learning:

- **`api_reference.md`** - Complete API reference covering LightningModule hooks, Trainer parameters, Callbacks, DataModules, Loggers, and common patterns
- **`distributed_training.md`** - In-depth guide for distributed training strategies (DDP, FSDP, DeepSpeed), multi-node setup, memory optimization, and troubleshooting

Load references when needing detailed information:
```python
# Example: Load distributed training reference
# See references/distributed_training.md for comprehensive distributed training guide
```

## Troubleshooting

### Common Issues

**Out of Memory:**
- Reduce batch size
- Use gradient accumulation
- Enable mixed precision (`precision='16-mixed'`)
- Use FSDP or DeepSpeed for large models
- Enable activation checkpointing

**Slow Training:**
- Use multiple DataLoader workers (`num_workers > 0`)
- Enable `pin_memory=True` and `persistent_workers=True`
- Enable `benchmark=True` for fixed input sizes
- Profile with `profiler='simple'`

**Validation Not Running:**
- Check `check_val_every_n_epoch` setting
- Ensure validation data provided
- Verify `validation_step()` implemented

**Checkpoints Not Saving:**
- Ensure `enable_checkpointing=True`
- Check `ModelCheckpoint` callback configuration
- Verify `monitor` metric exists in logs

## Additional Resources

- Official Documentation: https://lightning.ai/docs/pytorch/stable/
- GitHub: https://github.com/Lightning-AI/lightning
- Community: https://lightning.ai/community

When unclear about specific functionality, refer to `references/api_reference.md` for detailed API documentation or `references/distributed_training.md` for distributed training specifics.
