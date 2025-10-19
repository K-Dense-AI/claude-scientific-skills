# PyTorch Lightning API Reference

Comprehensive reference for PyTorch Lightning core APIs, hooks, and components.

## LightningModule

The LightningModule is the core abstraction for organizing PyTorch code in Lightning.

### Essential Hooks

#### `__init__(self, *args, **kwargs)`
Initialize the model, define layers, and save hyperparameters.

```python
def __init__(self, learning_rate=1e-3, hidden_dim=128):
    super().__init__()
    self.save_hyperparameters()  # Saves all args to self.hparams
    self.model = nn.Sequential(...)
```

#### `forward(self, x)`
Define the forward pass for inference. Called by `predict_step` by default.

```python
def forward(self, x):
    return self.model(x)
```

#### `training_step(self, batch, batch_idx)`
Define the training loop logic. Return loss for automatic optimization.

```python
def training_step(self, batch, batch_idx):
    x, y = batch
    y_hat = self(x)
    loss = F.cross_entropy(y_hat, y)
    self.log('train_loss', loss)
    return loss
```

#### `validation_step(self, batch, batch_idx)`
Define the validation loop logic. Model automatically in eval mode with no gradients.

```python
def validation_step(self, batch, batch_idx):
    x, y = batch
    y_hat = self(x)
    loss = F.cross_entropy(y_hat, y)
    self.log('val_loss', loss)
    return loss
```

#### `test_step(self, batch, batch_idx)`
Define the test loop logic. Only runs when `trainer.test()` is called.

```python
def test_step(self, batch, batch_idx):
    x, y = batch
    y_hat = self(x)
    loss = F.cross_entropy(y_hat, y)
    self.log('test_loss', loss)
    return loss
```

#### `predict_step(self, batch, batch_idx, dataloader_idx=0)`
Define prediction logic for inference. Defaults to calling `forward()`.

```python
def predict_step(self, batch, batch_idx, dataloader_idx=0):
    x, y = batch
    return self(x)
```

#### `configure_optimizers(self)`
Return optimizer(s) and optional learning rate scheduler(s).

```python
def configure_optimizers(self):
    optimizer = torch.optim.Adam(self.parameters(), lr=self.hparams.learning_rate)
    scheduler = ReduceLROnPlateau(optimizer, mode='min')
    return {
        "optimizer": optimizer,
        "lr_scheduler": {
            "scheduler": scheduler,
            "monitor": "val_loss",
            "interval": "epoch",
            "frequency": 1,
        }
    }
```

### Lifecycle Hooks

#### Epoch-Level Hooks
- `on_train_epoch_start()` - Called at the start of each training epoch
- `on_train_epoch_end()` - Called at the end of each training epoch
- `on_validation_epoch_start()` - Called at the start of validation epoch
- `on_validation_epoch_end()` - Called at the end of validation epoch
- `on_test_epoch_start()` - Called at the start of test epoch
- `on_test_epoch_end()` - Called at the end of test epoch

#### Batch-Level Hooks
- `on_train_batch_start(batch, batch_idx)` - Called before training batch
- `on_train_batch_end(outputs, batch, batch_idx)` - Called after training batch
- `on_validation_batch_start(batch, batch_idx)` - Called before validation batch
- `on_validation_batch_end(outputs, batch, batch_idx)` - Called after validation batch

#### Training Lifecycle
- `on_fit_start()` - Called at the start of fit
- `on_fit_end()` - Called at the end of fit
- `on_train_start()` - Called at the start of training
- `on_train_end()` - Called at the end of training

### Logging

#### `self.log(name, value, **kwargs)`
Log a metric to all configured loggers.

**Common Parameters:**
- `on_step` (bool) - Log at each batch step
- `on_epoch` (bool) - Log at the end of epoch (automatically aggregated)
- `prog_bar` (bool) - Display in progress bar
- `logger` (bool) - Send to logger
- `sync_dist` (bool) - Synchronize across all distributed processes
- `reduce_fx` (str) - Reduction function for distributed ("mean", "sum", etc.)

```python
self.log('train_loss', loss, on_step=True, on_epoch=True, prog_bar=True)
```

#### `self.log_dict(dictionary, **kwargs)`
Log multiple metrics at once.

```python
metrics = {'loss': loss, 'acc': acc, 'f1': f1}
self.log_dict(metrics, on_step=True, on_epoch=True)
```

### Device Management

- `self.device` - Current device (automatically managed)
- `self.to(device)` - Move model to device (usually handled automatically)

**Best Practice:** Create tensors on model's device:
```python
new_tensor = torch.zeros(10, device=self.device)
```

### Hyperparameter Management

#### `self.save_hyperparameters(*args, **kwargs)`
Automatically save init arguments to `self.hparams` and checkpoints.

```python
def __init__(self, learning_rate, hidden_dim):
    super().__init__()
    self.save_hyperparameters()  # Saves all args
    # Access via self.hparams.learning_rate, self.hparams.hidden_dim
```

---

## Trainer

The Trainer automates the training loop and engineering complexity.

### Core Parameters

#### Training Duration
- `max_epochs` (int) - Maximum number of epochs (default: 1000)
- `min_epochs` (int) - Minimum number of epochs
- `max_steps` (int) - Maximum number of optimizer steps
- `min_steps` (int) - Minimum number of optimizer steps
- `max_time` (str/dict) - Maximum training time ("DD:HH:MM:SS" or dict)

#### Hardware Configuration
- `accelerator` (str) - Hardware to use: "cpu", "gpu", "tpu", "auto"
- `devices` (int/list) - Number or specific device IDs: 1, 4, [0,2], "auto"
- `num_nodes` (int) - Number of GPU nodes for distributed training
- `strategy` (str) - Training strategy: "ddp", "fsdp", "deepspeed", etc.

#### Data Management
- `limit_train_batches` (int/float) - Limit training batches (0.0-1.0 for %, int for count)
- `limit_val_batches` (int/float) - Limit validation batches
- `limit_test_batches` (int/float) - Limit test batches
- `limit_predict_batches` (int/float) - Limit prediction batches

#### Validation
- `check_val_every_n_epoch` (int) - Run validation every N epochs
- `val_check_interval` (int/float) - Validate every N batches or fraction
- `num_sanity_val_steps` (int) - Validation steps before training (default: 2)

#### Optimization
- `gradient_clip_val` (float) - Clip gradients by value
- `gradient_clip_algorithm` (str) - "value" or "norm"
- `accumulate_grad_batches` (int) - Accumulate gradients over K batches
- `precision` (str) - Training precision: "32-true", "16-mixed", "bf16-mixed", "64-true"

#### Logging and Checkpointing
- `logger` (Logger/list) - Logger instance(s) or True/False
- `log_every_n_steps` (int) - Logging frequency
- `enable_checkpointing` (bool) - Enable automatic checkpointing
- `callbacks` (list) - List of callback instances
- `default_root_dir` (str) - Default path for logs and checkpoints

#### Debugging
- `fast_dev_run` (bool/int) - Run N batches for quick testing
- `overfit_batches` (int/float) - Overfit on limited data for debugging
- `detect_anomaly` (bool) - Enable PyTorch anomaly detection
- `profiler` (str/Profiler) - Profile training: "simple", "advanced", or custom

#### Performance
- `benchmark` (bool) - Enable cudnn.benchmark for performance
- `deterministic` (bool) - Enable deterministic training for reproducibility
- `sync_batchnorm` (bool) - Synchronize batch norm across GPUs

### Training Methods

#### `trainer.fit(model, train_dataloaders=None, val_dataloaders=None, datamodule=None, ckpt_path=None)`
Run the full training routine.

```python
trainer.fit(model, train_loader, val_loader)
# Or with DataModule
trainer.fit(model, datamodule=dm)
# Resume from checkpoint
trainer.fit(model, train_loader, val_loader, ckpt_path="path/to/checkpoint.ckpt")
```

#### `trainer.validate(model, dataloaders=None, datamodule=None, ckpt_path=None)`
Run validation independently.

```python
trainer.validate(model, val_loader)
```

#### `trainer.test(model, dataloaders=None, datamodule=None, ckpt_path=None)`
Run test evaluation.

```python
trainer.test(model, test_loader)
# Or load from checkpoint
trainer.test(ckpt_path="best_model.ckpt", datamodule=dm)
```

#### `trainer.predict(model, dataloaders=None, datamodule=None, ckpt_path=None)`
Run inference predictions.

```python
predictions = trainer.predict(model, predict_loader)
```

---

## LightningDataModule

Encapsulates all data processing logic in a reusable class.

### Core Methods

#### `prepare_data(self)`
Download and prepare data (called once on single process).
Do NOT set state here (no self.x = y).

```python
def prepare_data(self):
    # Download datasets
    datasets.MNIST(self.data_dir, train=True, download=True)
    datasets.MNIST(self.data_dir, train=False, download=True)
```

#### `setup(self, stage=None)`
Load data and create splits (called on every process/GPU).
Setting state is OK here.

**stage parameter:** "fit", "validate", "test", or "predict"

```python
def setup(self, stage=None):
    if stage == "fit" or stage is None:
        full_dataset = datasets.MNIST(self.data_dir, train=True)
        self.train_dataset, self.val_dataset = random_split(full_dataset, [55000, 5000])

    if stage == "test" or stage is None:
        self.test_dataset = datasets.MNIST(self.data_dir, train=False)
```

#### DataLoader Methods
- `train_dataloader(self)` - Return training DataLoader
- `val_dataloader(self)` - Return validation DataLoader
- `test_dataloader(self)` - Return test DataLoader
- `predict_dataloader(self)` - Return prediction DataLoader

```python
def train_dataloader(self):
    return DataLoader(self.train_dataset, batch_size=32, shuffle=True)
```

### Optional Methods
- `teardown(stage=None)` - Cleanup after training/testing
- `state_dict()` - Save state for checkpointing
- `load_state_dict(state_dict)` - Load state from checkpoint

---

## Callbacks

Extend training with modular, reusable functionality.

### Built-in Callbacks

#### ModelCheckpoint
Save model checkpoints based on monitored metrics.

```python
from lightning.pytorch.callbacks import ModelCheckpoint

checkpoint_callback = ModelCheckpoint(
    dirpath='checkpoints/',
    filename='{epoch}-{val_loss:.2f}',
    monitor='val_loss',
    mode='min',
    save_top_k=3,
    save_last=True,
    verbose=True,
)
```

**Key Parameters:**
- `monitor` - Metric to monitor
- `mode` - "min" or "max"
- `save_top_k` - Save top K models
- `save_last` - Always save last checkpoint
- `every_n_epochs` - Save every N epochs

#### EarlyStopping
Stop training when metric stops improving.

```python
from lightning.pytorch.callbacks import EarlyStopping

early_stop = EarlyStopping(
    monitor='val_loss',
    patience=10,
    mode='min',
    verbose=True,
)
```

#### LearningRateMonitor
Log learning rate values.

```python
from lightning.pytorch.callbacks import LearningRateMonitor

lr_monitor = LearningRateMonitor(logging_interval='epoch')
```

#### RichProgressBar
Display rich progress bar with metrics.

```python
from lightning.pytorch.callbacks import RichProgressBar

progress_bar = RichProgressBar()
```

### Custom Callbacks

Create custom callbacks by inheriting from `Callback`.

```python
from lightning.pytorch.callbacks import Callback

class MyCallback(Callback):
    def on_train_start(self, trainer, pl_module):
        print("Training starting!")

    def on_train_epoch_end(self, trainer, pl_module):
        print(f"Epoch {trainer.current_epoch} ended")

    def on_validation_end(self, trainer, pl_module):
        val_loss = trainer.callback_metrics.get('val_loss')
        print(f"Validation loss: {val_loss}")
```

**Common Hooks:**
- `on_train_start/end`
- `on_train_epoch_start/end`
- `on_validation_epoch_start/end`
- `on_test_epoch_start/end`
- `on_before_backward/on_after_backward`
- `on_before_optimizer_step`

---

## Loggers

Track experiments with various logging frameworks.

### TensorBoardLogger
```python
from lightning.pytorch.loggers import TensorBoardLogger

logger = TensorBoardLogger(save_dir='logs/', name='my_experiment')
trainer = Trainer(logger=logger)
```

### WandbLogger
```python
from lightning.pytorch.loggers import WandbLogger

logger = WandbLogger(project='my_project', name='experiment_1')
trainer = Trainer(logger=logger)
```

### MLFlowLogger
```python
from lightning.pytorch.loggers import MLFlowLogger

logger = MLFlowLogger(experiment_name='my_exp', tracking_uri='file:./ml-runs')
trainer = Trainer(logger=logger)
```

### CSVLogger
```python
from lightning.pytorch.loggers import CSVLogger

logger = CSVLogger(save_dir='logs/', name='my_experiment')
trainer = Trainer(logger=logger)
```

### Multiple Loggers
```python
loggers = [
    TensorBoardLogger('logs/'),
    CSVLogger('logs/'),
]
trainer = Trainer(logger=loggers)
```

---

## Common Patterns

### Reproducibility
```python
from lightning.pytorch import seed_everything

seed_everything(42, workers=True)
trainer = Trainer(deterministic=True)
```

### Mixed Precision Training
```python
trainer = Trainer(precision='16-mixed')  # or 'bf16-mixed'
```

### Multi-GPU Training
```python
# Data parallel (DDP)
trainer = Trainer(accelerator='gpu', devices=4, strategy='ddp')

# Model parallel (FSDP)
trainer = Trainer(accelerator='gpu', devices=4, strategy='fsdp')
```

### Gradient Accumulation
```python
trainer = Trainer(accumulate_grad_batches=4)  # Effective batch size = 4x
```

### Learning Rate Finding
```python
from lightning.pytorch.tuner import Tuner

trainer = Trainer()
tuner = Tuner(trainer)
lr_finder = tuner.lr_find(model, train_dataloader)
model.hparams.learning_rate = lr_finder.suggestion()
```

### Loading from Checkpoint
```python
# Load model
model = MyLightningModule.load_from_checkpoint('checkpoint.ckpt')

# Resume training
trainer.fit(model, ckpt_path='checkpoint.ckpt')
```
