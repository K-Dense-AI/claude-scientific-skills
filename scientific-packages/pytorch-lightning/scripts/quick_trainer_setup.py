"""
Helper script to quickly set up a PyTorch Lightning Trainer with common configurations.

This script provides preset configurations for different training scenarios
and makes it easy to create a Trainer with best practices.
"""

import lightning as L
from lightning.pytorch.callbacks import (
    ModelCheckpoint,
    EarlyStopping,
    LearningRateMonitor,
    RichProgressBar,
    ModelSummary,
)
from lightning.pytorch.loggers import TensorBoardLogger, CSVLogger


def create_trainer(
    preset: str = "default",
    max_epochs: int = 100,
    accelerator: str = "auto",
    devices: int = 1,
    log_dir: str = "./logs",
    experiment_name: str = "lightning_experiment",
    enable_checkpointing: bool = True,
    enable_early_stopping: bool = True,
    **kwargs
):
    """
    Create a Lightning Trainer with preset configurations.

    Args:
        preset: Configuration preset - "default", "fast_dev", "production", "distributed"
        max_epochs: Maximum number of training epochs
        accelerator: Device to use ("auto", "gpu", "cpu", "tpu")
        devices: Number of devices to use
        log_dir: Directory for logs and checkpoints
        experiment_name: Name for the experiment
        enable_checkpointing: Whether to enable model checkpointing
        enable_early_stopping: Whether to enable early stopping
        **kwargs: Additional arguments to pass to Trainer

    Returns:
        Configured Lightning Trainer instance
    """

    callbacks = []
    logger_list = []

    # Configure based on preset
    if preset == "fast_dev":
        # Fast development run - minimal epochs, quick debugging
        config = {
            "fast_dev_run": False,
            "max_epochs": 3,
            "limit_train_batches": 100,
            "limit_val_batches": 50,
            "log_every_n_steps": 10,
            "enable_progress_bar": True,
            "enable_model_summary": True,
        }

    elif preset == "production":
        # Production-ready configuration with all bells and whistles
        config = {
            "max_epochs": max_epochs,
            "precision": "16-mixed",
            "gradient_clip_val": 1.0,
            "log_every_n_steps": 50,
            "enable_progress_bar": True,
            "enable_model_summary": True,
            "deterministic": True,
            "benchmark": True,
        }

        # Add model checkpointing
        if enable_checkpointing:
            callbacks.append(
                ModelCheckpoint(
                    dirpath=f"{log_dir}/{experiment_name}/checkpoints",
                    filename="{epoch}-{val_loss:.2f}",
                    monitor="val_loss",
                    mode="min",
                    save_top_k=3,
                    save_last=True,
                    verbose=True,
                )
            )

        # Add early stopping
        if enable_early_stopping:
            callbacks.append(
                EarlyStopping(
                    monitor="val_loss",
                    patience=10,
                    mode="min",
                    verbose=True,
                )
            )

        # Add learning rate monitor
        callbacks.append(LearningRateMonitor(logging_interval="epoch"))

        # Add TensorBoard logger
        logger_list.append(
            TensorBoardLogger(
                save_dir=log_dir,
                name=experiment_name,
                version=None,
            )
        )

    elif preset == "distributed":
        # Distributed training configuration
        config = {
            "max_epochs": max_epochs,
            "strategy": "ddp",
            "precision": "16-mixed",
            "sync_batchnorm": True,
            "use_distributed_sampler": True,
            "log_every_n_steps": 50,
            "enable_progress_bar": True,
        }

        # Add model checkpointing
        if enable_checkpointing:
            callbacks.append(
                ModelCheckpoint(
                    dirpath=f"{log_dir}/{experiment_name}/checkpoints",
                    filename="{epoch}-{val_loss:.2f}",
                    monitor="val_loss",
                    mode="min",
                    save_top_k=3,
                    save_last=True,
                )
            )

    else:  # default
        # Default configuration - balanced for most use cases
        config = {
            "max_epochs": max_epochs,
            "log_every_n_steps": 50,
            "enable_progress_bar": True,
            "enable_model_summary": True,
        }

        # Add basic checkpointing
        if enable_checkpointing:
            callbacks.append(
                ModelCheckpoint(
                    dirpath=f"{log_dir}/{experiment_name}/checkpoints",
                    filename="{epoch}-{val_loss:.2f}",
                    monitor="val_loss",
                    save_last=True,
                )
            )

        # Add CSV logger
        logger_list.append(
            CSVLogger(
                save_dir=log_dir,
                name=experiment_name,
            )
        )

    # Add progress bar
    if config.get("enable_progress_bar", True):
        callbacks.append(RichProgressBar())

    # Merge with provided kwargs
    final_config = {
        **config,
        "accelerator": accelerator,
        "devices": devices,
        "callbacks": callbacks,
        "logger": logger_list if logger_list else True,
        **kwargs,
    }

    # Create and return trainer
    return L.Trainer(**final_config)


def create_debugging_trainer():
    """Create a trainer optimized for debugging."""
    return create_trainer(
        preset="fast_dev",
        max_epochs=1,
        limit_train_batches=10,
        limit_val_batches=5,
        num_sanity_val_steps=2,
    )


def create_gpu_trainer(num_gpus: int = 1, precision: str = "16-mixed"):
    """Create a trainer optimized for GPU training."""
    return create_trainer(
        preset="production",
        accelerator="gpu",
        devices=num_gpus,
        precision=precision,
    )


def create_distributed_trainer(num_gpus: int = 2, num_nodes: int = 1):
    """Create a trainer for distributed training across multiple GPUs."""
    return create_trainer(
        preset="distributed",
        accelerator="gpu",
        devices=num_gpus,
        num_nodes=num_nodes,
        strategy="ddp",
    )


# Example usage
if __name__ == "__main__":
    print("Creating different trainer configurations...\n")

    # 1. Default trainer
    print("1. Default trainer:")
    trainer_default = create_trainer(preset="default", max_epochs=50)
    print(f"   Max epochs: {trainer_default.max_epochs}")
    print(f"   Accelerator: {trainer_default.accelerator}")
    print(f"   Callbacks: {len(trainer_default.callbacks)}")
    print()

    # 2. Fast development trainer
    print("2. Fast development trainer:")
    trainer_dev = create_trainer(preset="fast_dev")
    print(f"   Max epochs: {trainer_dev.max_epochs}")
    print(f"   Train batches limit: {trainer_dev.limit_train_batches}")
    print()

    # 3. Production trainer
    print("3. Production trainer:")
    trainer_prod = create_trainer(
        preset="production",
        max_epochs=100,
        experiment_name="my_experiment"
    )
    print(f"   Max epochs: {trainer_prod.max_epochs}")
    print(f"   Precision: {trainer_prod.precision}")
    print(f"   Callbacks: {len(trainer_prod.callbacks)}")
    print()

    # 4. Debugging trainer
    print("4. Debugging trainer:")
    trainer_debug = create_debugging_trainer()
    print(f"   Max epochs: {trainer_debug.max_epochs}")
    print(f"   Train batches: {trainer_debug.limit_train_batches}")
    print()

    # 5. GPU trainer
    print("5. GPU trainer:")
    trainer_gpu = create_gpu_trainer(num_gpus=1)
    print(f"   Accelerator: {trainer_gpu.accelerator}")
    print(f"   Precision: {trainer_gpu.precision}")
    print()

    print("All trainer configurations created successfully!")
