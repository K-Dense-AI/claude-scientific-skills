"""
Template for creating a PyTorch Lightning LightningModule.

This template includes all common hooks and patterns for building
a Lightning model with best practices.
"""

import lightning as L
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.optim import Adam, SGD
from torch.optim.lr_scheduler import ReduceLROnPlateau, StepLR


class TemplateLightningModule(L.LightningModule):
    """Template LightningModule with all common hooks and patterns."""

    def __init__(
        self,
        # Model architecture parameters
        input_dim: int = 784,
        hidden_dim: int = 128,
        output_dim: int = 10,
        # Optimization parameters
        learning_rate: float = 1e-3,
        optimizer_type: str = "adam",
        scheduler_type: str = None,
        # Other hyperparameters
        dropout: float = 0.1,
    ):
        super().__init__()

        # Save hyperparameters for checkpointing and logging
        self.save_hyperparameters()

        # Define model architecture
        self.model = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, output_dim)
        )

        # Define loss function
        self.criterion = nn.CrossEntropyLoss()

        # For tracking validation outputs (optional)
        self.validation_step_outputs = []

    def forward(self, x):
        """Forward pass for inference."""
        return self.model(x)

    def training_step(self, batch, batch_idx):
        """Training step - called for each training batch."""
        x, y = batch

        # Forward pass
        logits = self(x)
        loss = self.criterion(logits, y)

        # Calculate accuracy
        preds = torch.argmax(logits, dim=1)
        acc = (preds == y).float().mean()

        # Log metrics
        self.log("train_loss", loss, on_step=True, on_epoch=True, prog_bar=True)
        self.log("train_acc", acc, on_step=True, on_epoch=True, prog_bar=True)

        return loss

    def validation_step(self, batch, batch_idx):
        """Validation step - called for each validation batch."""
        x, y = batch

        # Forward pass (model automatically in eval mode)
        logits = self(x)
        loss = self.criterion(logits, y)

        # Calculate accuracy
        preds = torch.argmax(logits, dim=1)
        acc = (preds == y).float().mean()

        # Log metrics
        self.log("val_loss", loss, on_step=False, on_epoch=True, prog_bar=True)
        self.log("val_acc", acc, on_step=False, on_epoch=True, prog_bar=True)

        # Optional: store outputs for epoch-level processing
        self.validation_step_outputs.append({"loss": loss, "acc": acc})

        return loss

    def on_validation_epoch_end(self):
        """Called at the end of validation epoch."""
        # Optional: process all validation outputs
        if self.validation_step_outputs:
            avg_loss = torch.stack([x["loss"] for x in self.validation_step_outputs]).mean()
            avg_acc = torch.stack([x["acc"] for x in self.validation_step_outputs]).mean()

            # Log epoch-level metrics if needed
            # self.log("val_epoch_loss", avg_loss)
            # self.log("val_epoch_acc", avg_acc)

            # Clear outputs
            self.validation_step_outputs.clear()

    def test_step(self, batch, batch_idx):
        """Test step - called for each test batch."""
        x, y = batch

        # Forward pass
        logits = self(x)
        loss = self.criterion(logits, y)

        # Calculate accuracy
        preds = torch.argmax(logits, dim=1)
        acc = (preds == y).float().mean()

        # Log metrics
        self.log("test_loss", loss, on_step=False, on_epoch=True)
        self.log("test_acc", acc, on_step=False, on_epoch=True)

        return loss

    def predict_step(self, batch, batch_idx, dataloader_idx=0):
        """Prediction step - called for each prediction batch."""
        x, y = batch
        logits = self(x)
        preds = torch.argmax(logits, dim=1)
        return preds

    def configure_optimizers(self):
        """Configure optimizer and learning rate scheduler."""
        # Create optimizer
        if self.hparams.optimizer_type.lower() == "adam":
            optimizer = Adam(self.parameters(), lr=self.hparams.learning_rate)
        elif self.hparams.optimizer_type.lower() == "sgd":
            optimizer = SGD(self.parameters(), lr=self.hparams.learning_rate, momentum=0.9)
        else:
            raise ValueError(f"Unknown optimizer: {self.hparams.optimizer_type}")

        # Configure with scheduler if specified
        if self.hparams.scheduler_type:
            if self.hparams.scheduler_type.lower() == "reduce_on_plateau":
                scheduler = ReduceLROnPlateau(optimizer, mode="min", factor=0.5, patience=5)
                return {
                    "optimizer": optimizer,
                    "lr_scheduler": {
                        "scheduler": scheduler,
                        "monitor": "val_loss",
                        "interval": "epoch",
                        "frequency": 1,
                    }
                }
            elif self.hparams.scheduler_type.lower() == "step":
                scheduler = StepLR(optimizer, step_size=10, gamma=0.1)
                return {
                    "optimizer": optimizer,
                    "lr_scheduler": {
                        "scheduler": scheduler,
                        "interval": "epoch",
                        "frequency": 1,
                    }
                }

        return optimizer

    # Optional: Additional hooks for custom behavior

    def on_train_start(self):
        """Called at the beginning of training."""
        pass

    def on_train_epoch_start(self):
        """Called at the beginning of each training epoch."""
        pass

    def on_train_epoch_end(self):
        """Called at the end of each training epoch."""
        pass

    def on_train_end(self):
        """Called at the end of training."""
        pass


# Example usage
if __name__ == "__main__":
    # Create model
    model = TemplateLightningModule(
        input_dim=784,
        hidden_dim=128,
        output_dim=10,
        learning_rate=1e-3,
        optimizer_type="adam",
        scheduler_type="reduce_on_plateau"
    )

    # Create trainer
    trainer = L.Trainer(
        max_epochs=10,
        accelerator="auto",
        devices=1,
        log_every_n_steps=50,
    )

    # Note: You would need to provide dataloaders
    # trainer.fit(model, train_dataloader, val_dataloader)

    print("Template LightningModule created successfully!")
    print(f"Model hyperparameters: {model.hparams}")
