"""
Template for creating a PyTorch Lightning DataModule.

This template includes all common hooks and patterns for organizing
data processing workflows with best practices.
"""

import lightning as L
from torch.utils.data import DataLoader, Dataset, random_split
import torch


class TemplateDataset(Dataset):
    """Example dataset - replace with your actual dataset."""

    def __init__(self, data, targets, transform=None):
        self.data = data
        self.targets = targets
        self.transform = transform

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        x = self.data[idx]
        y = self.targets[idx]

        if self.transform:
            x = self.transform(x)

        return x, y


class TemplateDataModule(L.LightningDataModule):
    """Template DataModule with all common hooks and patterns."""

    def __init__(
        self,
        data_dir: str = "./data",
        batch_size: int = 32,
        num_workers: int = 4,
        train_val_split: tuple = (0.8, 0.2),
        seed: int = 42,
        pin_memory: bool = True,
        persistent_workers: bool = True,
    ):
        super().__init__()

        # Save hyperparameters
        self.save_hyperparameters()

        # Initialize attributes
        self.data_dir = data_dir
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.train_val_split = train_val_split
        self.seed = seed
        self.pin_memory = pin_memory
        self.persistent_workers = persistent_workers

        # Placeholders for datasets
        self.train_dataset = None
        self.val_dataset = None
        self.test_dataset = None
        self.predict_dataset = None

        # Placeholder for transforms
        self.train_transform = None
        self.val_transform = None
        self.test_transform = None

    def prepare_data(self):
        """
        Download and prepare data (called only on 1 GPU/TPU in distributed settings).
        Use this for downloading, tokenizing, etc. Do NOT set state here (no self.x = y).
        """
        # Example: Download datasets
        # datasets.MNIST(self.data_dir, train=True, download=True)
        # datasets.MNIST(self.data_dir, train=False, download=True)
        pass

    def setup(self, stage: str = None):
        """
        Load data and create train/val/test splits (called on every GPU/TPU in distributed).
        Use this for splitting, creating datasets, etc. Setting state is OK here (self.x = y).

        Args:
            stage: Either 'fit', 'validate', 'test', or 'predict'
        """

        # Fit stage: setup training and validation datasets
        if stage == "fit" or stage is None:
            # Load full dataset
            # Example: full_dataset = datasets.MNIST(self.data_dir, train=True, transform=self.train_transform)

            # Create dummy data for template
            full_data = torch.randn(1000, 784)
            full_targets = torch.randint(0, 10, (1000,))
            full_dataset = TemplateDataset(full_data, full_targets, transform=self.train_transform)

            # Split into train and validation
            train_size = int(len(full_dataset) * self.train_val_split[0])
            val_size = len(full_dataset) - train_size

            self.train_dataset, self.val_dataset = random_split(
                full_dataset,
                [train_size, val_size],
                generator=torch.Generator().manual_seed(self.seed)
            )

            # Apply validation transform if different from train
            if self.val_transform:
                self.val_dataset.dataset.transform = self.val_transform

        # Test stage: setup test dataset
        if stage == "test" or stage is None:
            # Example: self.test_dataset = datasets.MNIST(
            #     self.data_dir, train=False, transform=self.test_transform
            # )

            # Create dummy test data for template
            test_data = torch.randn(200, 784)
            test_targets = torch.randint(0, 10, (200,))
            self.test_dataset = TemplateDataset(test_data, test_targets, transform=self.test_transform)

        # Predict stage: setup prediction dataset
        if stage == "predict" or stage is None:
            # Example: self.predict_dataset = YourCustomDataset(...)

            # Create dummy predict data for template
            predict_data = torch.randn(100, 784)
            predict_targets = torch.zeros(100, dtype=torch.long)
            self.predict_dataset = TemplateDataset(predict_data, predict_targets)

    def train_dataloader(self):
        """Return training dataloader."""
        return DataLoader(
            self.train_dataset,
            batch_size=self.batch_size,
            shuffle=True,
            num_workers=self.num_workers,
            pin_memory=self.pin_memory,
            persistent_workers=self.persistent_workers if self.num_workers > 0 else False,
        )

    def val_dataloader(self):
        """Return validation dataloader."""
        return DataLoader(
            self.val_dataset,
            batch_size=self.batch_size,
            shuffle=False,
            num_workers=self.num_workers,
            pin_memory=self.pin_memory,
            persistent_workers=self.persistent_workers if self.num_workers > 0 else False,
        )

    def test_dataloader(self):
        """Return test dataloader."""
        return DataLoader(
            self.test_dataset,
            batch_size=self.batch_size,
            shuffle=False,
            num_workers=self.num_workers,
            pin_memory=self.pin_memory,
            persistent_workers=self.persistent_workers if self.num_workers > 0 else False,
        )

    def predict_dataloader(self):
        """Return prediction dataloader."""
        return DataLoader(
            self.predict_dataset,
            batch_size=self.batch_size,
            shuffle=False,
            num_workers=self.num_workers,
            pin_memory=self.pin_memory,
            persistent_workers=self.persistent_workers if self.num_workers > 0 else False,
        )

    def teardown(self, stage: str = None):
        """Clean up after fit, validate, test, or predict."""
        # Example: close database connections, clear caches, etc.
        pass

    def state_dict(self):
        """Save state for checkpointing."""
        # Return anything you want to save in the checkpoint
        return {}

    def load_state_dict(self, state_dict):
        """Load state from checkpoint."""
        # Restore state from checkpoint
        pass


# Example usage
if __name__ == "__main__":
    # Create datamodule
    datamodule = TemplateDataModule(
        data_dir="./data",
        batch_size=32,
        num_workers=4,
        train_val_split=(0.8, 0.2),
    )

    # Prepare and setup data
    datamodule.prepare_data()
    datamodule.setup("fit")

    # Get dataloaders
    train_loader = datamodule.train_dataloader()
    val_loader = datamodule.val_dataloader()

    print("Template DataModule created successfully!")
    print(f"Train batches: {len(train_loader)}")
    print(f"Val batches: {len(val_loader)}")
    print(f"Batch size: {datamodule.batch_size}")

    # Test a batch
    batch = next(iter(train_loader))
    x, y = batch
    print(f"Batch shape: {x.shape}, {y.shape}")
