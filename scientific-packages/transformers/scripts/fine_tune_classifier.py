#!/usr/bin/env python3
"""
Fine-tune a transformer model for text classification.

This script demonstrates the complete workflow for fine-tuning a pre-trained
model on a classification task using the Trainer API.
"""

import numpy as np
from datasets import load_dataset
from transformers import (
    AutoTokenizer,
    AutoModelForSequenceClassification,
    TrainingArguments,
    Trainer,
    DataCollatorWithPadding,
)
import evaluate


def load_and_prepare_data(dataset_name="imdb", model_name="distilbert-base-uncased", max_samples=None):
    """
    Load dataset and tokenize.

    Args:
        dataset_name: Name of the dataset to load
        model_name: Name of the model/tokenizer to use
        max_samples: Limit number of samples (for quick testing)

    Returns:
        tokenized_datasets, tokenizer
    """
    print(f"Loading dataset: {dataset_name}")
    dataset = load_dataset(dataset_name)

    # Optionally limit samples for quick testing
    if max_samples:
        dataset["train"] = dataset["train"].select(range(max_samples))
        dataset["test"] = dataset["test"].select(range(min(max_samples, len(dataset["test"]))))

    print(f"Loading tokenizer: {model_name}")
    tokenizer = AutoTokenizer.from_pretrained(model_name)

    def tokenize_function(examples):
        return tokenizer(
            examples["text"],
            padding="max_length",
            truncation=True,
            max_length=512
        )

    print("Tokenizing dataset...")
    tokenized_datasets = dataset.map(tokenize_function, batched=True)

    return tokenized_datasets, tokenizer


def create_model(model_name, num_labels, id2label, label2id):
    """
    Create classification model.

    Args:
        model_name: Name of the pre-trained model
        num_labels: Number of classification labels
        id2label: Dictionary mapping label IDs to names
        label2id: Dictionary mapping label names to IDs

    Returns:
        model
    """
    print(f"Loading model: {model_name}")
    model = AutoModelForSequenceClassification.from_pretrained(
        model_name,
        num_labels=num_labels,
        id2label=id2label,
        label2id=label2id
    )
    return model


def define_compute_metrics(metric_name="accuracy"):
    """
    Define function to compute metrics during evaluation.

    Args:
        metric_name: Name of the metric to use

    Returns:
        compute_metrics function
    """
    metric = evaluate.load(metric_name)

    def compute_metrics(eval_pred):
        logits, labels = eval_pred
        predictions = np.argmax(logits, axis=-1)
        return metric.compute(predictions=predictions, references=labels)

    return compute_metrics


def train_model(model, tokenizer, train_dataset, eval_dataset, output_dir="./results"):
    """
    Train the model.

    Args:
        model: The model to train
        tokenizer: The tokenizer
        train_dataset: Training dataset
        eval_dataset: Evaluation dataset
        output_dir: Directory for checkpoints and logs

    Returns:
        trained model, trainer
    """
    # Define training arguments
    training_args = TrainingArguments(
        output_dir=output_dir,
        num_train_epochs=3,
        per_device_train_batch_size=16,
        per_device_eval_batch_size=64,
        learning_rate=2e-5,
        weight_decay=0.01,
        eval_strategy="epoch",
        save_strategy="epoch",
        load_best_model_at_end=True,
        metric_for_best_model="accuracy",
        logging_dir=f"{output_dir}/logs",
        logging_steps=100,
        save_total_limit=2,
        fp16=False,  # Set to True if using GPU with fp16 support
    )

    # Create data collator
    data_collator = DataCollatorWithPadding(tokenizer=tokenizer)

    # Create trainer
    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=train_dataset,
        eval_dataset=eval_dataset,
        data_collator=data_collator,
        compute_metrics=define_compute_metrics("accuracy"),
    )

    # Train
    print("\nStarting training...")
    trainer.train()

    # Evaluate
    print("\nEvaluating model...")
    eval_results = trainer.evaluate()
    print(f"Evaluation results: {eval_results}")

    return model, trainer


def test_inference(model, tokenizer, id2label):
    """
    Test the trained model with sample texts.

    Args:
        model: Trained model
        tokenizer: Tokenizer
        id2label: Dictionary mapping label IDs to names
    """
    print("\n=== Testing Inference ===")

    test_texts = [
        "This movie was absolutely fantastic! I loved every minute of it.",
        "Terrible film. Waste of time and money.",
        "It was okay, nothing special but not bad either."
    ]

    for text in test_texts:
        inputs = tokenizer(text, return_tensors="pt", truncation=True, max_length=512)
        outputs = model(**inputs)
        predictions = outputs.logits.argmax(-1)
        predicted_label = id2label[predictions.item()]
        confidence = outputs.logits.softmax(-1).max().item()

        print(f"\nText: {text}")
        print(f"Prediction: {predicted_label} (confidence: {confidence:.3f})")


def main():
    """Main training pipeline."""
    # Configuration
    DATASET_NAME = "imdb"
    MODEL_NAME = "distilbert-base-uncased"
    OUTPUT_DIR = "./results"
    MAX_SAMPLES = None  # Set to a small number (e.g., 1000) for quick testing

    # Label mapping
    id2label = {0: "negative", 1: "positive"}
    label2id = {"negative": 0, "positive": 1}
    num_labels = len(id2label)

    print("=" * 60)
    print("Fine-Tuning Text Classification Model")
    print("=" * 60)

    # Load and prepare data
    tokenized_datasets, tokenizer = load_and_prepare_data(
        dataset_name=DATASET_NAME,
        model_name=MODEL_NAME,
        max_samples=MAX_SAMPLES
    )

    # Create model
    model = create_model(
        model_name=MODEL_NAME,
        num_labels=num_labels,
        id2label=id2label,
        label2id=label2id
    )

    # Train model
    model, trainer = train_model(
        model=model,
        tokenizer=tokenizer,
        train_dataset=tokenized_datasets["train"],
        eval_dataset=tokenized_datasets["test"],
        output_dir=OUTPUT_DIR
    )

    # Save final model
    print(f"\nSaving model to {OUTPUT_DIR}/final_model")
    trainer.save_model(f"{OUTPUT_DIR}/final_model")
    tokenizer.save_pretrained(f"{OUTPUT_DIR}/final_model")

    # Test inference
    test_inference(model, tokenizer, id2label)

    print("\n" + "=" * 60)
    print("Training completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
