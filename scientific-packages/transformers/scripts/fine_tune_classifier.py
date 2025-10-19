#!/usr/bin/env python3
"""
Complete example for fine-tuning a text classification model.

This script demonstrates the full workflow:
1. Load dataset
2. Preprocess with tokenizer
3. Configure model
4. Train with Trainer
5. Evaluate and save

Usage:
    python fine_tune_classifier.py --model bert-base-uncased --dataset imdb --epochs 3
"""

import argparse
from datasets import load_dataset
from transformers import (
    AutoTokenizer,
    AutoModelForSequenceClassification,
    TrainingArguments,
    Trainer,
    DataCollatorWithPadding,
)
import evaluate
import numpy as np


def compute_metrics(eval_pred):
    """Compute accuracy and F1 score."""
    metric_accuracy = evaluate.load("accuracy")
    metric_f1 = evaluate.load("f1")

    logits, labels = eval_pred
    predictions = np.argmax(logits, axis=-1)

    accuracy = metric_accuracy.compute(predictions=predictions, references=labels)
    f1 = metric_f1.compute(predictions=predictions, references=labels)

    return {"accuracy": accuracy["accuracy"], "f1": f1["f1"]}


def main():
    parser = argparse.ArgumentParser(description="Fine-tune a text classification model")
    parser.add_argument(
        "--model",
        type=str,
        default="bert-base-uncased",
        help="Pretrained model name or path",
    )
    parser.add_argument(
        "--dataset",
        type=str,
        default="imdb",
        help="Dataset name from Hugging Face Hub",
    )
    parser.add_argument(
        "--max-samples",
        type=int,
        default=None,
        help="Maximum samples to use (for quick testing)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="./results",
        help="Output directory for checkpoints",
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=3,
        help="Number of training epochs",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=16,
        help="Batch size per device",
    )
    parser.add_argument(
        "--learning-rate",
        type=float,
        default=2e-5,
        help="Learning rate",
    )
    parser.add_argument(
        "--push-to-hub",
        action="store_true",
        help="Push model to Hugging Face Hub after training",
    )

    args = parser.parse_args()

    print("=" * 60)
    print("Text Classification Fine-Tuning")
    print("=" * 60)
    print(f"Model: {args.model}")
    print(f"Dataset: {args.dataset}")
    print(f"Epochs: {args.epochs}")
    print(f"Batch size: {args.batch_size}")
    print(f"Learning rate: {args.learning_rate}")
    print("=" * 60)

    # 1. Load dataset
    print("\n[1/5] Loading dataset...")
    dataset = load_dataset(args.dataset)

    if args.max_samples:
        dataset["train"] = dataset["train"].select(range(args.max_samples))
        dataset["test"] = dataset["test"].select(range(args.max_samples // 5))

    print(f"Train samples: {len(dataset['train'])}")
    print(f"Test samples: {len(dataset['test'])}")

    # 2. Preprocess
    print("\n[2/5] Preprocessing data...")
    tokenizer = AutoTokenizer.from_pretrained(args.model)

    def preprocess_function(examples):
        return tokenizer(examples["text"], truncation=True, max_length=512)

    tokenized_dataset = dataset.map(preprocess_function, batched=True)
    data_collator = DataCollatorWithPadding(tokenizer=tokenizer)

    # 3. Load model
    print("\n[3/5] Loading model...")

    # Determine number of labels
    num_labels = len(set(dataset["train"]["label"]))

    model = AutoModelForSequenceClassification.from_pretrained(
        args.model,
        num_labels=num_labels,
    )

    print(f"Number of labels: {num_labels}")
    print(f"Model parameters: {model.num_parameters():,}")

    # 4. Configure training
    print("\n[4/5] Configuring training...")
    training_args = TrainingArguments(
        output_dir=args.output_dir,
        learning_rate=args.learning_rate,
        per_device_train_batch_size=args.batch_size,
        per_device_eval_batch_size=args.batch_size,
        num_train_epochs=args.epochs,
        weight_decay=0.01,
        eval_strategy="epoch",
        save_strategy="epoch",
        load_best_model_at_end=True,
        push_to_hub=args.push_to_hub,
        logging_steps=100,
    )

    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=tokenized_dataset["train"],
        eval_dataset=tokenized_dataset["test"],
        tokenizer=tokenizer,
        data_collator=data_collator,
        compute_metrics=compute_metrics,
    )

    # 5. Train
    print("\n[5/5] Training...")
    print("-" * 60)
    trainer.train()

    # Evaluate
    print("\n" + "=" * 60)
    print("Final Evaluation")
    print("=" * 60)
    metrics = trainer.evaluate()

    print(f"Accuracy: {metrics['eval_accuracy']:.4f}")
    print(f"F1 Score: {metrics['eval_f1']:.4f}")
    print(f"Loss: {metrics['eval_loss']:.4f}")

    # Save
    print("\n" + "=" * 60)
    print(f"Saving model to {args.output_dir}")
    trainer.save_model(args.output_dir)
    tokenizer.save_pretrained(args.output_dir)

    if args.push_to_hub:
        print("Pushing to Hugging Face Hub...")
        trainer.push_to_hub()

    print("=" * 60)
    print("Training complete!")
    print("=" * 60)

    # Quick inference example
    print("\nQuick inference example:")
    from transformers import pipeline

    classifier = pipeline(
        "text-classification",
        model=args.output_dir,
        tokenizer=args.output_dir,
    )

    example_text = "This is a great example of how to use transformers!"
    result = classifier(example_text)
    print(f"Text: {example_text}")
    print(f"Prediction: {result[0]['label']} (score: {result[0]['score']:.4f})")


if __name__ == "__main__":
    main()
