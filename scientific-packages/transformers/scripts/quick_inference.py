#!/usr/bin/env python3
"""
Quick inference script using Transformers pipelines.

This script demonstrates how to use various pipeline tasks for quick inference
without manually managing models, tokenizers, or preprocessing.

Usage:
    python quick_inference.py --task text-generation --model gpt2 --input "Hello world"
    python quick_inference.py --task sentiment-analysis --input "I love this!"
"""

import argparse
from transformers import pipeline, infer_device


def main():
    parser = argparse.ArgumentParser(description="Quick inference with Transformers pipelines")
    parser.add_argument(
        "--task",
        type=str,
        required=True,
        help="Pipeline task (text-generation, sentiment-analysis, question-answering, etc.)",
    )
    parser.add_argument(
        "--model",
        type=str,
        default=None,
        help="Model name or path (default: use task default)",
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Input text for inference",
    )
    parser.add_argument(
        "--context",
        type=str,
        default=None,
        help="Context for question-answering tasks",
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=50,
        help="Maximum generation length",
    )
    parser.add_argument(
        "--device",
        type=str,
        default=None,
        help="Device (cuda, cpu, or auto-detect)",
    )

    args = parser.parse_args()

    # Auto-detect device if not specified
    if args.device is None:
        device = infer_device()
    else:
        device = args.device

    print(f"Using device: {device}")
    print(f"Task: {args.task}")
    print(f"Model: {args.model or 'default'}")
    print("-" * 50)

    # Create pipeline
    pipe = pipeline(
        args.task,
        model=args.model,
        device=device,
    )

    # Run inference based on task
    if args.task == "question-answering":
        if args.context is None:
            print("Error: --context required for question-answering")
            return
        result = pipe(question=args.input, context=args.context)
        print(f"Question: {args.input}")
        print(f"Context: {args.context}")
        print(f"\nAnswer: {result['answer']}")
        print(f"Score: {result['score']:.4f}")

    elif args.task == "text-generation":
        result = pipe(args.input, max_length=args.max_length)
        print(f"Prompt: {args.input}")
        print(f"\nGenerated: {result[0]['generated_text']}")

    elif args.task in ["sentiment-analysis", "text-classification"]:
        result = pipe(args.input)
        print(f"Text: {args.input}")
        print(f"\nLabel: {result[0]['label']}")
        print(f"Score: {result[0]['score']:.4f}")

    else:
        # Generic handling for other tasks
        result = pipe(args.input)
        print(f"Input: {args.input}")
        print(f"\nResult: {result}")


if __name__ == "__main__":
    main()
