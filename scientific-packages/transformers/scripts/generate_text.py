#!/usr/bin/env python3
"""
Text generation with various strategies.

This script demonstrates different generation strategies:
- Greedy decoding
- Beam search
- Sampling with temperature
- Top-k and top-p sampling

Usage:
    python generate_text.py --model gpt2 --prompt "The future of AI" --strategy sampling
"""

import argparse
import torch
from transformers import AutoModelForCausalLM, AutoTokenizer


def generate_with_greedy(model, tokenizer, prompt, max_length):
    """Greedy decoding (deterministic)."""
    print("\n" + "=" * 60)
    print("GREEDY DECODING")
    print("=" * 60)

    inputs = tokenizer(prompt, return_tensors="pt").to(model.device)

    outputs = model.generate(
        **inputs,
        max_new_tokens=max_length,
        pad_token_id=tokenizer.eos_token_id,
    )

    text = tokenizer.decode(outputs[0], skip_special_tokens=True)
    print(f"\nPrompt: {prompt}")
    print(f"\nGenerated:\n{text}")


def generate_with_beam_search(model, tokenizer, prompt, max_length, num_beams=5):
    """Beam search for higher quality."""
    print("\n" + "=" * 60)
    print(f"BEAM SEARCH (num_beams={num_beams})")
    print("=" * 60)

    inputs = tokenizer(prompt, return_tensors="pt").to(model.device)

    outputs = model.generate(
        **inputs,
        max_new_tokens=max_length,
        num_beams=num_beams,
        early_stopping=True,
        no_repeat_ngram_size=2,
        pad_token_id=tokenizer.eos_token_id,
    )

    text = tokenizer.decode(outputs[0], skip_special_tokens=True)
    print(f"\nPrompt: {prompt}")
    print(f"\nGenerated:\n{text}")


def generate_with_sampling(model, tokenizer, prompt, max_length, temperature=0.8):
    """Sampling with temperature."""
    print("\n" + "=" * 60)
    print(f"SAMPLING (temperature={temperature})")
    print("=" * 60)

    inputs = tokenizer(prompt, return_tensors="pt").to(model.device)

    outputs = model.generate(
        **inputs,
        max_new_tokens=max_length,
        do_sample=True,
        temperature=temperature,
        pad_token_id=tokenizer.eos_token_id,
    )

    text = tokenizer.decode(outputs[0], skip_special_tokens=True)
    print(f"\nPrompt: {prompt}")
    print(f"\nGenerated:\n{text}")


def generate_with_top_k_top_p(model, tokenizer, prompt, max_length, top_k=50, top_p=0.95, temperature=0.8):
    """Top-k and top-p (nucleus) sampling."""
    print("\n" + "=" * 60)
    print(f"TOP-K TOP-P SAMPLING (k={top_k}, p={top_p}, temp={temperature})")
    print("=" * 60)

    inputs = tokenizer(prompt, return_tensors="pt").to(model.device)

    outputs = model.generate(
        **inputs,
        max_new_tokens=max_length,
        do_sample=True,
        top_k=top_k,
        top_p=top_p,
        temperature=temperature,
        repetition_penalty=1.2,
        no_repeat_ngram_size=3,
        pad_token_id=tokenizer.eos_token_id,
    )

    text = tokenizer.decode(outputs[0], skip_special_tokens=True)
    print(f"\nPrompt: {prompt}")
    print(f"\nGenerated:\n{text}")


def generate_multiple(model, tokenizer, prompt, max_length, num_sequences=3):
    """Generate multiple diverse sequences."""
    print("\n" + "=" * 60)
    print(f"MULTIPLE SEQUENCES (n={num_sequences})")
    print("=" * 60)

    inputs = tokenizer(prompt, return_tensors="pt").to(model.device)

    outputs = model.generate(
        **inputs,
        max_new_tokens=max_length,
        do_sample=True,
        num_return_sequences=num_sequences,
        temperature=0.9,
        top_p=0.95,
        pad_token_id=tokenizer.eos_token_id,
    )

    print(f"\nPrompt: {prompt}\n")
    for i, output in enumerate(outputs, 1):
        text = tokenizer.decode(output, skip_special_tokens=True)
        print(f"\n--- Sequence {i} ---\n{text}\n")


def main():
    parser = argparse.ArgumentParser(description="Text generation with various strategies")
    parser.add_argument(
        "--model",
        type=str,
        default="gpt2",
        help="Model name or path",
    )
    parser.add_argument(
        "--prompt",
        type=str,
        required=True,
        help="Input prompt for generation",
    )
    parser.add_argument(
        "--strategy",
        type=str,
        default="all",
        choices=["greedy", "beam", "sampling", "top_k_top_p", "multiple", "all"],
        help="Generation strategy to use",
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=100,
        help="Maximum number of new tokens to generate",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="auto",
        help="Device (cuda, cpu, or auto)",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=0.8,
        help="Sampling temperature",
    )
    parser.add_argument(
        "--quantize",
        action="store_true",
        help="Use 8-bit quantization",
    )

    args = parser.parse_args()

    print("=" * 60)
    print("Text Generation Demo")
    print("=" * 60)
    print(f"Model: {args.model}")
    print(f"Strategy: {args.strategy}")
    print(f"Max length: {args.max_length}")
    print(f"Device: {args.device}")
    print("=" * 60)

    # Load model and tokenizer
    print("\nLoading model...")

    if args.device == "auto":
        device_map = "auto"
        device = None
    else:
        device_map = None
        device = args.device

    model_kwargs = {"device_map": device_map} if device_map else {}

    if args.quantize:
        print("Using 8-bit quantization...")
        model_kwargs["load_in_8bit"] = True

    model = AutoModelForCausalLM.from_pretrained(args.model, **model_kwargs)
    tokenizer = AutoTokenizer.from_pretrained(args.model)

    if device and not device_map:
        model = model.to(device)

    print(f"Model loaded on: {model.device if hasattr(model, 'device') else 'multiple devices'}")

    # Generate based on strategy
    strategies = {
        "greedy": lambda: generate_with_greedy(model, tokenizer, args.prompt, args.max_length),
        "beam": lambda: generate_with_beam_search(model, tokenizer, args.prompt, args.max_length),
        "sampling": lambda: generate_with_sampling(model, tokenizer, args.prompt, args.max_length, args.temperature),
        "top_k_top_p": lambda: generate_with_top_k_top_p(model, tokenizer, args.prompt, args.max_length),
        "multiple": lambda: generate_multiple(model, tokenizer, args.prompt, args.max_length),
    }

    if args.strategy == "all":
        for strategy_fn in strategies.values():
            strategy_fn()
    else:
        strategies[args.strategy]()

    print("\n" + "=" * 60)
    print("Generation complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
