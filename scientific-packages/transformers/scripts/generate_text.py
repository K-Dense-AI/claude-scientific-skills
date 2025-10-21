#!/usr/bin/env python3
"""
Text generation with different decoding strategies.

This script demonstrates various text generation approaches using
different sampling and decoding strategies.
"""

from transformers import AutoModelForCausalLM, AutoTokenizer, GenerationConfig


def load_model_and_tokenizer(model_name="gpt2"):
    """
    Load model and tokenizer.

    Args:
        model_name: Name of the model to load

    Returns:
        model, tokenizer
    """
    print(f"Loading model: {model_name}")
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModelForCausalLM.from_pretrained(model_name)

    # Set pad token if not already set
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token

    return model, tokenizer


def generate_with_greedy(model, tokenizer, prompt, max_new_tokens=50):
    """Greedy decoding - always picks highest probability token."""
    print("\n=== Greedy Decoding ===")
    print(f"Prompt: {prompt}")

    inputs = tokenizer(prompt, return_tensors="pt")
    outputs = model.generate(
        **inputs,
        max_new_tokens=max_new_tokens,
        do_sample=False,
        num_beams=1,
        pad_token_id=tokenizer.pad_token_id
    )

    generated_text = tokenizer.decode(outputs[0], skip_special_tokens=True)
    print(f"Generated: {generated_text}\n")


def generate_with_beam_search(model, tokenizer, prompt, max_new_tokens=50, num_beams=5):
    """Beam search - explores multiple hypotheses."""
    print("\n=== Beam Search ===")
    print(f"Prompt: {prompt}")

    inputs = tokenizer(prompt, return_tensors="pt")
    outputs = model.generate(
        **inputs,
        max_new_tokens=max_new_tokens,
        num_beams=num_beams,
        early_stopping=True,
        no_repeat_ngram_size=2,
        pad_token_id=tokenizer.pad_token_id
    )

    generated_text = tokenizer.decode(outputs[0], skip_special_tokens=True)
    print(f"Generated: {generated_text}\n")


def generate_with_sampling(model, tokenizer, prompt, max_new_tokens=50,
                           temperature=0.7, top_k=50, top_p=0.9):
    """Sampling with temperature, top-k, and nucleus (top-p) sampling."""
    print("\n=== Sampling (Temperature + Top-K + Top-P) ===")
    print(f"Prompt: {prompt}")
    print(f"Parameters: temperature={temperature}, top_k={top_k}, top_p={top_p}")

    inputs = tokenizer(prompt, return_tensors="pt")
    outputs = model.generate(
        **inputs,
        max_new_tokens=max_new_tokens,
        do_sample=True,
        temperature=temperature,
        top_k=top_k,
        top_p=top_p,
        pad_token_id=tokenizer.pad_token_id
    )

    generated_text = tokenizer.decode(outputs[0], skip_special_tokens=True)
    print(f"Generated: {generated_text}\n")


def generate_multiple_sequences(model, tokenizer, prompt, max_new_tokens=50,
                                 num_return_sequences=3):
    """Generate multiple diverse sequences."""
    print("\n=== Multiple Sequences (with Sampling) ===")
    print(f"Prompt: {prompt}")
    print(f"Generating {num_return_sequences} sequences...")

    inputs = tokenizer(prompt, return_tensors="pt")
    outputs = model.generate(
        **inputs,
        max_new_tokens=max_new_tokens,
        do_sample=True,
        temperature=0.8,
        top_p=0.95,
        num_return_sequences=num_return_sequences,
        pad_token_id=tokenizer.pad_token_id
    )

    for i, output in enumerate(outputs):
        generated_text = tokenizer.decode(output, skip_special_tokens=True)
        print(f"\nSequence {i+1}: {generated_text}")
    print()


def generate_with_config(model, tokenizer, prompt):
    """Use GenerationConfig for reusable configuration."""
    print("\n=== Using GenerationConfig ===")
    print(f"Prompt: {prompt}")

    # Create a generation config
    generation_config = GenerationConfig(
        max_new_tokens=50,
        do_sample=True,
        temperature=0.7,
        top_p=0.9,
        top_k=50,
        repetition_penalty=1.2,
        no_repeat_ngram_size=3,
        pad_token_id=tokenizer.pad_token_id
    )

    inputs = tokenizer(prompt, return_tensors="pt")
    outputs = model.generate(**inputs, generation_config=generation_config)

    generated_text = tokenizer.decode(outputs[0], skip_special_tokens=True)
    print(f"Generated: {generated_text}\n")


def compare_temperatures(model, tokenizer, prompt, max_new_tokens=50):
    """Compare different temperature settings."""
    print("\n=== Temperature Comparison ===")
    print(f"Prompt: {prompt}\n")

    temperatures = [0.3, 0.7, 1.0, 1.5]

    for temp in temperatures:
        inputs = tokenizer(prompt, return_tensors="pt")
        outputs = model.generate(
            **inputs,
            max_new_tokens=max_new_tokens,
            do_sample=True,
            temperature=temp,
            top_p=0.9,
            pad_token_id=tokenizer.pad_token_id
        )

        generated_text = tokenizer.decode(outputs[0], skip_special_tokens=True)
        print(f"Temperature {temp}: {generated_text}\n")


def main():
    """Run all generation examples."""
    print("=" * 70)
    print("Text Generation Examples")
    print("=" * 70)

    # Load model and tokenizer
    model, tokenizer = load_model_and_tokenizer("gpt2")

    # Example prompts
    story_prompt = "Once upon a time in a distant galaxy"
    factual_prompt = "The three branches of the US government are"

    # Demonstrate different strategies
    generate_with_greedy(model, tokenizer, story_prompt)
    generate_with_beam_search(model, tokenizer, factual_prompt)
    generate_with_sampling(model, tokenizer, story_prompt)
    generate_multiple_sequences(model, tokenizer, story_prompt, num_return_sequences=3)
    generate_with_config(model, tokenizer, story_prompt)
    compare_temperatures(model, tokenizer, story_prompt)

    print("=" * 70)
    print("All generation examples completed!")
    print("=" * 70)


if __name__ == "__main__":
    main()
