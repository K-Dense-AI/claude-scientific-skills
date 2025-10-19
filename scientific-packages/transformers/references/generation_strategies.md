# Text Generation Strategies

Comprehensive guide to text generation methods in Transformers for controlling output quality, creativity, and diversity.

## Overview

Text generation is the process of predicting tokens sequentially using a language model. The choice of generation strategy significantly impacts output quality, diversity, and computational cost.

**When to use each strategy:**
- **Greedy**: Fast, deterministic, good for short outputs or when consistency is critical
- **Beam Search**: Better quality for tasks with clear "correct" answers (translation, summarization)
- **Sampling**: Creative, diverse outputs for open-ended generation (stories, dialogue)
- **Top-k/Top-p**: Balanced creativity and coherence

## Basic Generation Methods

### Greedy Decoding

Selects the highest probability token at each step. Fast but prone to repetition and suboptimal sequences.

```python
from transformers import AutoModelForCausalLM, AutoTokenizer

model = AutoModelForCausalLM.from_pretrained("gpt2")
tokenizer = AutoTokenizer.from_pretrained("gpt2")

inputs = tokenizer("The future of AI", return_tensors="pt")

# Greedy decoding (default)
outputs = model.generate(**inputs, max_new_tokens=50)
print(tokenizer.decode(outputs[0]))
```

**Characteristics:**
- Deterministic (always same output for same input)
- Fast (single forward pass per token)
- Prone to repetition in longer sequences
- Best for: Short generations, deterministic applications

**Parameters:**
```python
outputs = model.generate(
    **inputs,
    max_new_tokens=50,              # Number of tokens to generate
    min_length=10,                  # Minimum total length
    pad_token_id=tokenizer.pad_token_id,
)
```

### Beam Search

Maintains multiple hypotheses (beams) and selects the sequence with highest overall probability.

```python
outputs = model.generate(
    **inputs,
    max_new_tokens=50,
    num_beams=5,                    # Number of beams
    early_stopping=True,            # Stop when all beams finish
    no_repeat_ngram_size=2,         # Prevent 2-gram repetition
)
```

**Characteristics:**
- Higher quality than greedy for tasks with "correct" answers
- Slower than greedy (num_beams forward passes per step)
- Still can suffer from repetition
- Best for: Translation, summarization, QA generation

**Advanced Parameters:**
```python
outputs = model.generate(
    **inputs,
    num_beams=5,
    num_beam_groups=1,              # Diverse beam search groups
    diversity_penalty=0.0,          # Penalty for similar beams
    length_penalty=1.0,             # >1: longer sequences, <1: shorter
    early_stopping=True,            # Stop when num_beams sequences finish
    no_repeat_ngram_size=2,         # Block repeating n-grams
    num_return_sequences=1,         # Return top-k sequences (≤ num_beams)
)
```

**Length Penalty:**
- `length_penalty > 1.0`: Favor longer sequences
- `length_penalty = 1.0`: No penalty
- `length_penalty < 1.0`: Favor shorter sequences

### Sampling (Multinomial)

Randomly sample tokens according to the probability distribution.

```python
outputs = model.generate(
    **inputs,
    max_new_tokens=50,
    do_sample=True,                 # Enable sampling
    temperature=1.0,                # Sampling temperature
    num_beams=1,                    # Must be 1 for sampling
)
```

**Characteristics:**
- Non-deterministic (different output each time)
- More diverse and creative than greedy/beam search
- Can produce incoherent output if not controlled
- Best for: Creative writing, dialogue, open-ended generation

**Temperature Parameter:**
```python
# Low temperature (0.1-0.7): More focused, less random
outputs = model.generate(**inputs, do_sample=True, temperature=0.5)

# Medium temperature (0.7-1.0): Balanced
outputs = model.generate(**inputs, do_sample=True, temperature=0.8)

# High temperature (1.0-2.0): More random, more creative
outputs = model.generate(**inputs, do_sample=True, temperature=1.5)
```

- `temperature → 0`: Approaches greedy decoding
- `temperature = 1.0`: Sample from original distribution
- `temperature > 1.0`: Flatter distribution, more random
- `temperature < 1.0`: Sharper distribution, more confident

## Advanced Sampling Methods

### Top-k Sampling

Sample from only the k most likely tokens.

```python
outputs = model.generate(
    **inputs,
    do_sample=True,
    max_new_tokens=50,
    top_k=50,                       # Consider top 50 tokens
    temperature=0.8,
)
```

**How it works:**
1. Filter to top-k most probable tokens
2. Renormalize probabilities
3. Sample from filtered distribution

**Choosing k:**
- `k=1`: Equivalent to greedy decoding
- `k=10-50`: More focused, coherent output
- `k=100-500`: More diverse output
- Too high k: Includes low-probability tokens (noise)
- Too low k: Less diverse, may miss good alternatives

### Top-p (Nucleus) Sampling

Sample from the smallest set of tokens whose cumulative probability ≥ p.

```python
outputs = model.generate(
    **inputs,
    do_sample=True,
    max_new_tokens=50,
    top_p=0.95,                     # Nucleus probability
    temperature=0.8,
)
```

**How it works:**
1. Sort tokens by probability
2. Find smallest set with cumulative probability ≥ p
3. Sample from this set

**Choosing p:**
- `p=0.9-0.95`: Good balance (recommended)
- `p=1.0`: Sample from full distribution
- Higher p: More diverse, might include unlikely tokens
- Lower p: More focused, like top-k with adaptive k

**Top-p vs Top-k:**
- Top-p adapts to probability distribution shape
- Top-k is fixed regardless of distribution
- Top-p generally better for variable-quality contexts
- Can combine: `top_k=50, top_p=0.95` (apply both filters)

### Combining Strategies

```python
# Recommended for high-quality open-ended generation
outputs = model.generate(
    **inputs,
    do_sample=True,
    max_new_tokens=100,
    temperature=0.8,                # Moderate temperature
    top_k=50,                       # Limit to top 50 tokens
    top_p=0.95,                     # Nucleus sampling
    repetition_penalty=1.2,         # Discourage repetition
    no_repeat_ngram_size=3,         # Block 3-gram repetition
)
```

## Controlling Generation Quality

### Repetition Control

Prevent models from repeating themselves:

```python
outputs = model.generate(
    **inputs,
    max_new_tokens=100,

    # Method 1: Repetition penalty
    repetition_penalty=1.2,         # Penalize repeated tokens (>1.0)

    # Method 2: Block n-gram repetition
    no_repeat_ngram_size=3,         # Never repeat 3-grams

    # Method 3: Encoder repetition penalty (for seq2seq)
    encoder_repetition_penalty=1.0, # Penalize input tokens
)
```

**Repetition Penalty Values:**
- `1.0`: No penalty
- `1.0-1.5`: Mild penalty (recommended: 1.1-1.3)
- `>1.5`: Strong penalty (may harm coherence)

### Length Control

```python
outputs = model.generate(
    **inputs,

    # Hard constraints
    min_length=20,                  # Minimum total length
    max_length=100,                 # Maximum total length
    max_new_tokens=50,              # Maximum new tokens (excluding input)

    # Soft constraints (with beam search)
    length_penalty=1.0,             # Encourage longer/shorter outputs

    # Early stopping
    early_stopping=True,            # Stop when condition met
)
```

### Bad Words and Forced Tokens

```python
# Prevent specific tokens
bad_words_ids = [
    tokenizer.encode("badword1", add_special_tokens=False),
    tokenizer.encode("badword2", add_special_tokens=False),
]

outputs = model.generate(
    **inputs,
    bad_words_ids=bad_words_ids,
)

# Force specific tokens
force_words_ids = [
    tokenizer.encode("important", add_special_tokens=False),
]

outputs = model.generate(
    **inputs,
    force_words_ids=force_words_ids,
)
```

## Streaming Generation

Generate and process tokens as they're produced:

```python
from transformers import TextStreamer, TextIteratorStreamer
from threading import Thread

# Simple streaming (prints to stdout)
streamer = TextStreamer(tokenizer, skip_prompt=True)
outputs = model.generate(**inputs, streamer=streamer, max_new_tokens=100)

# Iterator streaming (for custom processing)
streamer = TextIteratorStreamer(tokenizer, skip_prompt=True)

generation_kwargs = dict(**inputs, streamer=streamer, max_new_tokens=100)
thread = Thread(target=model.generate, kwargs=generation_kwargs)
thread.start()

for text in streamer:
    print(text, end="", flush=True)

thread.join()
```

## Advanced Techniques

### Contrastive Search

Balance coherence and diversity using contrastive objective:

```python
outputs = model.generate(
    **inputs,
    max_new_tokens=50,
    penalty_alpha=0.6,              # Contrastive penalty
    top_k=4,                        # Consider top-4 tokens
)
```

**When to use:**
- Open-ended text generation
- Reduces repetition without sacrificing coherence
- Good alternative to sampling

### Diverse Beam Search

Generate multiple diverse outputs:

```python
outputs = model.generate(
    **inputs,
    max_new_tokens=50,
    num_beams=10,
    num_beam_groups=5,              # 5 groups of 2 beams each
    diversity_penalty=1.0,          # Penalty for similar beams
    num_return_sequences=5,         # Return 5 diverse outputs
)
```

### Constrained Beam Search

Force output to include specific phrases:

```python
from transformers import PhrasalConstraint

constraints = [
    PhrasalConstraint(
        tokenizer("machine learning", add_special_tokens=False).input_ids
    ),
]

outputs = model.generate(
    **inputs,
    constraints=constraints,
    num_beams=10,                   # Requires beam search
)
```

## Speculative Decoding

Accelerate generation using a smaller draft model:

```python
from transformers import AutoModelForCausalLM

# Load main and assistant models
model = AutoModelForCausalLM.from_pretrained("large-model")
assistant_model = AutoModelForCausalLM.from_pretrained("small-model")

# Generate with speculative decoding
outputs = model.generate(
    **inputs,
    assistant_model=assistant_model,
    do_sample=True,
    temperature=0.8,
)
```

**Benefits:**
- 2-3x faster generation
- Identical output distribution to regular generation
- Works with sampling and greedy decoding

## Recipe: Recommended Settings by Task

### Creative Writing / Dialogue

```python
outputs = model.generate(
    **inputs,
    do_sample=True,
    max_new_tokens=200,
    temperature=0.9,
    top_p=0.95,
    top_k=50,
    repetition_penalty=1.2,
    no_repeat_ngram_size=3,
)
```

### Translation / Summarization

```python
outputs = model.generate(
    **inputs,
    num_beams=5,
    max_new_tokens=150,
    early_stopping=True,
    length_penalty=1.0,
    no_repeat_ngram_size=2,
)
```

### Code Generation

```python
outputs = model.generate(
    **inputs,
    max_new_tokens=300,
    temperature=0.2,                # Low temperature for correctness
    top_p=0.95,
    do_sample=True,
)
```

### Chatbot / Instruction Following

```python
outputs = model.generate(
    **inputs,
    do_sample=True,
    max_new_tokens=256,
    temperature=0.7,
    top_p=0.9,
    repetition_penalty=1.15,
)
```

### Factual QA / Information Extraction

```python
outputs = model.generate(
    **inputs,
    max_new_tokens=50,
    num_beams=3,
    early_stopping=True,
    # Or greedy for very short answers:
    # (no special parameters needed)
)
```

## Debugging Generation

### Check Token Probabilities

```python
outputs = model.generate(
    **inputs,
    max_new_tokens=20,
    output_scores=True,             # Return generation scores
    return_dict_in_generate=True,   # Return as dict
)

# Access generation scores
scores = outputs.scores  # Tuple of tensors (seq_len, vocab_size)

# Get token probabilities
import torch
probs = torch.softmax(scores[0], dim=-1)
```

### Monitor Generation Process

```python
from transformers import LogitsProcessor, LogitsProcessorList

class DebugLogitsProcessor(LogitsProcessor):
    def __call__(self, input_ids, scores):
        # Print top 5 tokens at each step
        top_tokens = scores[0].topk(5)
        print(f"Top 5 tokens: {top_tokens}")
        return scores

outputs = model.generate(
    **inputs,
    max_new_tokens=10,
    logits_processor=LogitsProcessorList([DebugLogitsProcessor()]),
)
```

## Common Issues and Solutions

**Issue: Repetitive output**
- Solution: Increase `repetition_penalty` (1.2-1.5), set `no_repeat_ngram_size=3`
- For sampling: Increase `temperature`, enable `top_p`

**Issue: Incoherent output**
- Solution: Lower `temperature` (0.5-0.8), use beam search
- Set `top_k=50` or `top_p=0.9` to filter unlikely tokens

**Issue: Too short output**
- Solution: Increase `min_length`, set `length_penalty > 1.0` (beam search)
- Check if EOS token is being generated early

**Issue: Too slow generation**
- Solution: Use greedy instead of beam search
- Reduce `num_beams`
- Try speculative decoding with assistant model
- Use smaller model variant

**Issue: Output doesn't follow format**
- Solution: Use constrained beam search
- Add format examples to prompt
- Use `bad_words_ids` to prevent format-breaking tokens

## Performance Optimization

```python
# Use half precision
model = AutoModelForCausalLM.from_pretrained(
    "model-name",
    torch_dtype=torch.float16,
    device_map="auto"
)

# Use KV cache optimization (default, but can be disabled)
outputs = model.generate(**inputs, use_cache=True)

# Batch generation
inputs = tokenizer(["Prompt 1", "Prompt 2"], return_tensors="pt", padding=True)
outputs = model.generate(**inputs, max_new_tokens=50)

# Static cache for longer sequences (if supported)
outputs = model.generate(**inputs, cache_implementation="static")
```

This guide covers the main generation strategies. For task-specific examples, see `task_patterns.md`.
