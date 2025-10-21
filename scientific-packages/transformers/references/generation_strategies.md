# Text Generation Strategies

Transformers provides flexible text generation capabilities through the `generate()` method, supporting multiple decoding strategies and configuration options.

## Basic Generation

```python
from transformers import AutoModelForCausalLM, AutoTokenizer

model = AutoModelForCausalLM.from_pretrained("gpt2")
tokenizer = AutoTokenizer.from_pretrained("gpt2")

inputs = tokenizer("Once upon a time", return_tensors="pt")
outputs = model.generate(**inputs, max_new_tokens=50)
generated_text = tokenizer.decode(outputs[0])
```

## Decoding Strategies

### 1. Greedy Decoding

Selects the token with highest probability at each step. Deterministic but can be repetitive.

```python
outputs = model.generate(
    **inputs,
    max_new_tokens=50,
    do_sample=False,
    num_beams=1  # Greedy is default when num_beams=1 and do_sample=False
)
```

### 2. Beam Search

Explores multiple hypotheses simultaneously, keeping top-k candidates at each step.

```python
outputs = model.generate(
    **inputs,
    max_new_tokens=50,
    num_beams=5,  # Number of beams
    early_stopping=True,  # Stop when all beams reach EOS
    no_repeat_ngram_size=2,  # Prevent repeating n-grams
)
```

**Key parameters:**
- `num_beams`: Number of beams (higher = more thorough but slower)
- `early_stopping`: Stop when all beams finish (True/False)
- `length_penalty`: Exponential penalty for length (>1.0 favors longer sequences)
- `no_repeat_ngram_size`: Prevent repeating n-grams

### 3. Sampling (Multinomial)

Samples from probability distribution, introducing randomness and diversity.

```python
outputs = model.generate(
    **inputs,
    max_new_tokens=50,
    do_sample=True,
    temperature=0.7,  # Controls randomness (lower = more focused)
    top_k=50,  # Consider only top-k tokens
    top_p=0.9,  # Nucleus sampling (cumulative probability threshold)
)
```

**Key parameters:**
- `temperature`: Scales logits before softmax (0.1-2.0 typical range)
  - Lower (0.1-0.7): More focused, deterministic
  - Higher (0.8-1.5): More creative, random
- `top_k`: Sample from top-k tokens only
- `top_p`: Nucleus sampling - sample from smallest set with cumulative probability > p

### 4. Beam Search with Sampling

Combines beam search with sampling for diverse but coherent outputs.

```python
outputs = model.generate(
    **inputs,
    max_new_tokens=50,
    num_beams=5,
    do_sample=True,
    temperature=0.8,
    top_k=50,
)
```

### 5. Contrastive Search

Balances coherence and diversity using contrastive objective.

```python
outputs = model.generate(
    **inputs,
    max_new_tokens=50,
    penalty_alpha=0.6,  # Contrastive penalty
    top_k=4,  # Consider top-k candidates
)
```

### 6. Assisted Decoding

Uses a smaller "assistant" model to speed up generation of larger model.

```python
from transformers import AutoModelForCausalLM

model = AutoModelForCausalLM.from_pretrained("gpt2-large")
assistant_model = AutoModelForCausalLM.from_pretrained("gpt2")

outputs = model.generate(
    **inputs,
    assistant_model=assistant_model,
    max_new_tokens=50,
)
```

## GenerationConfig

Configure generation parameters with `GenerationConfig` for reusability.

```python
from transformers import GenerationConfig

generation_config = GenerationConfig(
    max_new_tokens=100,
    do_sample=True,
    temperature=0.7,
    top_p=0.9,
    top_k=50,
    repetition_penalty=1.2,
    no_repeat_ngram_size=3,
)

# Use with model
outputs = model.generate(**inputs, generation_config=generation_config)

# Save and load
generation_config.save_pretrained("./config")
loaded_config = GenerationConfig.from_pretrained("./config")
```

## Key Parameters Reference

### Output Length Control

- `max_length`: Maximum total tokens (input + output)
- `max_new_tokens`: Maximum new tokens to generate (recommended over max_length)
- `min_length`: Minimum total tokens
- `min_new_tokens`: Minimum new tokens to generate

### Sampling Parameters

- `temperature`: Sampling temperature (0.1-2.0, default 1.0)
- `top_k`: Top-k sampling (1-100, typically 50)
- `top_p`: Nucleus sampling (0.0-1.0, typically 0.9)
- `do_sample`: Enable sampling (True/False)

### Beam Search Parameters

- `num_beams`: Number of beams (1-20, typically 5)
- `early_stopping`: Stop when beams finish (True/False)
- `length_penalty`: Length penalty (>1.0 favors longer, <1.0 favors shorter)
- `num_beam_groups`: Diverse beam search groups
- `diversity_penalty`: Penalty for similar beams

### Repetition Control

- `repetition_penalty`: Penalty for repeating tokens (1.0-2.0, default 1.0)
- `no_repeat_ngram_size`: Prevent repeating n-grams (2-5 typical)
- `encoder_repetition_penalty`: Penalty for repeating encoder tokens

### Special Tokens

- `bos_token_id`: Beginning of sequence token
- `eos_token_id`: End of sequence token (or list of tokens)
- `pad_token_id`: Padding token
- `forced_bos_token_id`: Force specific token at beginning
- `forced_eos_token_id`: Force specific token at end

### Multiple Sequences

- `num_return_sequences`: Number of sequences to return
- `num_beam_groups`: Number of diverse beam groups

## Advanced Generation Techniques

### Constrained Generation

Force generation to include specific tokens or follow patterns.

```python
from transformers import PhrasalConstraint

constraints = [
    PhrasalConstraint(tokenizer("New York", add_special_tokens=False).input_ids)
]

outputs = model.generate(
    **inputs,
    constraints=constraints,
    num_beams=5,
)
```

### Streaming Generation

Generate tokens one at a time for real-time display.

```python
from transformers import TextIteratorStreamer
from threading import Thread

streamer = TextIteratorStreamer(tokenizer, skip_special_tokens=True)

generation_kwargs = dict(
    **inputs,
    max_new_tokens=100,
    streamer=streamer,
)

thread = Thread(target=model.generate, kwargs=generation_kwargs)
thread.start()

for new_text in streamer:
    print(new_text, end="", flush=True)

thread.join()
```

### Logit Processors

Customize token selection with custom logit processors.

```python
from transformers import LogitsProcessor, LogitsProcessorList

class CustomLogitsProcessor(LogitsProcessor):
    def __call__(self, input_ids, scores):
        # Modify scores here
        return scores

logits_processor = LogitsProcessorList([CustomLogitsProcessor()])

outputs = model.generate(
    **inputs,
    logits_processor=logits_processor,
)
```

### Stopping Criteria

Define custom stopping conditions.

```python
from transformers import StoppingCriteria, StoppingCriteriaList

class CustomStoppingCriteria(StoppingCriteria):
    def __call__(self, input_ids, scores, **kwargs):
        # Return True to stop generation
        return False

stopping_criteria = StoppingCriteriaList([CustomStoppingCriteria()])

outputs = model.generate(
    **inputs,
    stopping_criteria=stopping_criteria,
)
```

## Best Practices

### For Creative Tasks (Stories, Dialogue)
```python
outputs = model.generate(
    **inputs,
    max_new_tokens=200,
    do_sample=True,
    temperature=0.8,
    top_p=0.95,
    repetition_penalty=1.2,
    no_repeat_ngram_size=3,
)
```

### For Factual Tasks (Summaries, QA)
```python
outputs = model.generate(
    **inputs,
    max_new_tokens=100,
    num_beams=4,
    early_stopping=True,
    no_repeat_ngram_size=2,
    length_penalty=1.0,
)
```

### For Chat/Instruction Following
```python
outputs = model.generate(
    **inputs,
    max_new_tokens=512,
    do_sample=True,
    temperature=0.7,
    top_p=0.9,
    repetition_penalty=1.1,
)
```

## Vision-Language Model Generation

For models like LLaVA, BLIP-2, etc.:

```python
from transformers import AutoProcessor, AutoModelForVision2Seq
from PIL import Image

model = AutoModelForVision2Seq.from_pretrained("llava-hf/llava-1.5-7b-hf")
processor = AutoProcessor.from_pretrained("llava-hf/llava-1.5-7b-hf")

image = Image.open("image.jpg")
inputs = processor(text="Describe this image", images=image, return_tensors="pt")

outputs = model.generate(
    **inputs,
    max_new_tokens=100,
    do_sample=True,
    temperature=0.7,
)

generated_text = processor.decode(outputs[0], skip_special_tokens=True)
```

## Performance Optimization

### Use KV Cache
```python
# KV cache is enabled by default
outputs = model.generate(**inputs, use_cache=True)
```

### Mixed Precision
```python
import torch

with torch.cuda.amp.autocast():
    outputs = model.generate(**inputs, max_new_tokens=100)
```

### Batch Generation
```python
texts = ["Prompt 1", "Prompt 2", "Prompt 3"]
inputs = tokenizer(texts, return_tensors="pt", padding=True)
outputs = model.generate(**inputs, max_new_tokens=50)
```

### Quantization
```python
from transformers import BitsAndBytesConfig

quantization_config = BitsAndBytesConfig(
    load_in_4bit=True,
    bnb_4bit_compute_dtype=torch.float16
)

model = AutoModelForCausalLM.from_pretrained(
    "meta-llama/Llama-2-7b-hf",
    quantization_config=quantization_config,
    device_map="auto"
)
```
