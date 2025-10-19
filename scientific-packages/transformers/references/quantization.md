# Model Quantization Guide

Comprehensive guide to reducing model memory footprint through quantization while maintaining accuracy.

## Overview

Quantization reduces memory requirements by storing model weights in lower precision formats (int8, int4) instead of full precision (float32). This enables:
- Running larger models on limited hardware
- Faster inference (reduced memory bandwidth)
- Lower deployment costs
- Enabling fine-tuning of models that wouldn't fit in memory

**Tradeoffs:**
- Slight accuracy loss (typically < 1-2%)
- Initial quantization overhead
- Some methods require calibration data

## Quick Comparison

| Method | Precision | Speed | Accuracy | Fine-tuning | Hardware | Setup |
|--------|-----------|-------|----------|-------------|----------|-------|
| **Bitsandbytes** | 4/8-bit | Fast | High | Yes (PEFT) | CUDA, CPU | Easy |
| **GPTQ** | 2-8-bit | Very Fast | High | Limited | CUDA, ROCm, Metal | Medium |
| **AWQ** | 4-bit | Very Fast | High | Yes (PEFT) | CUDA, ROCm | Medium |
| **GGUF** | 1-8-bit | Medium | Variable | No | CPU-optimized | Easy |
| **HQQ** | 1-8-bit | Fast | High | Yes | Multi-platform | Medium |

## Bitsandbytes (BnB)

On-the-fly quantization with excellent PEFT fine-tuning support.

### 8-bit Quantization

```python
from transformers import AutoModelForCausalLM, AutoTokenizer

model = AutoModelForCausalLM.from_pretrained(
    "meta-llama/Llama-2-7b-hf",
    load_in_8bit=True,              # Enable 8-bit quantization
    device_map="auto",              # Automatic device placement
)

tokenizer = AutoTokenizer.from_pretrained("meta-llama/Llama-2-7b-hf")

# Use normally
inputs = tokenizer("Hello, how are you?", return_tensors="pt").to("cuda")
outputs = model.generate(**inputs, max_new_tokens=50)
```

**Memory Savings:**
- 7B model: ~14GB → ~7GB (50% reduction)
- 13B model: ~26GB → ~13GB
- 70B model: ~140GB → ~70GB

**Characteristics:**
- Fast inference
- Minimal accuracy loss
- Works with PEFT (LoRA, QLoRA)
- Supports CPU and CUDA GPUs

### 4-bit Quantization (QLoRA)

```python
from transformers import AutoModelForCausalLM, BitsAndBytesConfig
import torch

# Configure 4-bit quantization
bnb_config = BitsAndBytesConfig(
    load_in_4bit=True,                      # Enable 4-bit quantization
    bnb_4bit_quant_type="nf4",              # Quantization type ("nf4" or "fp4")
    bnb_4bit_compute_dtype=torch.float16,   # Computation dtype
    bnb_4bit_use_double_quant=True,         # Nested quantization for more savings
)

model = AutoModelForCausalLM.from_pretrained(
    "meta-llama/Llama-2-7b-hf",
    quantization_config=bnb_config,
    device_map="auto",
)
```

**Memory Savings:**
- 7B model: ~14GB → ~4GB (70% reduction)
- 13B model: ~26GB → ~7GB
- 70B model: ~140GB → ~35GB

**Quantization Types:**
- `nf4`: Normal Float 4 (recommended, better quality)
- `fp4`: Float Point 4 (slightly more memory efficient)

**Compute Dtype:**
```python
# For better quality
bnb_4bit_compute_dtype=torch.float16

# For best performance on Ampere+ GPUs
bnb_4bit_compute_dtype=torch.bfloat16
```

**Double Quantization:**
```python
# Enable for additional ~0.4 bits/param savings
bnb_4bit_use_double_quant=True  # Quantize the quantization constants
```

### Fine-tuning with QLoRA

```python
from transformers import AutoModelForCausalLM, BitsAndBytesConfig, TrainingArguments, Trainer
from peft import LoraConfig, get_peft_model, prepare_model_for_kbit_training
import torch

# Load quantized model
bnb_config = BitsAndBytesConfig(
    load_in_4bit=True,
    bnb_4bit_quant_type="nf4",
    bnb_4bit_compute_dtype=torch.float16,
    bnb_4bit_use_double_quant=True,
)

model = AutoModelForCausalLM.from_pretrained(
    "meta-llama/Llama-2-7b-hf",
    quantization_config=bnb_config,
    device_map="auto",
)

# Prepare for training
model = prepare_model_for_kbit_training(model)

# Configure LoRA
lora_config = LoraConfig(
    r=16,
    lora_alpha=32,
    target_modules=["q_proj", "k_proj", "v_proj", "o_proj"],
    lora_dropout=0.05,
    bias="none",
    task_type="CAUSAL_LM"
)

model = get_peft_model(model, lora_config)

# Train normally
trainer = Trainer(model=model, args=training_args, ...)
trainer.train()
```

## GPTQ

Post-training quantization requiring calibration, optimized for inference speed.

### Loading GPTQ Models

```python
from transformers import AutoModelForCausalLM, AutoTokenizer, GPTQConfig

# Load pre-quantized GPTQ model
model = AutoModelForCausalLM.from_pretrained(
    "TheBloke/Llama-2-7B-GPTQ",         # Pre-quantized model
    device_map="auto",
    revision="gptq-4bit-32g-actorder_True",  # Specific quantization config
)

# Or quantize yourself
gptq_config = GPTQConfig(
    bits=4,                              # 2, 3, 4, 8 bits
    dataset="c4",                        # Calibration dataset
    tokenizer=tokenizer,
)

model = AutoModelForCausalLM.from_pretrained(
    "meta-llama/Llama-2-7b-hf",
    device_map="auto",
    quantization_config=gptq_config,
)

# Save quantized model
model.save_pretrained("llama-2-7b-gptq")
```

**Configuration Options:**
```python
gptq_config = GPTQConfig(
    bits=4,                              # Quantization bits
    group_size=128,                      # Group size for quantization (128, 32, -1)
    dataset="c4",                        # Calibration dataset
    desc_act=False,                      # Activation order (can improve accuracy)
    sym=True,                            # Symmetric quantization
    damp_percent=0.1,                    # Dampening factor
)
```

**Characteristics:**
- Fastest inference among quantization methods
- Requires one-time calibration (slow)
- Best when using pre-quantized models from Hub
- Limited fine-tuning support
- Excellent for production deployment

## AWQ (Activation-aware Weight Quantization)

Protects important weights for better quality.

### Loading AWQ Models

```python
from transformers import AutoModelForCausalLM, AwqConfig

# Load pre-quantized AWQ model
model = AutoModelForCausalLM.from_pretrained(
    "TheBloke/Llama-2-7B-AWQ",
    device_map="auto",
)

# Or quantize yourself
awq_config = AwqConfig(
    bits=4,                              # 4-bit quantization
    group_size=128,                      # Quantization group size
    zero_point=True,                     # Use zero-point quantization
    version="GEMM",                      # Quantization version
)

model = AutoModelForCausalLM.from_pretrained(
    "meta-llama/Llama-2-7b-hf",
    quantization_config=awq_config,
    device_map="auto",
)
```

**Characteristics:**
- Better accuracy than GPTQ at same bit width
- Excellent inference speed
- Supports PEFT fine-tuning
- Requires calibration data

### Fine-tuning AWQ Models

```python
from peft import LoraConfig, get_peft_model

# AWQ models support LoRA fine-tuning
lora_config = LoraConfig(
    r=16,
    lora_alpha=32,
    target_modules=["q_proj", "v_proj"],
    lora_dropout=0.05,
    task_type="CAUSAL_LM"
)

model = get_peft_model(model, lora_config)
trainer = Trainer(model=model, ...)
trainer.train()
```

## GGUF (GGML Format)

CPU-optimized quantization format, popular in llama.cpp ecosystem.

### Using GGUF Models

```python
from transformers import AutoModelForCausalLM, AutoTokenizer

# Load GGUF model
model = AutoModelForCausalLM.from_pretrained(
    "TheBloke/Llama-2-7B-GGUF",
    gguf_file="llama-2-7b.Q4_K_M.gguf",  # Specific quantization file
    device_map="auto",
)

tokenizer = AutoTokenizer.from_pretrained("TheBloke/Llama-2-7B-GGUF")
```

**GGUF Quantization Types:**
- `Q4_0`: 4-bit, smallest, lowest quality
- `Q4_K_M`: 4-bit, medium quality (recommended)
- `Q5_K_M`: 5-bit, good quality
- `Q6_K`: 6-bit, high quality
- `Q8_0`: 8-bit, very high quality

**Characteristics:**
- Optimized for CPU inference
- Wide range of bit depths (1-8)
- Good for Apple Silicon (M1/M2)
- No fine-tuning support
- Excellent for local/edge deployment

## HQQ (Half-Quadratic Quantization)

Flexible quantization with good accuracy retention.

### Using HQQ

```python
from transformers import AutoModelForCausalLM, HqqConfig

hqq_config = HqqConfig(
    nbits=4,                             # Quantization bits
    group_size=64,                       # Group size
    quant_zero=False,                    # Quantize zero point
    quant_scale=False,                   # Quantize scale
    axis=0,                              # Quantization axis
)

model = AutoModelForCausalLM.from_pretrained(
    "meta-llama/Llama-2-7b-hf",
    quantization_config=hqq_config,
    device_map="auto",
)
```

**Characteristics:**
- Very fast quantization
- No calibration data needed
- Support for 1-8 bits
- Can serialize/deserialize
- Good accuracy vs size tradeoff

## Choosing a Quantization Method

### Decision Tree

**For inference only:**
1. Need fastest inference? → **GPTQ or AWQ** (use pre-quantized models)
2. CPU-only deployment? → **GGUF**
3. Want easiest setup? → **Bitsandbytes 8-bit**
4. Need extreme compression? → **GGUF Q4_0 or HQQ 2-bit**

**For fine-tuning:**
1. Limited VRAM? → **QLoRA (BnB 4-bit + LoRA)**
2. Want best accuracy? → **Bitsandbytes 8-bit + LoRA**
3. Need very large models? → **QLoRA with double quantization**

**For production:**
1. Latency-critical? → **GPTQ or AWQ**
2. Cost-optimized? → **Bitsandbytes 8-bit**
3. CPU deployment? → **GGUF**

## Memory Requirements

Approximate memory for Llama-2 7B model:

| Method | Memory | vs FP16 |
|--------|--------|---------|
| FP32 | 28GB | 2x |
| FP16 / BF16 | 14GB | 1x |
| 8-bit (BnB) | 7GB | 0.5x |
| 4-bit (QLoRA) | 3.5GB | 0.25x |
| 4-bit Double Quant | 3GB | 0.21x |
| GPTQ 4-bit | 4GB | 0.29x |
| AWQ 4-bit | 4GB | 0.29x |

**Note:** Add ~1-2GB for inference activations, KV cache, and framework overhead.

## Best Practices

### For Training

```python
# QLoRA recommended configuration
bnb_config = BitsAndBytesConfig(
    load_in_4bit=True,
    bnb_4bit_quant_type="nf4",
    bnb_4bit_compute_dtype=torch.bfloat16,  # BF16 if available
    bnb_4bit_use_double_quant=True,
)

# LoRA configuration
lora_config = LoraConfig(
    r=16,                                    # Rank (8, 16, 32, 64)
    lora_alpha=32,                           # Scaling (typically 2*r)
    target_modules=["q_proj", "k_proj", "v_proj", "o_proj", "gate_proj", "up_proj", "down_proj"],
    lora_dropout=0.05,
    bias="none",
    task_type="CAUSAL_LM"
)
```

### For Inference

```python
# High-speed inference
model = AutoModelForCausalLM.from_pretrained(
    "TheBloke/Llama-2-7B-GPTQ",
    device_map="auto",
    torch_dtype=torch.float16,           # Use FP16 for activations
)

# Balanced quality/speed
model = AutoModelForCausalLM.from_pretrained(
    "meta-llama/Llama-2-7b-hf",
    load_in_8bit=True,
    device_map="auto",
)

# Maximum compression
model = AutoModelForCausalLM.from_pretrained(
    "meta-llama/Llama-2-7b-hf",
    quantization_config=BitsAndBytesConfig(
        load_in_4bit=True,
        bnb_4bit_quant_type="nf4",
        bnb_4bit_compute_dtype=torch.float16,
        bnb_4bit_use_double_quant=True,
    ),
    device_map="auto",
)
```

### Multi-GPU Setups

```python
# Automatically distribute across GPUs
model = AutoModelForCausalLM.from_pretrained(
    "meta-llama/Llama-2-70b-hf",
    load_in_4bit=True,
    device_map="auto",                   # Automatic distribution
    max_memory={0: "20GB", 1: "20GB"},   # Optional: limit per GPU
)

# Manual device map
device_map = {
    "model.embed_tokens": 0,
    "model.layers.0": 0,
    "model.layers.1": 0,
    # ... distribute layers ...
    "model.norm": 1,
    "lm_head": 1,
}

model = AutoModelForCausalLM.from_pretrained(
    "meta-llama/Llama-2-70b-hf",
    load_in_4bit=True,
    device_map=device_map,
)
```

## Troubleshooting

**Issue: OOM during quantization**
```python
# Solution: Use low_cpu_mem_usage
model = AutoModelForCausalLM.from_pretrained(
    "model-name",
    quantization_config=config,
    device_map="auto",
    low_cpu_mem_usage=True,              # Reduce CPU memory during loading
)
```

**Issue: Slow quantization**
```python
# GPTQ/AWQ take time to calibrate
# Solution: Use pre-quantized models from Hub
model = AutoModelForCausalLM.from_pretrained("TheBloke/Model-GPTQ")

# Or use BnB for instant quantization
model = AutoModelForCausalLM.from_pretrained("model-name", load_in_4bit=True)
```

**Issue: Poor quality after quantization**
```python
# Try different quantization types
bnb_config = BitsAndBytesConfig(
    load_in_4bit=True,
    bnb_4bit_quant_type="nf4",           # Try "nf4" instead of "fp4"
    bnb_4bit_compute_dtype=torch.bfloat16,  # Use BF16 if available
)

# Or use 8-bit instead of 4-bit
model = AutoModelForCausalLM.from_pretrained("model-name", load_in_8bit=True)
```

**Issue: Can't fine-tune quantized model**
```python
# Ensure using compatible quantization method
from peft import prepare_model_for_kbit_training

model = prepare_model_for_kbit_training(model)

# Only BnB and AWQ support PEFT fine-tuning
# GPTQ has limited support, GGUF doesn't support fine-tuning
```

## Performance Benchmarks

Approximate generation speed (tokens/sec) for Llama-2 7B on A100 40GB:

| Method | Speed | Memory |
|--------|-------|--------|
| FP16 | 100 tok/s | 14GB |
| 8-bit | 90 tok/s | 7GB |
| 4-bit QLoRA | 70 tok/s | 4GB |
| GPTQ 4-bit | 95 tok/s | 4GB |
| AWQ 4-bit | 95 tok/s | 4GB |

**Note:** Actual performance varies by hardware, sequence length, and batch size.

## Resources

- **Pre-quantized models:** Search "GPTQ" or "AWQ" on Hugging Face Hub
- **BnB documentation:** https://github.com/TimDettmers/bitsandbytes
- **PEFT library:** https://github.com/huggingface/peft
- **QLoRA paper:** https://arxiv.org/abs/2305.14314

For task-specific quantization examples, see `training_guide.md`.
