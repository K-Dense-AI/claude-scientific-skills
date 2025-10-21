---
name: transformers
description: "Essential toolkit for Hugging Face Transformers library enabling state-of-the-art machine learning across natural language processing, computer vision, audio processing, and multimodal applications. Use this skill for: loading and using pretrained transformer models (BERT, GPT, T5, RoBERTa, DistilBERT, BART, T5, ViT, CLIP, Whisper, Llama, Mistral), implementing text generation and completion, fine-tuning models for custom tasks, text classification and sentiment analysis, question answering and reading comprehension, named entity recognition and token classification, text summarization and translation, image classification and object detection, speech recognition and audio processing, multimodal tasks combining text and images, parameter-efficient fine-tuning with LoRA and adapters, model quantization and optimization, training custom transformer models, implementing chat interfaces and conversational AI, working with tokenizers and text preprocessing, handling model inference and deployment, managing GPU memory and device allocation, implementing custom training loops, using pipelines for quick inference, working with Hugging Face Hub for model sharing, and any machine learning task involving transformer architectures or attention mechanisms."
---

# Transformers

## Overview

Transformers is Hugging Face's flagship library providing unified access to over 1 million pretrained models for machine learning across text, vision, audio, and multimodal domains. The library serves as a standardized model-definition framework compatible with PyTorch, TensorFlow, and JAX, emphasizing ease of use through three core components:

- **Pipeline**: Simple, optimized inference API for common tasks
- **AutoClasses**: Automatic model/tokenizer selection from pretrained checkpoints
- **Trainer**: Full-featured training loop with distributed training, mixed precision, and optimization

The library prioritizes accessibility with pretrained models that reduce computational costs and carbon footprint while providing compatibility across major training frameworks (PyTorch-Lightning, DeepSpeed, vLLM, etc.).

## Quick Start with Pipelines

Use pipelines for simple, efficient inference without managing models, tokenizers, or preprocessing manually. Pipelines abstract complexity into a single function call.

### Basic Pipeline Usage

```python
from transformers import pipeline

# Text classification
classifier = pipeline("text-classification")
result = classifier("This restaurant is awesome")
# [{'label': 'POSITIVE', 'score': 0.9998}]

# Text generation
generator = pipeline("text-generation", model="meta-llama/Llama-2-7b-hf")
generator("The secret to baking a good cake is", max_length=50)

# Question answering
qa = pipeline("question-answering")
qa(question="What is extractive QA?", context="Extractive QA is...")

# Image classification
img_classifier = pipeline("image-classification")
img_classifier("path/to/image.jpg")

# Automatic speech recognition
transcriber = pipeline("automatic-speech-recognition")
transcriber("audio_file.mp3")
```

### Available Pipeline Tasks

**NLP Tasks:**
- `text-classification`, `token-classification`, `question-answering`
- `fill-mask`, `summarization`, `translation`
- `text-generation`, `conversational`
- `zero-shot-classification`, `sentiment-analysis`

**Vision Tasks:**
- `image-classification`, `image-segmentation`, `object-detection`
- `depth-estimation`, `image-to-image`, `zero-shot-image-classification`

**Audio Tasks:**
- `automatic-speech-recognition`, `audio-classification`
- `text-to-audio`, `zero-shot-audio-classification`

**Multimodal Tasks:**
- `visual-question-answering`, `document-question-answering`
- `image-to-text`, `zero-shot-object-detection`

### Pipeline Best Practices

**Device Management:**
```python
from transformers import pipeline, infer_device

device = infer_device()  # Auto-detect best device
pipe = pipeline("text-generation", model="...", device=device)
```

**Batch Processing:**
```python
# Process multiple inputs efficiently
results = classifier(["Text 1", "Text 2", "Text 3"])

# Use KeyDataset for large datasets
from transformers.pipelines.pt_utils import KeyDataset
from datasets import load_dataset

dataset = load_dataset("imdb", split="test")
for result in pipe(KeyDataset(dataset, "text")):
    print(result)
```

**Memory Optimization:**
```python
# Use half-precision for faster inference
pipe = pipeline("text-generation", model="...",
                torch_dtype=torch.float16, device="cuda")
```

## Core Components

### AutoClasses for Model Loading

AutoClasses automatically select the correct architecture based on pretrained checkpoints.

```python
from transformers import (
    AutoModel, AutoTokenizer, AutoConfig,
    AutoModelForCausalLM, AutoModelForSequenceClassification
)

# Load any model by checkpoint name
tokenizer = AutoTokenizer.from_pretrained("bert-base-uncased")
model = AutoModel.from_pretrained("bert-base-uncased")

# Task-specific model classes
causal_lm = AutoModelForCausalLM.from_pretrained("gpt2")
classifier = AutoModelForSequenceClassification.from_pretrained(
    "bert-base-uncased",
    num_labels=3
)

# Load with device and dtype optimization
model = AutoModelForCausalLM.from_pretrained(
    "meta-llama/Llama-2-7b-hf",
    device_map="auto",      # Automatically distribute across devices
    torch_dtype="auto"      # Use optimal dtype
)
```

**Key Parameters:**
- `device_map="auto"`: Optimal device allocation (CPU/GPU/multi-GPU)
- `torch_dtype`: Control precision (torch.float16, torch.bfloat16, "auto")
- `trust_remote_code`: Enable custom model code (use cautiously)
- `use_fast`: Enable Rust-backed fast tokenizers (default True)

### Tokenization

Tokenizers convert text to model-compatible tensor inputs.

```python
from transformers import AutoTokenizer

tokenizer = AutoTokenizer.from_pretrained("bert-base-uncased")

# Basic tokenization
tokens = tokenizer.tokenize("Hello, how are you?")
# ['hello', ',', 'how', 'are', 'you', '?']

# Encoding (text → token IDs)
encoded = tokenizer("Hello, how are you?", return_tensors="pt")
# {'input_ids': tensor([[...]], 'attention_mask': tensor([[...]])}

# Batch encoding with padding and truncation
batch = tokenizer(
    ["Short text", "This is a much longer text..."],
    padding=True,           # Pad to longest in batch
    truncation=True,        # Truncate to model's max length
    max_length=512,
    return_tensors="pt"
)

# Decoding (token IDs → text)
text = tokenizer.decode(encoded['input_ids'][0])
```

**Special Tokens:**
```python
# Access special tokens
tokenizer.pad_token      # Padding token
tokenizer.cls_token      # Classification token
tokenizer.sep_token      # Separator token
tokenizer.mask_token     # Mask token (for MLM)

# Add custom tokens
tokenizer.add_tokens(["[CUSTOM]"])
tokenizer.add_special_tokens({'additional_special_tokens': ['[NEW]']})

# Resize model embeddings to match new vocabulary
model.resize_token_embeddings(len(tokenizer))
```

### Image Processors

For vision tasks, use image processors instead of tokenizers.

```python
from transformers import AutoImageProcessor

processor = AutoImageProcessor.from_pretrained("google/vit-base-patch16-224")

# Process single image
from PIL import Image
image = Image.open("path/to/image.jpg")
inputs = processor(image, return_tensors="pt")
# Returns: {'pixel_values': tensor([[...]])}

# Batch processing
images = [Image.open(f"img{i}.jpg") for i in range(3)]
inputs = processor(images, return_tensors="pt")
```

### Processors for Multimodal Models

Multimodal models use processors that combine image and text processing.

```python
from transformers import AutoProcessor

processor = AutoProcessor.from_pretrained("microsoft/git-base")

# Process image + text caption
inputs = processor(
    images=image,
    text="A description of the image",
    return_tensors="pt",
    padding=True
)
```

## Model Inference

### Basic Inference Pattern

```python
from transformers import AutoModelForCausalLM, AutoTokenizer

# Load model and tokenizer
model = AutoModelForCausalLM.from_pretrained("gpt2")
tokenizer = AutoTokenizer.from_pretrained("gpt2")

# Tokenize input
inputs = tokenizer("The future of AI is", return_tensors="pt")

# Generate (for causal LM)
outputs = model.generate(**inputs, max_length=50)
text = tokenizer.decode(outputs[0])

# Or get model outputs directly
outputs = model(**inputs)
logits = outputs.logits  # Shape: (batch_size, seq_len, vocab_size)
```

### Text Generation Strategies

For generative models, control generation behavior with parameters:

```python
# Greedy decoding (default)
output = model.generate(inputs, max_length=50)

# Beam search (multiple hypothesis)
output = model.generate(
    inputs,
    max_length=50,
    num_beams=5,           # Keep top 5 beams
    early_stopping=True
)

# Sampling with temperature
output = model.generate(
    inputs,
    max_length=50,
    do_sample=True,
    temperature=0.7,       # Lower = more focused, higher = more random
    top_k=50,              # Sample from top 50 tokens
    top_p=0.95             # Nucleus sampling
)

# Streaming generation
from transformers import TextStreamer

streamer = TextStreamer(tokenizer)
model.generate(**inputs, streamer=streamer, max_length=100)
```

**Generation Parameters:**
- `max_length` / `max_new_tokens`: Control output length
- `num_beams`: Beam search width (1 = greedy)
- `temperature`: Randomness (0.7-1.0 typical)
- `top_k`: Sample from top k tokens
- `top_p`: Nucleus sampling threshold
- `repetition_penalty`: Discourage repetition (>1.0)

Refer to `references/generation_strategies.md` for detailed information on choosing appropriate strategies.

## Training and Fine-Tuning

### Training Workflow Overview

1. **Load dataset** → 2. **Preprocess** → 3. **Configure training** → 4. **Train** → 5. **Evaluate** → 6. **Save/Share**

### Text Classification Example

```python
from transformers import (
    AutoTokenizer, AutoModelForSequenceClassification,
    TrainingArguments, Trainer, DataCollatorWithPadding
)
from datasets import load_dataset

# 1. Load dataset
dataset = load_dataset("imdb")

# 2. Preprocess
tokenizer = AutoTokenizer.from_pretrained("bert-base-uncased")

def preprocess(examples):
    return tokenizer(examples["text"], truncation=True)

tokenized = dataset.map(preprocess, batched=True)
data_collator = DataCollatorWithPadding(tokenizer=tokenizer)

# 3. Load model
model = AutoModelForSequenceClassification.from_pretrained(
    "bert-base-uncased",
    num_labels=2,
    id2label={0: "negative", 1: "positive"},
    label2id={"negative": 0, "positive": 1}
)

# 4. Configure training
training_args = TrainingArguments(
    output_dir="./results",
    learning_rate=2e-5,
    per_device_train_batch_size=16,
    per_device_eval_batch_size=16,
    num_train_epochs=3,
    weight_decay=0.01,
    eval_strategy="epoch",
    save_strategy="epoch",
    load_best_model_at_end=True,
    push_to_hub=False,
)

# 5. Train
trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=tokenized["train"],
    eval_dataset=tokenized["test"],
    tokenizer=tokenizer,
    data_collator=data_collator,
)

trainer.train()

# 6. Evaluate and save
metrics = trainer.evaluate()
trainer.save_model("./my-finetuned-model")
trainer.push_to_hub()  # Share to Hugging Face Hub
```

### Vision Task Fine-Tuning

```python
from transformers import (
    AutoImageProcessor, AutoModelForImageClassification,
    TrainingArguments, Trainer
)
from datasets import load_dataset

# Load dataset
dataset = load_dataset("food101", split="train[:5000]")

# Image preprocessing
processor = AutoImageProcessor.from_pretrained("google/vit-base-patch16-224")

def transform(examples):
    examples["pixel_values"] = [
        processor(img.convert("RGB"), return_tensors="pt")["pixel_values"][0]
        for img in examples["image"]
    ]
    return examples

dataset = dataset.with_transform(transform)

# Load model
model = AutoModelForImageClassification.from_pretrained(
    "google/vit-base-patch16-224",
    num_labels=101,  # 101 food categories
    ignore_mismatched_sizes=True
)

# Training (similar pattern to text)
training_args = TrainingArguments(
    output_dir="./vit-food101",
    remove_unused_columns=False,  # Keep image data
    eval_strategy="epoch",
    save_strategy="epoch",
    learning_rate=5e-5,
    per_device_train_batch_size=32,
    num_train_epochs=3,
)

trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=dataset,
    tokenizer=processor,
)

trainer.train()
```

### Sequence-to-Sequence Tasks

For tasks like summarization, translation, use Seq2SeqTrainer:

```python
from transformers import (
    AutoTokenizer, AutoModelForSeq2SeqLM,
    Seq2SeqTrainingArguments, Seq2SeqTrainer,
    DataCollatorForSeq2Seq
)

tokenizer = AutoTokenizer.from_pretrained("t5-small")
model = AutoModelForSeq2SeqLM.from_pretrained("t5-small")

def preprocess(examples):
    # Prefix input for T5
    inputs = ["summarize: " + doc for doc in examples["text"]]
    model_inputs = tokenizer(inputs, max_length=1024, truncation=True)

    # Tokenize targets
    labels = tokenizer(
        examples["summary"],
        max_length=128,
        truncation=True
    )
    model_inputs["labels"] = labels["input_ids"]
    return model_inputs

tokenized_dataset = dataset.map(preprocess, batched=True)

training_args = Seq2SeqTrainingArguments(
    output_dir="./t5-summarization",
    eval_strategy="epoch",
    learning_rate=2e-5,
    per_device_train_batch_size=8,
    num_train_epochs=3,
    predict_with_generate=True,  # Important for seq2seq
)

trainer = Seq2SeqTrainer(
    model=model,
    args=training_args,
    train_dataset=tokenized_dataset["train"],
    eval_dataset=tokenized_dataset["test"],
    tokenizer=tokenizer,
    data_collator=DataCollatorForSeq2Seq(tokenizer=tokenizer),
)

trainer.train()
```

### Important TrainingArguments

```python
TrainingArguments(
    # Essential
    output_dir="./results",
    num_train_epochs=3,
    per_device_train_batch_size=8,
    learning_rate=2e-5,

    # Evaluation
    eval_strategy="epoch",        # or "steps"
    eval_steps=500,               # if eval_strategy="steps"

    # Checkpointing
    save_strategy="epoch",
    save_steps=500,
    save_total_limit=2,           # Keep only 2 best checkpoints
    load_best_model_at_end=True,
    metric_for_best_model="accuracy",

    # Optimization
    gradient_accumulation_steps=4,
    warmup_steps=500,
    weight_decay=0.01,
    max_grad_norm=1.0,

    # Mixed Precision
    fp16=True,                    # For Nvidia GPUs
    bf16=True,                    # For Ampere+ GPUs (better)

    # Logging
    logging_steps=100,
    report_to="tensorboard",      # or "wandb", "mlflow"

    # Memory Optimization
    gradient_checkpointing=True,
    optim="adamw_torch",          # or "adafactor" for memory

    # Distributed Training
    ddp_find_unused_parameters=False,
)
```

Refer to `references/training_guide.md` for comprehensive training patterns and optimization strategies.

## Performance Optimization

### Model Quantization

Reduce memory footprint while maintaining accuracy:

```python
from transformers import AutoModelForCausalLM, BitsAndBytesConfig

# 8-bit quantization
model = AutoModelForCausalLM.from_pretrained(
    "meta-llama/Llama-2-7b-hf",
    load_in_8bit=True,
    device_map="auto"
)

# 4-bit quantization (even smaller)
bnb_config = BitsAndBytesConfig(
    load_in_4bit=True,
    bnb_4bit_quant_type="nf4",
    bnb_4bit_compute_dtype=torch.float16,
    bnb_4bit_use_double_quant=True,
)

model = AutoModelForCausalLM.from_pretrained(
    "meta-llama/Llama-2-7b-hf",
    quantization_config=bnb_config,
    device_map="auto"
)
```

**Quantization Methods:**
- **Bitsandbytes**: 4/8-bit on-the-fly quantization, supports PEFT fine-tuning
- **GPTQ**: 2/3/4/8-bit, requires calibration, very fast inference
- **AWQ**: 4-bit activation-aware, balanced speed/accuracy

Refer to `references/quantization.md` for detailed comparison and usage patterns.

### Training Optimization

```python
# Gradient accumulation (simulate larger batch)
training_args = TrainingArguments(
    per_device_train_batch_size=4,
    gradient_accumulation_steps=8,  # Effective batch = 4 * 8 = 32
)

# Gradient checkpointing (reduce memory, slower)
training_args = TrainingArguments(
    gradient_checkpointing=True,
)

# Mixed precision training
training_args = TrainingArguments(
    bf16=True,  # or fp16=True
)

# Efficient optimizer
training_args = TrainingArguments(
    optim="adafactor",  # Lower memory than AdamW
)
```

**Key Strategies:**
- **Batch sizes**: Use powers of 2 (8, 16, 32, 64, 128)
- **Gradient accumulation**: Enables larger effective batch sizes
- **Gradient checkpointing**: Reduces memory ~60%, increases time ~20%
- **Mixed precision**: bf16 for Ampere+ GPUs, fp16 for older
- **torch.compile**: Optimize model graph (PyTorch 2.0+)

## Advanced Features

### Custom Training Loop

For maximum control, bypass Trainer:

```python
from torch.utils.data import DataLoader
from transformers import AdamW, get_scheduler

# Prepare data
train_dataloader = DataLoader(tokenized_dataset, batch_size=8, shuffle=True)

# Setup optimizer and scheduler
optimizer = AdamW(model.parameters(), lr=5e-5)
scheduler = get_scheduler(
    "linear",
    optimizer=optimizer,
    num_warmup_steps=0,
    num_training_steps=len(train_dataloader) * num_epochs
)

# Training loop
model.train()
for epoch in range(num_epochs):
    for batch in train_dataloader:
        batch = {k: v.to(device) for k, v in batch.items()}

        outputs = model(**batch)
        loss = outputs.loss
        loss.backward()

        optimizer.step()
        scheduler.step()
        optimizer.zero_grad()
```

### Parameter-Efficient Fine-Tuning (PEFT)

Use PEFT library with transformers for efficient fine-tuning:

```python
from peft import LoraConfig, get_peft_model

# Configure LoRA
lora_config = LoraConfig(
    r=16,                   # Low-rank dimension
    lora_alpha=32,
    target_modules=["q_proj", "v_proj"],  # Which layers to adapt
    lora_dropout=0.05,
    bias="none",
    task_type="CAUSAL_LM"
)

# Apply to model
model = AutoModelForCausalLM.from_pretrained("meta-llama/Llama-2-7b-hf")
model = get_peft_model(model, lora_config)

# Now train as usual - only LoRA parameters train
trainer = Trainer(model=model, ...)
trainer.train()
```

### Chat Templates

Apply chat templates for instruction-tuned models:

```python
tokenizer = AutoTokenizer.from_pretrained("meta-llama/Llama-2-7b-chat-hf")

messages = [
    {"role": "system", "content": "You are a helpful assistant."},
    {"role": "user", "content": "What is machine learning?"},
]

# Format according to model's chat template
formatted = tokenizer.apply_chat_template(
    messages,
    tokenize=False,
    add_generation_prompt=True
)

# Tokenize and generate
inputs = tokenizer(formatted, return_tensors="pt")
outputs = model.generate(**inputs, max_length=200)
response = tokenizer.decode(outputs[0])
```

### Multi-GPU Training

```python
# Automatic with Trainer - no code changes needed
# Just run with: accelerate launch train.py

# Or use PyTorch DDP explicitly
training_args = TrainingArguments(
    output_dir="./results",
    ddp_find_unused_parameters=False,
    # ... other args
)

# For larger models, use FSDP
training_args = TrainingArguments(
    output_dir="./results",
    fsdp="full_shard auto_wrap",
    fsdp_config={
        "fsdp_transformer_layer_cls_to_wrap": ["BertLayer"],
    },
)
```

## Task-Specific Patterns

### Question Answering (Extractive)

```python
from transformers import pipeline

qa = pipeline("question-answering", model="distilbert-base-cased-distilled-squad")

result = qa(
    question="What is extractive QA?",
    context="Extractive QA extracts the answer from the given context..."
)
# {'answer': 'extracts the answer from the given context', 'score': 0.97, ...}
```

### Named Entity Recognition

```python
ner = pipeline("token-classification", model="dslim/bert-base-NER")

result = ner("My name is John and I live in New York")
# [{'entity': 'B-PER', 'word': 'John', ...}, {'entity': 'B-LOC', 'word': 'New York', ...}]
```

### Image Captioning

```python
from transformers import AutoProcessor, AutoModelForCausalLM

processor = AutoProcessor.from_pretrained("microsoft/git-base")
model = AutoModelForCausalLM.from_pretrained("microsoft/git-base")

from PIL import Image
image = Image.open("image.jpg")

inputs = processor(images=image, return_tensors="pt")
outputs = model.generate(**inputs, max_length=50)
caption = processor.batch_decode(outputs, skip_special_tokens=True)[0]
```

### Speech Recognition

```python
transcriber = pipeline(
    "automatic-speech-recognition",
    model="openai/whisper-base"
)

result = transcriber("audio.mp3")
# {'text': 'This is the transcribed text...'}

# With timestamps
result = transcriber("audio.mp3", return_timestamps=True)
```

## Common Patterns and Best Practices

### Saving and Loading Models

```python
# Save entire model
model.save_pretrained("./my-model")
tokenizer.save_pretrained("./my-model")

# Load later
model = AutoModel.from_pretrained("./my-model")
tokenizer = AutoTokenizer.from_pretrained("./my-model")

# Push to Hugging Face Hub
model.push_to_hub("username/my-model")
tokenizer.push_to_hub("username/my-model")

# Load from Hub
model = AutoModel.from_pretrained("username/my-model")
```

### Error Handling

```python
from transformers import AutoModel
import torch

try:
    model = AutoModel.from_pretrained("model-name")
except OSError:
    print("Model not found - check internet connection or model name")
except torch.cuda.OutOfMemoryError:
    print("GPU memory exceeded - try quantization or smaller batch size")
```

### Device Management

```python
import torch

# Check device availability
device = "cuda" if torch.cuda.is_available() else "cpu"

# Move model to device
model = model.to(device)

# Or use device_map for automatic distribution
model = AutoModel.from_pretrained("model-name", device_map="auto")

# For inputs
inputs = tokenizer(text, return_tensors="pt").to(device)
```

### Memory Management

```python
import torch

# Clear CUDA cache
torch.cuda.empty_cache()

# Use context manager for inference
with torch.no_grad():
    outputs = model(**inputs)

# Delete unused models
del model
torch.cuda.empty_cache()
```

## Resources

This skill includes comprehensive reference documentation and example scripts:

### scripts/

- `quick_inference.py`: Ready-to-use script for running inference with pipelines
- `fine_tune_classifier.py`: Complete example for fine-tuning a text classifier
- `generate_text.py`: Text generation with various strategies

Execute scripts directly or read them as implementation templates.

### references/

- `api_reference.md`: Comprehensive API documentation for key classes
- `training_guide.md`: Detailed training patterns, optimization, and troubleshooting
- `generation_strategies.md`: In-depth guide to text generation methods
- `quantization.md`: Model quantization techniques comparison and usage
- `task_patterns.md`: Quick reference for common task implementations

Load reference files when you need detailed information on specific topics. References contain extensive examples, parameter explanations, and best practices.

## Troubleshooting

**Import errors:**
```bash
pip install transformers
pip install accelerate  # For device_map="auto"
pip install bitsandbytes  # For quantization
```

**CUDA out of memory:**
- Reduce batch size
- Enable gradient checkpointing
- Use gradient accumulation
- Try quantization (8-bit or 4-bit)
- Use smaller model variant

**Slow training:**
- Enable mixed precision (fp16/bf16)
- Increase batch size (if memory allows)
- Use torch.compile (PyTorch 2.0+)
- Check data loading isn't bottleneck

**Poor generation quality:**
- Adjust temperature (lower = more focused)
- Try different decoding strategies (beam search vs sampling)
- Increase max_length if outputs cut off
- Use repetition_penalty to reduce repetition

For task-specific guidance, consult the appropriate reference file in the `references/` directory.
