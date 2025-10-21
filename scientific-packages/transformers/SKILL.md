---
name: transformers
description: Work with state-of-the-art machine learning models for NLP, computer vision, audio, and multimodal tasks using HuggingFace Transformers. This skill should be used when fine-tuning pre-trained models, performing inference with pipelines, generating text, training sequence models, or working with BERT, GPT, T5, ViT, and other transformer architectures. Covers model loading, tokenization, training with Trainer API, text generation strategies, and task-specific patterns for classification, NER, QA, summarization, translation, and image tasks. (plugin:scientific-packages@claude-scientific-skills)
---

# Transformers

## Overview

The Transformers library provides state-of-the-art machine learning models for NLP, computer vision, audio, and multimodal tasks. Apply this skill for quick inference through pipelines, comprehensive training via the Trainer API, and flexible text generation with various decoding strategies.

## Core Capabilities

### 1. Quick Inference with Pipelines

For rapid inference without complex setup, use the `pipeline()` API. Pipelines abstract away tokenization, model invocation, and post-processing.

```python
from transformers import pipeline

# Text classification
classifier = pipeline("text-classification")
result = classifier("This product is amazing!")

# Named entity recognition
ner = pipeline("token-classification")
entities = ner("Sarah works at Microsoft in Seattle")

# Question answering
qa = pipeline("question-answering")
answer = qa(question="What is the capital?", context="Paris is the capital of France.")

# Text generation
generator = pipeline("text-generation", model="gpt2")
text = generator("Once upon a time", max_length=50)

# Image classification
image_classifier = pipeline("image-classification")
predictions = image_classifier("image.jpg")
```

**When to use pipelines:**
- Quick prototyping and testing
- Simple inference tasks without custom logic
- Demonstrations and examples
- Production inference for standard tasks

**Available pipeline tasks:**
- **NLP**: text-classification, token-classification, question-answering, summarization, translation, text-generation, fill-mask, zero-shot-classification
- **Vision**: image-classification, object-detection, image-segmentation, depth-estimation, zero-shot-image-classification
- **Audio**: automatic-speech-recognition, audio-classification, text-to-audio
- **Multimodal**: image-to-text, visual-question-answering, image-text-to-text

For comprehensive pipeline documentation, see `references/pipelines.md`.

### 2. Model Training and Fine-Tuning

Use the Trainer API for comprehensive model training with support for distributed training, mixed precision, and advanced optimization.

**Basic training workflow:**

```python
from transformers import (
    AutoTokenizer,
    AutoModelForSequenceClassification,
    TrainingArguments,
    Trainer
)
from datasets import load_dataset

# 1. Load and tokenize data
dataset = load_dataset("imdb")
tokenizer = AutoTokenizer.from_pretrained("bert-base-uncased")

def tokenize_function(examples):
    return tokenizer(examples["text"], padding="max_length", truncation=True)

tokenized_datasets = dataset.map(tokenize_function, batched=True)

# 2. Load model
model = AutoModelForSequenceClassification.from_pretrained(
    "bert-base-uncased",
    num_labels=2
)

# 3. Configure training
training_args = TrainingArguments(
    output_dir="./results",
    num_train_epochs=3,
    per_device_train_batch_size=16,
    eval_strategy="epoch",
    save_strategy="epoch",
    load_best_model_at_end=True,
)

# 4. Create trainer and train
trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=tokenized_datasets["train"],
    eval_dataset=tokenized_datasets["test"],
)

trainer.train()
```

**Key training features:**
- Mixed precision training (fp16/bf16)
- Distributed training (multi-GPU, multi-node)
- Gradient accumulation
- Learning rate scheduling with warmup
- Checkpoint management
- Hyperparameter search
- Push to Hugging Face Hub

For detailed training documentation, see `references/training.md`.

### 3. Text Generation

Generate text using various decoding strategies including greedy decoding, beam search, sampling, and more.

**Generation strategies:**

```python
from transformers import AutoModelForCausalLM, AutoTokenizer

model = AutoModelForCausalLM.from_pretrained("gpt2")
tokenizer = AutoTokenizer.from_pretrained("gpt2")
inputs = tokenizer("Once upon a time", return_tensors="pt")

# Greedy decoding (deterministic)
outputs = model.generate(**inputs, max_new_tokens=50)

# Beam search (explores multiple hypotheses)
outputs = model.generate(
    **inputs,
    max_new_tokens=50,
    num_beams=5,
    early_stopping=True
)

# Sampling (creative, diverse)
outputs = model.generate(
    **inputs,
    max_new_tokens=50,
    do_sample=True,
    temperature=0.7,
    top_p=0.9,
    top_k=50
)
```

**Generation parameters:**
- `temperature`: Controls randomness (0.1-2.0)
- `top_k`: Sample from top-k tokens
- `top_p`: Nucleus sampling threshold
- `num_beams`: Number of beams for beam search
- `repetition_penalty`: Discourage repetition
- `no_repeat_ngram_size`: Prevent repeating n-grams

For comprehensive generation documentation, see `references/generation_strategies.md`.

### 4. Task-Specific Patterns

Common task patterns with appropriate model classes:

**Text Classification:**
```python
from transformers import AutoModelForSequenceClassification

model = AutoModelForSequenceClassification.from_pretrained(
    "bert-base-uncased",
    num_labels=3,
    id2label={0: "negative", 1: "neutral", 2: "positive"}
)
```

**Named Entity Recognition (Token Classification):**
```python
from transformers import AutoModelForTokenClassification

model = AutoModelForTokenClassification.from_pretrained(
    "bert-base-uncased",
    num_labels=9  # Number of entity types
)
```

**Question Answering:**
```python
from transformers import AutoModelForQuestionAnswering

model = AutoModelForQuestionAnswering.from_pretrained("bert-base-uncased")
```

**Summarization and Translation (Seq2Seq):**
```python
from transformers import AutoModelForSeq2SeqLM

model = AutoModelForSeq2SeqLM.from_pretrained("t5-base")
```

**Image Classification:**
```python
from transformers import AutoModelForImageClassification

model = AutoModelForImageClassification.from_pretrained(
    "google/vit-base-patch16-224",
    num_labels=num_classes
)
```

For detailed task-specific workflows including data preprocessing, training, and evaluation, see `references/task_patterns.md`.

## Auto Classes

Use Auto classes for automatic architecture selection based on model checkpoints:

```python
from transformers import (
    AutoTokenizer,           # Tokenization
    AutoModel,               # Base model (hidden states)
    AutoModelForSequenceClassification,
    AutoModelForTokenClassification,
    AutoModelForQuestionAnswering,
    AutoModelForCausalLM,    # GPT-style
    AutoModelForMaskedLM,    # BERT-style
    AutoModelForSeq2SeqLM,   # T5, BART
    AutoProcessor,           # For multimodal models
    AutoImageProcessor,      # For vision models
)

# Load any model by name
tokenizer = AutoTokenizer.from_pretrained("bert-base-uncased")
model = AutoModelForSequenceClassification.from_pretrained("bert-base-uncased")
```

For comprehensive API documentation, see `references/api_reference.md`.

## Model Loading and Optimization

**Device placement:**
```python
model = AutoModel.from_pretrained("bert-base-uncased", device_map="auto")
```

**Mixed precision:**
```python
model = AutoModel.from_pretrained(
    "model-name",
    torch_dtype=torch.float16  # or torch.bfloat16
)
```

**Quantization:**
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

## Common Workflows

### Quick Inference Workflow
1. Choose appropriate pipeline for task
2. Load pipeline with optional model specification
3. Pass inputs and get results
4. For batch processing, pass list of inputs

**See:** `scripts/quick_inference.py` for comprehensive pipeline examples

### Training Workflow
1. Load and preprocess dataset using 🤗 Datasets
2. Tokenize data with appropriate tokenizer
3. Load pre-trained model for specific task
4. Configure TrainingArguments
5. Create Trainer with model, data, and compute_metrics
6. Train with `trainer.train()`
7. Evaluate with `trainer.evaluate()`
8. Save model and optionally push to Hub

**See:** `scripts/fine_tune_classifier.py` for complete training example

### Text Generation Workflow
1. Load causal or seq2seq language model
2. Load tokenizer and tokenize prompt
3. Choose generation strategy (greedy, beam search, sampling)
4. Configure generation parameters
5. Generate with `model.generate()`
6. Decode output tokens to text

**See:** `scripts/generate_text.py` for generation strategy examples

## Best Practices

1. **Use Auto classes** for flexibility across different model architectures
2. **Batch processing** for efficiency - process multiple inputs at once
3. **Device management** - use `device_map="auto"` for automatic placement
4. **Memory optimization** - enable fp16/bf16 or quantization for large models
5. **Checkpoint management** - save checkpoints regularly and load best model
6. **Pipeline for quick tasks** - use pipelines for standard inference tasks
7. **Custom metrics** - define compute_metrics for task-specific evaluation
8. **Gradient accumulation** - use for large effective batch sizes on limited memory
9. **Learning rate warmup** - typically 5-10% of total training steps
10. **Hub integration** - push trained models to Hub for sharing and versioning

## Resources

### scripts/
Executable Python scripts demonstrating common Transformers workflows:

- `quick_inference.py` - Pipeline examples for NLP, vision, audio, and multimodal tasks
- `fine_tune_classifier.py` - Complete fine-tuning workflow with Trainer API
- `generate_text.py` - Text generation with various decoding strategies

Run scripts directly to see examples in action:
```bash
python scripts/quick_inference.py
python scripts/fine_tune_classifier.py
python scripts/generate_text.py
```

### references/
Comprehensive reference documentation loaded into context as needed:

- `api_reference.md` - Core classes and APIs (Auto classes, Trainer, GenerationConfig, etc.)
- `pipelines.md` - All available pipelines organized by modality with examples
- `training.md` - Training patterns, TrainingArguments, distributed training, callbacks
- `generation_strategies.md` - Text generation methods, decoding strategies, parameters
- `task_patterns.md` - Complete workflows for common tasks (classification, NER, QA, summarization, etc.)

When working on specific tasks or features, load the relevant reference file for detailed guidance.

## Additional Information

- **Official Documentation**: https://huggingface.co/docs/transformers/index
- **Model Hub**: https://huggingface.co/models (1M+ pre-trained models)
- **Datasets Hub**: https://huggingface.co/datasets
- **Installation**: `pip install transformers datasets evaluate accelerate`
- **GPU Support**: Requires PyTorch or TensorFlow with CUDA
- **Framework Support**: PyTorch (primary), TensorFlow, JAX/Flax
