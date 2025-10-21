# Transformers API Reference

This reference covers the core classes and APIs in the Transformers library.

## Core Auto Classes

Auto classes provide a convenient way to automatically select the appropriate architecture based on model name or checkpoint.

### AutoTokenizer

```python
from transformers import AutoTokenizer

# Load tokenizer
tokenizer = AutoTokenizer.from_pretrained("bert-base-uncased")

# Tokenize single text
encoded = tokenizer("Hello, how are you?")
# Returns: {'input_ids': [...], 'attention_mask': [...]}

# Tokenize with options
encoded = tokenizer(
    "Hello, how are you?",
    padding="max_length",
    truncation=True,
    max_length=512,
    return_tensors="pt"  # "pt" for PyTorch, "tf" for TensorFlow
)

# Tokenize pairs (for classification, QA, etc.)
encoded = tokenizer(
    "Question or sentence A",
    "Context or sentence B",
    padding=True,
    truncation=True
)

# Batch tokenization
texts = ["Text 1", "Text 2", "Text 3"]
encoded = tokenizer(texts, padding=True, truncation=True)

# Decode tokens back to text
text = tokenizer.decode(token_ids, skip_special_tokens=True)

# Batch decode
texts = tokenizer.batch_decode(batch_token_ids, skip_special_tokens=True)
```

**Key Parameters:**
- `padding`: "max_length", "longest", or True (pad to max in batch)
- `truncation`: True or strategy ("longest_first", "only_first", "only_second")
- `max_length`: Maximum sequence length
- `return_tensors`: "pt" (PyTorch), "tf" (TensorFlow), "np" (NumPy)
- `return_attention_mask`: Return attention masks (default True)
- `return_token_type_ids`: Return token type IDs for pairs (default True)
- `add_special_tokens`: Add special tokens like [CLS], [SEP] (default True)

**Special Properties:**
- `tokenizer.vocab_size`: Size of vocabulary
- `tokenizer.pad_token_id`: ID of padding token
- `tokenizer.eos_token_id`: ID of end-of-sequence token
- `tokenizer.bos_token_id`: ID of beginning-of-sequence token
- `tokenizer.unk_token_id`: ID of unknown token

### AutoModel

Base model class that outputs hidden states.

```python
from transformers import AutoModel

model = AutoModel.from_pretrained("bert-base-uncased")

# Forward pass
outputs = model(**inputs)

# Access hidden states
last_hidden_state = outputs.last_hidden_state  # [batch_size, seq_length, hidden_size]
pooler_output = outputs.pooler_output  # [batch_size, hidden_size]

# Get all hidden states
model = AutoModel.from_pretrained("bert-base-uncased", output_hidden_states=True)
outputs = model(**inputs)
all_hidden_states = outputs.hidden_states  # Tuple of tensors
```

### Task-Specific Auto Classes

```python
from transformers import (
    AutoModelForSequenceClassification,
    AutoModelForTokenClassification,
    AutoModelForQuestionAnswering,
    AutoModelForCausalLM,
    AutoModelForMaskedLM,
    AutoModelForSeq2SeqLM,
    AutoModelForImageClassification,
    AutoModelForObjectDetection,
    AutoModelForVision2Seq,
)

# Sequence classification (sentiment, topic, etc.)
model = AutoModelForSequenceClassification.from_pretrained(
    "bert-base-uncased",
    num_labels=3,
    id2label={0: "negative", 1: "neutral", 2: "positive"},
    label2id={"negative": 0, "neutral": 1, "positive": 2}
)

# Token classification (NER, POS tagging)
model = AutoModelForTokenClassification.from_pretrained(
    "bert-base-uncased",
    num_labels=9  # Number of entity types
)

# Question answering
model = AutoModelForQuestionAnswering.from_pretrained("bert-base-uncased")

# Causal language modeling (GPT-style)
model = AutoModelForCausalLM.from_pretrained("gpt2")

# Masked language modeling (BERT-style)
model = AutoModelForMaskedLM.from_pretrained("bert-base-uncased")

# Sequence-to-sequence (T5, BART)
model = AutoModelForSeq2SeqLM.from_pretrained("t5-base")

# Image classification
model = AutoModelForImageClassification.from_pretrained("google/vit-base-patch16-224")

# Object detection
model = AutoModelForObjectDetection.from_pretrained("facebook/detr-resnet-50")

# Vision-to-text (image captioning, VQA)
model = AutoModelForVision2Seq.from_pretrained("microsoft/git-base")
```

### AutoProcessor

For multimodal models that need both text and image processing.

```python
from transformers import AutoProcessor

# For vision-language models
processor = AutoProcessor.from_pretrained("microsoft/git-base")

# Process image and text
from PIL import Image
image = Image.open("image.jpg")
inputs = processor(images=image, text="caption", return_tensors="pt")

# For audio models
processor = AutoProcessor.from_pretrained("openai/whisper-base")
inputs = processor(audio, sampling_rate=16000, return_tensors="pt")
```

### AutoImageProcessor

For vision-only models.

```python
from transformers import AutoImageProcessor

processor = AutoImageProcessor.from_pretrained("google/vit-base-patch16-224")

# Process single image
from PIL import Image
image = Image.open("image.jpg")
inputs = processor(image, return_tensors="pt")

# Batch processing
images = [Image.open(f"image{i}.jpg") for i in range(10)]
inputs = processor(images, return_tensors="pt")
```

## Model Loading Options

### from_pretrained Parameters

```python
model = AutoModel.from_pretrained(
    "model-name",
    # Device and precision
    device_map="auto",  # Automatic device placement
    torch_dtype=torch.float16,  # Use fp16
    low_cpu_mem_usage=True,  # Reduce CPU memory during loading

    # Quantization
    load_in_8bit=True,  # 8-bit quantization
    load_in_4bit=True,  # 4-bit quantization

    # Model configuration
    num_labels=3,  # For classification
    id2label={...},  # Label mapping
    label2id={...},

    # Outputs
    output_hidden_states=True,
    output_attentions=True,

    # Trust remote code
    trust_remote_code=True,  # For custom models

    # Caching
    cache_dir="./cache",
    force_download=False,
    resume_download=True,
)
```

### Quantization with BitsAndBytes

```python
from transformers import BitsAndBytesConfig, AutoModelForCausalLM

# 4-bit quantization
quantization_config = BitsAndBytesConfig(
    load_in_4bit=True,
    bnb_4bit_compute_dtype=torch.float16,
    bnb_4bit_use_double_quant=True,
    bnb_4bit_quant_type="nf4"
)

model = AutoModelForCausalLM.from_pretrained(
    "meta-llama/Llama-2-7b-hf",
    quantization_config=quantization_config,
    device_map="auto"
)
```

## Training Components

### TrainingArguments

See `training.md` for comprehensive coverage. Key parameters:

```python
from transformers import TrainingArguments

args = TrainingArguments(
    output_dir="./results",
    num_train_epochs=3,
    per_device_train_batch_size=16,
    per_device_eval_batch_size=64,
    learning_rate=2e-5,
    weight_decay=0.01,
    eval_strategy="epoch",
    save_strategy="epoch",
    load_best_model_at_end=True,
    metric_for_best_model="accuracy",
    fp16=True,
    logging_steps=100,
    save_total_limit=2,
)
```

### Trainer

```python
from transformers import Trainer

trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=train_dataset,
    eval_dataset=eval_dataset,
    compute_metrics=compute_metrics,
    data_collator=data_collator,
    callbacks=[callback1, callback2],
)

# Train
trainer.train()

# Resume from checkpoint
trainer.train(resume_from_checkpoint=True)

# Evaluate
metrics = trainer.evaluate()

# Predict
predictions = trainer.predict(test_dataset)

# Hyperparameter search
best_trial = trainer.hyperparameter_search(
    direction="maximize",
    backend="optuna",
    n_trials=10,
)

# Save model
trainer.save_model("./final_model")

# Push to Hub
trainer.push_to_hub(commit_message="Training complete")
```

### Data Collators

```python
from transformers import (
    DataCollatorWithPadding,
    DataCollatorForTokenClassification,
    DataCollatorForSeq2Seq,
    DataCollatorForLanguageModeling,
    DefaultDataCollator,
)

# For classification/regression with dynamic padding
data_collator = DataCollatorWithPadding(tokenizer=tokenizer)

# For token classification (NER)
data_collator = DataCollatorForTokenClassification(tokenizer=tokenizer)

# For seq2seq tasks
data_collator = DataCollatorForSeq2Seq(tokenizer=tokenizer, model=model)

# For language modeling
data_collator = DataCollatorForLanguageModeling(
    tokenizer=tokenizer,
    mlm=True,  # True for masked LM, False for causal LM
    mlm_probability=0.15
)

# Default (no special handling)
data_collator = DefaultDataCollator()
```

## Generation Components

### GenerationConfig

See `generation_strategies.md` for comprehensive coverage.

```python
from transformers import GenerationConfig

config = GenerationConfig(
    max_new_tokens=100,
    do_sample=True,
    temperature=0.7,
    top_p=0.9,
    top_k=50,
    num_beams=5,
    repetition_penalty=1.2,
    no_repeat_ngram_size=3,
)

# Use with model
outputs = model.generate(**inputs, generation_config=config)
```

### generate() Method

```python
outputs = model.generate(
    input_ids=inputs.input_ids,
    attention_mask=inputs.attention_mask,
    max_new_tokens=100,
    do_sample=True,
    temperature=0.7,
    top_p=0.9,
    num_return_sequences=3,
    eos_token_id=tokenizer.eos_token_id,
    pad_token_id=tokenizer.pad_token_id,
)
```

## Pipeline API

See `pipelines.md` for comprehensive coverage.

```python
from transformers import pipeline

# Basic usage
pipe = pipeline("task-name", model="model-name", device=0)
results = pipe(inputs)

# With custom model
from transformers import AutoModelForSequenceClassification, AutoTokenizer

model = AutoModelForSequenceClassification.from_pretrained("model-name")
tokenizer = AutoTokenizer.from_pretrained("model-name")
pipe = pipeline("text-classification", model=model, tokenizer=tokenizer)
```

## Configuration Classes

### Model Configuration

```python
from transformers import AutoConfig

# Load configuration
config = AutoConfig.from_pretrained("bert-base-uncased")

# Access configuration
print(config.hidden_size)
print(config.num_attention_heads)
print(config.num_hidden_layers)

# Modify configuration
config.num_labels = 5
config.output_hidden_states = True

# Create model with config
model = AutoModel.from_config(config)

# Save configuration
config.save_pretrained("./config")
```

## Utilities

### Hub Utilities

```python
from huggingface_hub import login, snapshot_download

# Login
login(token="hf_...")

# Download model
snapshot_download(repo_id="model-name", cache_dir="./cache")

# Push to Hub
model.push_to_hub("username/model-name", commit_message="Initial commit")
tokenizer.push_to_hub("username/model-name")
```

### Evaluation Metrics

```python
import evaluate

# Load metric
metric = evaluate.load("accuracy")

# Compute metric
results = metric.compute(predictions=predictions, references=labels)

# Common metrics
accuracy = evaluate.load("accuracy")
precision = evaluate.load("precision")
recall = evaluate.load("recall")
f1 = evaluate.load("f1")
bleu = evaluate.load("bleu")
rouge = evaluate.load("rouge")
```

## Model Outputs

All models return dataclass objects with named attributes:

```python
# Sequence classification output
outputs = model(**inputs)
logits = outputs.logits  # [batch_size, num_labels]
loss = outputs.loss  # If labels provided

# Causal LM output
outputs = model(**inputs)
logits = outputs.logits  # [batch_size, seq_length, vocab_size]
past_key_values = outputs.past_key_values  # KV cache

# Seq2Seq output
outputs = model(**inputs, labels=labels)
loss = outputs.loss
logits = outputs.logits
encoder_last_hidden_state = outputs.encoder_last_hidden_state

# Access as dict
outputs_dict = outputs.to_tuple()  # or dict(outputs)
```

## Best Practices

1. **Use Auto classes**: AutoModel, AutoTokenizer for flexibility
2. **Device management**: Use `device_map="auto"` for multi-GPU
3. **Memory optimization**: Use `torch_dtype=torch.float16` and quantization
4. **Caching**: Set `cache_dir` to avoid re-downloading
5. **Batch processing**: Process multiple inputs at once for efficiency
6. **Trust remote code**: Only set `trust_remote_code=True` for trusted sources
