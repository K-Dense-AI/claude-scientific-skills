# Transformers API Reference

This document provides comprehensive API reference for the most commonly used classes and methods in the Transformers library.

## Core Model Classes

### PreTrainedModel

Base class for all models. Handles loading, saving, and common model operations.

**Key Methods:**

```python
from transformers import PreTrainedModel

# Load pretrained model
model = ModelClass.from_pretrained(
    pretrained_model_name_or_path,
    config=None,                   # Custom config
    cache_dir=None,                # Custom cache location
    force_download=False,          # Force re-download
    resume_download=False,         # Resume interrupted download
    proxies=None,                  # HTTP proxies
    local_files_only=False,        # Only use cached files
    token=None,                    # HF auth token
    revision="main",               # Git branch/tag
    trust_remote_code=False,       # Allow custom model code
    device_map=None,               # Device allocation ("auto", "cpu", "cuda:0", etc.)
    torch_dtype=None,              # Model dtype (torch.float16, "auto", etc.)
    low_cpu_mem_usage=False,       # Reduce CPU memory during loading
    **model_kwargs
)

# Save model
model.save_pretrained(
    save_directory,
    save_config=True,              # Save config.json
    state_dict=None,               # Custom state dict
    save_function=torch.save,      # Custom save function
    push_to_hub=False,             # Upload to Hub
    max_shard_size="5GB",          # Max checkpoint size
    safe_serialization=True,       # Use SafeTensors format
    variant=None,                  # Model variant name
)

# Generate text (for generative models)
outputs = model.generate(
    inputs=None,                   # Input token IDs
    max_length=20,                 # Max total length
    max_new_tokens=None,           # Max new tokens to generate
    min_length=0,                  # Minimum length
    do_sample=False,               # Enable sampling
    early_stopping=False,          # Stop when num_beams finish
    num_beams=1,                   # Beam search width
    temperature=1.0,               # Sampling temperature
    top_k=50,                      # Top-k sampling
    top_p=1.0,                     # Nucleus sampling
    repetition_penalty=1.0,        # Penalize repetition
    length_penalty=1.0,            # Beam search length penalty
    no_repeat_ngram_size=0,        # Block repeated n-grams
    num_return_sequences=1,        # Number of sequences to return
    **model_kwargs
)

# Resize token embeddings (after adding tokens)
new_embeddings = model.resize_token_embeddings(
    new_num_tokens,
    pad_to_multiple_of=None
)

# Utility methods
num_params = model.num_parameters(only_trainable=False)
model.gradient_checkpointing_enable()  # Enable gradient checkpointing
model.enable_input_require_grads()     # For PEFT with frozen models
```

### AutoModel Classes

Automatically instantiate the correct model architecture.

**Available Classes:**

- `AutoModel`: Base model (returns hidden states)
- `AutoModelForCausalLM`: Causal language modeling (GPT-style)
- `AutoModelForMaskedLM`: Masked language modeling (BERT-style)
- `AutoModelForSeq2SeqLM`: Sequence-to-sequence (T5, BART)
- `AutoModelForSequenceClassification`: Text classification
- `AutoModelForTokenClassification`: Token classification (NER)
- `AutoModelForQuestionAnswering`: Extractive QA
- `AutoModelForImageClassification`: Image classification
- `AutoModelForObjectDetection`: Object detection
- `AutoModelForSemanticSegmentation`: Semantic segmentation
- `AutoModelForAudioClassification`: Audio classification
- `AutoModelForSpeechSeq2Seq`: Speech-to-text
- `AutoModelForVision2Seq`: Image captioning, VQA

**Usage:**

```python
from transformers import AutoModel, AutoConfig

# Load with default configuration
model = AutoModel.from_pretrained("bert-base-uncased")

# Load with custom configuration
config = AutoConfig.from_pretrained("bert-base-uncased")
config.hidden_dropout_prob = 0.2
model = AutoModel.from_pretrained("bert-base-uncased", config=config)

# Register custom models
from transformers import AutoConfig, AutoModel

AutoConfig.register("my-model", MyModelConfig)
AutoModel.register(MyModelConfig, MyModel)
```

## Tokenizer Classes

### PreTrainedTokenizer / PreTrainedTokenizerFast

Convert text to token IDs and vice versa.

**Key Methods:**

```python
from transformers import AutoTokenizer

tokenizer = AutoTokenizer.from_pretrained(
    pretrained_model_name_or_path,
    use_fast=True,                 # Use fast (Rust) tokenizer if available
    revision="main",
    **kwargs
)

# Encoding (text → token IDs)
encoded = tokenizer(
    text,                          # String or List[str]
    text_pair=None,                # Second sequence for pairs
    add_special_tokens=True,       # Add [CLS], [SEP], etc.
    padding=False,                 # True, False, "longest", "max_length"
    truncation=False,              # True, False, "longest_first", "only_first", "only_second"
    max_length=None,               # Max sequence length
    stride=0,                      # Overlap for split sequences
    return_tensors=None,           # "pt" (PyTorch), "tf" (TensorFlow), "np" (NumPy)
    return_token_type_ids=None,    # Return token type IDs
    return_attention_mask=None,    # Return attention mask
    return_overflowing_tokens=False,  # Return overflowing tokens
    return_special_tokens_mask=False, # Return special token mask
    return_offsets_mapping=False,  # Return char-level offsets (fast only)
    return_length=False,           # Return sequence lengths
    **kwargs
)

# Decoding (token IDs → text)
text = tokenizer.decode(
    token_ids,
    skip_special_tokens=False,     # Remove special tokens
    clean_up_tokenization_spaces=True,  # Clean up spacing
)

# Batch decoding
texts = tokenizer.batch_decode(
    sequences,
    skip_special_tokens=False,
    clean_up_tokenization_spaces=True,
)

# Tokenization (text → tokens)
tokens = tokenizer.tokenize(text, **kwargs)

# Convert tokens to IDs
ids = tokenizer.convert_tokens_to_ids(tokens)

# Convert IDs to tokens
tokens = tokenizer.convert_ids_to_tokens(ids)

# Add new tokens
num_added = tokenizer.add_tokens(["[NEW_TOKEN1]", "[NEW_TOKEN2]"])

# Add special tokens
tokenizer.add_special_tokens({
    "bos_token": "[BOS]",
    "eos_token": "[EOS]",
    "unk_token": "[UNK]",
    "sep_token": "[SEP]",
    "pad_token": "[PAD]",
    "cls_token": "[CLS]",
    "mask_token": "[MASK]",
    "additional_special_tokens": ["[SPECIAL1]", "[SPECIAL2]"],
})

# Chat template formatting
formatted = tokenizer.apply_chat_template(
    conversation,                  # List[Dict[str, str]] with "role" and "content"
    chat_template=None,            # Custom template
    add_generation_prompt=False,   # Add prompt for model to continue
    tokenize=True,                 # Return token IDs
    padding=False,
    truncation=False,
    max_length=None,
    return_tensors=None,
    return_dict=True,
)

# Save tokenizer
tokenizer.save_pretrained(save_directory)

# Get vocab size
vocab_size = len(tokenizer)

# Get special tokens
pad_token = tokenizer.pad_token
pad_token_id = tokenizer.pad_token_id
# Similar for: bos, eos, unk, sep, cls, mask
```

**Special Token Attributes:**

```python
tokenizer.bos_token         # Beginning of sequence
tokenizer.eos_token         # End of sequence
tokenizer.unk_token         # Unknown token
tokenizer.sep_token         # Separator token
tokenizer.pad_token         # Padding token
tokenizer.cls_token         # Classification token
tokenizer.mask_token        # Mask token

# Corresponding IDs
tokenizer.bos_token_id
tokenizer.eos_token_id
# ... etc
```

## Image Processors

### AutoImageProcessor

Preprocess images for vision models.

**Key Methods:**

```python
from transformers import AutoImageProcessor

processor = AutoImageProcessor.from_pretrained("google/vit-base-patch16-224")

# Process images
inputs = processor(
    images,                        # PIL Image, np.array, torch.Tensor, or List
    return_tensors="pt",           # "pt", "tf", "np", None
    do_resize=True,                # Resize to model size
    size=None,                     # Target size dict
    resample=None,                 # Resampling method
    do_rescale=True,               # Rescale pixel values
    do_normalize=True,             # Normalize with mean/std
    image_mean=None,               # Custom mean
    image_std=None,                # Custom std
    do_center_crop=False,          # Center crop
    crop_size=None,                # Crop size
    **kwargs
)

# Returns: BatchFeature with 'pixel_values' key
```

## Training Components

### TrainingArguments

Configuration for the Trainer class.

**Essential Arguments:**

```python
from transformers import TrainingArguments

args = TrainingArguments(
    # ===== Output & Logging =====
    output_dir="./results",              # REQUIRED: Output directory
    overwrite_output_dir=False,          # Overwrite output directory

    # ===== Training Parameters =====
    num_train_epochs=3.0,                # Number of epochs
    max_steps=-1,                        # Max training steps (overrides epochs)
    per_device_train_batch_size=8,       # Train batch size per device
    per_device_eval_batch_size=8,        # Eval batch size per device
    gradient_accumulation_steps=1,       # Accumulation steps

    # ===== Learning Rate & Optimization =====
    learning_rate=5e-5,                  # Initial learning rate
    weight_decay=0.0,                    # Weight decay
    adam_beta1=0.9,                      # Adam beta1
    adam_beta2=0.999,                    # Adam beta2
    adam_epsilon=1e-8,                   # Adam epsilon
    max_grad_norm=1.0,                   # Gradient clipping
    optim="adamw_torch",                 # Optimizer ("adamw_torch", "adafactor", "adamw_8bit")

    # ===== Learning Rate Scheduler =====
    lr_scheduler_type="linear",          # Scheduler type
    warmup_steps=0,                      # Warmup steps
    warmup_ratio=0.0,                    # Warmup ratio (alternative to steps)

    # ===== Evaluation =====
    eval_strategy="no",                  # "no", "steps", "epoch"
    eval_steps=None,                     # Eval every N steps
    eval_delay=0,                        # Delay first eval
    eval_accumulation_steps=None,        # Accumulate eval outputs

    # ===== Checkpointing =====
    save_strategy="steps",               # "no", "steps", "epoch"
    save_steps=500,                      # Save every N steps
    save_total_limit=None,               # Max checkpoints to keep
    save_safetensors=True,               # Save as SafeTensors
    save_on_each_node=False,             # Save on each node (distributed)

    # ===== Best Model Selection =====
    load_best_model_at_end=False,        # Load best checkpoint at end
    metric_for_best_model=None,          # Metric to use
    greater_is_better=None,              # True if higher is better

    # ===== Logging =====
    logging_dir=None,                    # TensorBoard log directory
    logging_strategy="steps",            # "no", "steps", "epoch"
    logging_steps=500,                   # Log every N steps
    logging_first_step=False,            # Log first step
    logging_nan_inf_filter=True,         # Filter NaN/Inf

    # ===== Mixed Precision =====
    fp16=False,                          # Use fp16 training
    fp16_opt_level="O1",                 # Apex AMP optimization level
    fp16_backend="auto",                 # "auto", "apex", "cpu_amp"
    bf16=False,                          # Use bfloat16 training
    bf16_full_eval=False,                # Use bf16 for evaluation
    tf32=None,                           # Use TF32 (Ampere+ GPUs)

    # ===== Memory Optimization =====
    gradient_checkpointing=False,        # Enable gradient checkpointing
    gradient_checkpointing_kwargs=None,  # Kwargs for gradient checkpointing
    torch_empty_cache_steps=None,        # Clear cache every N steps

    # ===== Distributed Training =====
    local_rank=-1,                       # Local rank for distributed
    ddp_backend=None,                    # "nccl", "gloo", "mpi", "ccl"
    ddp_find_unused_parameters=None,     # Find unused parameters
    ddp_bucket_cap_mb=None,              # DDP bucket size
    fsdp="",                             # FSDP configuration
    fsdp_config=None,                    # FSDP config dict
    deepspeed=None,                      # DeepSpeed config

    # ===== Hub Integration =====
    push_to_hub=False,                   # Push to Hugging Face Hub
    hub_model_id=None,                   # Hub model ID
    hub_strategy="every_save",           # "every_save", "checkpoint", "end"
    hub_token=None,                      # Hub authentication token
    hub_private_repo=False,              # Make repo private

    # ===== Data Handling =====
    dataloader_num_workers=0,            # DataLoader workers
    dataloader_pin_memory=True,          # Pin memory
    dataloader_drop_last=False,          # Drop last incomplete batch
    dataloader_prefetch_factor=None,     # Prefetch factor
    remove_unused_columns=True,          # Remove unused dataset columns
    label_names=None,                    # Label column names

    # ===== Other =====
    seed=42,                             # Random seed
    data_seed=None,                      # Data sampling seed
    jit_mode_eval=False,                 # Use PyTorch JIT for eval
    use_ipex=False,                      # Use Intel Extension for PyTorch
    torch_compile=False,                 # Use torch.compile()
    torch_compile_backend=None,          # Compile backend
    torch_compile_mode=None,             # Compile mode
    include_inputs_for_metrics=False,    # Pass inputs to compute_metrics
    skip_memory_metrics=True,            # Skip memory profiling
)
```

### Trainer

Main training class with full training loop.

**Key Methods:**

```python
from transformers import Trainer

trainer = Trainer(
    model=None,                          # Model to train
    args=None,                           # TrainingArguments
    data_collator=None,                  # Data collator
    train_dataset=None,                  # Training dataset
    eval_dataset=None,                   # Evaluation dataset
    tokenizer=None,                      # Tokenizer
    model_init=None,                     # Function to instantiate model
    compute_metrics=None,                # Function to compute metrics
    callbacks=None,                      # List of callbacks
    optimizers=(None, None),             # (optimizer, scheduler) tuple
    preprocess_logits_for_metrics=None,  # Preprocess logits before metrics
)

# Train model
train_result = trainer.train(
    resume_from_checkpoint=None,         # Resume from checkpoint
    trial=None,                          # Optuna/Ray trial
    ignore_keys_for_eval=None,           # Keys to ignore in eval
)

# Evaluate model
eval_result = trainer.evaluate(
    eval_dataset=None,                   # Eval dataset (default: self.eval_dataset)
    ignore_keys=None,                    # Keys to ignore
    metric_key_prefix="eval",            # Prefix for metric names
)

# Make predictions
predictions = trainer.predict(
    test_dataset,                        # Test dataset
    ignore_keys=None,                    # Keys to ignore
    metric_key_prefix="test",            # Metric prefix
)
# Returns: PredictionOutput(predictions, label_ids, metrics)

# Save model
trainer.save_model(output_dir=None)

# Push to Hub
trainer.push_to_hub(
    commit_message="End of training",
    blocking=True,
    **kwargs
)

# Hyperparameter search
best_trial = trainer.hyperparameter_search(
    hp_space=None,                       # Hyperparameter search space
    compute_objective=None,              # Objective function
    n_trials=20,                         # Number of trials
    direction="minimize",                # "minimize" or "maximize"
    backend=None,                        # "optuna", "ray", "sigopt"
    **kwargs
)

# Create optimizer
optimizer = trainer.create_optimizer()

# Create scheduler
scheduler = trainer.create_scheduler(
    num_training_steps,
    optimizer=None
)

# Log metrics
trainer.log_metrics(split, metrics)
trainer.save_metrics(split, metrics)

# Save checkpoint
trainer.save_state()

# Access current step/epoch
current_step = trainer.state.global_step
current_epoch = trainer.state.epoch

# Access training logs
logs = trainer.state.log_history
```

### Seq2SeqTrainer

Specialized trainer for sequence-to-sequence models.

```python
from transformers import Seq2SeqTrainer, Seq2SeqTrainingArguments

# Use Seq2SeqTrainingArguments with additional parameters
training_args = Seq2SeqTrainingArguments(
    output_dir="./results",
    predict_with_generate=True,          # Use generate() for evaluation
    generation_max_length=None,          # Max length for generation
    generation_num_beams=None,           # Num beams for generation
    **other_training_arguments
)

# Trainer usage is identical to Trainer
trainer = Seq2SeqTrainer(
    model=model,
    args=training_args,
    train_dataset=train_dataset,
    eval_dataset=eval_dataset,
    tokenizer=tokenizer,
    data_collator=data_collator,
    compute_metrics=compute_metrics,
)
```

## Pipeline Classes

### pipeline()

Unified inference API for all tasks.

```python
from transformers import pipeline

pipe = pipeline(
    task=None,                           # Task name (required)
    model=None,                          # Model name/path or model object
    config=None,                         # Model config
    tokenizer=None,                      # Tokenizer
    feature_extractor=None,              # Feature extractor
    image_processor=None,                # Image processor
    framework=None,                      # "pt" or "tf"
    revision=None,                       # Model revision
    use_fast=True,                       # Use fast tokenizer
    token=None,                          # HF token
    device=None,                         # Device (-1 for CPU, 0+ for GPU)
    device_map=None,                     # Device map for multi-GPU
    torch_dtype=None,                    # Model dtype
    trust_remote_code=False,             # Allow custom code
    model_kwargs=None,                   # Additional model kwargs
    pipeline_class=None,                 # Custom pipeline class
    **kwargs
)

# Use pipeline
results = pipe(
    inputs,                              # Input data
    **task_specific_parameters
)
```

## Data Collators

Batch and pad data for training.

```python
from transformers import (
    DataCollatorWithPadding,             # Dynamic padding for classification
    DataCollatorForTokenClassification,  # Padding for token classification
    DataCollatorForSeq2Seq,              # Padding for seq2seq
    DataCollatorForLanguageModeling,     # MLM/CLM data collation
    default_data_collator,               # Simple collator (no padding)
)

# Text classification
data_collator = DataCollatorWithPadding(
    tokenizer=tokenizer,
    padding=True,
    max_length=None,
    pad_to_multiple_of=None,
)

# Token classification
data_collator = DataCollatorForTokenClassification(
    tokenizer=tokenizer,
    padding=True,
    max_length=None,
    pad_to_multiple_of=None,
    label_pad_token_id=-100,
)

# Seq2Seq
data_collator = DataCollatorForSeq2Seq(
    tokenizer=tokenizer,
    model=None,
    padding=True,
    max_length=None,
    pad_to_multiple_of=None,
    label_pad_token_id=-100,
)

# Language modeling
data_collator = DataCollatorForLanguageModeling(
    tokenizer=tokenizer,
    mlm=True,                            # Masked LM (False for causal LM)
    mlm_probability=0.15,                # Mask probability
    pad_to_multiple_of=None,
)
```

## Optimization & Scheduling

```python
from transformers import (
    AdamW,                               # AdamW optimizer
    Adafactor,                           # Adafactor optimizer
    get_scheduler,                       # Get LR scheduler
    get_linear_schedule_with_warmup,
    get_cosine_schedule_with_warmup,
    get_polynomial_decay_schedule_with_warmup,
)

# Create optimizer
optimizer = AdamW(
    model.parameters(),
    lr=5e-5,
    betas=(0.9, 0.999),
    eps=1e-8,
    weight_decay=0.01,
)

# Create scheduler
scheduler = get_scheduler(
    name="linear",                       # "linear", "cosine", "polynomial", "constant"
    optimizer=optimizer,
    num_warmup_steps=0,
    num_training_steps=total_steps,
)

# Or use specific schedulers
scheduler = get_linear_schedule_with_warmup(
    optimizer,
    num_warmup_steps=warmup_steps,
    num_training_steps=total_steps,
)

scheduler = get_cosine_schedule_with_warmup(
    optimizer,
    num_warmup_steps=warmup_steps,
    num_training_steps=total_steps,
    num_cycles=0.5,
)
```

## Configuration Classes

```python
from transformers import AutoConfig

# Load configuration
config = AutoConfig.from_pretrained(
    pretrained_model_name_or_path,
    **kwargs
)

# Common configuration attributes
config.vocab_size                        # Vocabulary size
config.hidden_size                       # Hidden layer size
config.num_hidden_layers                 # Number of layers
config.num_attention_heads               # Attention heads
config.intermediate_size                 # FFN intermediate size
config.hidden_dropout_prob               # Dropout probability
config.attention_probs_dropout_prob      # Attention dropout
config.max_position_embeddings           # Max sequence length

# Save configuration
config.save_pretrained(save_directory)

# Create model from config
from transformers import AutoModel
model = AutoModel.from_config(config)
```

## Utility Functions

```python
from transformers import (
    set_seed,                            # Set random seed
    logging,                             # Logging utilities
)

# Set seed for reproducibility
set_seed(42)

# Configure logging
logging.set_verbosity_info()
logging.set_verbosity_warning()
logging.set_verbosity_error()
logging.set_verbosity_debug()

# Get logger
logger = logging.get_logger(__name__)
```

## Model Outputs

All models return model-specific output classes (subclasses of `ModelOutput`):

```python
# Common output attributes
outputs.loss                             # Loss (if labels provided)
outputs.logits                           # Model logits
outputs.hidden_states                    # All hidden states (if output_hidden_states=True)
outputs.attentions                       # Attention weights (if output_attentions=True)

# Seq2Seq specific
outputs.encoder_last_hidden_state
outputs.encoder_hidden_states
outputs.encoder_attentions
outputs.decoder_hidden_states
outputs.decoder_attentions
outputs.cross_attentions

# Access as dict or tuple
logits = outputs.logits
logits = outputs["logits"]
loss, logits = outputs.to_tuple()[:2]
```

This reference covers the most commonly used API components. For complete documentation, refer to https://huggingface.co/docs/transformers.
