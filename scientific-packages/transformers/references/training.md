# Training with Transformers

Transformers provides comprehensive training capabilities through the `Trainer` API, supporting distributed training, mixed precision, and advanced optimization techniques.

## Basic Training Workflow

```python
from transformers import (
    AutoModelForSequenceClassification,
    AutoTokenizer,
    Trainer,
    TrainingArguments
)
from datasets import load_dataset

# 1. Load and preprocess data
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

# 3. Define training arguments
training_args = TrainingArguments(
    output_dir="./results",
    num_train_epochs=3,
    per_device_train_batch_size=16,
    per_device_eval_batch_size=64,
    learning_rate=2e-5,
    eval_strategy="epoch",
    save_strategy="epoch",
    load_best_model_at_end=True,
)

# 4. Create trainer
trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=tokenized_datasets["train"],
    eval_dataset=tokenized_datasets["test"],
)

# 5. Train
trainer.train()

# 6. Evaluate
trainer.evaluate()

# 7. Save model
trainer.save_model("./final_model")
```

## TrainingArguments Configuration

### Essential Parameters

**Output and Logging:**
- `output_dir`: Directory for checkpoints and outputs (required)
- `logging_dir`: TensorBoard log directory (default: `{output_dir}/runs`)
- `logging_steps`: Log every N steps (default: 500)
- `logging_strategy`: "steps" or "epoch"

**Training Duration:**
- `num_train_epochs`: Number of epochs (default: 3.0)
- `max_steps`: Max training steps (overrides num_train_epochs if set)

**Batch Size and Gradient Accumulation:**
- `per_device_train_batch_size`: Batch size per device (default: 8)
- `per_device_eval_batch_size`: Eval batch size per device (default: 8)
- `gradient_accumulation_steps`: Accumulate gradients over N steps (default: 1)
- Effective batch size = `per_device_train_batch_size * gradient_accumulation_steps * num_gpus`

**Learning Rate:**
- `learning_rate`: Peak learning rate (default: 5e-5)
- `lr_scheduler_type`: Scheduler type ("linear", "cosine", "constant", etc.)
- `warmup_steps`: Warmup steps (default: 0)
- `warmup_ratio`: Warmup as fraction of total steps

**Evaluation:**
- `eval_strategy`: "no", "steps", or "epoch" (default: "no")
- `eval_steps`: Evaluate every N steps (if eval_strategy="steps")
- `eval_delay`: Delay evaluation until N steps

**Checkpointing:**
- `save_strategy`: "no", "steps", or "epoch" (default: "steps")
- `save_steps`: Save checkpoint every N steps (default: 500)
- `save_total_limit`: Keep only N most recent checkpoints
- `load_best_model_at_end`: Load best checkpoint at end (default: False)
- `metric_for_best_model`: Metric to determine best model

**Optimization:**
- `optim`: Optimizer ("adamw_torch", "adamw_hf", "sgd", etc.)
- `weight_decay`: Weight decay coefficient (default: 0.0)
- `adam_beta1`, `adam_beta2`: Adam optimizer betas
- `adam_epsilon`: Epsilon for Adam (default: 1e-8)
- `max_grad_norm`: Max gradient norm for clipping (default: 1.0)

### Mixed Precision Training

```python
training_args = TrainingArguments(
    output_dir="./results",
    fp16=True,  # Use fp16 on NVIDIA GPUs
    fp16_opt_level="O1",  # O0, O1, O2, O3 (Apex levels)
    # or
    bf16=True,  # Use bf16 on Ampere+ GPUs (better than fp16)
)
```

### Distributed Training

**DataParallel (single-node multi-GPU):**
```python
# Automatic with multiple GPUs
training_args = TrainingArguments(
    output_dir="./results",
    per_device_train_batch_size=16,  # Per GPU
)
# Run: python script.py
```

**DistributedDataParallel (multi-node or multi-GPU):**
```bash
# Single node, multiple GPUs
python -m torch.distributed.launch --nproc_per_node=4 script.py

# Or use accelerate
accelerate config
accelerate launch script.py
```

**DeepSpeed Integration:**
```python
training_args = TrainingArguments(
    output_dir="./results",
    deepspeed="ds_config.json",  # DeepSpeed config file
)
```

### Advanced Features

**Gradient Checkpointing (reduce memory):**
```python
training_args = TrainingArguments(
    output_dir="./results",
    gradient_checkpointing=True,
)
```

**Compilation with torch.compile:**
```python
training_args = TrainingArguments(
    output_dir="./results",
    torch_compile=True,
    torch_compile_backend="inductor",  # or "cudagraphs"
)
```

**Push to Hub:**
```python
training_args = TrainingArguments(
    output_dir="./results",
    push_to_hub=True,
    hub_model_id="username/model-name",
    hub_strategy="every_save",  # or "end"
)
```

## Custom Training Components

### Custom Metrics

```python
import evaluate
import numpy as np

metric = evaluate.load("accuracy")

def compute_metrics(eval_pred):
    logits, labels = eval_pred
    predictions = np.argmax(logits, axis=-1)
    return metric.compute(predictions=predictions, references=labels)

trainer = Trainer(
    model=model,
    args=training_args,
    compute_metrics=compute_metrics,
)
```

### Custom Loss Function

```python
class CustomTrainer(Trainer):
    def compute_loss(self, model, inputs, return_outputs=False):
        labels = inputs.pop("labels")
        outputs = model(**inputs)
        logits = outputs.logits

        # Custom loss calculation
        loss_fct = torch.nn.CrossEntropyLoss(weight=class_weights)
        loss = loss_fct(logits.view(-1, self.model.config.num_labels), labels.view(-1))

        return (loss, outputs) if return_outputs else loss
```

### Data Collator

```python
from transformers import DataCollatorWithPadding

data_collator = DataCollatorWithPadding(tokenizer=tokenizer)

trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=train_dataset,
    data_collator=data_collator,
)
```

### Callbacks

```python
from transformers import TrainerCallback

class CustomCallback(TrainerCallback):
    def on_epoch_end(self, args, state, control, **kwargs):
        print(f"Epoch {state.epoch} completed!")
        return control

trainer = Trainer(
    model=model,
    args=training_args,
    callbacks=[CustomCallback],
)
```

## Hyperparameter Search

```python
def model_init():
    return AutoModelForSequenceClassification.from_pretrained(
        "bert-base-uncased",
        num_labels=2
    )

trainer = Trainer(
    model_init=model_init,
    args=training_args,
    train_dataset=train_dataset,
    eval_dataset=eval_dataset,
    compute_metrics=compute_metrics,
)

# Optuna-based search
best_trial = trainer.hyperparameter_search(
    direction="maximize",
    backend="optuna",
    n_trials=10,
    hp_space=lambda trial: {
        "learning_rate": trial.suggest_float("learning_rate", 1e-5, 5e-5, log=True),
        "per_device_train_batch_size": trial.suggest_categorical("per_device_train_batch_size", [8, 16, 32]),
        "num_train_epochs": trial.suggest_int("num_train_epochs", 2, 5),
    }
)
```

## Training Best Practices

1. **Start with small learning rates**: 2e-5 to 5e-5 for fine-tuning
2. **Use warmup**: 5-10% of total steps for learning rate warmup
3. **Monitor training**: Use eval_strategy="epoch" or "steps" to track progress
4. **Save checkpoints**: Set save_strategy and save_total_limit
5. **Use mixed precision**: Enable fp16 or bf16 for faster training
6. **Gradient accumulation**: For large effective batch sizes on limited memory
7. **Load best model**: Set load_best_model_at_end=True to avoid overfitting
8. **Push to Hub**: Enable push_to_hub for easy model sharing and versioning

## Common Training Patterns

### Classification
```python
model = AutoModelForSequenceClassification.from_pretrained(
    "bert-base-uncased",
    num_labels=num_classes,
    id2label=id2label,
    label2id=label2id
)
```

### Question Answering
```python
model = AutoModelForQuestionAnswering.from_pretrained("bert-base-uncased")
```

### Token Classification (NER)
```python
model = AutoModelForTokenClassification.from_pretrained(
    "bert-base-uncased",
    num_labels=num_tags,
    id2label=id2label,
    label2id=label2id
)
```

### Sequence-to-Sequence
```python
model = AutoModelForSeq2SeqLM.from_pretrained("t5-base")
```

### Causal Language Modeling
```python
model = AutoModelForCausalLM.from_pretrained("gpt2")
```

### Masked Language Modeling
```python
model = AutoModelForMaskedLM.from_pretrained("bert-base-uncased")
```
