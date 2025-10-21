# Common Task Patterns

This document provides common patterns and workflows for typical tasks using Transformers.

## Text Classification

### Binary or Multi-class Classification

```python
from transformers import (
    AutoTokenizer,
    AutoModelForSequenceClassification,
    TrainingArguments,
    Trainer
)
from datasets import load_dataset
import evaluate
import numpy as np

# Load dataset
dataset = load_dataset("imdb")

# Tokenize
tokenizer = AutoTokenizer.from_pretrained("bert-base-uncased")

def tokenize_function(examples):
    return tokenizer(examples["text"], padding="max_length", truncation=True, max_length=512)

tokenized_datasets = dataset.map(tokenize_function, batched=True)

# Load model
id2label = {0: "negative", 1: "positive"}
label2id = {"negative": 0, "positive": 1}

model = AutoModelForSequenceClassification.from_pretrained(
    "bert-base-uncased",
    num_labels=2,
    id2label=id2label,
    label2id=label2id
)

# Metrics
metric = evaluate.load("accuracy")

def compute_metrics(eval_pred):
    logits, labels = eval_pred
    predictions = np.argmax(logits, axis=-1)
    return metric.compute(predictions=predictions, references=labels)

# Train
training_args = TrainingArguments(
    output_dir="./results",
    eval_strategy="epoch",
    learning_rate=2e-5,
    per_device_train_batch_size=16,
    per_device_eval_batch_size=64,
    num_train_epochs=3,
    weight_decay=0.01,
    load_best_model_at_end=True,
)

trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=tokenized_datasets["train"],
    eval_dataset=tokenized_datasets["test"],
    compute_metrics=compute_metrics,
)

trainer.train()

# Inference
text = "This movie was fantastic!"
inputs = tokenizer(text, return_tensors="pt")
outputs = model(**inputs)
predictions = outputs.logits.argmax(-1)
print(id2label[predictions.item()])
```

## Named Entity Recognition (Token Classification)

```python
from transformers import (
    AutoTokenizer,
    AutoModelForTokenClassification,
    TrainingArguments,
    Trainer,
    DataCollatorForTokenClassification
)
from datasets import load_dataset
import evaluate
import numpy as np

# Load dataset
dataset = load_dataset("conll2003")

# Tokenize (align labels with tokenized words)
tokenizer = AutoTokenizer.from_pretrained("bert-base-uncased")

def tokenize_and_align_labels(examples):
    tokenized_inputs = tokenizer(
        examples["tokens"],
        truncation=True,
        is_split_into_words=True
    )

    labels = []
    for i, label in enumerate(examples["ner_tags"]):
        word_ids = tokenized_inputs.word_ids(batch_index=i)
        label_ids = []
        previous_word_idx = None
        for word_idx in word_ids:
            if word_idx is None:
                label_ids.append(-100)
            elif word_idx != previous_word_idx:
                label_ids.append(label[word_idx])
            else:
                label_ids.append(-100)
            previous_word_idx = word_idx
        labels.append(label_ids)

    tokenized_inputs["labels"] = labels
    return tokenized_inputs

tokenized_datasets = dataset.map(tokenize_and_align_labels, batched=True)

# Model
label_list = dataset["train"].features["ner_tags"].feature.names
model = AutoModelForTokenClassification.from_pretrained(
    "bert-base-uncased",
    num_labels=len(label_list)
)

# Data collator
data_collator = DataCollatorForTokenClassification(tokenizer)

# Metrics
metric = evaluate.load("seqeval")

def compute_metrics(eval_pred):
    predictions, labels = eval_pred
    predictions = np.argmax(predictions, axis=2)

    true_labels = [[label_list[l] for l in label if l != -100] for label in labels]
    true_predictions = [
        [label_list[p] for (p, l) in zip(prediction, label) if l != -100]
        for prediction, label in zip(predictions, labels)
    ]

    return metric.compute(predictions=true_predictions, references=true_labels)

# Train
training_args = TrainingArguments(
    output_dir="./results",
    eval_strategy="epoch",
    learning_rate=2e-5,
    per_device_train_batch_size=16,
    num_train_epochs=3,
    weight_decay=0.01,
)

trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=tokenized_datasets["train"],
    eval_dataset=tokenized_datasets["validation"],
    data_collator=data_collator,
    compute_metrics=compute_metrics,
)

trainer.train()
```

## Question Answering

```python
from transformers import (
    AutoTokenizer,
    AutoModelForQuestionAnswering,
    TrainingArguments,
    Trainer,
    DefaultDataCollator
)
from datasets import load_dataset

# Load dataset
dataset = load_dataset("squad")

# Tokenize
tokenizer = AutoTokenizer.from_pretrained("bert-base-uncased")

def preprocess_function(examples):
    questions = [q.strip() for q in examples["question"]]
    inputs = tokenizer(
        questions,
        examples["context"],
        max_length=384,
        truncation="only_second",
        return_offsets_mapping=True,
        padding="max_length",
    )

    offset_mapping = inputs.pop("offset_mapping")
    answers = examples["answers"]
    start_positions = []
    end_positions = []

    for i, offset in enumerate(offset_mapping):
        answer = answers[i]
        start_char = answer["answer_start"][0]
        end_char = start_char + len(answer["text"][0])

        # Find start and end token positions
        sequence_ids = inputs.sequence_ids(i)
        context_start = sequence_ids.index(1)
        context_end = len(sequence_ids) - 1 - sequence_ids[::-1].index(1)

        if offset[context_start][0] > end_char or offset[context_end][1] < start_char:
            start_positions.append(0)
            end_positions.append(0)
        else:
            idx = context_start
            while idx <= context_end and offset[idx][0] <= start_char:
                idx += 1
            start_positions.append(idx - 1)

            idx = context_end
            while idx >= context_start and offset[idx][1] >= end_char:
                idx -= 1
            end_positions.append(idx + 1)

    inputs["start_positions"] = start_positions
    inputs["end_positions"] = end_positions
    return inputs

tokenized_datasets = dataset.map(preprocess_function, batched=True, remove_columns=dataset["train"].column_names)

# Model
model = AutoModelForQuestionAnswering.from_pretrained("bert-base-uncased")

# Train
training_args = TrainingArguments(
    output_dir="./results",
    eval_strategy="epoch",
    learning_rate=2e-5,
    per_device_train_batch_size=16,
    num_train_epochs=3,
)

trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=tokenized_datasets["train"],
    eval_dataset=tokenized_datasets["validation"],
    data_collator=DefaultDataCollator(),
)

trainer.train()

# Inference
question = "What is the capital of France?"
context = "Paris is the capital and most populous city of France."
inputs = tokenizer(question, context, return_tensors="pt")
outputs = model(**inputs)

start_pos = outputs.start_logits.argmax()
end_pos = outputs.end_logits.argmax()
answer_tokens = inputs.input_ids[0][start_pos:end_pos+1]
answer = tokenizer.decode(answer_tokens)
```

## Text Summarization

```python
from transformers import (
    AutoTokenizer,
    AutoModelForSeq2SeqLM,
    TrainingArguments,
    Trainer,
    DataCollatorForSeq2Seq
)
from datasets import load_dataset
import evaluate

# Load dataset
dataset = load_dataset("cnn_dailymail", "3.0.0")

# Tokenize
tokenizer = AutoTokenizer.from_pretrained("t5-small")

def preprocess_function(examples):
    inputs = ["summarize: " + doc for doc in examples["article"]]
    model_inputs = tokenizer(inputs, max_length=1024, truncation=True)

    labels = tokenizer(
        text_target=examples["highlights"],
        max_length=128,
        truncation=True
    )

    model_inputs["labels"] = labels["input_ids"]
    return model_inputs

tokenized_datasets = dataset.map(preprocess_function, batched=True)

# Model
model = AutoModelForSeq2SeqLM.from_pretrained("t5-small")

# Data collator
data_collator = DataCollatorForSeq2Seq(tokenizer, model=model)

# Metrics
rouge = evaluate.load("rouge")

def compute_metrics(eval_pred):
    predictions, labels = eval_pred
    decoded_preds = tokenizer.batch_decode(predictions, skip_special_tokens=True)
    labels = np.where(labels != -100, labels, tokenizer.pad_token_id)
    decoded_labels = tokenizer.batch_decode(labels, skip_special_tokens=True)

    result = rouge.compute(
        predictions=decoded_preds,
        references=decoded_labels,
        use_stemmer=True
    )

    return {k: round(v, 4) for k, v in result.items()}

# Train
training_args = TrainingArguments(
    output_dir="./results",
    eval_strategy="epoch",
    learning_rate=2e-5,
    per_device_train_batch_size=8,
    per_device_eval_batch_size=8,
    num_train_epochs=3,
    predict_with_generate=True,
)

trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=tokenized_datasets["train"],
    eval_dataset=tokenized_datasets["validation"],
    data_collator=data_collator,
    compute_metrics=compute_metrics,
)

trainer.train()

# Inference
text = "Long article text..."
inputs = tokenizer("summarize: " + text, return_tensors="pt", max_length=1024, truncation=True)
outputs = model.generate(**inputs, max_length=128, num_beams=4)
summary = tokenizer.decode(outputs[0], skip_special_tokens=True)
```

## Translation

```python
from transformers import (
    AutoTokenizer,
    AutoModelForSeq2SeqLM,
    TrainingArguments,
    Trainer,
    DataCollatorForSeq2Seq
)
from datasets import load_dataset

# Load dataset
dataset = load_dataset("wmt16", "de-en")

# Tokenize
tokenizer = AutoTokenizer.from_pretrained("t5-small")

def preprocess_function(examples):
    inputs = [f"translate German to English: {de}" for de in examples["de"]]
    model_inputs = tokenizer(inputs, max_length=128, truncation=True)

    labels = tokenizer(
        text_target=examples["en"],
        max_length=128,
        truncation=True
    )

    model_inputs["labels"] = labels["input_ids"]
    return model_inputs

tokenized_datasets = dataset.map(preprocess_function, batched=True)

# Model and training (similar to summarization)
model = AutoModelForSeq2SeqLM.from_pretrained("t5-small")

# Inference
text = "Guten Tag, wie geht es Ihnen?"
inputs = tokenizer(f"translate German to English: {text}", return_tensors="pt")
outputs = model.generate(**inputs, max_length=128)
translation = tokenizer.decode(outputs[0], skip_special_tokens=True)
```

## Causal Language Modeling (Training from Scratch or Fine-tuning)

```python
from transformers import (
    AutoTokenizer,
    AutoModelForCausalLM,
    TrainingArguments,
    Trainer,
    DataCollatorForLanguageModeling
)
from datasets import load_dataset

# Load dataset
dataset = load_dataset("wikitext", "wikitext-2-raw-v1")

# Tokenize
tokenizer = AutoTokenizer.from_pretrained("gpt2")
tokenizer.pad_token = tokenizer.eos_token

def tokenize_function(examples):
    return tokenizer(examples["text"], truncation=True, max_length=512)

tokenized_datasets = dataset.map(tokenize_function, batched=True, remove_columns=["text"])

# Group texts into chunks
block_size = 128

def group_texts(examples):
    concatenated_examples = {k: sum(examples[k], []) for k in examples.keys()}
    total_length = len(concatenated_examples[list(examples.keys())[0]])
    total_length = (total_length // block_size) * block_size
    result = {
        k: [t[i:i + block_size] for i in range(0, total_length, block_size)]
        for k, t in concatenated_examples.items()
    }
    result["labels"] = result["input_ids"].copy()
    return result

lm_datasets = tokenized_datasets.map(group_texts, batched=True)

# Model
model = AutoModelForCausalLM.from_pretrained("gpt2")

# Data collator
data_collator = DataCollatorForLanguageModeling(tokenizer=tokenizer, mlm=False)

# Train
training_args = TrainingArguments(
    output_dir="./results",
    eval_strategy="epoch",
    learning_rate=2e-5,
    per_device_train_batch_size=8,
    num_train_epochs=3,
)

trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=lm_datasets["train"],
    eval_dataset=lm_datasets["validation"],
    data_collator=data_collator,
)

trainer.train()
```

## Image Classification

```python
from transformers import (
    AutoImageProcessor,
    AutoModelForImageClassification,
    TrainingArguments,
    Trainer
)
from datasets import load_dataset
from torchvision.transforms import Compose, Resize, ToTensor, Normalize
import numpy as np
import evaluate

# Load dataset
dataset = load_dataset("food101", split="train[:5000]")

# Prepare image transforms
image_processor = AutoImageProcessor.from_pretrained("google/vit-base-patch16-224")

normalize = Normalize(mean=image_processor.image_mean, std=image_processor.image_std)
size = image_processor.size["height"]

transforms = Compose([
    Resize((size, size)),
    ToTensor(),
    normalize,
])

def preprocess_function(examples):
    examples["pixel_values"] = [transforms(img.convert("RGB")) for img in examples["image"]]
    return examples

dataset = dataset.with_transform(preprocess_function)

# Model
model = AutoModelForImageClassification.from_pretrained(
    "google/vit-base-patch16-224",
    num_labels=len(dataset["train"].features["label"].names),
    ignore_mismatched_sizes=True
)

# Metrics
metric = evaluate.load("accuracy")

def compute_metrics(eval_pred):
    predictions = np.argmax(eval_pred.predictions, axis=1)
    return metric.compute(predictions=predictions, references=eval_pred.label_ids)

# Train
training_args = TrainingArguments(
    output_dir="./results",
    eval_strategy="epoch",
    learning_rate=5e-5,
    per_device_train_batch_size=16,
    num_train_epochs=3,
)

trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=dataset["train"],
    eval_dataset=dataset["validation"],
    compute_metrics=compute_metrics,
)

trainer.train()
```

## Vision-Language Tasks (Image Captioning)

```python
from transformers import (
    AutoProcessor,
    AutoModelForVision2Seq,
    TrainingArguments,
    Trainer
)
from datasets import load_dataset
from PIL import Image

# Load dataset
dataset = load_dataset("ybelkada/football-dataset")

# Processor
processor = AutoProcessor.from_pretrained("microsoft/git-base")

def preprocess_function(examples):
    images = [Image.open(img).convert("RGB") for img in examples["image"]]
    texts = examples["caption"]

    inputs = processor(images=images, text=texts, padding="max_length", truncation=True)
    inputs["labels"] = inputs["input_ids"]
    return inputs

dataset = dataset.map(preprocess_function, batched=True)

# Model
model = AutoModelForVision2Seq.from_pretrained("microsoft/git-base")

# Train
training_args = TrainingArguments(
    output_dir="./results",
    eval_strategy="epoch",
    learning_rate=5e-5,
    per_device_train_batch_size=8,
    num_train_epochs=3,
)

trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=dataset["train"],
    eval_dataset=dataset["test"],
)

trainer.train()

# Inference
image = Image.open("image.jpg")
inputs = processor(images=image, return_tensors="pt")
outputs = model.generate(**inputs)
caption = processor.decode(outputs[0], skip_special_tokens=True)
```

## Best Practices Summary

1. **Use appropriate Auto* classes**: AutoTokenizer, AutoModel, etc. for model loading
2. **Proper preprocessing**: Tokenize, align labels, handle special cases
3. **Data collators**: Use appropriate collators for dynamic padding
4. **Metrics**: Load and compute relevant metrics for evaluation
5. **Training arguments**: Configure properly for task and hardware
6. **Inference**: Use pipeline() for quick inference, or manual tokenization for custom needs
