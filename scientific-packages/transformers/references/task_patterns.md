# Task-Specific Patterns

Quick reference for implementing common tasks with Transformers. Each pattern includes the complete workflow from data loading to inference.

## Text Classification

Classify text into predefined categories (sentiment, topic, intent, etc.).

```python
from transformers import (
    AutoTokenizer, AutoModelForSequenceClassification,
    TrainingArguments, Trainer, DataCollatorWithPadding
)
from datasets import load_dataset

# 1. Load data
dataset = load_dataset("imdb")

# 2. Preprocess
tokenizer = AutoTokenizer.from_pretrained("bert-base-uncased")

def preprocess(examples):
    return tokenizer(examples["text"], truncation=True, max_length=512)

tokenized = dataset.map(preprocess, batched=True)

# 3. Model
model = AutoModelForSequenceClassification.from_pretrained(
    "bert-base-uncased",
    num_labels=2,
    id2label={0: "negative", 1: "positive"},
    label2id={"negative": 0, "positive": 1}
)

# 4. Train
training_args = TrainingArguments(
    output_dir="./results",
    learning_rate=2e-5,
    per_device_train_batch_size=16,
    num_train_epochs=3,
    eval_strategy="epoch",
)

trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=tokenized["train"],
    eval_dataset=tokenized["test"],
    tokenizer=tokenizer,
    data_collator=DataCollatorWithPadding(tokenizer=tokenizer),
)

trainer.train()

# 5. Inference
text = "This movie was fantastic!"
inputs = tokenizer(text, return_tensors="pt")
outputs = model(**inputs)
predictions = outputs.logits.argmax(-1)
print(model.config.id2label[predictions.item()])  # "positive"
```

## Token Classification (NER)

Label each token in text (named entities, POS tags, etc.).

```python
from transformers import AutoTokenizer, AutoModelForTokenClassification
from datasets import load_dataset

# Load data (tokens and NER tags)
dataset = load_dataset("conll2003")

tokenizer = AutoTokenizer.from_pretrained("bert-base-cased")

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
                label_ids.append(-100)  # Special tokens
            elif word_idx != previous_word_idx:
                label_ids.append(label[word_idx])
            else:
                label_ids.append(-100)  # Subword tokens
            previous_word_idx = word_idx
        labels.append(label_ids)

    tokenized_inputs["labels"] = labels
    return tokenized_inputs

tokenized_dataset = dataset.map(tokenize_and_align_labels, batched=True)

# Model
label_list = dataset["train"].features["ner_tags"].feature.names
model = AutoModelForTokenClassification.from_pretrained(
    "bert-base-cased",
    num_labels=len(label_list),
    id2label={i: label for i, label in enumerate(label_list)},
    label2id={label: i for i, label in enumerate(label_list)}
)

# Training similar to classification
# ... (use Trainer with DataCollatorForTokenClassification)
```

## Question Answering (Extractive)

Extract answer spans from context.

```python
from transformers import AutoTokenizer, AutoModelForQuestionAnswering

tokenizer = AutoTokenizer.from_pretrained("distilbert-base-cased-distilled-squad")
model = AutoModelForQuestionAnswering.from_pretrained("distilbert-base-cased-distilled-squad")

question = "What is the capital of France?"
context = "Paris is the capital and most populous city of France."

inputs = tokenizer(question, context, return_tensors="pt")
outputs = model(**inputs)

# Get answer span
answer_start = outputs.start_logits.argmax()
answer_end = outputs.end_logits.argmax() + 1
answer = tokenizer.decode(inputs["input_ids"][0][answer_start:answer_end])
print(answer)  # "Paris"
```

## Text Generation

Generate text continuations.

```python
from transformers import AutoModelForCausalLM, AutoTokenizer

model = AutoModelForCausalLM.from_pretrained("gpt2")
tokenizer = AutoTokenizer.from_pretrained("gpt2")

prompt = "In the future, artificial intelligence will"
inputs = tokenizer(prompt, return_tensors="pt")

outputs = model.generate(
    **inputs,
    max_new_tokens=100,
    do_sample=True,
    temperature=0.8,
    top_p=0.95,
    repetition_penalty=1.2,
)

generated_text = tokenizer.decode(outputs[0], skip_special_tokens=True)
print(generated_text)
```

## Summarization

Condense long text into summaries.

```python
from transformers import (
    AutoTokenizer, AutoModelForSeq2SeqLM,
    Seq2SeqTrainingArguments, Seq2SeqTrainer,
    DataCollatorForSeq2Seq
)

tokenizer = AutoTokenizer.from_pretrained("t5-small")
model = AutoModelForSeq2SeqLM.from_pretrained("t5-small")

def preprocess(examples):
    inputs = ["summarize: " + doc for doc in examples["document"]]
    model_inputs = tokenizer(inputs, max_length=1024, truncation=True)

    labels = tokenizer(
        examples["summary"],
        max_length=128,
        truncation=True
    )
    model_inputs["labels"] = labels["input_ids"]
    return model_inputs

tokenized_dataset = dataset.map(preprocess, batched=True)

# Training
training_args = Seq2SeqTrainingArguments(
    output_dir="./results",
    predict_with_generate=True,  # Important for seq2seq
    eval_strategy="epoch",
    learning_rate=2e-5,
    per_device_train_batch_size=8,
    num_train_epochs=3,
)

trainer = Seq2SeqTrainer(
    model=model,
    args=training_args,
    train_dataset=tokenized_dataset["train"],
    eval_dataset=tokenized_dataset["validation"],
    tokenizer=tokenizer,
    data_collator=DataCollatorForSeq2Seq(tokenizer=tokenizer),
)

trainer.train()

# Inference
text = "Long article text here..."
inputs = tokenizer("summarize: " + text, return_tensors="pt", max_length=1024, truncation=True)
outputs = model.generate(**inputs, max_length=150, num_beams=4, early_stopping=True)
summary = tokenizer.decode(outputs[0], skip_special_tokens=True)
```

## Translation

Translate text between languages.

```python
from transformers import pipeline

translator = pipeline("translation_en_to_fr", model="Helsinki-NLP/opus-mt-en-fr")
result = translator("Hello, how are you?")
print(result[0]["translation_text"])  # "Bonjour, comment allez-vous?"

# For fine-tuning, similar to summarization with Seq2SeqTrainer
```

## Image Classification

Classify images into categories.

```python
from transformers import (
    AutoImageProcessor, AutoModelForImageClassification,
    TrainingArguments, Trainer
)
from datasets import load_dataset
from PIL import Image

# Load data
dataset = load_dataset("food101", split="train[:1000]")

# Preprocess
processor = AutoImageProcessor.from_pretrained("google/vit-base-patch16-224")

def transform(examples):
    examples["pixel_values"] = [
        processor(img.convert("RGB"), return_tensors="pt")["pixel_values"][0]
        for img in examples["image"]
    ]
    return examples

dataset = dataset.with_transform(transform)

# Model
model = AutoModelForImageClassification.from_pretrained(
    "google/vit-base-patch16-224",
    num_labels=101,
    ignore_mismatched_sizes=True
)

# Training
training_args = TrainingArguments(
    output_dir="./results",
    remove_unused_columns=False,  # Keep image data
    eval_strategy="epoch",
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

# Inference
image = Image.open("food.jpg")
inputs = processor(image, return_tensors="pt")
outputs = model(**inputs)
predicted_class = outputs.logits.argmax(-1).item()
print(model.config.id2label[predicted_class])
```

## Object Detection

Detect and localize objects in images.

```python
from transformers import pipeline
from PIL import Image

detector = pipeline("object-detection", model="facebook/detr-resnet-50")

image = Image.open("street.jpg")
results = detector(image)

for result in results:
    print(f"{result['label']}: {result['score']:.2f} at {result['box']}")
    # car: 0.98 at {'xmin': 123, 'ymin': 456, 'xmax': 789, 'ymax': 1011}
```

## Image Segmentation

Segment images into regions.

```python
from transformers import pipeline

segmenter = pipeline("image-segmentation", model="facebook/detr-resnet-50-panoptic")

image = "path/to/image.jpg"
segments = segmenter(image)

for segment in segments:
    print(f"{segment['label']}: {segment['score']:.2f}")
    # Access mask: segment['mask']
```

## Image Captioning

Generate textual descriptions of images.

```python
from transformers import AutoProcessor, AutoModelForCausalLM
from PIL import Image

processor = AutoProcessor.from_pretrained("microsoft/git-base")
model = AutoModelForCausalLM.from_pretrained("microsoft/git-base")

image = Image.open("photo.jpg")
inputs = processor(images=image, return_tensors="pt")

outputs = model.generate(**inputs, max_length=50)
caption = processor.batch_decode(outputs, skip_special_tokens=True)[0]
print(caption)  # "a dog sitting on grass"
```

## Speech Recognition (ASR)

Transcribe speech to text.

```python
from transformers import pipeline

transcriber = pipeline(
    "automatic-speech-recognition",
    model="openai/whisper-base"
)

result = transcriber("audio.mp3")
print(result["text"])  # "Hello, this is a test."

# With timestamps
result = transcriber("audio.mp3", return_timestamps=True)
for chunk in result["chunks"]:
    print(f"[{chunk['timestamp'][0]:.1f}s - {chunk['timestamp'][1]:.1f}s]: {chunk['text']}")
```

## Text-to-Speech

Generate speech from text.

```python
from transformers import pipeline

synthesizer = pipeline("text-to-speech", model="microsoft/speecht5_tts")

result = synthesizer("Hello, how are you today?")
# result["audio"] contains the waveform
# result["sampling_rate"] contains the sample rate

# Save audio
import scipy
scipy.io.wavfile.write("output.wav", result["sampling_rate"], result["audio"][0])
```

## Visual Question Answering

Answer questions about images.

```python
from transformers import pipeline
from PIL import Image

vqa = pipeline("visual-question-answering", model="dandelin/vilt-b32-finetuned-vqa")

image = Image.open("photo.jpg")
question = "What color is the car?"

result = vqa(image=image, question=question)
print(result[0]["answer"])  # "red"
```

## Document Question Answering

Extract information from documents (PDFs, images with text).

```python
from transformers import pipeline

doc_qa = pipeline("document-question-answering", model="impira/layoutlm-document-qa")

result = doc_qa(
    image="invoice.png",
    question="What is the total amount?"
)

print(result["answer"])  # "$1,234.56"
```

## Zero-Shot Classification

Classify without training data.

```python
from transformers import pipeline

classifier = pipeline("zero-shot-classification", model="facebook/bart-large-mnli")

text = "This is a delicious Italian restaurant with great pasta."
candidate_labels = ["food", "travel", "technology", "sports"]

result = classifier(text, candidate_labels)
print(result["labels"][0])  # "food"
print(result["scores"][0])  # 0.95
```

## Few-Shot Learning with LLMs

Use large language models for few-shot tasks.

```python
from transformers import AutoModelForCausalLM, AutoTokenizer

model = AutoModelForCausalLM.from_pretrained("meta-llama/Llama-2-7b-hf")
tokenizer = AutoTokenizer.from_pretrained("meta-llama/Llama-2-7b-hf")

# Few-shot prompt
prompt = """
Classify the sentiment: positive, negative, or neutral.

Text: "I love this product!"
Sentiment: positive

Text: "This is terrible."
Sentiment: negative

Text: "It's okay, nothing special."
Sentiment: neutral

Text: "Best purchase ever!"
Sentiment:"""

inputs = tokenizer(prompt, return_tensors="pt")
outputs = model.generate(**inputs, max_new_tokens=5, temperature=0.1)
response = tokenizer.decode(outputs[0], skip_special_tokens=True)
print(response.split("Sentiment:")[-1].strip())  # "positive"
```

## Instruction-Following / Chat

Use instruction-tuned models.

```python
from transformers import AutoModelForCausalLM, AutoTokenizer

tokenizer = AutoTokenizer.from_pretrained("meta-llama/Llama-2-7b-chat-hf")
model = AutoModelForCausalLM.from_pretrained("meta-llama/Llama-2-7b-chat-hf")

messages = [
    {"role": "system", "content": "You are a helpful assistant."},
    {"role": "user", "content": "What is machine learning?"},
]

formatted = tokenizer.apply_chat_template(
    messages,
    tokenize=False,
    add_generation_prompt=True
)

inputs = tokenizer(formatted, return_tensors="pt")
outputs = model.generate(**inputs, max_new_tokens=200, temperature=0.7)
response = tokenizer.decode(outputs[0], skip_special_tokens=True)

# Extract assistant response
assistant_response = response.split("[/INST]")[-1].strip()
print(assistant_response)
```

## Embeddings / Semantic Search

Generate embeddings for semantic similarity.

```python
from transformers import AutoTokenizer, AutoModel
import torch

tokenizer = AutoTokenizer.from_pretrained("sentence-transformers/all-MiniLM-L6-v2")
model = AutoModel.from_pretrained("sentence-transformers/all-MiniLM-L6-v2")

def get_embedding(text):
    inputs = tokenizer(text, return_tensors="pt", padding=True, truncation=True)
    with torch.no_grad():
        outputs = model(**inputs)

    # Mean pooling
    embeddings = outputs.last_hidden_state.mean(dim=1)
    return embeddings

# Get embeddings
text1 = "Machine learning is a subset of AI"
text2 = "AI includes machine learning"

emb1 = get_embedding(text1)
emb2 = get_embedding(text2)

# Compute similarity
similarity = torch.nn.functional.cosine_similarity(emb1, emb2)
print(f"Similarity: {similarity.item():.4f}")  # ~0.85
```

## Multimodal Understanding (CLIP)

Connect vision and language.

```python
from transformers import CLIPProcessor, CLIPModel
from PIL import Image

model = CLIPModel.from_pretrained("openai/clip-vit-base-patch32")
processor = CLIPProcessor.from_pretrained("openai/clip-vit-base-patch32")

image = Image.open("photo.jpg")
texts = ["a dog", "a cat", "a car", "a house"]

inputs = processor(text=texts, images=image, return_tensors="pt", padding=True)
outputs = model(**inputs)

# Get similarity scores
logits_per_image = outputs.logits_per_image
probs = logits_per_image.softmax(dim=1)

for text, prob in zip(texts, probs[0]):
    print(f"{text}: {prob.item():.4f}")
```

## Common Evaluation Metrics

```python
from datasets import load_metric

# Accuracy (classification)
metric = load_metric("accuracy")
predictions = [0, 1, 1, 0]
references = [0, 1, 0, 0]
result = metric.compute(predictions=predictions, references=references)

# F1 Score (classification, NER)
metric = load_metric("f1")
result = metric.compute(predictions=predictions, references=references)

# BLEU (translation)
metric = load_metric("bleu")
predictions = ["hello there general kenobi"]
references = [["hello there general kenobi", "hello there!"]]
result = metric.compute(predictions=predictions, references=references)

# ROUGE (summarization)
metric = load_metric("rouge")
predictions = ["summary text"]
references = ["reference summary"]
result = metric.compute(predictions=predictions, references=references)
```

## Common Data Collators

```python
from transformers import (
    DataCollatorWithPadding,
    DataCollatorForTokenClassification,
    DataCollatorForSeq2Seq,
    DataCollatorForLanguageModeling,
)

# Classification: dynamic padding
DataCollatorWithPadding(tokenizer=tokenizer)

# NER: pad labels too
DataCollatorForTokenClassification(tokenizer=tokenizer)

# Seq2Seq: pad inputs and labels
DataCollatorForSeq2Seq(tokenizer=tokenizer, model=model)

# Language modeling: create MLM masks
DataCollatorForLanguageModeling(tokenizer=tokenizer, mlm=True, mlm_probability=0.15)
```

This covers the most common task patterns. For detailed parameter tuning, see `api_reference.md` and `generation_strategies.md`.
