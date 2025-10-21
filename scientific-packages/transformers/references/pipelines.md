# Transformers Pipelines

Pipelines provide a simple and optimized interface for inference across many machine learning tasks. They abstract away the complexity of tokenization, model invocation, and post-processing.

## Usage Pattern

```python
from transformers import pipeline

# Basic usage
classifier = pipeline("text-classification")
result = classifier("This movie was amazing!")

# With specific model
classifier = pipeline("text-classification", model="distilbert-base-uncased-finetuned-sst-2-english")
result = classifier("This movie was amazing!")
```

## Natural Language Processing Pipelines

### Text Classification
```python
classifier = pipeline("text-classification")
classifier("I love this product!")
# [{'label': 'POSITIVE', 'score': 0.9998}]
```

### Zero-Shot Classification
```python
classifier = pipeline("zero-shot-classification")
classifier("This is about climate change", candidate_labels=["politics", "science", "sports"])
```

### Token Classification (NER)
```python
ner = pipeline("token-classification")
ner("My name is Sarah and I work at Microsoft in Seattle")
```

### Question Answering
```python
qa = pipeline("question-answering")
qa(question="What is the capital?", context="The capital of France is Paris.")
```

### Text Generation
```python
generator = pipeline("text-generation")
generator("Once upon a time", max_length=50)
```

### Text2Text Generation
```python
generator = pipeline("text2text-generation", model="t5-base")
generator("translate English to French: Hello")
```

### Summarization
```python
summarizer = pipeline("summarization")
summarizer("Long article text here...", max_length=130, min_length=30)
```

### Translation
```python
translator = pipeline("translation_en_to_fr")
translator("Hello, how are you?")
```

### Fill Mask
```python
unmasker = pipeline("fill-mask")
unmasker("Paris is the [MASK] of France.")
```

### Feature Extraction
```python
extractor = pipeline("feature-extraction")
embeddings = extractor("This is a sentence")
```

### Document Question Answering
```python
doc_qa = pipeline("document-question-answering")
doc_qa(image="document.png", question="What is the invoice number?")
```

### Table Question Answering
```python
table_qa = pipeline("table-question-answering")
table_qa(table=data, query="How many employees?")
```

## Computer Vision Pipelines

### Image Classification
```python
classifier = pipeline("image-classification")
classifier("cat.jpg")
```

### Zero-Shot Image Classification
```python
classifier = pipeline("zero-shot-image-classification")
classifier("cat.jpg", candidate_labels=["cat", "dog", "bird"])
```

### Object Detection
```python
detector = pipeline("object-detection")
detector("street.jpg")
```

### Image Segmentation
```python
segmenter = pipeline("image-segmentation")
segmenter("image.jpg")
```

### Image-to-Image
```python
img2img = pipeline("image-to-image", model="lllyasviel/sd-controlnet-canny")
img2img("input.jpg")
```

### Depth Estimation
```python
depth = pipeline("depth-estimation")
depth("image.jpg")
```

### Video Classification
```python
classifier = pipeline("video-classification")
classifier("video.mp4")
```

### Keypoint Matching
```python
matcher = pipeline("keypoint-matching")
matcher(image1="img1.jpg", image2="img2.jpg")
```

## Audio Pipelines

### Automatic Speech Recognition
```python
asr = pipeline("automatic-speech-recognition")
asr("audio.wav")
```

### Audio Classification
```python
classifier = pipeline("audio-classification")
classifier("audio.wav")
```

### Zero-Shot Audio Classification
```python
classifier = pipeline("zero-shot-audio-classification")
classifier("audio.wav", candidate_labels=["speech", "music", "noise"])
```

### Text-to-Audio/Text-to-Speech
```python
synthesizer = pipeline("text-to-audio")
audio = synthesizer("Hello, how are you today?")
```

## Multimodal Pipelines

### Image-to-Text (Image Captioning)
```python
captioner = pipeline("image-to-text")
captioner("image.jpg")
```

### Visual Question Answering
```python
vqa = pipeline("visual-question-answering")
vqa(image="image.jpg", question="What color is the car?")
```

### Image-Text-to-Text (VLMs)
```python
vlm = pipeline("image-text-to-text")
vlm(images="image.jpg", text="Describe this image in detail")
```

### Zero-Shot Object Detection
```python
detector = pipeline("zero-shot-object-detection")
detector("image.jpg", candidate_labels=["car", "person", "tree"])
```

## Pipeline Configuration

### Common Parameters

- `model`: Specify model identifier or path
- `device`: Set device (0 for GPU, -1 for CPU, or "cuda:0")
- `batch_size`: Process multiple inputs at once
- `torch_dtype`: Set precision (torch.float16, torch.bfloat16)

```python
# GPU with half precision
pipe = pipeline("text-generation", model="gpt2", device=0, torch_dtype=torch.float16)

# Batch processing
pipe(["text 1", "text 2", "text 3"], batch_size=8)
```

### Task-Specific Parameters

Each pipeline accepts task-specific parameters in the call:

```python
# Text generation
generator("prompt", max_length=100, temperature=0.7, top_p=0.9, num_return_sequences=3)

# Summarization
summarizer("text", max_length=130, min_length=30, do_sample=False)

# Translation
translator("text", max_length=512, num_beams=4)
```

## Best Practices

1. **Reuse pipelines**: Create once, use multiple times for efficiency
2. **Batch processing**: Use batches for multiple inputs to maximize throughput
3. **GPU acceleration**: Set `device=0` for GPU when available
4. **Model selection**: Choose task-specific models for best results
5. **Memory management**: Use `torch_dtype=torch.float16` for large models
