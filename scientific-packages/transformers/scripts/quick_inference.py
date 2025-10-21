#!/usr/bin/env python3
"""
Quick inference using Transformers pipelines.

This script demonstrates how to quickly use pre-trained models for inference
across various tasks using the pipeline API.
"""

from transformers import pipeline


def text_classification_example():
    """Sentiment analysis example."""
    print("=== Text Classification ===")
    classifier = pipeline("text-classification")
    result = classifier("I love using Transformers! It makes NLP so easy.")
    print(f"Result: {result}\n")


def named_entity_recognition_example():
    """Named Entity Recognition example."""
    print("=== Named Entity Recognition ===")
    ner = pipeline("token-classification", aggregation_strategy="simple")
    text = "My name is Sarah and I work at Microsoft in Seattle"
    entities = ner(text)
    for entity in entities:
        print(f"{entity['word']}: {entity['entity_group']} (score: {entity['score']:.3f})")
    print()


def question_answering_example():
    """Question Answering example."""
    print("=== Question Answering ===")
    qa = pipeline("question-answering")
    context = "Paris is the capital and most populous city of France. It is located in northern France."
    question = "What is the capital of France?"
    answer = qa(question=question, context=context)
    print(f"Question: {question}")
    print(f"Answer: {answer['answer']} (score: {answer['score']:.3f})\n")


def text_generation_example():
    """Text generation example."""
    print("=== Text Generation ===")
    generator = pipeline("text-generation", model="gpt2")
    prompt = "Once upon a time in a land far away"
    generated = generator(prompt, max_length=50, num_return_sequences=1)
    print(f"Prompt: {prompt}")
    print(f"Generated: {generated[0]['generated_text']}\n")


def summarization_example():
    """Text summarization example."""
    print("=== Summarization ===")
    summarizer = pipeline("summarization")
    article = """
    The Transformers library provides thousands of pretrained models to perform tasks
    on texts such as classification, information extraction, question answering,
    summarization, translation, text generation, etc in over 100 languages. Its aim
    is to make cutting-edge NLP easier to use for everyone. The library provides APIs
    to quickly download and use pretrained models on a given text, fine-tune them on
    your own datasets then share them with the community on the model hub.
    """
    summary = summarizer(article, max_length=50, min_length=25, do_sample=False)
    print(f"Summary: {summary[0]['summary_text']}\n")


def translation_example():
    """Translation example."""
    print("=== Translation ===")
    translator = pipeline("translation_en_to_fr")
    text = "Hello, how are you today?"
    translation = translator(text)
    print(f"English: {text}")
    print(f"French: {translation[0]['translation_text']}\n")


def zero_shot_classification_example():
    """Zero-shot classification example."""
    print("=== Zero-Shot Classification ===")
    classifier = pipeline("zero-shot-classification")
    text = "This is a breaking news story about a major earthquake."
    candidate_labels = ["politics", "sports", "science", "breaking news"]
    result = classifier(text, candidate_labels)
    print(f"Text: {text}")
    print("Predictions:")
    for label, score in zip(result['labels'], result['scores']):
        print(f"  {label}: {score:.3f}")
    print()


def image_classification_example():
    """Image classification example (requires PIL)."""
    print("=== Image Classification ===")
    try:
        from PIL import Image
        import requests

        classifier = pipeline("image-classification")
        url = "http://images.cocodataset.org/val2017/000000039769.jpg"
        image = Image.open(requests.get(url, stream=True).raw)
        predictions = classifier(image)
        print("Top predictions:")
        for pred in predictions[:3]:
            print(f"  {pred['label']}: {pred['score']:.3f}")
        print()
    except ImportError:
        print("PIL not installed. Skipping image classification example.\n")


def main():
    """Run all examples."""
    print("Transformers Quick Inference Examples")
    print("=" * 50 + "\n")

    # Text tasks
    text_classification_example()
    named_entity_recognition_example()
    question_answering_example()
    text_generation_example()
    summarization_example()
    translation_example()
    zero_shot_classification_example()

    # Vision task (optional)
    image_classification_example()

    print("=" * 50)
    print("All examples completed!")


if __name__ == "__main__":
    main()
