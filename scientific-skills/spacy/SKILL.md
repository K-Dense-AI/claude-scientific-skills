---
name: spacy
description: Industrial-strength NLP with spaCy. Use for tokenization, part-of-speech tagging, dependency parsing, named entity recognition (NER), entity linking, rule-based matching, text classification, sentence segmentation, and biomedical/scientific text processing with SciSpaCy. Supports custom NER model training for domain-specific entity extraction.
license: MIT license
metadata:
    skill-author: Adem Usta
---

# spaCy: Natural Language Processing

## Overview

spaCy is the leading production-grade NLP library in Python. It provides fast, accurate pipelines for tokenization, POS tagging, dependency parsing, named entity recognition, text classification, and more. Combined with SciSpaCy, it becomes a powerful tool for biomedical and scientific text mining — extracting genes, chemicals, diseases, and other domain entities from research literature.

## Installation

```bash
# Install spaCy using uv
uv pip install spacy

# Download English model (medium, with word vectors)
python -m spacy download en_core_web_md

# For best accuracy (transformer-based)
uv pip install "spacy[transformers]"
python -m spacy download en_core_web_trf

# For biomedical/scientific text
uv pip install scispacy
pip install https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.4/en_core_sci_md-0.5.4.tar.gz

# Commonly used with
uv pip install pandas numpy
```

### Available Models

| Model | Size | Accuracy | Best For |
|-------|------|----------|----------|
| `en_core_web_sm` | 12 MB | Good | Quick prototyping |
| `en_core_web_md` | 40 MB | Better | General NLP with word vectors |
| `en_core_web_lg` | 560 MB | Better | Word vector similarity tasks |
| `en_core_web_trf` | 438 MB | Best | Highest accuracy (transformer) |
| `en_core_sci_sm` | 17 MB | Good | Biomedical text (fast) |
| `en_core_sci_md` | 83 MB | Better | Biomedical text (recommended) |
| `en_core_sci_lg` | 382 MB | Best | Biomedical text (with vectors) |

## When to Use This Skill

Use the spaCy skill when:

- Extracting named entities (persons, organizations, genes, chemicals, diseases)
- Processing scientific or biomedical text for NER
- Tokenizing and parsing text with linguistic annotations
- Building custom NER models for domain-specific entity types
- Performing rule-based text matching (patterns, phrases)
- Classifying text into categories
- Analyzing syntactic structure (dependency trees, noun chunks)
- Preprocessing text for downstream ML models
- Linking entities to knowledge bases (UMLS, MeSH)

Consider alternatives when:
- You need large language model generation → Use **transformers**
- You need topic modeling → Use **scikit-learn** with TF-IDF
- You need sentiment analysis at scale → Use **transformers**
- You need only tokenization → Use simple regex or `nltk`

## Quick Start

### Basic NLP Pipeline

```python
import spacy

# Load model
nlp = spacy.load("en_core_web_md")

# Process text
doc = nlp("BRCA1 mutations increase the risk of breast cancer. "
          "Trastuzumab is used for HER2-positive patients at "
          "Memorial Sloan Kettering Cancer Center in New York.")

# Access tokens
for token in doc:
    print(f"{token.text:15s} {token.pos_:6s} {token.dep_:10s} {token.lemma_}")

# Named entities
for ent in doc.ents:
    print(f"{ent.text:30s} {ent.label_:15s} ({ent.start_char}-{ent.end_char})")

# Sentences
for sent in doc.sents:
    print(f"→ {sent.text}")
```

### Biomedical NER with SciSpaCy

```python
import spacy
import scispacy

# Load biomedical model
nlp = spacy.load("en_core_sci_md")

text = """
Metformin activates AMPK signaling pathway and reduces hepatic glucose
production. In a randomized controlled trial, patients with type 2 diabetes
receiving 500mg metformin twice daily showed significant reduction in HbA1c
levels (p < 0.001) compared to placebo.
"""

doc = nlp(text)

# Extract biomedical entities
for ent in doc.ents:
    print(f"{ent.text:40s} → {ent.label_}")

# Output examples:
# Metformin                                → CHEMICAL
# AMPK                                     → GENE_OR_GENE_PRODUCT
# hepatic glucose production               → BIOLOGICAL_PROCESS
# type 2 diabetes                          → DISEASE
# HbA1c                                    → CHEMICAL
```

## Core Capabilities

### 1. Named Entity Recognition (NER)

Out-of-the-box entity types for English models:

| Label | Description | Example |
|-------|-------------|---------|
| `PERSON` | People, characters | "Marie Curie" |
| `ORG` | Companies, institutions | "MIT", "WHO" |
| `GPE` | Countries, cities, states | "New York", "Germany" |
| `DATE` | Dates and periods | "January 2024" |
| `MONEY` | Monetary values | "$50 million" |
| `CARDINAL` | Numerals | "three", "42" |
| `PRODUCT` | Objects, vehicles | "iPhone" |
| `EVENT` | Named events | "World Cup" |
| `LAW` | Named laws or documents | "GDPR" |
| `WORK_OF_ART` | Titles of works | "Nature" |

```python
# Extract entities with context
doc = nlp("Albert Einstein published his theory of general relativity in 1915.")
for ent in doc.ents:
    print(f"{ent.text} ({ent.label_}): "
          f"'{doc[max(0, ent.start-2):ent.end+2].text}'")
```

### 2. Rule-Based Matching

Pattern matching beyond simple regex — using linguistic features.

```python
from spacy.matcher import Matcher, PhraseMatcher

nlp = spacy.load("en_core_web_md")

# Token-based pattern matching
matcher = Matcher(nlp.vocab)

# Match drug dosage patterns: "500mg", "250 mg", "1.5 g"
pattern = [
    {"LIKE_NUM": True},
    {"LOWER": {"IN": ["mg", "g", "ml", "µg", "mcg", "iu"]}}
]
matcher.add("DOSAGE", [pattern])

doc = nlp("Administer 500 mg of amoxicillin every 8 hours")
matches = matcher(doc)
for match_id, start, end in matches:
    span = doc[start:end]
    print(f"Found dosage: {span.text}")

# Phrase matching (fast exact matching)
phrase_matcher = PhraseMatcher(nlp.vocab, attr="LOWER")
gene_names = ["brca1", "tp53", "egfr", "her2", "kras", "pten", "akt1"]
patterns = [nlp.make_doc(name) for name in gene_names]
phrase_matcher.add("GENE", patterns)

doc = nlp("Mutations in BRCA1 and TP53 are common in breast cancer")
matches = phrase_matcher(doc)
for match_id, start, end in matches:
    print(f"Gene: {doc[start:end].text}")
```

### 3. Custom NER Training

Train a custom NER model for domain-specific entities.

```python
import spacy
from spacy.tokens import DocBin
from spacy.training import Example
import json

# Step 1: Prepare training data
TRAINING_DATA = [
    ("Aspirin 100mg was prescribed for chest pain", {
        "entities": [(0, 7, "DRUG"), (8, 13, "DOSAGE"), (34, 44, "SYMPTOM")]
    }),
    ("Patient received metformin 500mg for diabetes", {
        "entities": [(17, 26, "DRUG"), (27, 32, "DOSAGE"), (37, 45, "DISEASE")]
    }),
]

# Step 2: Create config (recommended: use CLI)
# python -m spacy init config config.cfg --lang en --pipeline ner

# Step 3: Convert to DocBin format
nlp = spacy.blank("en")
db = DocBin()
for text, annotations in TRAINING_DATA:
    doc = nlp.make_doc(text)
    example = Example.from_dict(doc, annotations)
    db.add(example.reference)
db.to_disk("./train.spacy")

# Step 4: Train via CLI
# python -m spacy train config.cfg --output ./model --paths.train ./train.spacy
```

### 4. Dependency Parsing and Syntax

```python
doc = nlp("The KRAS mutation drives tumor growth through the MAPK pathway.")

# Noun chunks (base noun phrases)
for chunk in doc.noun_chunks:
    print(f"{chunk.text:30s} root={chunk.root.text}, dep={chunk.root.dep_}")

# Dependency tree
for token in doc:
    print(f"{token.text:15s} ──{token.dep_:10s}──▶ {token.head.text}")

# Find subject-verb-object triples
def extract_svo(doc):
    triples = []
    for token in doc:
        if token.dep_ == "nsubj":
            subject = token.text
            verb = token.head.text
            objects = [child.text for child in token.head.children
                       if child.dep_ in ("dobj", "attr", "prep")]
            triples.append((subject, verb, " ".join(objects)))
    return triples

triples = extract_svo(doc)
for s, v, o in triples:
    print(f"{s} → {v} → {o}")
```

### 5. Text Classification

```python
import spacy
from spacy.training import Example

# Add text categorizer to pipeline
nlp = spacy.blank("en")
textcat = nlp.add_pipe("textcat_multilabel")
textcat.add_label("METHODS")
textcat.add_label("RESULTS")
textcat.add_label("DISCUSSION")

# Train with examples
train_data = [
    ("We used PCR to amplify the target sequence", {"cats": {"METHODS": 1.0, "RESULTS": 0.0, "DISCUSSION": 0.0}}),
    ("Expression levels were significantly higher", {"cats": {"METHODS": 0.0, "RESULTS": 1.0, "DISCUSSION": 0.0}}),
]

# Training loop
optimizer = nlp.begin_training()
for epoch in range(20):
    losses = {}
    for text, annotations in train_data:
        doc = nlp.make_doc(text)
        example = Example.from_dict(doc, annotations)
        nlp.update([example], sgd=optimizer, losses=losses)
    print(f"Epoch {epoch}: {losses}")
```

### 6. Word Vectors and Similarity

```python
# Requires a model with vectors (md or lg)
nlp = spacy.load("en_core_web_md")

# Document similarity
doc1 = nlp("CRISPR gene editing technology for cancer therapy")
doc2 = nlp("Genome modification tools for treating tumors")
doc3 = nlp("Stock market analysis for portfolio optimization")

print(f"doc1 ↔ doc2: {doc1.similarity(doc2):.3f}")  # High similarity
print(f"doc1 ↔ doc3: {doc1.similarity(doc3):.3f}")  # Low similarity

# Word similarity
token1 = nlp("gene")[0]
token2 = nlp("genome")[0]
print(f"gene ↔ genome: {token1.similarity(token2):.3f}")

# Find most similar words
import numpy as np
target = nlp("mutation").vector
words = ["variation", "deletion", "expression", "weather", "finance"]
for word in words:
    sim = np.dot(target, nlp(word).vector) / (
        np.linalg.norm(target) * np.linalg.norm(nlp(word).vector))
    print(f"mutation ↔ {word}: {sim:.3f}")
```

## Common Scientific Workflows

### Extract Gene-Disease Relations from Abstracts

```python
import spacy
from spacy.matcher import Matcher

nlp = spacy.load("en_core_sci_md")

abstracts = [
    "Mutations in BRCA1 are associated with hereditary breast and ovarian cancer.",
    "TP53 loss-of-function variants drive tumorigenesis in Li-Fraumeni syndrome.",
    "EGFR amplification is a hallmark of glioblastoma multiforme.",
]

for text in abstracts:
    doc = nlp(text)
    entities = [(ent.text, ent.label_) for ent in doc.ents]
    print(f"\nText: {text[:60]}...")
    for name, label in entities:
        print(f"  {name:30s} → {label}")
```

### Scientific Literature Processing Pipeline

```python
import spacy
import pandas as pd

nlp = spacy.load("en_core_sci_md")

# Process abstracts in batch (more efficient)
abstracts = pd.read_csv("pubmed_abstracts.csv")["abstract"].tolist()

# Batch processing
results = []
for doc in nlp.pipe(abstracts, batch_size=50, n_process=2):
    entities = [(ent.text, ent.label_, ent.start_char, ent.end_char)
                for ent in doc.ents]
    results.append({
        "n_entities": len(entities),
        "entities": entities,
        "n_sentences": len(list(doc.sents)),
    })

df_results = pd.DataFrame(results)
print(f"Processed {len(df_results)} abstracts")
print(f"Total entities: {df_results['n_entities'].sum()}")
```

### UMLS Entity Linking (SciSpaCy)

```python
import spacy
import scispacy
from scispacy.linking import EntityLinker

nlp = spacy.load("en_core_sci_md")

# Add UMLS entity linker
nlp.add_pipe("scispacy_linker", config={
    "resolve_abbreviations": True,
    "linker_name": "umls"
})

doc = nlp("Aspirin inhibits cyclooxygenase-2 enzyme activity, "
          "reducing prostaglandin synthesis.")

# Get linked entities
linker = nlp.get_pipe("scispacy_linker")
for ent in doc.ents:
    print(f"\n{ent.text} ({ent.label_})")
    for umls_ent in ent._.kb_ents:
        cui, score = umls_ent
        concept = linker.kb.cui_to_entity[cui]
        print(f"  CUI: {cui}, Score: {score:.3f}")
        print(f"  Name: {concept.canonical_name}")
        print(f"  Types: {concept.types}")
```

## Key Parameters to Adjust

### Pipeline Configuration
- `nlp.max_length`: Maximum doc length (default 1M chars)
- `nlp.pipe(texts, batch_size=X)`: Batch size for processing
- `nlp.pipe(texts, n_process=N)`: Parallel processing workers

### Matcher Configuration
- `attr`: Token attribute to match on ('TEXT', 'LOWER', 'LEMMA', 'POS')
- Pattern operators: `{"OP": "?"}` (optional), `{"OP": "+"}` (one or more), `{"OP": "*"}` (zero or more)

### Training
- `max_epochs`: Number of training passes
- `learn_rate`: Learning rate (default 0.001)
- `batch_size`: Training batch size
- `dropout`: Dropout rate (default 0.1)

## Common Pitfalls and Best Practices

1. **Use `nlp.pipe()` for batch processing**: Much faster than calling `nlp()` in a loop
2. **Choose the right model**: `sm` for speed, `trf` for accuracy, `sci_*` for biomedical
3. **Disable unused components**: `nlp.pipe(texts, disable=["parser", "tagger"])` for speed
4. **Handle long text**: Increase `nlp.max_length` or split documents into chunks
5. **SciSpaCy for science**: Always prefer `en_core_sci_*` models for scientific text
6. **Custom NER via CLI**: Use `spacy train` instead of manual training loops for production
7. **Validate entity spans**: Check for overlapping or misaligned entity annotations
8. **Save/load custom models**: Use `nlp.to_disk("model_dir")` and `spacy.load("model_dir")`

## Integration with Other Tools

- **SciSpaCy**: Biomedical NER, abbreviation detection, UMLS/MeSH linking
- **pandas**: Process DataFrames with `nlp.pipe(df['text'])`
- **transformers (Hugging Face)**: Use `spacy-transformers` for transformer-based pipelines
- **prodigy**: Active learning annotation tool (by spaCy team)
- **pubmed-database**: Extract entities from PubMed abstracts
- **biorxiv-database**: Process bioRxiv preprints

## Additional Resources

- **Official Documentation**: https://spacy.io/
- **spaCy 101**: https://spacy.io/usage/spacy-101
- **Training Guide**: https://spacy.io/usage/training
- **SciSpaCy**: https://allenai.github.io/scispacy/
- **Models**: https://spacy.io/models
- **Universe (plugins)**: https://spacy.io/universe

## Reference Documentation

- **API Reference**: See [references/api_reference.md](references/api_reference.md) for pipeline component and API listing
