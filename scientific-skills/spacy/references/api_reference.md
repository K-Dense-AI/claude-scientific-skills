# spaCy API Quick Reference

## Loading Models

| Command | Description |
|---------|-------------|
| `spacy.load("en_core_web_sm")` | Small English model |
| `spacy.load("en_core_web_md")` | Medium English model (w/ vectors) |
| `spacy.load("en_core_web_lg")` | Large English model (w/ vectors) |
| `spacy.load("en_core_web_trf")` | Transformer English model |
| `spacy.load("en_core_sci_md")` | SciSpaCy biomedical model |
| `spacy.blank("en")` | Blank English model |

## Doc Object

| Attribute/Method | Returns | Description |
|------------------|---------|-------------|
| `doc.text` | str | Original text |
| `doc.ents` | tuple[Span] | Named entities |
| `doc.sents` | generator[Span] | Sentences |
| `doc.noun_chunks` | generator[Span] | Base noun phrases |
| `doc.vector` | ndarray | Document vector |
| `doc.similarity(other)` | float | Cosine similarity |
| `doc[i]` | Token | Token by index |
| `doc[i:j]` | Span | Slice of tokens |
| `len(doc)` | int | Number of tokens |

## Token Object

| Attribute | Type | Description |
|-----------|------|-------------|
| `token.text` | str | Original text |
| `token.lemma_` | str | Lemma (base form) |
| `token.pos_` | str | Coarse POS tag |
| `token.tag_` | str | Fine-grained POS tag |
| `token.dep_` | str | Dependency label |
| `token.head` | Token | Syntactic parent |
| `token.children` | generator | Syntactic children |
| `token.ent_type_` | str | Entity type (if entity) |
| `token.ent_iob_` | str | IOB entity tag (B/I/O) |
| `token.is_alpha` | bool | Alphabetic characters only |
| `token.is_stop` | bool | Stop word |
| `token.is_punct` | bool | Punctuation |
| `token.is_digit` | bool | Digit characters only |
| `token.like_num` | bool | Resembles a number |
| `token.shape_` | str | Word shape (e.g., "Xxxxx") |
| `token.vector` | ndarray | Word vector |
| `token.similarity(other)` | float | Cosine similarity |

## Span Object (Entity / Sentence)

| Attribute | Type | Description |
|-----------|------|-------------|
| `span.text` | str | Span text |
| `span.label_` | str | Entity label |
| `span.start` | int | Start token index |
| `span.end` | int | End token index (exclusive) |
| `span.start_char` | int | Start character offset |
| `span.end_char` | int | End character offset |
| `span.root` | Token | Root token of span |
| `span.vector` | ndarray | Span vector |
| `span.similarity(other)` | float | Cosine similarity |

## POS Tags (Universal Dependencies)

| Tag | Description | Example |
|-----|-------------|---------|
| `ADJ` | Adjective | "significant", "novel" |
| `ADP` | Adposition | "in", "of", "with" |
| `ADV` | Adverb | "significantly", "highly" |
| `NOUN` | Noun | "gene", "mutation", "study" |
| `VERB` | Verb | "inhibits", "activates" |
| `PROPN` | Proper noun | "BRCA1", "Harvard" |
| `DET` | Determiner | "the", "a", "this" |
| `PRON` | Pronoun | "it", "we", "they" |
| `NUM` | Numeral | "42", "three" |
| `PUNCT` | Punctuation | ".", ",", "(" |

## Dependency Labels (Common)

| Label | Description | Example |
|-------|-------------|---------|
| `nsubj` | Nominal subject | "Mutations → cause" |
| `dobj` | Direct object | "inhibits → enzyme" |
| `amod` | Adjectival modifier | "significant → reduction" |
| `prep` | Prepositional modifier | "role → in" |
| `pobj` | Object of preposition | "in → cancer" |
| `compound` | Compound modifier | "breast → cancer" |
| `ROOT` | Root of sentence | main verb |

## Matcher Pattern Syntax

### Token Attributes for Matching

| Key | Type | Description |
|-----|------|-------------|
| `TEXT` | str | Exact text |
| `LOWER` | str | Lowercase text |
| `LEMMA` | str | Token lemma |
| `POS` | str | Coarse POS |
| `TAG` | str | Fine-grained POS |
| `DEP` | str | Dependency label |
| `SHAPE` | str | Word shape |
| `IS_ALPHA` | bool | Alphabetic |
| `IS_DIGIT` | bool | Digit |
| `IS_STOP` | bool | Stop word |
| `LIKE_NUM` | bool | Number-like |
| `LIKE_EMAIL` | bool | Email-like |
| `LIKE_URL` | bool | URL-like |

### Pattern Operators

| Operator | Meaning | Example |
|----------|---------|---------|
| `{"OP": "!"}` | Negation (0 times) | Token must not match |
| `{"OP": "?"}` | Optional (0 or 1) | `[{"LOWER": "very", "OP": "?"}]` |
| `{"OP": "+"}` | One or more | `[{"POS": "ADJ", "OP": "+"}]` |
| `{"OP": "*"}` | Zero or more | `[{"IS_PUNCT": True, "OP": "*"}]` |

### Pattern Examples

```python
# Gene mutation pattern: "BRCA1 mutation" or "p53 variant"
[{"POS": "PROPN"}, {"LOWER": {"IN": ["mutation", "variant", "deletion"]}}]

# Dosage: "500 mg" or "1.5 g"
[{"LIKE_NUM": True}, {"LOWER": {"IN": ["mg", "g", "ml", "µg"]}}]

# p-value: "p < 0.05" or "p = 0.001"
[{"LOWER": "p"}, {"LOWER": {"IN": ["<", ">", "=", "≤", "≥"]}}, {"LIKE_NUM": True}]
```

## Pipeline Components

| Component | Description | Creates |
|-----------|-------------|---------|
| `tok2vec` | Token-to-vector encoding | `token.tensor` |
| `tagger` | Part-of-speech tagger | `token.pos_`, `token.tag_` |
| `parser` | Dependency parser | `token.dep_`, `token.head` |
| `ner` | Named entity recognizer | `doc.ents` |
| `textcat` | Text classifier | `doc.cats` |
| `lemmatizer` | Lemmatizer | `token.lemma_` |
| `sentencizer` | Sentence segmenter | `doc.sents` |
| `entity_linker` | Entity linking | `ent.kb_id_` |
| `entity_ruler` | Rule-based NER | `doc.ents` |

## Batch Processing

```python
# Efficient batch processing
docs = nlp.pipe(texts, batch_size=50, n_process=2)

# Disable unused components for speed
docs = nlp.pipe(texts, disable=["parser", "tagger"])

# As tuples (with context)
data = [("text1", {"id": 1}), ("text2", {"id": 2})]
for doc, context in nlp.pipe(data, as_tuples=True):
    print(context["id"], len(doc.ents))
```

## CLI Commands

| Command | Purpose |
|---------|---------|
| `python -m spacy download <model>` | Download model |
| `python -m spacy info` | System info |
| `python -m spacy validate` | Validate models |
| `python -m spacy init config config.cfg` | Generate config |
| `python -m spacy train config.cfg --output ./model` | Train model |
| `python -m spacy evaluate ./model ./test.spacy` | Evaluate model |
| `python -m spacy convert data.json ./` | Convert data |
| `python -m spacy package ./model ./packages` | Package model |

## SciSpaCy Extensions

| Component | Description |
|-----------|-------------|
| `abbreviation_detector` | Detects abbreviations and their definitions |
| `scispacy_linker` | Links entities to UMLS/MeSH/GO/HPO |

### Linker Names
`umls`, `mesh`, `go`, `hpo`, `rxnorm`, `drugbank`

### Abbreviation Detection

```python
from scispacy.abbreviation import AbbreviationDetector
nlp.add_pipe("abbreviation_detector")
doc = nlp("Polymerase Chain Reaction (PCR) is used for amplification.")
for abrv in doc._.abbreviations:
    print(f"{abrv} ← {abrv._.long_form}")
```
