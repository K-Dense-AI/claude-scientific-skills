.PHONY: help install retrieve import dedup screen prisma extract meta clean test

help:
	@echo "POD Systematic Review Pipeline"
	@echo ""
	@echo "Available commands:"
	@echo "  make install      - Install dependencies"
	@echo "  make retrieve     - Retrieve citations from PubMed"
	@echo "  make import       - Import citations from file"
	@echo "  make dedup        - Deduplicate citations"
	@echo "  make screen       - Launch screening application"
	@echo "  make prisma       - Generate PRISMA flow diagram"
	@echo "  make extract      - Create extraction templates"
	@echo "  make meta         - Run meta-analysis"
	@echo "  make test         - Run tests"
	@echo "  make clean        - Clean generated files"

install:
	pip install -e ".[dev]"

retrieve:
	pod-review retrieve --max-results 10000

import:
	@echo "Usage: pod-review import-records <file_path> [--format ris|bibtex|csv]"

dedup:
	pod-review dedup --input-dir data --output-dir data/processed

load-screening:
	pod-review load-screening --input-file data/processed/deduplicated.jsonl

screen:
	streamlit run src/pod_review/screening/app.py

prisma:
	pod-review prisma --output-dir outputs/prisma

extract:
	pod-review create-templates --output-dir outputs/extraction

meta:
	@echo "Usage: pod-review meta-analyze <extraction_file> <risk_factor> [--effect-type OR|RR|HR]"
	@echo "Example: pod-review meta-analyze outputs/extraction/extraction_template.csv 'Age' --effect-type OR"

test:
	pytest tests/ -v --cov=pod_review --cov-report=term-missing

clean:
	rm -rf data/raw/* data/processed/* data/imported/* outputs/*
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
	find . -type d -name "*.egg-info" -exec rm -rf {} +
	rm -f .coverage

.env:
	@if [ ! -f .env ]; then \
		echo "Creating .env from .env.example..."; \
		cp .env.example .env; \
		echo "Please edit .env and add your API credentials"; \
	fi
