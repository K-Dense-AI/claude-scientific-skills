"""Configuration for AI Code Execution System."""

import os
from dotenv import load_dotenv

load_dotenv()

# AI backend: "gemini" | "claude" | "openai"
AI_BACKEND = os.getenv("AI_BACKEND", "groq")

# Retry on execution error
MAX_RETRIES = 3

# Warn if .md is too long (approximate tokens)
MAX_MD_TOKENS = 100000

# E2B sandbox timeout in seconds
SANDBOX_TIMEOUT = 120

# E2B file paths
DATA_DIR = "/home/user/data/"
OUTPUT_DIR = "/home/user/output/"

# Model defaults per backend
MODEL_DEFAULTS = {
    "gemini": "gemini-2.0-flash",
    "claude": "claude-haiku-4-5",
    "openai": "gpt-4o-mini",
    "groq": "llama-3.3-70b-versatile",
}

# API key env var names
API_KEY_VARS = {
    "gemini": "GEMINI_API_KEY",
    "claude": "ANTHROPIC_API_KEY",
    "openai": "OPENAI_API_KEY",
    "groq": "GROQ_API_KEY",
}
