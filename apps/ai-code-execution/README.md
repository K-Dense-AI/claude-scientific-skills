# AI-Powered Code Execution

Upload a `.md` instruction file (and optionally data files), and an AI generates Python code and executes it in a cloud sandbox (E2B). No coding or local setup required.

## Setup

1. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

2. Copy `.env.example` to `.env` and fill in your API keys:
   ```bash
   cp .env.example .env
   ```

3. Run:
   ```bash
   streamlit run app.py
   ```

## API Keys

| Service | Where to get |
|---------|-------------|
| Gemini | https://aistudio.google.com/app/apikey |
| Claude | https://console.anthropic.com/ |
| OpenAI | https://platform.openai.com/api-keys |
| E2B | https://e2b.dev/dashboard |

## Instruction Modes

- **NO CODE**: Describe what you want in plain text. The AI writes all the code.
- **WITH CODE**: Include ` ```python ` blocks in your `.md`. The AI completes/fixes/runs them.

The system auto-detects the mode based on whether Python code blocks exist in your `.md` file.

## Configuration

Edit `config.py` or set environment variables:
- `AI_BACKEND`: `"gemini"` (default), `"claude"`, or `"openai"`
- `MAX_RETRIES`: Number of retry attempts on execution failure (default: 3)
- `SANDBOX_TIMEOUT`: E2B sandbox timeout in seconds (default: 120)
