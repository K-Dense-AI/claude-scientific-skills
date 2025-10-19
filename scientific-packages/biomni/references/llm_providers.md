# LLM Provider Configuration Guide

This document provides comprehensive configuration instructions for all LLM providers supported by Biomni.

## Overview

Biomni supports multiple LLM providers through a unified interface. Configure providers using:
- Environment variables
- `.env` files
- Runtime configuration via `default_config`

## Quick Reference Table

| Provider | Recommended For | API Key Required | Cost | Setup Complexity |
|----------|----------------|------------------|------|------------------|
| Anthropic Claude | Most biomedical tasks | Yes | Medium | Easy |
| OpenAI | General tasks | Yes | Medium-High | Easy |
| Azure OpenAI | Enterprise deployment | Yes | Varies | Medium |
| Google Gemini | Multimodal tasks | Yes | Medium | Easy |
| Groq | Fast inference | Yes | Low | Easy |
| Ollama | Local/offline use | No | Free | Medium |
| AWS Bedrock | AWS ecosystem | Yes | Varies | Hard |
| Biomni-R0 | Complex biological reasoning | No | Free | Hard |

## Anthropic Claude (Recommended)

### Overview

Claude models from Anthropic provide excellent biological reasoning capabilities and are the recommended choice for most Biomni tasks.

### Setup

1. **Obtain API Key:**
   - Sign up at https://console.anthropic.com/
   - Navigate to API Keys section
   - Generate a new key

2. **Configure Environment:**

   **Option A: Environment Variable**
   ```bash
   export ANTHROPIC_API_KEY="sk-ant-api03-..."
   ```

   **Option B: .env File**
   ```bash
   # .env file in project root
   ANTHROPIC_API_KEY=sk-ant-api03-...
   ```

3. **Set Model in Code:**
   ```python
   from biomni.config import default_config

   # Claude Sonnet 4 (Recommended)
   default_config.llm = "claude-sonnet-4-20250514"

   # Claude Opus 4 (Most capable)
   default_config.llm = "claude-opus-4-20250514"

   # Claude 3.5 Sonnet (Previous version)
   default_config.llm = "claude-3-5-sonnet-20241022"
   ```

### Available Models

| Model | Context Window | Strengths | Best For |
|-------|---------------|-----------|----------|
| `claude-sonnet-4-20250514` | 200K tokens | Balanced performance, cost-effective | Most biomedical tasks |
| `claude-opus-4-20250514` | 200K tokens | Highest capability, complex reasoning | Difficult multi-step analyses |
| `claude-3-5-sonnet-20241022` | 200K tokens | Fast, reliable | Standard workflows |
| `claude-3-opus-20240229` | 200K tokens | Strong reasoning | Legacy support |

### Advanced Configuration

```python
from biomni.config import default_config

# Use Claude with custom parameters
default_config.llm = "claude-sonnet-4-20250514"
default_config.timeout_seconds = 1800

# Optional: Custom API endpoint (for proxy/enterprise)
default_config.api_base = "https://your-proxy.com/v1"
```

### Cost Estimation

Approximate costs per 1M tokens (as of January 2025):
- Input: $3-15 depending on model
- Output: $15-75 depending on model

For a typical biomedical analysis (~50K tokens total): $0.50-$2.00

## OpenAI

### Overview

OpenAI's GPT models provide strong general capabilities suitable for diverse biomedical tasks.

### Setup

1. **Obtain API Key:**
   - Sign up at https://platform.openai.com/
   - Navigate to API Keys
   - Create new secret key

2. **Configure Environment:**

   ```bash
   export OPENAI_API_KEY="sk-proj-..."
   ```

   Or in `.env`:
   ```
   OPENAI_API_KEY=sk-proj-...
   ```

3. **Set Model:**
   ```python
   from biomni.config import default_config

   default_config.llm = "gpt-4o"          # Recommended
   # default_config.llm = "gpt-4"         # Previous flagship
   # default_config.llm = "gpt-4-turbo"   # Fast variant
   # default_config.llm = "gpt-3.5-turbo" # Budget option
   ```

### Available Models

| Model | Context Window | Strengths | Cost |
|-------|---------------|-----------|------|
| `gpt-4o` | 128K tokens | Fast, multimodal | Medium |
| `gpt-4-turbo` | 128K tokens | Fast inference | Medium |
| `gpt-4` | 8K tokens | Reliable | High |
| `gpt-3.5-turbo` | 16K tokens | Fast, cheap | Low |

### Cost Optimization

```python
# For exploratory analysis (budget-conscious)
default_config.llm = "gpt-3.5-turbo"

# For production analysis (quality-focused)
default_config.llm = "gpt-4o"
```

## Azure OpenAI

### Overview

Azure-hosted OpenAI models for enterprise users requiring data residency and compliance.

### Setup

1. **Azure Prerequisites:**
   - Active Azure subscription
   - Azure OpenAI resource created
   - Model deployment configured

2. **Environment Variables:**
   ```bash
   export AZURE_OPENAI_API_KEY="your-key"
   export AZURE_OPENAI_ENDPOINT="https://your-resource.openai.azure.com/"
   export AZURE_OPENAI_API_VERSION="2024-02-15-preview"
   ```

3. **Configuration:**
   ```python
   from biomni.config import default_config

   # Option 1: Use deployment name
   default_config.llm = "azure/your-deployment-name"

   # Option 2: Specify endpoint explicitly
   default_config.llm = "azure/gpt-4"
   default_config.api_base = "https://your-resource.openai.azure.com/"
   ```

### Deployment Setup

Azure OpenAI requires explicit model deployments:

1. Navigate to Azure OpenAI Studio
2. Create deployment for desired model (e.g., GPT-4)
3. Note the deployment name
4. Use deployment name in Biomni configuration

### Example Configuration

```python
from biomni.config import default_config
import os

# Set Azure credentials
os.environ['AZURE_OPENAI_API_KEY'] = 'your-key'
os.environ['AZURE_OPENAI_ENDPOINT'] = 'https://your-resource.openai.azure.com/'

# Configure Biomni to use Azure deployment
default_config.llm = "azure/gpt-4-biomni"  # Your deployment name
default_config.api_base = os.environ['AZURE_OPENAI_ENDPOINT']
```

## Google Gemini

### Overview

Google's Gemini models offer multimodal capabilities and competitive performance.

### Setup

1. **Obtain API Key:**
   - Visit https://makersuite.google.com/app/apikey
   - Create new API key

2. **Environment Configuration:**
   ```bash
   export GEMINI_API_KEY="your-key"
   ```

3. **Set Model:**
   ```python
   from biomni.config import default_config

   default_config.llm = "gemini/gemini-1.5-pro"
   # Or: default_config.llm = "gemini/gemini-pro"
   ```

### Available Models

| Model | Context Window | Strengths |
|-------|---------------|-----------|
| `gemini/gemini-1.5-pro` | 1M tokens | Very large context, multimodal |
| `gemini/gemini-pro` | 32K tokens | Balanced performance |

### Use Cases

Gemini excels at:
- Tasks requiring very large context windows
- Multimodal analysis (when incorporating images)
- Cost-effective alternative to GPT-4

```python
# For tasks with large context requirements
default_config.llm = "gemini/gemini-1.5-pro"
default_config.timeout_seconds = 2400  # May need longer timeout
```

## Groq

### Overview

Groq provides ultra-fast inference with open-source models, ideal for rapid iteration.

### Setup

1. **Get API Key:**
   - Sign up at https://console.groq.com/
   - Generate API key

2. **Configure:**
   ```bash
   export GROQ_API_KEY="gsk_..."
   ```

3. **Set Model:**
   ```python
   from biomni.config import default_config

   default_config.llm = "groq/llama-3.1-70b-versatile"
   # Or: default_config.llm = "groq/mixtral-8x7b-32768"
   ```

### Available Models

| Model | Context Window | Speed | Quality |
|-------|---------------|-------|---------|
| `groq/llama-3.1-70b-versatile` | 32K tokens | Very Fast | Good |
| `groq/mixtral-8x7b-32768` | 32K tokens | Very Fast | Good |
| `groq/llama-3-70b-8192` | 8K tokens | Ultra Fast | Moderate |

### Best Practices

```python
# For rapid prototyping and testing
default_config.llm = "groq/llama-3.1-70b-versatile"
default_config.timeout_seconds = 600  # Groq is fast

# Note: Quality may be lower than GPT-4/Claude for complex tasks
# Recommended for: QC, simple analyses, testing workflows
```

## Ollama (Local Deployment)

### Overview

Run LLMs entirely locally for offline use, data privacy, or cost savings.

### Setup

1. **Install Ollama:**
   ```bash
   # macOS/Linux
   curl -fsSL https://ollama.com/install.sh | sh

   # Or download from https://ollama.com/download
   ```

2. **Pull Models:**
   ```bash
   ollama pull llama3       # Meta Llama 3 (8B)
   ollama pull mixtral      # Mixtral (47B)
   ollama pull codellama    # Code-specialized
   ollama pull medllama     # Medical domain (if available)
   ```

3. **Start Ollama Server:**
   ```bash
   ollama serve  # Runs on http://localhost:11434
   ```

4. **Configure Biomni:**
   ```python
   from biomni.config import default_config

   default_config.llm = "ollama/llama3"
   default_config.api_base = "http://localhost:11434"
   ```

### Hardware Requirements

Minimum recommendations:
- **8B models:** 16GB RAM, CPU inference acceptable
- **70B models:** 64GB RAM, GPU highly recommended
- **Storage:** 5-50GB per model

### Model Selection

```python
# Fast, local, good for testing
default_config.llm = "ollama/llama3"

# Better quality (requires more resources)
default_config.llm = "ollama/mixtral"

# Code generation tasks
default_config.llm = "ollama/codellama"
```

### Advantages & Limitations

**Advantages:**
- Complete data privacy
- No API costs
- Offline operation
- Unlimited usage

**Limitations:**
- Lower quality than GPT-4/Claude for complex tasks
- Requires significant hardware
- Slower inference (especially on CPU)
- May struggle with specialized biomedical knowledge

## AWS Bedrock

### Overview

AWS-managed LLM service offering multiple model providers.

### Setup

1. **AWS Prerequisites:**
   - AWS account with Bedrock access
   - Model access enabled in Bedrock console
   - AWS credentials configured

2. **Configure AWS Credentials:**
   ```bash
   # Option 1: AWS CLI
   aws configure

   # Option 2: Environment variables
   export AWS_ACCESS_KEY_ID="your-key"
   export AWS_SECRET_ACCESS_KEY="your-secret"
   export AWS_REGION="us-east-1"
   ```

3. **Enable Model Access:**
   - Navigate to AWS Bedrock console
   - Request access to desired models
   - Wait for approval (may take hours/days)

4. **Configure Biomni:**
   ```python
   from biomni.config import default_config

   default_config.llm = "bedrock/anthropic.claude-3-sonnet"
   # Or: default_config.llm = "bedrock/anthropic.claude-v2"
   ```

### Available Models

Bedrock provides access to:
- Anthropic Claude models
- Amazon Titan models
- AI21 Jurassic models
- Cohere Command models
- Meta Llama models

### IAM Permissions

Required IAM policy:
```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "bedrock:InvokeModel",
        "bedrock:InvokeModelWithResponseStream"
      ],
      "Resource": "arn:aws:bedrock:*::foundation-model/*"
    }
  ]
}
```

### Example Configuration

```python
from biomni.config import default_config
import boto3

# Verify AWS credentials
session = boto3.Session()
credentials = session.get_credentials()
print(f"AWS Access Key: {credentials.access_key[:8]}...")

# Configure Biomni
default_config.llm = "bedrock/anthropic.claude-3-sonnet"
default_config.timeout_seconds = 1800
```

## Biomni-R0 (Local Specialized Model)

### Overview

Biomni-R0 is a 32B parameter reasoning model specifically trained for biological problem-solving. Provides the highest quality for complex biomedical reasoning but requires local deployment.

### Setup

1. **Hardware Requirements:**
   - GPU with 48GB+ VRAM (e.g., A100, H100)
   - Or multi-GPU setup (2x 24GB)
   - 100GB+ storage for model weights

2. **Install Dependencies:**
   ```bash
   pip install "sglang[all]"
   pip install flashinfer  # Optional but recommended
   ```

3. **Deploy Model:**
   ```bash
   python -m sglang.launch_server \
       --model-path snap-stanford/biomni-r0 \
       --host 0.0.0.0 \
       --port 30000 \
       --trust-remote-code \
       --mem-fraction-static 0.8
   ```

   For multi-GPU:
   ```bash
   python -m sglang.launch_server \
       --model-path snap-stanford/biomni-r0 \
       --host 0.0.0.0 \
       --port 30000 \
       --trust-remote-code \
       --tp 2  # Tensor parallelism across 2 GPUs
   ```

4. **Configure Biomni:**
   ```python
   from biomni.config import default_config

   default_config.llm = "openai/biomni-r0"
   default_config.api_base = "http://localhost:30000/v1"
   default_config.timeout_seconds = 2400  # Longer for complex reasoning
   ```

### When to Use Biomni-R0

Biomni-R0 excels at:
- Multi-step biological reasoning
- Complex experimental design
- Hypothesis generation and evaluation
- Literature-informed analysis
- Tasks requiring deep biological knowledge

```python
# For complex biological reasoning tasks
default_config.llm = "openai/biomni-r0"

agent.go("""
Design a comprehensive CRISPR screening experiment to identify synthetic
lethal interactions with TP53 mutations in cancer cells, including:
1. Rationale and hypothesis
2. Guide RNA library design strategy
3. Experimental controls
4. Statistical analysis plan
5. Expected outcomes and validation approach
""")
```

### Performance Comparison

| Model | Speed | Biological Reasoning | Code Quality | Cost |
|-------|-------|---------------------|--------------|------|
| GPT-4 | Fast | Good | Excellent | Medium |
| Claude Sonnet 4 | Fast | Excellent | Excellent | Medium |
| Biomni-R0 | Moderate | Outstanding | Good | Free (local) |

## Multi-Provider Strategy

### Intelligent Model Selection

Use different models for different task types:

```python
from biomni.agent import A1
from biomni.config import default_config

# Strategy 1: Task-based selection
def get_agent_for_task(task_complexity):
    if task_complexity == "simple":
        default_config.llm = "gpt-3.5-turbo"
        default_config.timeout_seconds = 300
    elif task_complexity == "medium":
        default_config.llm = "claude-sonnet-4-20250514"
        default_config.timeout_seconds = 1200
    else:  # complex
        default_config.llm = "openai/biomni-r0"
        default_config.timeout_seconds = 2400

    return A1(path='./data')

# Strategy 2: Fallback on failure
def execute_with_fallback(task):
    models = [
        "claude-sonnet-4-20250514",
        "gpt-4o",
        "claude-opus-4-20250514"
    ]

    for model in models:
        try:
            default_config.llm = model
            agent = A1(path='./data')
            agent.go(task)
            return
        except Exception as e:
            print(f"Failed with {model}: {e}, trying next...")

    raise Exception("All models failed")
```

### Cost Optimization Strategy

```python
# Phase 1: Rapid prototyping with cheap models
default_config.llm = "gpt-3.5-turbo"
agent.go("Quick exploratory analysis of dataset structure")

# Phase 2: Detailed analysis with high-quality models
default_config.llm = "claude-sonnet-4-20250514"
agent.go("Comprehensive differential expression analysis with pathway enrichment")

# Phase 3: Complex reasoning with specialized models
default_config.llm = "openai/biomni-r0"
agent.go("Generate biological hypotheses based on multi-omics integration")
```

## Troubleshooting

### Common Issues

**Issue: "API key not found"**
- Verify environment variable is set: `echo $ANTHROPIC_API_KEY`
- Check `.env` file exists and is in correct location
- Try setting key programmatically: `os.environ['ANTHROPIC_API_KEY'] = 'key'`

**Issue: "Rate limit exceeded"**
- Implement exponential backoff and retry
- Upgrade API tier if available
- Switch to alternative provider temporarily

**Issue: "Model not found"**
- Verify model identifier is correct
- Check API key has access to requested model
- For Azure: ensure deployment exists with exact name

**Issue: "Timeout errors"**
- Increase `default_config.timeout_seconds`
- Break complex tasks into smaller steps
- Consider using faster model for initial phases

**Issue: "Connection refused (Ollama/Biomni-R0)"**
- Verify local server is running
- Check port is not blocked by firewall
- Confirm `api_base` URL is correct

### Testing Configuration

```python
from biomni.utils import list_available_models, validate_environment

# Check environment setup
status = validate_environment()
print("Environment Status:", status)

# List available models based on configured keys
models = list_available_models()
print("Available Models:", models)

# Test specific model
try:
    from biomni.agent import A1
    agent = A1(path='./data', llm='claude-sonnet-4-20250514')
    agent.go("Print 'Configuration successful!'")
except Exception as e:
    print(f"Configuration test failed: {e}")
```

## Best Practices Summary

1. **For most users:** Start with Claude Sonnet 4 or GPT-4o
2. **For cost sensitivity:** Use GPT-3.5-turbo for exploration, Claude Sonnet 4 for production
3. **For privacy/offline:** Deploy Ollama locally
4. **For complex reasoning:** Use Biomni-R0 if hardware available
5. **For enterprise:** Consider Azure OpenAI or AWS Bedrock
6. **For speed:** Use Groq for rapid iteration

7. **Always:**
   - Set appropriate timeouts
   - Implement error handling and retries
   - Log model and configuration for reproducibility
   - Test configuration before production use
