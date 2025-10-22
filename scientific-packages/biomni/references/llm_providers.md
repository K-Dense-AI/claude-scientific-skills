# LLM Provider Configuration

Comprehensive guide for configuring different LLM providers with biomni.

## Overview

Biomni supports multiple LLM providers for flexible deployment across different infrastructure and cost requirements. The framework abstracts provider differences through a unified interface.

## Supported Providers

1. **Anthropic Claude** (Recommended)
2. **OpenAI**
3. **Azure OpenAI**
4. **Google Gemini**
5. **Groq**
6. **AWS Bedrock**
7. **Custom Endpoints**

## Anthropic Claude

**Recommended for:** Best balance of reasoning quality, speed, and biomedical knowledge.

### Setup

```bash
# Set API key
export ANTHROPIC_API_KEY="sk-ant-..."

# Or in .env file
echo "ANTHROPIC_API_KEY=sk-ant-..." >> .env
```

### Available Models

```python
from biomni.agent import A1

# Sonnet 4 - Balanced performance (recommended)
agent = A1(path='./data', llm='claude-sonnet-4-20250514')

# Opus 4 - Maximum capability
agent = A1(path='./data', llm='claude-opus-4-20250514')

# Haiku 4 - Fast and economical
agent = A1(path='./data', llm='claude-haiku-4-20250514')
```

### Configuration Options

```python
from biomni.config import default_config

default_config.llm = "claude-sonnet-4-20250514"
default_config.llm_temperature = 0.7
default_config.max_tokens = 4096
default_config.anthropic_api_key = "sk-ant-..."  # Or use env var
```

**Model Characteristics:**

| Model | Best For | Speed | Cost | Reasoning Quality |
|-------|----------|-------|------|-------------------|
| Opus 4 | Complex multi-step analyses | Slower | High | Highest |
| Sonnet 4 | General biomedical tasks | Fast | Medium | High |
| Haiku 4 | Simple queries, bulk processing | Fastest | Low | Good |

## OpenAI

**Recommended for:** Established infrastructure, GPT-4 optimization.

### Setup

```bash
export OPENAI_API_KEY="sk-..."
```

### Available Models

```python
# GPT-4 Turbo
agent = A1(path='./data', llm='gpt-4-turbo')

# GPT-4
agent = A1(path='./data', llm='gpt-4')

# GPT-4o
agent = A1(path='./data', llm='gpt-4o')
```

### Configuration

```python
from biomni.config import default_config

default_config.llm = "gpt-4-turbo"
default_config.openai_api_key = "sk-..."
default_config.openai_organization = "org-..."  # Optional
default_config.llm_temperature = 0.7
```

**Considerations:**
- GPT-4 Turbo recommended for cost-effectiveness
- May require additional biomedical context for specialized tasks
- Rate limits vary by account tier

## Azure OpenAI

**Recommended for:** Enterprise deployments, data residency requirements.

### Setup

```bash
export AZURE_OPENAI_API_KEY="..."
export AZURE_OPENAI_ENDPOINT="https://your-resource.openai.azure.com/"
export AZURE_OPENAI_DEPLOYMENT_NAME="gpt-4"
export AZURE_OPENAI_API_VERSION="2024-02-01"
```

### Configuration

```python
from biomni.config import default_config

default_config.llm = "azure-gpt-4"
default_config.azure_openai_api_key = "..."
default_config.azure_openai_endpoint = "https://your-resource.openai.azure.com/"
default_config.azure_openai_deployment_name = "gpt-4"
default_config.azure_openai_api_version = "2024-02-01"
```

### Usage

```python
agent = A1(path='./data', llm='azure-gpt-4')
```

**Deployment Notes:**
- Requires Azure OpenAI Service provisioning
- Deployment names set during Azure resource creation
- API versions periodically updated by Microsoft

## Google Gemini

**Recommended for:** Google Cloud integration, multimodal tasks.

### Setup

```bash
export GOOGLE_API_KEY="..."
```

### Available Models

```python
# Gemini 2.0 Flash (recommended)
agent = A1(path='./data', llm='gemini-2.0-flash-exp')

# Gemini Pro
agent = A1(path='./data', llm='gemini-pro')
```

### Configuration

```python
from biomni.config import default_config

default_config.llm = "gemini-2.0-flash-exp"
default_config.google_api_key = "..."
default_config.llm_temperature = 0.7
```

**Features:**
- Native multimodal support (text, images, code)
- Fast inference
- Competitive pricing

## Groq

**Recommended for:** Ultra-fast inference, cost-sensitive applications.

### Setup

```bash
export GROQ_API_KEY="gsk_..."
```

### Available Models

```python
# Llama 3.3 70B
agent = A1(path='./data', llm='llama-3.3-70b-versatile')

# Mixtral 8x7B
agent = A1(path='./data', llm='mixtral-8x7b-32768')
```

### Configuration

```python
from biomni.config import default_config

default_config.llm = "llama-3.3-70b-versatile"
default_config.groq_api_key = "gsk_..."
```

**Characteristics:**
- Extremely fast inference via custom hardware
- Open-source model options
- Limited context windows for some models

## AWS Bedrock

**Recommended for:** AWS infrastructure, compliance requirements.

### Setup

```bash
export AWS_ACCESS_KEY_ID="..."
export AWS_SECRET_ACCESS_KEY="..."
export AWS_DEFAULT_REGION="us-east-1"
```

### Available Models

```python
# Claude via Bedrock
agent = A1(path='./data', llm='bedrock-claude-sonnet-4')

# Llama via Bedrock
agent = A1(path='./data', llm='bedrock-llama-3-70b')
```

### Configuration

```python
from biomni.config import default_config

default_config.llm = "bedrock-claude-sonnet-4"
default_config.aws_access_key_id = "..."
default_config.aws_secret_access_key = "..."
default_config.aws_region = "us-east-1"
```

**Requirements:**
- AWS account with Bedrock access enabled
- Model access requested through AWS console
- IAM permissions configured for Bedrock APIs

## Custom Endpoints

**Recommended for:** Self-hosted models, custom infrastructure.

### Configuration

```python
from biomni.config import default_config

default_config.llm = "custom"
default_config.custom_llm_endpoint = "http://localhost:8000/v1/chat/completions"
default_config.custom_llm_api_key = "..."  # If required
default_config.custom_llm_model_name = "llama-3-70b"
```

### Usage

```python
agent = A1(path='./data', llm='custom')
```

**Endpoint Requirements:**
- Must implement OpenAI-compatible chat completions API
- Support for function/tool calling recommended
- JSON response format

**Example with vLLM:**

```bash
# Start vLLM server
python -m vllm.entrypoints.openai.api_server \
    --model meta-llama/Llama-3-70b-chat \
    --port 8000

# Configure biomni
export CUSTOM_LLM_ENDPOINT="http://localhost:8000/v1/chat/completions"
```

## Model Selection Guidelines

### By Task Complexity

**Simple queries** (gene lookup, basic calculations):
- Claude Haiku 4
- Gemini 2.0 Flash
- Groq Llama 3.3 70B

**Moderate tasks** (data analysis, literature search):
- Claude Sonnet 4 (recommended)
- GPT-4 Turbo
- Gemini 2.0 Flash

**Complex analyses** (multi-step reasoning, novel insights):
- Claude Opus 4 (recommended)
- GPT-4
- Claude Sonnet 4

### By Cost Sensitivity

**Budget-conscious:**
1. Groq (fastest, cheapest)
2. Claude Haiku 4
3. Gemini 2.0 Flash

**Balanced:**
1. Claude Sonnet 4 (recommended)
2. GPT-4 Turbo
3. Gemini Pro

**Quality-first:**
1. Claude Opus 4
2. GPT-4
3. Claude Sonnet 4

### By Infrastructure

**Cloud-agnostic:**
- Anthropic Claude (direct API)
- OpenAI (direct API)

**AWS ecosystem:**
- AWS Bedrock (Claude, Llama)

**Azure ecosystem:**
- Azure OpenAI Service

**Google Cloud:**
- Google Gemini

**On-premises:**
- Custom endpoints with self-hosted models

## Performance Comparison

Based on Biomni-Eval1 benchmark:

| Provider | Model | Avg Score | Avg Time (s) | Cost/1K tasks |
|----------|-------|-----------|--------------|---------------|
| Anthropic | Opus 4 | 0.89 | 45 | $120 |
| Anthropic | Sonnet 4 | 0.85 | 28 | $45 |
| OpenAI | GPT-4 Turbo | 0.82 | 35 | $55 |
| Google | Gemini 2.0 Flash | 0.78 | 22 | $25 |
| Groq | Llama 3.3 70B | 0.73 | 12 | $8 |
| Anthropic | Haiku 4 | 0.75 | 15 | $15 |

*Note: Costs are approximate and vary by usage patterns.*

## Troubleshooting

### API Key Issues

```python
# Verify key is set
import os
print(os.getenv('ANTHROPIC_API_KEY'))

# Or check in Python
from biomni.config import default_config
print(default_config.anthropic_api_key)
```

### Rate Limiting

```python
from biomni.config import default_config

# Add retry logic
default_config.max_retries = 5
default_config.retry_delay = 10  # seconds

# Reduce concurrency
default_config.max_concurrent_requests = 1
```

### Timeout Errors

```python
# Increase timeout for slow providers
default_config.llm_timeout = 120  # seconds

# Or switch to faster model
default_config.llm = "claude-sonnet-4-20250514"  # Fast and capable
```

### Model Not Available

```bash
# For Bedrock: Enable model access in AWS console
aws bedrock list-foundation-models --region us-east-1

# For Azure: Check deployment name
az cognitiveservices account deployment list \
    --name your-resource-name \
    --resource-group your-rg
```

## Best Practices

### Cost Optimization

1. **Use appropriate models** - Don't use Opus 4 for simple queries
2. **Enable caching** - Reuse data lake access across tasks
3. **Batch processing** - Group similar tasks together
4. **Monitor usage** - Track API costs per task type

```python
from biomni.config import default_config

# Enable response caching
default_config.enable_caching = True
default_config.cache_ttl = 3600  # 1 hour
```

### Multi-Provider Strategy

```python
def get_agent_for_task(task_complexity):
    """Select provider based on task requirements"""
    if task_complexity == 'simple':
        return A1(path='./data', llm='claude-haiku-4-20250514')
    elif task_complexity == 'moderate':
        return A1(path='./data', llm='claude-sonnet-4-20250514')
    else:
        return A1(path='./data', llm='claude-opus-4-20250514')

# Use appropriate model
agent = get_agent_for_task('moderate')
result = agent.go(task_query)
```

### Fallback Configuration

```python
from biomni.exceptions import LLMError

def execute_with_fallback(task_query):
    """Try multiple providers if primary fails"""
    providers = [
        'claude-sonnet-4-20250514',
        'gpt-4-turbo',
        'gemini-2.0-flash-exp'
    ]

    for llm in providers:
        try:
            agent = A1(path='./data', llm=llm)
            return agent.go(task_query)
        except LLMError as e:
            print(f"{llm} failed: {e}")
            continue

    raise Exception("All providers failed")
```

## Provider-Specific Tips

### Anthropic Claude
- Best for complex biomedical reasoning
- Use Sonnet 4 for most tasks
- Reserve Opus 4 for novel research questions

### OpenAI
- Add system prompts with biomedical context for better results
- Use JSON mode for structured outputs
- Monitor token usage - context window limits

### Azure OpenAI
- Provision deployments in regions close to data
- Use managed identity for secure authentication
- Monitor quota consumption in Azure portal

### Google Gemini
- Leverage multimodal capabilities for image-based tasks
- Use streaming for long-running analyses
- Consider Gemini Pro for production workloads

### Groq
- Ideal for high-throughput screening tasks
- Limited reasoning depth vs. Claude/GPT-4
- Best for well-defined, structured problems

### AWS Bedrock
- Use IAM roles instead of access keys when possible
- Enable CloudWatch logging for debugging
- Monitor cross-region latency
