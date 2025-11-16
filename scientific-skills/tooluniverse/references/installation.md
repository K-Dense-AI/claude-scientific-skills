# ToolUniverse Installation and Setup

## Installation

```bash
uv pip install tooluniverse
```

## Basic Setup

### Python SDK
```python
from tooluniverse import ToolUniverse

# Initialize ToolUniverse
tu = ToolUniverse()

# Load all available tools (600+ scientific tools)
tu.load_tools()
```

## Model Context Protocol (MCP) Setup

ToolUniverse provides native MCP support for integration with Claude Desktop, Claude Code, and other MCP-compatible systems.

### Starting MCP Server
```bash
tooluniverse-smcp
```

This launches an MCP server that exposes ToolUniverse's 600+ tools through the Model Context Protocol.

### Claude Desktop Integration

Add to Claude Desktop configuration (~/.config/Claude/claude_desktop_config.json):
```json
{
  "mcpServers": {
    "tooluniverse": {
      "command": "tooluniverse-smcp"
    }
  }
}
```

### Claude Code Integration

ToolUniverse MCP server works natively with Claude Code through the MCP protocol.

## Integration with Other Platforms

### OpenRouter Integration
ToolUniverse integrates with OpenRouter for access to 100+ LLMs through a single API:
- GPT-5, Claude, Gemini
- Qwen, Deepseek
- Open-source models

### Supported LLM Platforms
- Claude Desktop and Claude Code
- Gemini CLI
- Qwen Code
- ChatGPT API
- GPT Codex CLI

## Requirements

- Python 3.8+
- For Tool_Finder (embedding-based search): GPU recommended
- For Tool_Finder_LLM: No GPU required (uses LLM-based search)

## Verification

Test installation:
```python
from tooluniverse import ToolUniverse

tu = ToolUniverse()
tu.load_tools()

# List first 5 tools to verify setup
tools = tu.list_tools(limit=5)
print(f"Loaded {len(tools)} tools successfully")
```
