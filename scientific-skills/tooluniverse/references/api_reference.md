# ToolUniverse Python API Reference

## Core Classes

### ToolUniverse

Main class for interacting with the ToolUniverse ecosystem.

```python
from tooluniverse import ToolUniverse

tu = ToolUniverse()
```

#### Methods

##### `load_tools()`
Load all available tools into the ToolUniverse instance.

```python
tu.load_tools()
```

**Returns:** None

**Side effects:** Loads 600+ tools into memory for discovery and execution.

---

##### `run(tool_config)`
Execute a tool with specified arguments.

**Parameters:**
- `tool_config` (dict): Configuration dictionary with keys:
  - `name` (str): Tool name to execute
  - `arguments` (dict): Tool-specific arguments

**Returns:** Tool-specific output (dict, list, str, or other types)

**Example:**
```python
result = tu.run({
    "name": "OpenTargets_get_associated_targets_by_disease_efoId",
    "arguments": {
        "efoId": "EFO_0000537"
    }
})
```

---

##### `list_tools(limit=None)`
List all available tools or a subset.

**Parameters:**
- `limit` (int, optional): Maximum number of tools to return. If None, returns all tools.

**Returns:** List of tool dictionaries

**Example:**
```python
# List all tools
all_tools = tu.list_tools()

# List first 20 tools
tools = tu.list_tools(limit=20)
```

---

##### `get_tool_info(tool_name)`
Get detailed information about a specific tool.

**Parameters:**
- `tool_name` (str): Name of the tool

**Returns:** Dictionary containing tool metadata, parameters, and documentation

**Example:**
```python
info = tu.get_tool_info("AlphaFold_get_structure")
print(info['description'])
print(info['parameters'])
```

---

## Built-in Discovery Tools

These are special tools that help find other tools in the ecosystem.

### Tool_Finder

Embedding-based semantic search for tools. Requires GPU.

```python
tools = tu.run({
    "name": "Tool_Finder",
    "arguments": {
        "description": "protein structure prediction",
        "limit": 10
    }
})
```

**Parameters:**
- `description` (str): Natural language description of desired functionality
- `limit` (int): Maximum number of tools to return

**Returns:** List of relevant tools with similarity scores

---

### Tool_Finder_LLM

LLM-based semantic search for tools. No GPU required.

```python
tools = tu.run({
    "name": "Tool_Finder_LLM",
    "arguments": {
        "description": "Find tools for RNA sequencing analysis",
        "limit": 10
    }
})
```

**Parameters:**
- `description` (str): Natural language query
- `limit` (int): Maximum number of tools to return

**Returns:** List of relevant tools

---

### Tool_Finder_Keyword

Fast keyword-based search through tool names and descriptions.

```python
tools = tu.run({
    "name": "Tool_Finder_Keyword",
    "arguments": {
        "description": "pathway enrichment",
        "limit": 10
    }
})
```

**Parameters:**
- `description` (str): Keywords to search for
- `limit` (int): Maximum number of tools to return

**Returns:** List of matching tools

---

## Tool Output Hooks

Post-processing hooks for tool results.

### Summarization Hook
```python
result = tu.run({
    "name": "some_tool",
    "arguments": {"param": "value"}
},
hooks={
    "summarize": {
        "format": "brief"  # or "detailed"
    }
})
```

### File Saving Hook
```python
result = tu.run({
    "name": "some_tool",
    "arguments": {"param": "value"}
},
hooks={
    "save_to_file": {
        "filename": "output.json",
        "format": "json"  # or "csv", "txt"
    }
})
```

---

## Model Context Protocol (MCP)

### Starting MCP Server

Command-line interface:
```bash
tooluniverse-smcp
```

This launches an MCP server that exposes all ToolUniverse tools through the Model Context Protocol.

**Configuration:**
- Default port: Automatically assigned
- Protocol: MCP standard
- Authentication: None required for local use

---

## Integration Modules

### OpenRouter Integration

Access 100+ LLMs through OpenRouter API:

```python
from tooluniverse import OpenRouterClient

client = OpenRouterClient(api_key="your_key")
response = client.chat("Analyze this protein sequence", model="anthropic/claude-3-5-sonnet")
```

---

## AI-Tool Interaction Protocol

ToolUniverse uses a standardized protocol for LLM-tool communication:

**Request Format:**
```json
{
  "name": "tool_name",
  "arguments": {
    "param1": "value1",
    "param2": "value2"
  }
}
```

**Response Format:**
```json
{
  "status": "success",
  "data": { ... },
  "metadata": {
    "execution_time": 1.23,
    "tool_version": "1.0.0"
  }
}
```

---

## Error Handling

```python
try:
    result = tu.run({
        "name": "some_tool",
        "arguments": {"param": "value"}
    })
except ToolNotFoundError as e:
    print(f"Tool not found: {e}")
except InvalidArgumentError as e:
    print(f"Invalid arguments: {e}")
except ToolExecutionError as e:
    print(f"Execution failed: {e}")
```

---

## Type Hints

```python
from typing import Dict, List, Any, Optional

def run_tool(
    tu: ToolUniverse,
    tool_name: str,
    arguments: Dict[str, Any]
) -> Any:
    """Execute a tool with type-safe arguments."""
    return tu.run({
        "name": tool_name,
        "arguments": arguments
    })
```

---

## Best Practices

1. **Initialize Once**: Create a single ToolUniverse instance and reuse it
2. **Load Tools Early**: Call `load_tools()` once at startup
3. **Cache Tool Info**: Store frequently used tool information
4. **Error Handling**: Always wrap tool execution in try-except blocks
5. **Type Validation**: Validate argument types before execution
6. **Resource Management**: Consider rate limits for remote APIs
7. **Logging**: Enable logging for production environments
