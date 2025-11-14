#!/usr/bin/env python3
"""
Example script demonstrating tool discovery in ToolUniverse.

This script shows how to search for tools using different methods:
- Embedding-based search (Tool_Finder)
- LLM-based search (Tool_Finder_LLM)
- Keyword-based search (Tool_Finder_Keyword)
"""

from tooluniverse import ToolUniverse


def main():
    # Initialize ToolUniverse
    print("Initializing ToolUniverse...")
    tu = ToolUniverse()
    tu.load_tools()
    print(f"Loaded {len(tu.list_tools())} tools\n")

    # Example 1: Keyword-based search (fastest)
    print("=" * 60)
    print("Example 1: Keyword Search for Disease-Target Tools")
    print("=" * 60)

    tools = tu.run({
        "name": "Tool_Finder_Keyword",
        "arguments": {
            "description": "disease target associations",
            "limit": 5
        }
    })

    print(f"Found {len(tools)} tools:")
    for idx, tool in enumerate(tools, 1):
        print(f"\n{idx}. {tool['name']}")
        print(f"   Description: {tool['description']}")

    # Example 2: LLM-based search (no GPU required)
    print("\n" + "=" * 60)
    print("Example 2: LLM Search for Protein Structure Tools")
    print("=" * 60)

    tools = tu.run({
        "name": "Tool_Finder_LLM",
        "arguments": {
            "description": "Find tools for predicting protein structures from sequences",
            "limit": 5
        }
    })

    print(f"Found {len(tools)} tools:")
    for idx, tool in enumerate(tools, 1):
        print(f"\n{idx}. {tool['name']}")
        print(f"   Description: {tool['description']}")

    # Example 3: Search by specific domain
    print("\n" + "=" * 60)
    print("Example 3: Search for Cheminformatics Tools")
    print("=" * 60)

    tools = tu.run({
        "name": "Tool_Finder_Keyword",
        "arguments": {
            "description": "molecular docking SMILES compound",
            "limit": 5
        }
    })

    print(f"Found {len(tools)} tools:")
    for idx, tool in enumerate(tools, 1):
        print(f"\n{idx}. {tool['name']}")
        print(f"   Description: {tool['description']}")

    # Example 4: Get detailed tool information
    print("\n" + "=" * 60)
    print("Example 4: Get Tool Details")
    print("=" * 60)

    if tools:
        tool_name = tools[0]['name']
        print(f"Getting details for: {tool_name}")

        tool_info = tu.get_tool_info(tool_name)
        print(f"\nTool: {tool_info['name']}")
        print(f"Description: {tool_info['description']}")
        print(f"Parameters: {tool_info.get('parameters', 'No parameters listed')}")


if __name__ == "__main__":
    main()
