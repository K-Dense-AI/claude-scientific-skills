#!/usr/bin/env python3
"""
Research Lookup Tool for Claude Code
Performs research queries using Perplexity Sonar Pro Search via OpenRouter.
"""

import os
import sys
import json
from typing import Dict, List, Optional

# Import the main research lookup class
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scripts'))
from research_lookup import ResearchLookup


def format_response(result: Dict) -> str:
    """Format the research result for display."""
    if not result["success"]:
        return f"❌ Research lookup failed: {result['error']}"

    response = result["response"]
    citations = result["citations"]

    # Format the output for Claude Code
    output = f"""🔍 **Research Results**

**Query:** {result['query']}
**Model:** {result['model']}
**Timestamp:** {result['timestamp']}

---

{response}

"""

    if citations:
        output += f"\n**Extracted Citations ({len(citations)}):**\n"
        for i, citation in enumerate(citations, 1):
            if citation.get("doi"):
                output += f"{i}. DOI: {citation['doi']}\n"
            elif citation.get("authors") and citation.get("year"):
                output += f"{i}. {citation['authors']} ({citation['year']})\n"
            else:
                output += f"{i}. {citation}\n"

    if result.get("usage"):
        usage = result["usage"]
        output += f"\n**Usage:** {usage.get('total_tokens', 'N/A')} tokens"

    return output


def main():
    """Main entry point for Claude Code tool."""
    # Check for API key
    if not os.getenv("OPENROUTER_API_KEY"):
        print("❌ Error: OPENROUTER_API_KEY environment variable not set")
        print("Please set it in your .env file or export it:")
        print("  export OPENROUTER_API_KEY='your_openrouter_api_key'")
        return 1

    # Get query from command line arguments
    if len(sys.argv) < 2:
        print("❌ Error: No query provided")
        print("Usage: python lookup.py 'your research query here'")
        return 1

    query = " ".join(sys.argv[1:])

    try:
        # Initialize research tool
        research = ResearchLookup()

        # Perform lookup
        print(f"🔍 Researching: {query}")
        result = research.lookup(query)

        # Format and output result
        formatted_output = format_response(result)
        print(formatted_output)

        # Return success code
        return 0 if result["success"] else 1

    except Exception as e:
        print(f"❌ Error: {str(e)}")
        return 1


if __name__ == "__main__":
    exit(main())
