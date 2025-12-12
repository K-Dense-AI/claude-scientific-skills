#!/usr/bin/env python3
"""
Research Information Lookup Tool
Uses Perplexity's Sonar Pro Search model through OpenRouter for academic research queries.
"""

import os
import json
import requests
import time
from datetime import datetime
from typing import Dict, List, Optional, Any
from urllib.parse import quote


class ResearchLookup:
    """Research information lookup using Perplexity Sonar models via OpenRouter."""

    # Available models
    MODELS = {
        "pro": "perplexity/sonar-pro",  # Fast lookup, cost-effective
        "reasoning": "perplexity/sonar-reasoning-pro",  # Deep analysis with reasoning
    }

    # Keywords that indicate complex queries requiring reasoning model
    REASONING_KEYWORDS = [
        "compare", "contrast", "analyze", "analysis", "evaluate", "critique",
        "versus", "vs", "vs.", "compared to", "differences between", "similarities",
        "meta-analysis", "systematic review", "synthesis", "integrate",
        "mechanism", "why", "how does", "how do", "explain", "relationship",
        "theoretical framework", "implications", "interpret", "reasoning",
        "controversy", "conflicting", "paradox", "debate", "reconcile",
        "pros and cons", "advantages and disadvantages", "trade-off", "tradeoff",
    ]

    def __init__(self, force_model: Optional[str] = None):
        """
        Initialize the research lookup tool.
        
        Args:
            force_model: Optional model override ('pro' or 'reasoning'). 
                        If None, model is auto-selected based on query complexity.
        """
        self.api_key = os.getenv("OPENROUTER_API_KEY")
        if not self.api_key:
            raise ValueError("OPENROUTER_API_KEY environment variable not set")

        self.base_url = "https://openrouter.ai/api/v1"
        self.force_model = force_model
        self.headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json",
            "HTTP-Referer": "https://scientific-writer.local",
            "X-Title": "Scientific Writer Research Tool"
        }

    def _select_model(self, query: str) -> str:
        """
        Select the appropriate model based on query complexity.
        
        Args:
            query: The research query
            
        Returns:
            Model identifier string
        """
        if self.force_model:
            return self.MODELS.get(self.force_model, self.MODELS["reasoning"])
        
        # Check for reasoning keywords (case-insensitive)
        query_lower = query.lower()
        for keyword in self.REASONING_KEYWORDS:
            if keyword in query_lower:
                return self.MODELS["reasoning"]
        
        # Check for multiple questions or complex structure
        question_count = query.count("?")
        if question_count >= 2:
            return self.MODELS["reasoning"]
        
        # Check for very long queries (likely complex)
        if len(query) > 200:
            return self.MODELS["reasoning"]
        
        # Default to pro for simple lookups
        return self.MODELS["pro"]

    def _make_request(self, messages: List[Dict[str, str]], model: str, **kwargs) -> Dict[str, Any]:
        """Make a request to the OpenRouter API with academic search mode."""
        data = {
            "model": model,
            "messages": messages,
            "max_tokens": 8000,
            "temperature": 0.1,  # Low temperature for factual research
            # Perplexity-specific parameters for academic search
            "search_mode": "academic",  # Prioritize scholarly sources (peer-reviewed papers, journals)
            "search_context_size": "high",  # Always use high context for deeper research
            **kwargs
        }

        try:
            response = requests.post(
                f"{self.base_url}/chat/completions",
                headers=self.headers,
                json=data,
                timeout=90  # Increased timeout for academic search
            )
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            raise Exception(f"API request failed: {str(e)}")

    def _format_research_prompt(self, query: str) -> str:
        """Format the query for optimal research results."""
        return f"""You are an expert research assistant. Please provide comprehensive, accurate research information for the following query: "{query}"

IMPORTANT INSTRUCTIONS:
1. Focus on ACADEMIC and SCIENTIFIC sources (peer-reviewed papers, reputable journals, institutional research)
2. Include RECENT information (prioritize 2020-2026 publications)
3. Provide COMPLETE citations with authors, title, journal/conference, year, and DOI when available
4. Structure your response with clear sections and proper attribution
5. Be comprehensive but concise - aim for 800-1200 words
6. Include key findings, methodologies, and implications when relevant
7. Note any controversies, limitations, or conflicting evidence

RESPONSE FORMAT:
- Start with a brief summary (2-3 sentences)
- Present key findings and studies in organized sections
- End with future directions or research gaps if applicable
- Include 5-8 high-quality citations at the end

Remember: This is for academic research purposes. Prioritize accuracy, completeness, and proper attribution."""

    def lookup(self, query: str) -> Dict[str, Any]:
        """Perform a research lookup for the given query."""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        # Select model based on query complexity
        model = self._select_model(query)

        # Format the research prompt
        research_prompt = self._format_research_prompt(query)

        # Prepare messages for the API with system message for academic mode
        messages = [
            {
                "role": "system", 
                "content": "You are an academic research assistant. Focus exclusively on scholarly sources: peer-reviewed journals, academic papers, research institutions, and reputable scientific publications. Prioritize recent academic literature (2020-2026) and provide complete citations with DOIs. Use academic/scholarly search mode."
            },
            {"role": "user", "content": research_prompt}
        ]

        try:
            # Make the API request
            response = self._make_request(messages, model)

            # Extract the response content
            if "choices" in response and len(response["choices"]) > 0:
                choice = response["choices"][0]
                if "message" in choice and "content" in choice["message"]:
                    content = choice["message"]["content"]

                    # Extract citations from API response (Perplexity provides these)
                    api_citations = self._extract_api_citations(response, choice)
                    
                    # Also extract citations from text as fallback
                    text_citations = self._extract_citations_from_text(content)
                    
                    # Combine: prioritize API citations, add text citations if no duplicates
                    citations = api_citations + text_citations

                    return {
                        "success": True,
                        "query": query,
                        "response": content,
                        "citations": citations,
                        "sources": api_citations,  # Separate field for API-provided sources
                        "timestamp": timestamp,
                        "model": model,
                        "usage": response.get("usage", {})
                    }
                else:
                    raise Exception("Invalid response format from API")
            else:
                raise Exception("No response choices received from API")

        except Exception as e:
            return {
                "success": False,
                "query": query,
                "error": str(e),
                "timestamp": timestamp,
                "model": model
            }

    def _extract_api_citations(self, response: Dict[str, Any], choice: Dict[str, Any]) -> List[Dict[str, str]]:
        """Extract citations from Perplexity API response fields."""
        citations = []
        
        # Perplexity returns citations in search_results field (new format)
        # Check multiple possible locations where OpenRouter might place them
        search_results = (
            response.get("search_results") or 
            choice.get("search_results") or
            choice.get("message", {}).get("search_results") or
            []
        )
        
        for result in search_results:
            citation = {
                "type": "source",
                "title": result.get("title", ""),
                "url": result.get("url", ""),
                "date": result.get("date", ""),
            }
            # Add snippet if available (newer API feature)
            if result.get("snippet"):
                citation["snippet"] = result.get("snippet")
            citations.append(citation)
        
        # Also check for legacy citations field (backward compatibility)
        legacy_citations = (
            response.get("citations") or
            choice.get("citations") or
            choice.get("message", {}).get("citations") or
            []
        )
        
        for url in legacy_citations:
            if isinstance(url, str):
                # Legacy format was just URLs
                citations.append({
                    "type": "source",
                    "url": url,
                    "title": "",
                    "date": ""
                })
            elif isinstance(url, dict):
                citations.append({
                    "type": "source",
                    "url": url.get("url", ""),
                    "title": url.get("title", ""),
                    "date": url.get("date", "")
                })
        
        return citations

    def _extract_citations_from_text(self, text: str) -> List[Dict[str, str]]:
        """Extract potential citations from the response text as fallback."""
        import re
        citations = []

        # Look for DOI patterns first (most reliable)
        # Matches: doi:10.xxx, DOI: 10.xxx, https://doi.org/10.xxx
        doi_pattern = r'(?:doi[:\s]*|https?://(?:dx\.)?doi\.org/)(10\.[0-9]{4,}/[^\s\)\]\,\[\<\>]+)'
        doi_matches = re.findall(doi_pattern, text, re.IGNORECASE)
        seen_dois = set()

        for doi in doi_matches:
            # Clean up DOI - remove trailing punctuation and brackets
            doi_clean = doi.strip().rstrip('.,;:)]')
            if doi_clean and doi_clean not in seen_dois:
                seen_dois.add(doi_clean)
                citations.append({
                    "type": "doi",
                    "doi": doi_clean,
                    "url": f"https://doi.org/{doi_clean}"
                })

        # Look for URLs that might be sources
        url_pattern = r'https?://[^\s\)\]\,\<\>\"\']+(?:arxiv\.org|pubmed|ncbi\.nlm\.nih\.gov|nature\.com|science\.org|wiley\.com|springer\.com|ieee\.org|acm\.org)[^\s\)\]\,\<\>\"\']*'
        url_matches = re.findall(url_pattern, text, re.IGNORECASE)
        seen_urls = set()
        
        for url in url_matches:
            url_clean = url.rstrip('.')
            if url_clean not in seen_urls:
                seen_urls.add(url_clean)
                citations.append({
                    "type": "url",
                    "url": url_clean
                })

        return citations

    def batch_lookup(self, queries: List[str], delay: float = 1.0) -> List[Dict[str, Any]]:
        """Perform multiple research lookups with optional delay between requests."""
        results = []

        for i, query in enumerate(queries):
            if i > 0 and delay > 0:
                time.sleep(delay)  # Rate limiting

            result = self.lookup(query)
            results.append(result)

            # Print progress
            print(f"[Research] Completed query {i+1}/{len(queries)}: {query[:50]}...")

        return results

    def get_model_info(self) -> Dict[str, Any]:
        """Get information about available models from OpenRouter."""
        try:
            response = requests.get(
                f"{self.base_url}/models",
                headers=self.headers,
                timeout=30
            )
            response.raise_for_status()
            return response.json()
        except Exception as e:
            return {"error": str(e)}


def main():
    """Command-line interface for testing the research lookup tool."""
    import argparse

    parser = argparse.ArgumentParser(description="Research Information Lookup Tool")
    parser.add_argument("query", nargs="?", help="Research query to look up")
    parser.add_argument("--model-info", action="store_true", help="Show available models")
    parser.add_argument("--batch", nargs="+", help="Run multiple queries")
    parser.add_argument("--force-model", choices=["pro", "reasoning"], 
                        help="Force specific model: 'pro' for fast lookup, 'reasoning' for deep analysis")

    args = parser.parse_args()

    # Check for API key
    if not os.getenv("OPENROUTER_API_KEY"):
        print("Error: OPENROUTER_API_KEY environment variable not set")
        print("Please set it in your .env file or export it:")
        print("  export OPENROUTER_API_KEY='your_openrouter_api_key'")
        return 1

    try:
        research = ResearchLookup(force_model=args.force_model)

        if args.model_info:
            print("Available models from OpenRouter:")
            models = research.get_model_info()
            if "data" in models:
                for model in models["data"]:
                    if "perplexity" in model["id"].lower():
                        print(f"  - {model['id']}: {model.get('name', 'N/A')}")
            return 0

        if not args.query and not args.batch:
            print("Error: No query provided. Use --model-info to see available models.")
            return 1

        if args.batch:
            print(f"Running batch research for {len(args.batch)} queries...")
            results = research.batch_lookup(args.batch)
        else:
            print(f"Researching: {args.query}")
            results = [research.lookup(args.query)]

        # Display results
        for i, result in enumerate(results):
            if result["success"]:
                print(f"\n{'='*80}")
                print(f"Query {i+1}: {result['query']}")
                print(f"Timestamp: {result['timestamp']}")
                print(f"Model: {result['model']}")
                print(f"{'='*80}")
                print(result["response"])

                # Display API-provided sources first (most reliable)
                sources = result.get("sources", [])
                if sources:
                    print(f"\n📚 Sources ({len(sources)}):")
                    for j, source in enumerate(sources):
                        title = source.get("title", "Untitled")
                        url = source.get("url", "")
                        date = source.get("date", "")
                        date_str = f" ({date})" if date else ""
                        print(f"  [{j+1}] {title}{date_str}")
                        if url:
                            print(f"      {url}")

                # Display additional text-extracted citations
                citations = result.get("citations", [])
                text_citations = [c for c in citations if c.get("type") in ("doi", "url")]
                if text_citations:
                    print(f"\n🔗 Additional References ({len(text_citations)}):")
                    for j, citation in enumerate(text_citations):
                        if citation.get("type") == "doi":
                            print(f"  [{j+1}] DOI: {citation.get('doi', '')} - {citation.get('url', '')}")
                        elif citation.get("type") == "url":
                            print(f"  [{j+1}] {citation.get('url', '')}")

                if result.get("usage"):
                    print(f"\nUsage: {result['usage']}")
            else:
                print(f"\nError in query {i+1}: {result['error']}")

        return 0

    except Exception as e:
        print(f"Error: {str(e)}")
        return 1


if __name__ == "__main__":
    exit(main())
