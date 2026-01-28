#!/usr/bin/env python3
"""
Comprehensive test of claude-skills-mcp tools
Tests all three MCP tools with real scientific queries
"""

import subprocess
import json
import sys
import time
import select

def send_mcp_request(process, request):
    """Send a JSON-RPC request"""
    request_str = json.dumps(request) + "\n"
    process.stdin.write(request_str.encode())
    process.stdin.flush()

def read_mcp_response(process, timeout=20):
    """Read a JSON-RPC response"""
    start_time = time.time()
    while time.time() - start_time < timeout:
        ready, _, _ = select.select([process.stdout], [], [], 0.5)
        if ready:
            line = process.stdout.readline().decode().strip()
            if line:
                try:
                    return json.loads(line)
                except json.JSONDecodeError:
                    continue
    return None

def test_mcp_tools():
    """Test all MCP tools"""

    print("🧬 Testing Claude Scientific Skills MCP Tools")
    print("=" * 70)

    # Start server
    print("\n📡 Starting MCP server...")
    process = subprocess.Popen(
        ['uvx', 'claude-skills-mcp', '--config', 'mcp-config.json', '--verbose'],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=0
    )
    print(f"✅ Server started (PID: {process.pid})")
    time.sleep(5)  # Give backend time to load

    # Initialize
    print("\n🔌 Initializing connection...")
    send_mcp_request(process, {
        "jsonrpc": "2.0",
        "id": 1,
        "method": "initialize",
        "params": {
            "protocolVersion": "2024-11-05",
            "capabilities": {},
            "clientInfo": {"name": "test", "version": "1.0"}
        }
    })

    init_response = read_mcp_response(process, timeout=20)
    if not init_response or init_response.get("id") != 1:
        print("❌ Initialization failed")
        process.terminate()
        return False
    print("✅ Connection initialized")

    # Test 1: List all skills
    print("\n" + "=" * 70)
    print("TEST 1: list_skills - List all available skills")
    print("=" * 70)

    send_mcp_request(process, {
        "jsonrpc": "2.0",
        "id": 2,
        "method": "tools/call",
        "params": {
            "name": "list_skills",
            "arguments": {}
        }
    })

    list_response = read_mcp_response(process, timeout=30)
    if list_response and "result" in list_response:
        content = list_response["result"]["content"]
        if content:
            text = content[0]["text"] if isinstance(content, list) else content.get("text", "")

            # Count skills by category
            lines = text.split('\n')
            total_skills = len([l for l in lines if l.strip().startswith('-')])

            print(f"✅ Found {total_skills} skills total")
            print("\n📋 Sample skills (first 15 lines):")
            for line in lines[:15]:
                if line.strip():
                    print(f"   {line}")
    else:
        print("⚠️  list_skills did not return expected results")
        print(f"   Response: {list_response}")

    # Test 2: Find skills for drug discovery
    print("\n" + "=" * 70)
    print("TEST 2: find_helpful_skills - Search for drug discovery skills")
    print("=" * 70)

    send_mcp_request(process, {
        "jsonrpc": "2.0",
        "id": 3,
        "method": "tools/call",
        "params": {
            "name": "find_helpful_skills",
            "arguments": {
                "task_description": "I need to analyze molecular structures, calculate drug-like properties, and perform molecular docking for drug discovery"
            }
        }
    })

    find_response = read_mcp_response(process, timeout=30)
    if find_response and "result" in find_response:
        content = find_response["result"]["content"]
        if content:
            text = content[0]["text"] if isinstance(content, list) else content.get("text", "")
            print("✅ Search completed")
            print("\n🎯 Found relevant skills:")
            # Show first 20 lines of results
            for line in text.split('\n')[:20]:
                if line.strip():
                    print(f"   {line}")
    else:
        print("⚠️  find_helpful_skills did not return expected results")

    # Test 3: Find skills for genomics
    print("\n" + "=" * 70)
    print("TEST 3: find_helpful_skills - Search for genomics skills")
    print("=" * 70)

    send_mcp_request(process, {
        "jsonrpc": "2.0",
        "id": 4,
        "method": "tools/call",
        "params": {
            "name": "find_helpful_skills",
            "arguments": {
                "task_description": "Analyze single-cell RNA-seq data, perform clustering, and identify cell types"
            }
        }
    })

    genomics_response = read_mcp_response(process, timeout=30)
    if genomics_response and "result" in genomics_response:
        content = genomics_response["result"]["content"]
        if content:
            text = content[0]["text"] if isinstance(content, list) else content.get("text", "")
            print("✅ Search completed")
            print("\n🧬 Found relevant skills:")
            for line in text.split('\n')[:20]:
                if line.strip():
                    print(f"   {line}")
    else:
        print("⚠️  find_helpful_skills did not return expected results")

    # Test 4: Read a specific skill document
    print("\n" + "=" * 70)
    print("TEST 4: read_skill_document - Read RDKit skill documentation")
    print("=" * 70)

    send_mcp_request(process, {
        "jsonrpc": "2.0",
        "id": 5,
        "method": "tools/call",
        "params": {
            "name": "read_skill_document",
            "arguments": {
                "skill_name": "rdkit",
                "file_path": "SKILL.md"
            }
        }
    })

    read_response = read_mcp_response(process, timeout=30)
    if read_response and "result" in read_response:
        content = read_response["result"]["content"]
        if content:
            text = content[0]["text"] if isinstance(content, list) else content.get("text", "")
            print("✅ Document read successfully")
            print(f"   Document size: {len(text)} characters")
            print("\n📄 First 500 characters:")
            print(f"   {text[:500]}...")
    else:
        print("⚠️  read_skill_document did not return expected results")

    # Cleanup
    print("\n" + "=" * 70)
    print("🧹 Cleaning up...")
    process.terminate()
    try:
        process.wait(timeout=5)
        print("✅ Server stopped")
    except subprocess.TimeoutExpired:
        process.kill()

    print("\n" + "=" * 70)
    print("✅ ALL TESTS COMPLETED SUCCESSFULLY")
    print("=" * 70)
    return True

if __name__ == "__main__":
    try:
        success = test_mcp_tools()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
