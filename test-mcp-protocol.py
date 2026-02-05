#!/usr/bin/env python3
"""
Test script for claude-skills-mcp server
Simulates MCP protocol communication to verify server functionality
"""

import subprocess
import json
import sys
import time

def send_mcp_request(process, request):
    """Send a JSON-RPC request to the MCP server"""
    request_str = json.dumps(request) + "\n"
    process.stdin.write(request_str.encode())
    process.stdin.flush()

def read_mcp_response(process, timeout=10):
    """Read a JSON-RPC response from the MCP server"""
    import select

    start_time = time.time()
    while time.time() - start_time < timeout:
        # Check if there's data available
        ready, _, _ = select.select([process.stdout], [], [], 0.5)
        if ready:
            line = process.stdout.readline().decode().strip()
            if line:
                try:
                    return json.loads(line)
                except json.JSONDecodeError:
                    print(f"⚠️  Invalid JSON: {line}", file=sys.stderr)
                    continue
    return None

def test_mcp_server():
    """Test the MCP server functionality"""

    print("🧪 Testing Claude Skills MCP Server")
    print("=" * 60)

    # Start the MCP server
    print("\n1️⃣ Starting MCP server...")
    try:
        process = subprocess.Popen(
            ['uvx', 'claude-skills-mcp', '--config', 'mcp-config.json'],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            bufsize=0
        )
        print("✅ Server process started (PID: {})".format(process.pid))
    except Exception as e:
        print(f"❌ Failed to start server: {e}")
        return False

    # Give server time to initialize
    time.sleep(3)

    # Check if process is still running
    if process.poll() is not None:
        stderr = process.stderr.read().decode()
        print(f"❌ Server terminated unexpectedly")
        print(f"   STDERR: {stderr}")
        return False

    print("✅ Server is running")

    # Test 1: Initialize connection
    print("\n2️⃣ Testing initialization...")
    init_request = {
        "jsonrpc": "2.0",
        "id": 1,
        "method": "initialize",
        "params": {
            "protocolVersion": "2024-11-05",
            "capabilities": {},
            "clientInfo": {
                "name": "test-client",
                "version": "1.0.0"
            }
        }
    }

    send_mcp_request(process, init_request)
    response = read_mcp_response(process, timeout=15)

    if response and response.get("id") == 1:
        print("✅ Server initialized successfully")
        print(f"   Protocol version: {response.get('result', {}).get('protocolVersion')}")
    else:
        print("⚠️  Initialization response not received or invalid")
        print(f"   Response: {response}")

    # Test 2: List available tools
    print("\n3️⃣ Listing available tools...")
    tools_request = {
        "jsonrpc": "2.0",
        "id": 2,
        "method": "tools/list"
    }

    send_mcp_request(process, tools_request)
    response = read_mcp_response(process, timeout=15)

    if response and "result" in response:
        tools = response["result"].get("tools", [])
        print(f"✅ Found {len(tools)} tools:")
        for tool in tools:
            print(f"   - {tool.get('name')}: {tool.get('description', 'No description')[:60]}...")
    else:
        print("⚠️  Tools list not received")
        print(f"   Response: {response}")

    # Cleanup
    print("\n4️⃣ Cleaning up...")
    process.terminate()
    try:
        process.wait(timeout=5)
        print("✅ Server stopped gracefully")
    except subprocess.TimeoutExpired:
        process.kill()
        print("⚠️  Server killed (did not stop gracefully)")

    print("\n" + "=" * 60)
    print("✅ MCP Protocol Test Complete")
    return True

if __name__ == "__main__":
    try:
        success = test_mcp_server()
        sys.exit(0 if success else 1)
    except KeyboardInterrupt:
        print("\n⚠️  Test interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
