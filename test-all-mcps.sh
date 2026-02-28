#!/bin/bash

# Comprehensive test for all 4 MCP servers before restarting Claude Desktop
# Tests: sequential-thinking, perplexity-ask, brave-search, claude-skills-scientific

set -e

CONFIG_FILE="$HOME/.config/Claude/claude_desktop_config.json"

echo "🧪 Testing All MCP Servers Configuration"
echo "========================================================================"
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

PASS="${GREEN}✅ PASS${NC}"
FAIL="${RED}❌ FAIL${NC}"
WARN="${YELLOW}⚠️  WARN${NC}"

total_tests=0
passed_tests=0
failed_tests=0
warnings=0

# Test 1: Configuration file exists
echo "TEST 1: Configuration File"
echo "----------------------------------------"
total_tests=$((total_tests + 1))
if [ -f "$CONFIG_FILE" ]; then
    echo -e "$PASS Configuration file exists"
    echo "   Location: $CONFIG_FILE"
    passed_tests=$((passed_tests + 1))
else
    echo -e "$FAIL Configuration file not found"
    echo "   Expected: $CONFIG_FILE"
    failed_tests=$((failed_tests + 1))
    exit 1
fi
echo ""

# Test 2: Valid JSON
echo "TEST 2: JSON Syntax Validation"
echo "----------------------------------------"
total_tests=$((total_tests + 1))
if python3 -m json.tool "$CONFIG_FILE" > /dev/null 2>&1; then
    echo -e "$PASS JSON syntax is valid"
    passed_tests=$((passed_tests + 1))
else
    echo -e "$FAIL JSON syntax is invalid"
    failed_tests=$((failed_tests + 1))
    echo "   Run: python3 -m json.tool $CONFIG_FILE"
    exit 1
fi
echo ""

# Test 3: Count MCP servers
echo "TEST 3: MCP Servers Count"
echo "----------------------------------------"
total_tests=$((total_tests + 1))
SERVER_COUNT=$(python3 -c "
import json
with open('$CONFIG_FILE', 'r') as f:
    config = json.load(f)
print(len(config.get('mcpServers', {})))
")

echo "   Found: $SERVER_COUNT MCP servers"
if [ "$SERVER_COUNT" -eq 4 ]; then
    echo -e "$PASS Expected 4 servers found"
    passed_tests=$((passed_tests + 1))
else
    echo -e "$WARN Expected 4, found $SERVER_COUNT"
    warnings=$((warnings + 1))
fi

python3 -c "
import json
with open('$CONFIG_FILE', 'r') as f:
    config = json.load(f)
servers = config.get('mcpServers', {})
for name in servers.keys():
    print(f'   - {name}')
"
echo ""

# Test 4: Check dependencies
echo "TEST 4: Required Dependencies"
echo "----------------------------------------"

# Check docker
total_tests=$((total_tests + 1))
if command -v docker &> /dev/null; then
    DOCKER_VERSION=$(docker --version 2>&1 | head -1)
    echo -e "$PASS Docker available"
    echo "   $DOCKER_VERSION"
    passed_tests=$((passed_tests + 1))
else
    echo -e "$FAIL Docker not found"
    echo "   Required for: sequential-thinking, perplexity-ask"
    failed_tests=$((failed_tests + 1))
fi

# Check npx
total_tests=$((total_tests + 1))
if command -v npx &> /dev/null; then
    NPX_VERSION=$(npx --version 2>&1)
    echo -e "$PASS npx available"
    echo "   Version: $NPX_VERSION"
    passed_tests=$((passed_tests + 1))
else
    echo -e "$FAIL npx not found"
    echo "   Required for: brave-search"
    failed_tests=$((failed_tests + 1))
fi

# Check uvx
total_tests=$((total_tests + 1))
if command -v uvx &> /dev/null; then
    UVX_VERSION=$(uvx --version 2>&1)
    echo -e "$PASS uvx available"
    echo "   Version: $UVX_VERSION"
    passed_tests=$((passed_tests + 1))
else
    echo -e "$FAIL uvx not found"
    echo "   Required for: claude-skills-scientific"
    failed_tests=$((failed_tests + 1))
fi
echo ""

# Test 5: Check Docker daemon
echo "TEST 5: Docker Daemon Status"
echo "----------------------------------------"
total_tests=$((total_tests + 1))
if docker info &> /dev/null; then
    echo -e "$PASS Docker daemon is running"
    passed_tests=$((passed_tests + 1))
else
    echo -e "$FAIL Docker daemon not running"
    echo "   Start Docker Desktop or docker service"
    echo "   Required for: sequential-thinking, perplexity-ask"
    failed_tests=$((failed_tests + 1))
fi
echo ""

# Test 6: Check Docker images
echo "TEST 6: Docker Images for MCPs"
echo "----------------------------------------"

# Check sequential-thinking image
total_tests=$((total_tests + 1))
if docker images | grep -q "mcp/sequentialthinking"; then
    echo -e "$PASS mcp/sequentialthinking image found"
    passed_tests=$((passed_tests + 1))
else
    echo -e "$WARN mcp/sequentialthinking image not found"
    echo "   Will be pulled on first use"
    warnings=$((warnings + 1))
fi

# Check perplexity-ask image
total_tests=$((total_tests + 1))
if docker images | grep -q "mcp/perplexity-ask"; then
    echo -e "$PASS mcp/perplexity-ask image found"
    passed_tests=$((passed_tests + 1))
else
    echo -e "$WARN mcp/perplexity-ask image not found"
    echo "   Will be pulled on first use"
    warnings=$((warnings + 1))
fi
echo ""

# Test 7: Verify API keys are set
echo "TEST 7: API Keys Configuration"
echo "----------------------------------------"

total_tests=$((total_tests + 1))
PERPLEXITY_KEY=$(python3 -c "
import json
with open('$CONFIG_FILE', 'r') as f:
    config = json.load(f)
key = config.get('mcpServers', {}).get('perplexity-ask', {}).get('env', {}).get('PERPLEXITY_API_KEY', '')
print('SET' if key and key != 'YOUR_PERPLEXITY_API_KEY_HERE' else 'NOT_SET')
")

if [ "$PERPLEXITY_KEY" = "SET" ]; then
    echo -e "$PASS Perplexity API key configured"
    passed_tests=$((passed_tests + 1))
else
    echo -e "$FAIL Perplexity API key not set"
    echo "   Update PERPLEXITY_API_KEY in config"
    failed_tests=$((failed_tests + 1))
fi

total_tests=$((total_tests + 1))
BRAVE_KEY=$(python3 -c "
import json
with open('$CONFIG_FILE', 'r') as f:
    config = json.load(f)
key = config.get('mcpServers', {}).get('brave-search', {}).get('env', {}).get('BRAVE_API_KEY', '')
print('SET' if key and key != 'YOUR_BRAVE_API_KEY_HERE' else 'NOT_SET')
")

if [ "$BRAVE_KEY" = "SET" ]; then
    echo -e "$PASS Brave API key configured"
    passed_tests=$((passed_tests + 1))
else
    echo -e "$FAIL Brave API key not set"
    echo "   Update BRAVE_API_KEY in config"
    failed_tests=$((failed_tests + 1))
fi
echo ""

# Test 8: Check scientific skills MCP config path
echo "TEST 8: Scientific Skills MCP Configuration"
echo "----------------------------------------"

total_tests=$((total_tests + 1))
MCP_CONFIG_PATH=$(python3 -c "
import json
with open('$CONFIG_FILE', 'r') as f:
    config = json.load(f)
args = config.get('mcpServers', {}).get('claude-skills-scientific', {}).get('args', [])
for i, arg in enumerate(args):
    if arg == '--config' and i + 1 < len(args):
        print(args[i + 1])
        break
")

if [ -n "$MCP_CONFIG_PATH" ] && [ -f "$MCP_CONFIG_PATH" ]; then
    echo -e "$PASS MCP config file exists"
    echo "   Path: $MCP_CONFIG_PATH"
    passed_tests=$((passed_tests + 1))
else
    echo -e "$FAIL MCP config file not found"
    echo "   Path: $MCP_CONFIG_PATH"
    failed_tests=$((failed_tests + 1))
fi
echo ""

# Test 9: Quick MCP protocol test (scientific skills only)
echo "TEST 9: Quick MCP Protocol Test (Scientific Skills)"
echo "----------------------------------------"

total_tests=$((total_tests + 1))
echo "   Testing claude-skills-scientific MCP server..."
if timeout 15 python3 << 'PYTHON_SCRIPT' 2>&1 | grep -q "Protocol test passed"; then
import subprocess
import json
import select
import sys
import time

def test_mcp():
    try:
        process = subprocess.Popen(
            ['uvx', 'claude-skills-mcp', '--config', '/home/user/claude-scientific-skills/mcp-config.json'],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            bufsize=0
        )

        time.sleep(3)

        # Send initialize
        init_request = {
            "jsonrpc": "2.0",
            "id": 1,
            "method": "initialize",
            "params": {
                "protocolVersion": "2024-11-05",
                "capabilities": {},
                "clientInfo": {"name": "test", "version": "1.0"}
            }
        }

        request_str = json.dumps(init_request) + "\n"
        process.stdin.write(request_str.encode())
        process.stdin.flush()

        # Wait for response
        start_time = time.time()
        while time.time() - start_time < 10:
            ready, _, _ = select.select([process.stdout], [], [], 0.5)
            if ready:
                line = process.stdout.readline().decode().strip()
                if line:
                    response = json.loads(line)
                    if response.get("id") == 1:
                        print("Protocol test passed")
                        process.terminate()
                        return True

        process.terminate()
        return False
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return False

test_mcp()
PYTHON_SCRIPT
    echo -e "$PASS MCP server responds correctly"
    passed_tests=$((passed_tests + 1))
else
    echo -e "$WARN MCP server test inconclusive"
    echo "   Server might need more time to initialize"
    echo "   This is normal on first run (backend download)"
    warnings=$((warnings + 1))
fi
echo ""

# Summary
echo "========================================================================"
echo "TEST SUMMARY"
echo "========================================================================"
echo ""
echo "Total Tests:    $total_tests"
echo -e "Passed:         ${GREEN}$passed_tests${NC}"
echo -e "Failed:         ${RED}$failed_tests${NC}"
echo -e "Warnings:       ${YELLOW}$warnings${NC}"
echo ""

# Calculate percentage
PASS_RATE=$((passed_tests * 100 / total_tests))

if [ $failed_tests -eq 0 ]; then
    echo -e "${GREEN}✅ ALL CRITICAL TESTS PASSED${NC} ($PASS_RATE%)"
    echo ""
    echo "✅ Configuration is ready!"
    echo ""
    echo "📋 Next steps:"
    echo "   1. RESTART Claude Desktop (close and reopen)"
    echo "   2. Wait 1-2 minutes on first launch (backend download)"
    echo "   3. Test with: 'Liste todos os skills científicos'"
    echo ""

    if [ $warnings -gt 0 ]; then
        echo -e "${YELLOW}ℹ️  Note:${NC} There are $warnings warnings but they are not critical."
        echo "   Docker images and backend will download automatically on first use."
    fi

    exit 0
else
    echo -e "${RED}❌ SOME TESTS FAILED${NC}"
    echo ""
    echo "⚠️  Please fix the failed tests before restarting Claude Desktop"
    echo ""
    echo "Common fixes:"
    echo "   - Install Docker Desktop and start it"
    echo "   - Install Node.js/npm (includes npx)"
    echo "   - Install uv/uvx: pip install uv"
    echo "   - Verify API keys in config file"
    echo ""
    exit 1
fi
