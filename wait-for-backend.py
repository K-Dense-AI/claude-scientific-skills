#!/usr/bin/env python3
"""
Wait for MCP backend to download and initialize completely
"""

import subprocess
import json
import sys
import time
import select

def send_request(process, request):
    request_str = json.dumps(request) + "\n"
    process.stdin.write(request_str.encode())
    process.stdin.flush()

def read_response(process, timeout=10):
    start_time = time.time()
    while time.time() - start_time < timeout:
        ready, _, _ = select.select([process.stdout], [], [], 0.5)
        if ready:
            line = process.stdout.readline().decode().strip()
            if line:
                try:
                    return json.loads(line)
                except:
                    continue
    return None

def check_backend_ready():
    """Check if backend is ready by testing list_skills"""

    process = subprocess.Popen(
        ['uvx', 'claude-skills-mcp', '--config', 'mcp-config.json'],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=0
    )

    time.sleep(3)

    # Initialize
    send_request(process, {
        "jsonrpc": "2.0",
        "id": 1,
        "method": "initialize",
        "params": {
            "protocolVersion": "2024-11-05",
            "capabilities": {},
            "clientInfo": {"name": "test", "version": "1.0"}
        }
    })

    init_resp = read_response(process, timeout=15)
    if not init_resp:
        process.terminate()
        return False

    # Try list_skills
    send_request(process, {
        "jsonrpc": "2.0",
        "id": 2,
        "method": "tools/call",
        "params": {
            "name": "list_skills",
            "arguments": {}
        }
    })

    list_resp = read_response(process, timeout=20)
    process.terminate()
    process.wait(timeout=5)

    if list_resp and "result" in list_resp:
        content = list_resp["result"]["content"]
        if content:
            text = content[0]["text"] if isinstance(content, list) else content.get("text", "")
            return "BACKEND LOADING" not in text

    return False

print("⏳ Waiting for MCP backend to download and initialize...")
print("   This happens only once and takes 30-120 seconds")
print("   Checking every 15 seconds...")
print()

max_attempts = 12  # 3 minutes total
attempt = 0

while attempt < max_attempts:
    attempt += 1
    print(f"   Attempt {attempt}/{max_attempts}...", end=" ", flush=True)

    if check_backend_ready():
        print("✅")
        print()
        print("🎉 Backend is ready!")
        print("   Skills are now loaded and searchable")
        sys.exit(0)
    else:
        print("⏳ (still loading)")
        if attempt < max_attempts:
            time.sleep(15)

print()
print("⚠️  Backend did not complete within expected time")
print("   This might indicate:")
print("   - Slow network connection")
print("   - Backend still downloading in background")
print("   - Try running the test again in a few minutes")
sys.exit(1)
