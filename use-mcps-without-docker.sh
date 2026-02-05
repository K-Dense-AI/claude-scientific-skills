#!/bin/bash

# Use only MCPs that don't require Docker
# Keeps: brave-search, claude-skills-scientific
# Removes: sequential-thinking, perplexity-ask

CONFIG_FILE="$HOME/.config/Claude/claude_desktop_config.json"

echo "🔧 Configurando MCPs sem Docker"
echo "================================"
echo ""

# Backup
BACKUP_FILE="${CONFIG_FILE}.backup-$(date +%Y%m%d-%H%M%S)"
if [ -f "$CONFIG_FILE" ]; then
    cp "$CONFIG_FILE" "$BACKUP_FILE"
    echo "💾 Backup criado: $BACKUP_FILE"
fi

# Create config with your Brave API key
cat > "$CONFIG_FILE" << 'EOF'
{
  "mcpServers": {
    "brave-search": {
      "command": "npx",
      "args": [
        "-y",
        "@modelcontextprotocol/server-brave-search"
      ],
      "env": {
        "BRAVE_API_KEY": "BSAIXjL1gRG1xvMGtHiy3i6OnoBKhkr"
      }
    },
    "claude-skills-scientific": {
      "command": "uvx",
      "args": [
        "claude-skills-mcp",
        "--config",
        "/home/user/claude-scientific-skills/mcp-config.json"
      ]
    }
  }
}
EOF

echo ""
echo "✅ Configuração atualizada!"
echo ""
echo "📋 MCPs ativos (sem Docker):"
echo "   ✓ brave-search (com sua API key)"
echo "   ✓ claude-skills-scientific"
echo ""
echo "📍 Local: $CONFIG_FILE"
echo ""
echo "📋 Próximo passo:"
echo "   REINICIE o Claude Desktop"
echo ""
