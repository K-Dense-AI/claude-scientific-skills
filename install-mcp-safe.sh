#!/bin/bash

# Script seguro para adicionar MCP ao Claude Desktop
# Preserva configurações existentes e faz backup

set -e  # Para em caso de erro

CONFIG_DIR="$HOME/.config/Claude"
CONFIG_FILE="$CONFIG_DIR/claude_desktop_config.json"
NEW_MCP_FILE=".claude-desktop-mcp.json"

echo "🔧 Instalação Segura de MCP no Claude Desktop"
echo "=============================================="
echo ""

# 1. Criar diretório se não existir
if [ ! -d "$CONFIG_DIR" ]; then
    echo "📁 Criando diretório de configuração..."
    mkdir -p "$CONFIG_DIR"
fi

# 2. Fazer backup se arquivo existir
if [ -f "$CONFIG_FILE" ]; then
    BACKUP_FILE="${CONFIG_FILE}.backup-$(date +%Y%m%d-%H%M%S)"
    echo "💾 Fazendo backup da configuração existente..."
    cp "$CONFIG_FILE" "$BACKUP_FILE"
    echo "   ✅ Backup salvo em: $BACKUP_FILE"
    echo ""

    # 3. Verificar se há outros MCPs configurados
    echo "🔍 Verificando MCPs existentes..."
    EXISTING_MCPS=$(python3 -c "
import json
try:
    with open('$CONFIG_FILE', 'r') as f:
        config = json.load(f)
    servers = config.get('mcpServers', {})
    if servers:
        print('Encontrados:')
        for name in servers.keys():
            print(f'  - {name}')
    else:
        print('Nenhum MCP encontrado')
except:
    print('Nenhum MCP encontrado')
" 2>/dev/null || echo "Nenhum MCP encontrado")

    echo "$EXISTING_MCPS"
    echo ""

    # 4. Mesclar configurações
    echo "🔀 Mesclando configurações..."
    python3 << 'PYTHON_SCRIPT'
import json
import sys

# Ler configuração existente
try:
    with open('$CONFIG_FILE', 'r') as f:
        existing_config = json.load(f)
except:
    existing_config = {"mcpServers": {}}

# Ler nova configuração
try:
    with open('$NEW_MCP_FILE', 'r') as f:
        new_config = json.load(f)
except Exception as e:
    print(f"❌ Erro ao ler nova configuração: {e}")
    sys.exit(1)

# Mesclar (novo MCP não sobrescreve existentes)
if 'mcpServers' not in existing_config:
    existing_config['mcpServers'] = {}

# Adicionar ou atualizar apenas o MCP científico
new_servers = new_config.get('mcpServers', {})
for server_name, server_config in new_servers.items():
    if server_name in existing_config['mcpServers']:
        print(f"⚠️  MCP '{server_name}' já existe, será atualizado")
    existing_config['mcpServers'][server_name] = server_config

# Salvar configuração mesclada
with open('$CONFIG_FILE', 'w') as f:
    json.dump(existing_config, f, indent=2)

print("✅ Configurações mescladas com sucesso!")
PYTHON_SCRIPT

else
    echo "📝 Nenhuma configuração existente encontrada"
    echo "   Criando nova configuração..."
    cp "$NEW_MCP_FILE" "$CONFIG_FILE"
    echo "   ✅ Configuração criada"
fi

echo ""
echo "=============================================="
echo "✅ INSTALAÇÃO CONCLUÍDA COM SUCESSO!"
echo "=============================================="
echo ""
echo "📋 Configuração final:"
python3 -c "
import json
with open('$CONFIG_FILE', 'r') as f:
    config = json.load(f)
servers = config.get('mcpServers', {})
print(f'Total de MCPs configurados: {len(servers)}')
for name in servers.keys():
    print(f'  ✓ {name}')
"
echo ""
echo "📍 Local: $CONFIG_FILE"
echo ""
echo "📋 Próximos passos:"
echo "   1. REINICIE o Claude Desktop (feche e abra novamente)"
echo "   2. Na primeira vez, aguarde 1-2 minutos (download do backend)"
echo "   3. Teste com: 'Liste todos os skills científicos disponíveis'"
echo ""

if [ -f "${CONFIG_FILE}.backup-"* ]; then
    echo "💡 Dica: Se algo der errado, você pode restaurar o backup:"
    echo "   cp ${CONFIG_FILE}.backup-* $CONFIG_FILE"
    echo ""
fi
