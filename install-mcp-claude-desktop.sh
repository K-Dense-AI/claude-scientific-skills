#!/bin/bash

echo "🔧 Configurando MCP no Claude Desktop..."
echo ""

# Criar diretório se não existir
mkdir -p ~/.config/Claude/

# Copiar configuração
cp .claude-desktop-mcp.json ~/.config/Claude/claude_desktop_config.json

# Verificar
if [ -f ~/.config/Claude/claude_desktop_config.json ]; then
    echo "✅ Configuração MCP instalada com sucesso!"
    echo ""
    echo "📍 Local: ~/.config/Claude/claude_desktop_config.json"
    echo ""
    echo "📋 Próximo passo:"
    echo "   1. REINICIE o Claude Desktop (feche e abra novamente)"
    echo "   2. Aguarde 1-2 minutos (download do backend na primeira vez)"
    echo "   3. Teste com: 'Liste todos os skills científicos disponíveis'"
    echo ""
    echo "🎉 Após reiniciar, você terá acesso a 83+ skills científicos!"
else
    echo "❌ Erro ao copiar configuração"
fi
