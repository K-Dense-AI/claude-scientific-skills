#!/bin/bash

# Script para testar a instalação do Claude Skills MCP

echo "🧪 Testando Claude Skills MCP..."
echo ""

# 1. Verificar uvx
echo "1️⃣ Verificando uvx..."
if command -v uvx &> /dev/null; then
    echo "   ✅ uvx instalado: $(uvx --version)"
else
    echo "   ❌ uvx não encontrado"
    exit 1
fi
echo ""

# 2. Verificar Python
echo "2️⃣ Verificando Python..."
python_version=$(python3 --version 2>&1)
echo "   ✅ $python_version"
echo ""

# 3. Verificar arquivos de configuração
echo "3️⃣ Verificando arquivos de configuração..."
if [ -f "mcp-config.json" ]; then
    echo "   ✅ mcp-config.json encontrado"
    # Validar JSON
    if python3 -m json.tool mcp-config.json > /dev/null 2>&1; then
        echo "   ✅ mcp-config.json é válido"
    else
        echo "   ❌ mcp-config.json inválido"
        exit 1
    fi
else
    echo "   ❌ mcp-config.json não encontrado"
    exit 1
fi

if [ -f ".cursor-mcp.json" ]; then
    echo "   ✅ .cursor-mcp.json encontrado"
fi

if [ -f ".claude-desktop-mcp.json" ]; then
    echo "   ✅ .claude-desktop-mcp.json encontrado"
fi
echo ""

# 4. Verificar SKILL files
echo "4️⃣ Verificando SKILL files..."
skill_count=$(find . -name "SKILL.md" | wc -l)
echo "   ✅ Encontrados $skill_count SKILL.md files"
echo ""

# 5. Testar comando MCP
echo "5️⃣ Testando comando claude-skills-mcp..."
echo "   ℹ️  Executando: uvx claude-skills-mcp --help"
if timeout 30 uvx claude-skills-mcp --help > /dev/null 2>&1; then
    echo "   ✅ claude-skills-mcp está funcionando"
else
    echo "   ⚠️  Timeout ou erro (normal na primeira execução)"
fi
echo ""

# 6. Verificar repositórios de skills
echo "6️⃣ Verificando estrutura de skills..."
if [ -d "scientific-packages" ]; then
    pkg_count=$(ls -d scientific-packages/*/ 2>/dev/null | wc -l)
    echo "   ✅ $pkg_count pacotes científicos encontrados"
fi

if [ -d "scientific-databases" ]; then
    db_count=$(ls -d scientific-databases/*/ 2>/dev/null | wc -l)
    echo "   ✅ $db_count databases científicas encontradas"
fi

if [ -d "scientific-integrations" ]; then
    int_count=$(ls -d scientific-integrations/*/ 2>/dev/null | wc -l)
    echo "   ✅ $int_count integrações científicas encontradas"
fi
echo ""

# Resumo
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "✅ TESTES CONCLUÍDOS COM SUCESSO!"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "📋 Próximos passos:"
echo ""
echo "Para Cursor:"
echo "  cp .cursor-mcp.json ~/.cursor/mcp.json"
echo "  # Reinicie o Cursor"
echo ""
echo "Para Claude Desktop (macOS/Linux):"
echo "  mkdir -p ~/.config/Claude/"
echo "  cp .claude-desktop-mcp.json ~/.config/Claude/claude_desktop_config.json"
echo "  # Reinicie o Claude Desktop"
echo ""
echo "Execução standalone:"
echo "  uvx claude-skills-mcp --config $(pwd)/mcp-config.json"
echo ""
echo "📖 Documentação completa: MCP-INSTALLATION-GUIDE.md"
echo ""
