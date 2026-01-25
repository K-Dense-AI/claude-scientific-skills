# Claude Scientific Skills - Guia de Instalação MCP

## 📋 Visão Geral

Este guia mostra como integrar os **83+ skills científicos** em qualquer cliente compatível com MCP (Model Context Protocol), incluindo:
- 🖱️ **Cursor** - IDE com AI integrado
- 💬 **Claude Desktop** - Aplicativo desktop oficial
- 🔧 **ChatGPT** - Via extensões MCP
- 🚀 **Outros clientes MCP** - Google ADK, OpenAI Agent SDK, etc.

## ✅ Pré-requisitos

- Python 3.11+ (✅ instalado: Python 3.11.14)
- uvx 0.8.17+ (✅ instalado: uvx 0.8.17)
- Cliente MCP (Cursor, Claude Desktop, etc.)

## 🚀 Instalação Rápida

### Opção 1: Configuração Automática para Cursor

1. **Copie a configuração para o Cursor:**
   ```bash
   cp .cursor-mcp.json ~/.cursor/mcp.json
   ```

2. **Reinicie o Cursor**

3. **Verifique a instalação:**
   - Abra o Cursor
   - O MCP server baixará o backend (~250 MB) automaticamente na primeira execução
   - Após 5 segundos, os 83+ skills científicos estarão disponíveis

### Opção 2: Configuração Manual para Claude Desktop

1. **Localize o arquivo de configuração:**
   - **macOS**: `~/Library/Application Support/Claude/claude_desktop_config.json`
   - **Windows**: `%APPDATA%\Claude\claude_desktop_config.json`
   - **Linux**: `~/.config/Claude/claude_desktop_config.json`

2. **Adicione a configuração:**
   ```bash
   # Para usuários macOS/Linux:
   mkdir -p ~/.config/Claude/
   cp .claude-desktop-mcp.json ~/.config/Claude/claude_desktop_config.json

   # Ou edite manualmente e adicione:
   ```

   ```json
   {
     "mcpServers": {
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
   ```

3. **Reinicie o Claude Desktop**

### Opção 3: Execução Standalone (para testes)

Execute diretamente via linha de comando:

```bash
# Com configuração padrão
uvx claude-skills-mcp

# Com configuração customizada
uvx claude-skills-mcp --config /home/user/claude-scientific-skills/mcp-config.json

# Com logging verbose
uvx claude-skills-mcp --config /home/user/claude-scientific-skills/mcp-config.json --verbose
```

## 🔧 Arquivos de Configuração

### 1. `mcp-config.json` (Configuração Principal)

```json
{
  "skill_sources": [
    {
      "type": "github",
      "owner": "anthropics",
      "repo": "anthropic-skills",
      "description": "Official Anthropic Skills"
    },
    {
      "type": "github",
      "owner": "K-Dense-AI",
      "repo": "claude-scientific-skills",
      "description": "K-Dense Scientific Skills Collection"
    },
    {
      "type": "local",
      "path": "/home/user/claude-scientific-skills",
      "description": "Local Scientific Skills Repository"
    }
  ],
  "embedding": {
    "model": "text-embedding-3-small",
    "dimensions": 1536
  },
  "content": {
    "max_file_size_kb": 500,
    "allowed_extensions": [".md", ".txt", ".py", ".json", ".yaml", ".yml"]
  },
  "server": {
    "host": "127.0.0.1",
    "port": 8765
  }
}
```

### 2. `.cursor-mcp.json` (Para Cursor)

Configuração específica para o Cursor IDE.

### 3. `.claude-desktop-mcp.json` (Para Claude Desktop)

Configuração específica para o Claude Desktop.

## 🎯 Como Usar

Após a instalação, você terá acesso a 3 ferramentas MCP:

### 1. **find_helpful_skills**
Busca semântica por skills relevantes:

```
Exemplo: "Preciso analisar dados de single-cell RNA-seq"
→ Retorna: Scanpy, AnnData, scvi-tools, PyDESeq2, etc.
```

### 2. **read_skill_document**
Lê documentação específica de um skill:

```
Exemplo: "Mostre-me como usar o RDKit para calcular propriedades moleculares"
→ Retorna: Documentação completa do RDKit SKILL
```

### 3. **list_skills**
Lista todos os skills disponíveis:

```
Retorna: Lista completa dos 83+ skills científicos disponíveis
```

## 📚 Skills Disponíveis

### Databases (25)
- PubMed, ChEMBL, UniProt, AlphaFold DB, PubChem, COSMIC, ClinVar, etc.

### Packages (50)
- BioPython, RDKit, Scanpy, PyTorch, DeepChem, DiffDock, Matplotlib, etc.

### Integrations (6)
- Benchling, DNAnexus, Opentrons, LabArchives, LatchBio, OMERO

### Methodologies
- Exploratory Data Analysis, Scientific Writing, Peer Review, etc.

## 🔍 Exemplos de Uso

### Drug Discovery
```
"Encontre inibidores de EGFR no ChEMBL com IC50 < 50nM,
analise suas relações estrutura-atividade com RDKit,
e faça docking virtual com DiffDock"
```

### Genomics Analysis
```
"Carregue este dataset 10X, faça análise single-cell com Scanpy,
identifique populações celulares, e compare com dados do
Cellxgene Census"
```

### Clinical Research
```
"Analise este VCF, anote todas as variantes usando Ensembl,
verifique significância clínica no ClinVar, e gere um relatório"
```

## 🐛 Troubleshooting

### Problema: Backend não baixa

**Solução:**
```bash
# Force o download do backend
uvx claude-skills-mcp --verbose
```

### Problema: Skills não aparecem

**Solução:**
1. Verifique se o MCP server está rodando:
   ```bash
   ps aux | grep claude-skills-mcp
   ```

2. Verifique logs do cliente (Cursor/Claude Desktop)

3. Reinicie o cliente

### Problema: Configuração não funciona

**Solução:**
```bash
# Valide o JSON de configuração
cat mcp-config.json | python -m json.tool

# Verifique se o caminho está correto
ls -la /home/user/claude-scientific-skills/mcp-config.json
```

### Problema: Permissões negadas

**Solução:**
```bash
# Dê permissão de execução
chmod +x ~/.local/bin/uvx

# Ou reinstale uvx
pip install --user --upgrade uv
```

## 🔄 Atualização

Para atualizar o claude-skills-mcp:

```bash
# Atualizar o pacote
uvx --reinstall claude-skills-mcp

# Atualizar skills do repositório local
cd /home/user/claude-scientific-skills
git pull origin main
```

## 📖 Recursos Adicionais

- **Repositório**: https://github.com/K-Dense-AI/claude-scientific-skills
- **MCP Server**: https://github.com/K-Dense-AI/claude-skills-mcp
- **Documentação MCP**: https://modelcontextprotocol.io/
- **K-Dense Enterprise**: https://k-dense.ai/

## 🎓 Próximos Passos

1. ✅ Instale seguindo as instruções acima
2. 🧪 Teste com exemplos simples
3. 📚 Explore os skills disponíveis com `list_skills`
4. 🚀 Use `find_helpful_skills` para encontrar o skill certo para sua tarefa
5. 🔬 Comece a usar para suas pesquisas científicas!

## 💡 Dicas

- Use busca semântica para encontrar skills relevantes rapidamente
- Os skills científicos já incluem exemplos práticos e best practices
- Configure o logging verbose durante desenvolvimento para debug
- O MCP server baixa automaticamente updates dos repositórios GitHub

---

**Versão**: 1.55.0
**Última Atualização**: 2026-01-25
**Suporte**: https://github.com/K-Dense-AI/claude-scientific-skills/issues
