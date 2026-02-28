# Configuração MCP - Resumo Final

Este repositório contém a configuração completa para integração MCP com Claude Desktop, incluindo **83+ scientific skills**.

## ✅ Configuração Atual (5 MCPs)

A configuração final em `~/.config/Claude/claude_desktop_config.json` inclui:

### 1. **sequential-thinking** (Docker)
- Raciocínio estruturado passo a passo
- Comando: `docker run --rm -i mcp/sequentialthinking`

### 2. **perplexity-ask** (Docker + API Key)
- Busca avançada com Perplexity AI
- Comando: `docker run -i --rm -e PERPLEXITY_API_KEY mcp/perplexity-ask`
- Requer: API key do Perplexity

### 3. **brave-search** (npx + API Key)
- Busca na web com Brave Search
- Comando: `npx -y @modelcontextprotocol/server-brave-search`
- Requer: API key do Brave

### 4. **mcp-reasoner** (Node.js)
- Análise e raciocínio lógico
- Comando: `node /Users/renatopanelli/mcp-reasoner/dist/index.js`
- Repositório: https://github.com/Jacck/mcp-reasoner

### 5. **claude-skills-scientific** (uvx)
- **83+ skills científicos** incluindo:
  - 25 databases (PubMed, ChEMBL, UniProt, AlphaFold DB...)
  - 50 packages (BioPython, RDKit, Scanpy, PyTorch...)
  - 6 integrações (Benchling, DNAnexus, Opentrons...)
- Comando: `uvx claude-skills-mcp --config mcp-config.json`

## 🚀 Como Usar

### Primeira Vez

1. **Reinicie o Claude Desktop**
   - Feche completamente (Quit/Sair)
   - Aguarde 5 segundos
   - Abra novamente

2. **Aguarde o Download Inicial** (2-3 minutos)
   - Docker images: `mcp/sequentialthinking`, `mcp/perplexity-ask`
   - Backend MCP: ~250MB (claude-skills-scientific)

3. **Teste os MCPs**
   ```
   "Liste todos os skills científicos disponíveis"
   "Use sequential thinking para resolver este problema"
   "Use Perplexity para buscar sobre CRISPR"
   "Use Brave Search para pesquisar AlphaFold"
   "Use mcp-reasoner para analisar este argumento"
   ```

## 📁 Arquivos Importantes

### Configuração
- `.claude-desktop-mcp-complete-with-reasoner.json` - Template completo (sem credenciais)
- `mcp-config.json` - Configuração do MCP científico

### Scripts (com credenciais, gitignored)
- `configure-all-5-mcps.sh` - Aplica configuração completa
- `restore-mcp-config.sh` - Restaura configuração com credenciais
- `add-mcp-reasoner.sh` - Adiciona apenas mcp-reasoner

### Testes
- `test-all-mcps.sh` - Testa todos os MCPs antes de reiniciar
- `test-mcp-protocol.py` - Testa protocolo MCP
- `test-mcp-tools.py` - Testa ferramentas MCP

### Documentação
- `MCP-INSTALLATION-GUIDE.md` - Guia completo de instalação
- `MCP-TECHNICAL-ASSESSMENT.md` - Avaliação técnica detalhada

## 🔒 Segurança

⚠️ **Importante:**
- API keys estão **apenas** em `~/.config/Claude/claude_desktop_config.json` (local)
- Scripts com credenciais estão no `.gitignore`
- Templates no repositório **não contêm** credenciais reais

## 📊 Status de Testes

Última execução: 2026-01-25

| Componente | Status |
|------------|--------|
| Configuração JSON | ✅ Válida |
| Total de MCPs | ✅ 5 detectados |
| Docker (host) | ✅ Instalado |
| npx | ✅ v10.9.4 |
| uvx | ✅ v0.8.17 |
| Node.js | ✅ Disponível |
| API Keys | ✅ Configuradas |
| mcp-reasoner path | ✅ Validado |

## 🎓 Exemplos de Uso

### Drug Discovery
```
"Encontre inibidores de EGFR no ChEMBL com IC50 < 50nM,
analise com RDKit suas propriedades ADMET,
e faça docking com DiffDock contra estrutura do AlphaFold"
```

### Genomics Analysis
```
"Carregue este dataset 10X, faça análise single-cell com Scanpy,
identifique populações celulares, e compare com dados do
Cellxgene Census"
```

### Research + Reasoning
```
"Use Perplexity para buscar os últimos papers sobre CRISPR,
depois use sequential thinking para analisar as metodologias,
e finalmente use mcp-reasoner para avaliar a consistência
das conclusões"
```

## 📖 Recursos

- **Repositório Principal:** https://github.com/K-Dense-AI/claude-scientific-skills
- **MCP Server:** https://github.com/K-Dense-AI/claude-skills-mcp
- **mcp-reasoner:** https://github.com/Jacck/mcp-reasoner
- **K-Dense Enterprise:** https://k-dense.ai/

## ✅ Checklist de Instalação

- [x] Docker instalado e rodando
- [x] Node.js/npm instalado (para npx)
- [x] uvx instalado
- [x] mcp-reasoner clonado e buildado
- [x] API keys do Perplexity configuradas
- [x] API keys do Brave configuradas
- [x] Configuração aplicada em `~/.config/Claude/claude_desktop_config.json`
- [ ] Claude Desktop reiniciado
- [ ] Todos os 5 MCPs testados

---

**Última Atualização:** 2026-01-25
**Branch:** `claude/add-scientific-skills-plugin-011CUg6mgwVqqKYPV1pQSSTx`
**Status:** ✅ Pronto para uso
