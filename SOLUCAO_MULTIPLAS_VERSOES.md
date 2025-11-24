# Solução: Suporte a Múltiplas Versões do Quantum ESPRESSO

## Problema Original (Traduzido)

A questão levantada foi:

> "A mesma máquina pode ter instalada mais de uma versão do quantum espresso, como o codes e o machines lidam com isso? E como o user vai pode escolher nesse caso? Outra coisa, o user pode escolher diferentes versões para diferentes tasks no mesmo script, e se ele escolher rodar na mesma máquina, não deveria haver uma quebra de conexão na hora de submeter os jobs, isso está contemplado?"

**Tradução:**
> "The same machine can have more than one version of Quantum ESPRESSO installed. How do codes and machines handle this? And how can the user choose in this case? Another thing, the user can choose different versions for different tasks in the same script, and if they choose to run on the same machine, shouldn't there be a connection break when submitting jobs? Is this covered?"

## Problema Identificado

1. **Múltiplas versões na mesma máquina**: Não havia suporte para gerenciar múltiplas instalações do QE
2. **Seleção de versão por tarefa**: Usuário não podia escolher qual versão usar para cada cálculo
3. **Persistência de conexão**: Preocupação sobre quebra de conexão SSH ao alternar entre versões

## Solução Implementada

### 1. Estrutura de Dados Estendida

A classe `CodesConfig` foi estendida para suportar múltiplas versões:

```python
@dataclass
class CodesConfig:
    machine_name: str
    codes: Dict[str, Code] = field(default_factory=dict)  # Códigos principais
    qe_version: Optional[str] = None  # Versão padrão
    versions: Optional[Dict[str, Dict]] = None  # NOVO: Versões específicas
```

**Formato multi-versão:**
```json
{
  "machine_name": "cluster1",
  "qe_version": "7.2",
  "versions": {
    "7.2": {
      "qe_prefix": "/opt/qe-7.2/bin",
      "modules": ["quantum-espresso/7.2"],
      "codes": {"pw": {...}, "hp": {...}}
    },
    "6.8": {
      "qe_prefix": "/opt/qe-6.8/bin",
      "modules": ["quantum-espresso/6.8"],
      "codes": {"pw": {...}, "hp": {...}}
    }
  }
}
```

### 2. API Estendida

#### Adicionar códigos a versões específicas:
```python
config = CodesConfig(machine_name="cluster1", qe_version="7.2")

# Adicionar QE 7.2
config.add_code(
    Code(name="pw", path="/opt/qe-7.2/bin/pw.x"),
    version="7.2"
)

# Adicionar QE 6.8
config.add_code(
    Code(name="pw", path="/opt/qe-6.8/bin/pw.x"),
    version="6.8"
)
```

#### Obter códigos de versões específicas:
```python
codes = load_codes_config("cluster1")

# Usar QE 7.2
pw_72 = codes.get_code("pw", version="7.2")

# Usar QE 6.8
pw_68 = codes.get_code("pw", version="6.8")
```

#### Listar versões disponíveis:
```python
versions = codes.list_versions()  # ['7.2', '6.8']
codes_72 = codes.list_codes(version="7.2")  # ['pw', 'hp', 'dos']
```

#### Adicionar versões programaticamente:
```python
from xespresso.codes import add_version_to_config

config = add_version_to_config(
    machine_name="cluster1",
    version="7.3",
    qe_prefix="/opt/qe-7.3/bin",
    modules=["quantum-espresso/7.3"]
)
```

### 3. Persistência de Conexão (Já Implementado!)

**Descoberta importante**: A persistência de conexão JÁ estava implementada no `RemoteExecutionMixin`!

O mixin já cacheia conexões por `(host, username)`, então:
- Alternar entre versões do QE **não quebra a conexão**
- Apenas o caminho do executável e os módulos mudam
- A conexão SSH é reutilizada automaticamente

**Como funciona:**
```python
# Primeira calculação com QE 7.2
pw_72 = codes.get_code("pw", version="7.2")
# ... executa cálculo ... (estabelece conexão SSH)

# Segunda calculação com QE 6.8
pw_68 = codes.get_code("pw", version="6.8")
# ... executa cálculo ... (REUTILIZA a conexão!)

# Terceira calculação com QE 7.2 novamente
pw_72_again = codes.get_code("pw", version="7.2")
# ... executa cálculo ... (AINDA usando mesma conexão!)
```

### 4. Compatibilidade Retroativa

Configurações antigas continuam funcionando:

```python
# Configuração antiga (ainda funciona!)
config = CodesConfig(machine_name="cluster1", qe_version="7.2")
config.add_code(Code(name="pw", path="/usr/bin/pw.x"))

# Acessar sem parâmetro de versão
pw = config.get_code("pw")  # Funciona como antes
```

## Testes Implementados

21 testes passando, incluindo 6 novos testes multi-versão:

1. `test_add_code_with_version` - Adicionar códigos a versões específicas
2. `test_list_versions` - Listar versões disponíveis
3. `test_list_codes_by_version` - Listar códigos por versão
4. `test_version_config_serialization` - Serialização de configurações multi-versão
5. `test_backward_compatibility` - Compatibilidade com configurações antigas
6. `test_get_code_fallback` - Fallback para códigos principais

```bash
$ python -m pytest tests/test_codes.py -v
...
============================== 21 passed in 0.45s ==============================
```

## Exemplos Criados

### 1. `examples/multi_version_example.py`
Demonstra setup básico de múltiplas versões:
- Criar configuração com QE 7.2
- Adicionar QE 6.8 à mesma máquina
- Carregar e inspecionar configuração
- Explicar persistência de conexão

### 2. `examples/workflow_multi_version.py`
Exemplo prático de workflow:
- SCF com QE 7.2 (recursos novos)
- DOS com QE 6.8 (estabilidade)
- Bands com QE 6.8 (consistência)
- Volta ao QE 7.2 (outro cálculo)
- Conexão mantida durante todo o processo

## Documentação Criada

### 1. `docs/MULTIPLE_VERSIONS.md`
- Overview completo do recurso
- Problemas resolvidos
- Exemplos de uso
- Casos de uso práticos
- Benefícios e detalhes de implementação

### 2. `docs/CODES_CONFIGURATION.md` (Atualizado)
- Formato de arquivo multi-versão
- API estendida com parâmetros de versão
- Exemplos de uso avançado
- Integração com machines

### 3. `NEW_FEATURES.md` (Atualizado)
- Seção sobre suporte multi-versão
- Quick start com múltiplas versões
- Links para documentação completa

## Resposta às Questões

### ❓ "Como o codes e o machines lidam com múltiplas versões?"

✅ **Resposta**: O módulo `codes` agora tem:
- Campo `versions` na configuração
- Métodos que aceitam parâmetro `version`
- Função `add_version_to_config()` para fácil gerenciamento

### ❓ "Como o user vai poder escolher?"

✅ **Resposta**: API simples e clara:
```python
codes = load_codes_config("cluster1")
pw_72 = codes.get_code("pw", version="7.2")  # Escolha explícita
pw_68 = codes.get_code("pw", version="6.8")  # Escolha explícita
```

### ❓ "Não deveria haver quebra de conexão?"

✅ **Resposta**: Correto! E a boa notícia é que **já estava implementado**:
- `RemoteExecutionMixin` cacheia conexões por `(host, username)`
- Alternar versões NÃO quebra a conexão
- Apenas caminho do executável e módulos mudam
- Performance mantida, overhead zero

## Benefícios da Solução

1. ✅ **Flexibilidade**: Use versões diferentes para tarefas diferentes
2. ✅ **Performance**: Conexão SSH persiste, sem overhead
3. ✅ **Confiabilidade**: Sem quebra de conexão ao alternar versões
4. ✅ **Organização**: Versionamento claro no código
5. ✅ **Compatibilidade**: Configs antigas continuam funcionando
6. ✅ **Bem testado**: 21 testes passando
7. ✅ **Bem documentado**: 3 documentos + 2 exemplos

## Casos de Uso Reais

### 1. Requisitos de Features
```python
# Usar QE 7.2 para novo formato de parâmetros Hubbard
pw_72 = codes.get_code("pw", version="7.2")
hp_72 = codes.get_code("hp", version="7.2")

# Usar QE 6.8 para compatibilidade legada
pw_68 = codes.get_code("pw", version="6.8")
```

### 2. Testes e Validação
```python
# Comparar resultados entre versões
for version in ["7.2", "6.8"]:
    pw = codes.get_code("pw", version=version)
    # ... executar mesmo cálculo com versões diferentes ...
```

### 3. Workflows Mistos
```python
# SCF com versão mais nova
pw_72 = codes.get_code("pw", version="7.2")
# ... executar SCF ...

# Pós-processamento com versão estável
dos_68 = codes.get_code("dos", version="6.8")
bands_68 = codes.get_code("bands", version="6.8")
# ... executar DOS e bands ...
```

## Resumo Técnico

### Arquivos Modificados
1. `xespresso/codes/config.py` - Classe `CodesConfig` estendida
2. `xespresso/codes/manager.py` - Funções `load_codes_config()` e `add_version_to_config()`
3. `xespresso/codes/__init__.py` - Export de nova função

### Arquivos Criados
1. `tests/test_codes.py` - 6 novos testes multi-versão
2. `examples/multi_version_example.py` - Exemplo de setup
3. `examples/workflow_multi_version.py` - Exemplo prático
4. `docs/MULTIPLE_VERSIONS.md` - Documentação completa

### Arquivos Atualizados
1. `docs/CODES_CONFIGURATION.md` - API estendida
2. `NEW_FEATURES.md` - Seção multi-versão

### Conexão Remota
**Nenhuma modificação necessária!** O `RemoteExecutionMixin` já implementava persistência corretamente.

## Conclusão

A solução implementada:

1. ✅ **Resolve** o problema de múltiplas versões na mesma máquina
2. ✅ **Fornece** API clara para seleção de versão
3. ✅ **Confirma** que conexão persiste (já estava implementado)
4. ✅ **Mantém** compatibilidade com código existente
5. ✅ **Adiciona** testes abrangentes (21 testes)
6. ✅ **Documenta** completamente com exemplos práticos

A implementação é **minimal**, **eficiente** e **bem testada**, resolvendo exatamente o problema descrito sem introduzir complexidade desnecessária.
