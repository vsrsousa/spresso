# Pseudopotentials Dict Fix - Replicando Lógica de PseudoConfig

## Problema
Quando usuário passa `pseudopotentials` como dict:
```python
wf = ConvergenceWorkflow.optimize_parameters(
    atoms=atoms,
    pseudopotentials={"Gd": "Gd.pbe-spfn-rrkjus_psl.1.0.0.UPF"},  # ← Dict
    machine="snake5",
    ...
)
```

O fluxo não era consistente com o `pseudopotentials_config`. Resultava em erro: 
```
RuntimeError: First batch of calculations FAILED (ALL 0 jobs failed)
```

## Causa
O dict era processado diferentemente de um config:

### Antes (INCORRETO):
1. `ConvergenceWorkflow.__init__`:
   - Chamava `discover_pseudopotential_directory` → obtinha caminhos **absolutos**
   - Armazenava caminhos absolutos em `self.pseudopotentials`
   - Não passava `pseudopotentials_base_path` para `CalculationWorkflow`

2. `CalculationWorkflow.__init__`:
   - Recebia `pseudopotentials` com caminhos absolutos
   - Não recebia `pseudopotentials_base_path`
   - Não setava `os.environ['ESPRESSO_PSEUDO']` apropriadamente

3. `Espresso` + `_transfer_pseudopotentials`:
   - Buscava pelos arquivos usando `search_dirs` (env var, home, etc)
   - Não encontrava porque paths eram absolutos ou estavam em formato errado

### Depois (CORRETO):
Agora replica **exatamente** a lógica do `pseudopotentials_config`:

1. **`ConvergenceWorkflow.__init__`** (quando dict):
   ```python
   resolved_pseudos, self.pseudopotentials_base_path = discover_pseudopotential_directory(pseudopotentials)
   
   # Extract FILENAMES from resolved absolute paths (same as pseudoconfig logic)
   self.pseudopotentials = {}
   for element, full_path in resolved_pseudos.items():
       filename = os.path.basename(full_path)
       self.pseudopotentials[element] = filename  # ← Armazena FILENAMES, não paths
   ```

2. **Pass `pseudopotentials_base_path` to `CalculationWorkflow`** (Fase 1 e 2):
   ```python
   wf_kwargs = {
       'atoms': self.atoms,
       'pseudopotentials': self.pseudopotentials,  # ← FILENAMES apenas
       'pseudopotentials_base_path': self.pseudopotentials_base_path,  # ← NOVO
       ...
   }
   wf1 = CalculationWorkflow(**wf_kwargs)
   ```

3. **`CalculationWorkflow.__init__`** (quando dict com base_path):
   ```python
   # If pseudopotentials_base_path is provided (from ConvergenceWorkflow), set env var
   if pseudopotentials_base_path:
       os.environ['ESPRESSO_PSEUDO'] = pseudopotentials_base_path
       logger.info(f"Set ESPRESSO_PSEUDO={pseudopotentials_base_path}")
   ```

4. **`CalculationWorkflow.submit_scf_batch_multiple`**:
   ```python
   temp_workflow = CalculationWorkflow(
       self.atoms,
       pseudopotentials=self.pseudopotentials,
       pseudopotentials_base_path=getattr(self, 'pseudopotentials_base_path', None),  # ← Pass along
       ...
   )
   ```

## Fluxo Completo Agora

```
User Input: pseudopotentials={"Gd": ".../Gd.pbe-spfn.UPF"}
    ↓
[ConvergenceWorkflow.__init__]
    ↓
discover_pseudopotential_directory(pseudos)
    ↓
resolved_pseudos = {
    "Gd": "/full/path/to/Gd.pbe-spfn.UPF"
}
base_path = "/full/path/to"
    ↓
Extract filenames:
self.pseudopotentials = {"Gd": "Gd.pbe-spfn.UPF"}
self.pseudopotentials_base_path = "/full/path/to"
    ↓
Set ESPRESSO_PSEUDO env var (Phase 1):
os.environ['ESPRESSO_PSEUDO'] = "/full/path/to"
    ↓
Create CalculationWorkflow(
    pseudopotentials={"Gd": "Gd.pbe-spfn.UPF"},
    pseudopotentials_base_path="/full/path/to"
)
    ↓
[CalculationWorkflow.__init__]
    ↓
Set ESPRESSO_PSEUDO again (garantido):
os.environ['ESPRESSO_PSEUDO'] = "/full/path/to"
    ↓
Espresso receives:
parameters["pseudopotentials"] = {"Gd": "Gd.pbe-spfn.UPF"}
ESPRESSO_PSEUDO env var = "/full/path/to"
    ↓
[Remote Execution]
    ↓
_transfer_pseudopotentials():
    pseudopotentials = {"Gd": "Gd.pbe-spfn.UPF"}
    search_dirs = [ESPRESSO_PSEUDO, env vars, home dirs]
    → Encontra: "/full/path/to/Gd.pbe-spfn.UPF"
    → Transfere pro remote
```

## Mudanças de Código

### 1. `xespresso/utils/pseudo_utils.py`
- **Fix**: Quando caminho é absoluto, agora extrai o diretório pai como `base_dir`
- Antes: `base_dir` era `None` para caminhos absolutos
- Depois: `base_dir = os.path.dirname(path)` quando é absoluto

### 2. `xespresso/workflow/convergence_workflow.py`
- **Change**: Quando usar dict, extrai FILENAMES (não caminhos absolutos) igual ao pseudoconfig
- **Change**: Passa `pseudopotentials_base_path` ao criar `CalculationWorkflow` (Fase 1 e 2)

### 3. `xespresso/workflow/calculation_workflow.py`
- **Change**: Recebe `pseudopotentials_base_path` como parâmetro
- **Change**: Seta `os.environ['ESPRESSO_PSEUDO']` quando base_path é fornecido
- **Change**: Passa `pseudopotentials_base_path` ao criar temp workflows em `submit_scf_batch_multiple`

## Benefícios

✅ **Consistency**: Dict e config agora seguem exatamente a mesma lógica  
✅ **Simplicity**: Não há duplicação de chamadas a `discover_pseudopotential_directory`  
✅ **Reliability**: ESPRESSO_PSEUDO é setado em todos os pontos necessários  
✅ **Transparency**: O mesmo caminho (`base_path`) é propagado por todo o workflow  

## Testando

```python
# Before: ❌ Erro "First batch FAILED (ALL 0 jobs failed)"
wf = ConvergenceWorkflow.optimize_parameters(
    atoms=atoms,
    pseudopotentials={"Gd": "Gd.pbe-spfn-rrkjus_psl.1.0.0.UPF"},
    precision='low',
    machine="snake5",
    code_version="7.4",
)

# After: ✅ Funciona!
```
