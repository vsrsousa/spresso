# Análise Completa do Fluxo de Pseudopotenciais

## PROBLEMA RELATADO
```
RuntimeError: ❌ CRITICAL ERROR: First batch of calculations FAILED
   Iteration 1: ALL 0 jobs failed
```

Isso significa que `len(completion) == 0` - nenhum job foi criado, não é que os jobs falharam.

## FLUXO COM PSEUDOPOTENTIALS_CONFIG (FUNCIONA)

### 1. ConvergenceWorkflow.__init__
```python
if pseudopotentials_config is not None:
    cfg = load_pseudopotentials_config(pseudopotentials_config, verbose=False)
    self.pseudopotentials_base_path = cfg.base_path
    for el, pseudo in cfg.pseudopotentials.items():
        if el in required_elements:
            filename = pseudo.filename  # ← FILENAME ONLY (não full path!)
            self.pseudopotentials[el] = filename
    self.ecutrho_ratio = get_ecutrho_ratio(required_elements, cfg)
```

Armazena:
- `self.pseudopotentials`: `{'Gd': 'Gd.pbe-spfn-rrkjus_psl.1.0.0.UPF'}`  (FILENAME ONLY)
- `self.pseudopotentials_base_path`: `/home/vinicius/scratch/projects/spresso/pseudo`

### 2. ConvergenceWorkflow.run_convergence_independent
```python
if self._pseudo_config_name:
    wf_kwargs['pseudopotentials_config'] = self._pseudo_config_name
else:
    wf_kwargs['pseudopotentials'] = self.pseudopotentials
    # pseudopotentials_base_path already set in __init__ via os.environ
```

Passa:
- `pseudopotentials_config='default'` (string, nome do config)

### 3. CalculationWorkflow.__init__
```python
if pseudopotentials_config is not None:
    pseudopotentials = self._load_pseudopotentials_from_config(pseudopotentials_config)
    # Armazena:
    self._pseudo_config = config
    self.pseudopotentials_base_path = config.base_path
```

Retorna:
- `pseudopotentials`: `{'Gd': 'Gd.pbe-spfn-rrkjus_psl.1.0.0.UPF'}`  (FILENAME ONLY)
- `self.pseudopotentials_base_path`: `/home/vinicius/scratch/projects/spresso/pseudo`

### 4. CalculationWorkflow.submit_scf_batch
```python
if self.pseudopotentials_base_path:
    os.environ['ESPRESSO_PSEUDO'] = self.pseudopotentials_base_path
    
params = {
    'pseudopotentials': self.pseudopotentials,  # {'Gd': 'Gd.pbe...UPF'}
    ...
}
# Cria Espresso calc COM pseudopotentials filenames
calc = Espresso(**params)
```

### 5. Remote Transfer em _transfer_pseudopotentials
```python
pseudopotentials = self.calc.parameters.get("pseudopotentials", {})
# {'Gd': 'Gd.pbe-spfn-rrkjus_psl.1.0.0.UPF'}

search_dirs = []
search_dirs.append(os.environ["ESPRESSO_PSEUDO"])  # /home/vinicius/.../pseudo
search_dirs.append(...)

for symbol, pseudo_file in pseudopotentials.items():  # pseudo_file = 'Gd.pbe...UPF'
    for pseudo_dir in search_dirs:
        local_path = os.path.join(pseudo_dir, pseudo_file)
        # /home/vinicius/.../pseudo/Gd.pbe...UPF ✓ FOUND!
        if os.path.exists(local_path):
            # Transfer
```

✅ **FUNCIONA**: Filenames + base_path via ESPRESSO_PSEUDO → encontra os arquivos

---

## FLUXO COM DICT PSEUDOPOTENCIAIS (QUEBRA)

### 1. ConvergenceWorkflow.__init__
```python
else:  # pseudopotentials_config is None, recebeu dict
    self.pseudopotentials, self.pseudopotentials_base_path = discover_pseudopotential_directory(pseudopotentials)
    # Armazena:
    # self.pseudopotentials: {'Gd': '/home/vinicius/.../Gd.pbe-spfn-rrkjus_psl.1.0.0.UPF'}  (FULL PATH!)
    # self.pseudopotentials_base_path: '/home/vinicius/scratch/projects/spresso/pseudo'
    
    self.ecutrho_ratio = 4.0  # Default
```

### 2. ConvergenceWorkflow.run_convergence_independent
```python
if self._pseudo_config_name:
    wf_kwargs['pseudopotentials_config'] = self._pseudo_config_name
else:
    wf_kwargs['pseudopotentials'] = self.pseudopotentials
    # pseudopotentials_base_path already set in __init__ via os.environ
```

Passa:
- `pseudopotentials`: `{'Gd': '/full/path/to/Gd.pbe...UPF'}`  (FULL PATHS)
- NÃO passa `pseudopotentials_base_path` como parâmetro

### 3. CalculationWorkflow.__init__
```python
elif pseudopotentials is None:
    raise ValueError(...)

self.original_pseudopotentials = pseudopotentials
# Armazena full paths direto, sem processar!
# self.pseudopotentials_base_path = None inicialmente
```

NÃO carrega de config, apenas armazena o dict como está.

### 4. CalculationWorkflow.submit_scf_batch_multiple
```python
temp_workflow = CalculationWorkflow(
    self.atoms,
    protocol=self.protocol,
    pseudopotentials=self.pseudopotentials,  # Full paths
    kspacing=params.get('kspacing', self.preset.get('kspacing')),
    input_data=input_data_override,
    queue=self.queue,
    **self.extra_kwargs
    # ← NÃO PASSA pseudopotentials_base_path como parâmetro!
)

# Depois tenta copiar:
if hasattr(self, 'pseudopotentials_base_path'):
    temp_workflow.pseudopotentials_base_path = self.pseudopotentials_base_path
```

⚠️ **PROBLEMA 1**: Não passa `pseudopotentials_base_path` como parâmetro ao `__init__`

### 5. CalculationWorkflow.submit_scf_batch
```python
if self.pseudopotentials_base_path:
    os.environ['ESPRESSO_PSEUDO'] = self.pseudopotentials_base_path
    # Este bloco PODE não executar se não foi copiado corretamente

params = {
    'pseudopotentials': self.pseudopotentials,  # Full paths
    # {'Gd': '/full/path/Gd.pbe...UPF'}
    ...
}

# Se pseudo_dir não foi setado:
if self.pseudopotentials_base_path and 'pseudo_dir' not in params['input_data']:
    params['input_data']['pseudo_dir'] = './pseudo'
```

### 6. Remote Transfer em _transfer_pseudopotentials
```python
pseudopotentials = self.calc.parameters.get("pseudopotentials", {})
# {'Gd': '/full/path/to/Gd.pbe-spfn-rrkjus_psl.1.0.0.UPF'}

search_dirs = []
if "ESPRESSO_PSEUDO" in os.environ:
    search_dirs.append(os.environ["ESPRESSO_PSEUDO"])  # Pode estar vazio!

for symbol, pseudo_file in pseudopotentials.items():
    # pseudo_file = '/full/path/Gd.pbe...UPF'
    for pseudo_dir in search_dirs:
        local_path = os.path.join(pseudo_dir, pseudo_file)
        # os.path.join(mixed_dir, '/absolute/path') → '/absolute/path'
        # Então tenta: os.path.exists('/full/path/Gd.pbe...UPF')
        # SE esse arquivo NÃO existe aqui, FALHA
```

❌ **PROBLEMA CRÍTICO**: 
- `pseudopotentials` tem **full paths** 
- Mas `_transfer_pseudopotentials` assume **filenames only**
- Quando faz `os.path.join(pseudo_dir, pseudo_file)` com um absolute path, retorna só o absolute path
- Se `ESPRESSO_PSEUDO` não está setado, `search_dirs` fica vazio ou com diretórios errados
- Resultado: arquivo não encontrado

---

## SOLUÇÃO

Precisa haver consistência:

### OPÇÃO A: Sempre armazenar SOB FILENAMES, deixar busca para transfer
```python
# Em ConvergenceWorkflow.__init__ com dict:
pseudopotentials, base_path = discover_pseudopotential_directory(pseudopotentials)
# Depois extrai FILENAMES e armazena base_path:
self.pseudopotentials = {el: os.path.basename(p) for el, p in pseudopotentials.items()}
self.pseudopotentials_base_path = base_path
os.environ['ESPRESSO_PSEUDO'] = base_path
```

### OPÇÃO B: Sempre armazenar FULL PATHS, deixar transfer reconhecer
```python
# Em _transfer_pseudopotentials:
for symbol, pseudo_file in pseudopotentials.items():
    if os.path.isabs(pseudo_file):
        # Já é full path, usa direto
        if os.path.exists(pseudo_file):
            remote_path = os.path.join(remote_pseudo_dir, os.path.basename(pseudo_file))
            self.remote.send_file(pseudo_file, remote_path)
    else:
        # É filename, busca nos search_dirs
        for pseudo_dir in search_dirs:
            local_path = os.path.join(pseudo_dir, pseudo_file)
            if os.path.exists(local_path):
                # Transfer
```

### OPÇÃO C: Normalizar para pseudo_config (mais limpo)
Passa SEMPRE como config (se foi dict, cria um config em memória temporário)
