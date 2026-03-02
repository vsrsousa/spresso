# Análise: Independência de Arquivos no Convergence Batch Mode

## Pergunta
Os arquivos em lote do `ecutwfc` são independentes do `kspacing`?

## Resposta: **NÃO COMPLETAMENTE**

## Análise Detalhada

### 1. O que é Independente (✅)

Os **arquivos de INPUT** (`.pwi`) são **independentes** em relação a:
- Diferentes `kspacing` → diferentes k-points grid → **ARQUIVOS DIFERENTES**
- Diferentes `ecutwfc` → mesmo diretor → **ARQUIVOS DIFERENTES**

Cada combinação `(ecutwfc, kspacing)` gera um arquivo `.pwi` **único com parâmetros diferentes**.

### 2. O que NÃO é Independente (❌)

Os **PSEUDOPOTENCIAIS** são **REENVIADOS múltiplas vezes**:

#### Fluxo no Batch Mode:

```
convergence_workflow.py:
  ├─ _run_convergence_study_batch()
  │  ├─ Para cada ecutwfc:
  │  │  ├─ Prepare batch_params[] com todos os kspacing
  │  │  └─ workflow.submit_scf_batch_multiple(batch_params)
  │  │     │
  │  │     └─ calculation_workflow.py:
  │  │        └─ submit_scf_batch_multiple()
  │  │           └─ Para cada parâmetro em batch_params:
  │  │              ├─ temp_workflow = CalculationWorkflow(kspacing=params.get('kspacing'))
  │  │              │  # NOVO WORKFLOW COM KSPACING DIFERENTE
  │  │              │
  │  │              └─ temp_workflow.submit_scf_batch()
  │  │                 ├─ calc.write_input()     # arquivo .pwi único
  │  │                 └─ calc.execute()         # chama run()
  │  │                    └─ remote_mixin.py:run()
  │  │                       ├─ _transfer_pseudopotentials()  # ❌ REENVIADO!
  │  │                       ├─ send_file(.pwi)
  │  │                       ├─ send_file(job_file)
  │  │                       └─ sbatch
```

### 3. O Problema Real

Para um convergence com:
- **6 valores de ecutwfc**: [40, 50, 60, 70, 80, 90]
- **4 valores de kspacing**: [0.5, 0.4, 0.3, 0.2]

**Total de submissões: 6 × 4 = 24 jobs**

No modo batch atual:
- Ecutwfc 40: submete 4 jobs (ksp=0.5, 0.4, 0.3, 0.2)
  - **Pseudo transferidos 4 vezes** ❌
- Ecutwfc 50: submete 4 jobs
  - **Pseudo transferidos 4 vezes** ❌
- ... (repetido para cada ecutwfc)

**Total: 24 transferências de pseudo!** 🔴

### 4. Impacto no Tempo

Se cada pseudopotencial leva ~1-2 segundos para transferir:
- **Com otimização**: 6 ecutwfc × 1 transfer by grupo ≈ 6-12 seg
- **Sem otimização (atual)**: 24 jobs × 1 transfer ≈ 24-48 seg

**Overhead: ~2-4x mais lento com `kspacing` variável!** ⏱️

## Verificação no Código

### convergence_workflow.py líneas 965-1000:
```python
# Para CADA ecutwfc
while current_ecutwfc <= max_ecutwfc:
    # Prepare batch_params para TODOS os kspacing
    batch_params = []  # Com kspacing diferente em cada item
    
    # AQUI É O PROBLEMA:
    batch_results = workflow.submit_scf_batch_multiple(batch_params)
    #                         ↓
    # calculation_workflow.py línea 1052
```

### calculation_workflow.py línea 1052:
```python
def submit_scf_batch_multiple(self, parameter_sets: List[Dict]):
    # Para CADA parâmetro em parameter_sets:
    for i, params in enumerate(parameter_sets):
        # NOVO WORKFLOW COM KSPACING DIFERENTE!
        temp_workflow = CalculationWorkflow(
            self.atoms,
            protocol=self.protocol,
            pseudopotentials=self.pseudopotentials,
            kspacing=params.get('kspacing'),  # ← MUDA CADA VEZ
            input_data=input_data_override,
            queue=self.queue,
            **self.extra_kwargs)
        
        # CADA WORKFLOW CHAMA EXECUTE() SEPARADAMENTE:
        result = temp_workflow.submit_scf_batch(label=label, wait_for_completion=False)
        #        ↓
        #        Chama: calc.execute() → run() → _transfer_pseudopotentials() ❌
```

### remote_mixin.py línea 155:
```python
def run(self):
    # ...
    self._transfer_pseudopotentials()  # ← REENVIADO AQUI
    self.remote.send_file(local_input, ...)
    self.remote.send_file(local_job, ...)
    self.submit_command()
```

## Recomendação

### Opção 1: Optimizar o Batch (RECOMENDADO)
Modificar `submit_scf_batch_multiple()` para:
1. Transferir pseudopotenciais **UMA VEZ** no início
2. Depois submeter todos os jobs rapidamente

```python
def submit_scf_batch_multiple(self, parameter_sets):
    # 1. Transferir pseudo UMA VEZ
    self._transfer_pseudopotentials()  
    
    # 2. Para cada parameter (sem reenviar pseudo):
    for params in parameter_sets:
        temp_workflow = CalculationWorkflow(...)
        # Desabilitar _transfer_pseudopotentials() aqui
        result = temp_workflow.submit_scf_batch_fast(label, skip_pseudo_transfer=True)
```

**Ganho**: 4-6x mais rápido para batch com múltiplos kspacing

### Opção 2: Usar Sequential Mode
Se a velocidade é crítica, usar modo sequencial onde os pseudopotenciais 
são reutilizados dentro da mesma conexão remota.

```python
return self._run_convergence_study_sequential(...)  # Reutiliza conexão
```

## Status Atual

❌ **NÃO otimizado**: Pseudopotenciais reenviados múltiplas vezes
✅ **Funcionalmente correto**: Todos os cálculos são executados
⚠️ **Performance**: Lento para grandes convergence studies com múltiplos kspacing

## Recomendação para Próximo Commit

1. Otimizar `submit_scf_batch_multiple()` para reuso de pseudo
2. Adicionar flag `--skip-pseudo-transfer` do Espresso/scheduler
3. Adicionar logging para mostrar quantas vezes pseudo foi transferido
4. Documentar na docstring do convergence_workflow sobre essa limitação

