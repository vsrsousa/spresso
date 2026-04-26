# Exemplo de uso: run_slab_relax() com nosym parameter

Este é o código correto que você pode usar agora:

```python
from xespresso.workflow.slab_workflow import SlabWorkflow

# Exemplo 1: Relaxação com parâmetros convergidos (RECOMENDADO)
relax_results = slab_wf.run_slab_relax(
    relax_type='relax',          # 'relax' (átomos) ou 'vc-relax' (célula)
    surfaces=[(1, 1, 1)],        # Relaxar apenas Au(111) para exemplo
    nlayers_test=[3, 4, 6, 8],   # ← Phase 4b: testar convergência de camadas
    use_primitive_cell=True,
    fmax=0.05,                   # Força máxima < 0.05 eV/Å
    machine='medusa',
    code_version='7.4.1',
    job_timeout=7200,            # 2 horas (relaxações são mais lentas)
    protocol=wf.protocol,
    # nosym não precisa ser especificado - usa True (padrão recomendado)
)
```

## Parâmetros:

### ✅ Novos parâmetros adicionados:

| Parâmetro | Tipo | Default | Descrição |
|-----------|------|---------|-----------|
| `nosym` | `bool` \| `None` | `None` → `True` | Desabilitar simetria durante relaxação |

### nosym behavior:

- **`nosym=None` (padrão)** → Usa `nosym=.true.` no QE
  - ✅ RECOMENDADO para geometry optimization
  - Desabilita operações de simetria para evitar violações
  - Safe default para qualquer estrutura

- **`nosym=True`** → Usa `nosym=.true.` no QE
  - Explicitly disable symmetry
  - Use quando sabe que relaxação pode quebrar simetrias

- **`nosym=False`** → Usa `nosym=.false.` no QE
  - Mantém simetrias ativas
  - ⚠️ Mais rápido mas arriscado - pode quebrar durante relaxação
  - Use apenas se confiante que estrutura mantém simetria

## Exemplos de uso:

### 1️⃣ Padrão (recomendado - sem especificar nosym):
```python
relax_results = slab_wf.run_slab_relax(
    relax_type='relax',
    surfaces=[(1, 1, 1)],
    nlayers_test=[3, 4, 6, 8],
    fmax=0.05,
    machine='medusa',
    # nosym não especificado → usa default True (seguro)
)
```

### 2️⃣ Explicitamente desabilitar simetria (conservador):
```python
relax_results = slab_wf.run_slab_relax(
    relax_type='relax',
    surfaces=[(1, 1, 1)],
    nlayers_test=[3, 4, 6, 8],
    nosym=True,  # ← Explicitamente desabilita
    fmax=0.05,
    machine='medusa',
    code_version='7.4.1',
    job_timeout=7200,
    protocol=wf.protocol,
    walltime='2:00:00',  # 2 horas no scheduler
)
```

### 3️⃣ Manter simetria ativada (risky, apenas se souber o que faz):
```python
relax_results = slab_wf.run_slab_relax(
    relax_type='relax',
    surfaces=[(1, 1, 1)],
    nlayers_test=[3, 4, 6, 8],
    nosym=False,  # ← Mantém simetrias (pode ser mais rápido)
    fmax=0.05,
)

# ⚠️ Aviso: pode quebrar se relaxação viola simetrias
```

### 4️⃣ Phase 4 (relaxação simples, sem nlayers_test):
```python
relax_results = slab_wf.run_slab_relax(
    relax_type='relax',
    surfaces=[(1, 1, 1)],
    # nlayers_test NÃO especificado → Phase 4 (relaxação única)
    fmax=0.05,
    nosym=True,  # Controla simetria também em Phase 4
    machine='medusa',
)
```

## Integração com workflow completo:

```python
# Phase 1: Convergência do bulk
bulk_conv = wf.run_bulk_convergence(
    max_ecutwfc=200.0,
    ecutwfc_step=10.0,
)

# Phase 3: Convergência do slab (vacuum)
slab_conv = slab_wf.run_slab_convergence(
    surface_index=(1, 1, 1),
    vacuum_test=[10, 12, 15, 18, 20, 25, 30],
    nlayers_test=[3, 4, 5, 6, 7],
    test_nlayers=False,  # ← Phase 3: só vacuum (SCF rápido)
    machine='medusa',
    protocol=wf.protocol,
)

# Phase 4b: Relaxação com convergência de camadas
relax_results = slab_wf.run_slab_relax(
    relax_type='relax',
    surfaces=[(1, 1, 1)],
    nlayers_test=[3, 4, 6, 8],  # ← Testa camadas com relaxação
    fmax=0.05,
    nosym=True,  # ← Desabilita simetria durante geometry opt
    machine='medusa',
    code_version='7.4.1',
    job_timeout=7200,
    protocol=wf.protocol,
    walltime='2:00:00',  # SLURM, PBS, etc - scheduler agnostic
)
```

## Resumo da correção:

| Antes | Depois |
|-------|--------|
| ❌ `nosym=False` → TypeError | ✅ `nosym` parâmetro adicionado |
| N/A | ✅ Default: `nosym=None` → `True` (seguro) |
| N/A | ✅ Documentação completa no docstring |
| N/A | ✅ Funciona em Phase 4 e Phase 4b |
| N/A | ✅ Parametrizável: `True`, `False`, ou `None` |

## Quando usar cada valor:

### Use `nosym=True` (ou deixar padrão None):
- ✅ Geometry optimization de superfícies
- ✅ Quando estrutura pode quebrar simetrias
- ✅ Quando quer resultado conservador/seguro
- ✅ **Recomendado para 99% dos casos**

### Use `nosym=False`:
- ⚠️ Quando estrutura é altamente simétrica
- ⚠️ Quando quer cálculos mais rápidos
- ⚠️ Quando confiante na estabilidade da simetria
- ⚠️ **Não recomendado para inicialização de relaxações**

## Verificação:

```python
import inspect
from xespresso.workflow.slab_workflow import SlabWorkflow

# Verificar que nosym está presente
sig = inspect.signature(SlabWorkflow.run_slab_relax)
print('nosym' in sig.parameters)  # True ✅
```
