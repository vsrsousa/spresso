# Exemplos Práticos: Ajustar Estrutura QE no run_slab_convergence

Exemplos diretos e copy-paste para ajustar a estrutura QE gerada pelo PyMatGen.

---

## 🎯 Exemplo 1: Teste Rápido com Vácuos Pequenos

**Use case**: Testes iniciais, verificar se tudo está funcionando

```python
from ase.build import bulk
from xespresso.workflow.slab_workflow import SlabWorkflow

bulk_au = bulk('Au', 'fcc', a=4.08)
slab_wf = SlabWorkflow(
    bulk_atoms=bulk_au,
    surface_indices=[(1, 1, 1)],
    protocol='moderate',
)
slab_wf.generate_slabs()

# ✅ EXEMPLO 1: Teste rápido
slab_conv = slab_wf.run_slab_convergence(
    surface_index=(1, 1, 1),
    use_primitive_cell=True,              # 1×1 = rápido
    nlayers_for_vacuum=4,                 # Slab fina para teste
    vacuum_test=[5, 7, 9, 11],            # Vácuos pequenos
    protocol='fast',                      # ecutwfc=30 Ry (bem rápido)
    convergence_tol=0.01,                 # Menos rigoroso (10 meV/atom)
    skip_calculations=False,
    label_prefix='au111_quick',
    machine='medusa',
)

print(f"Vácuo ótimo: {slab_conv['optimal_vacuum']:.1f} Å")
```

---

## 🎯 Exemplo 2: Convergência Realista com Múltiplos Vácuos

**Use case**: Produção real, vácuos mais realistas

```python
# ✅ EXEMPLO 2: Convergência realista
slab_conv = slab_wf.run_slab_convergence(
    surface_index=(1, 1, 1),
    use_primitive_cell=True,              # Primitiva (rápido)
    nlayers_for_vacuum=6,                 # Slab realista
    
    # Vácuos realistas para Au
    vacuum_test=[10, 12, 15, 18, 20, 25, 30],
    
    protocol='moderate',                  # ecutwfc=40 Ry
    convergence_tol=0.001,                # Critério rigoroso (1 meV/atom)
    
    skip_calculations=False,              # Fazer SCF real
    label_prefix='au111_conv',
    
    code_version="7.4.1",
    machine='medusa',
    walltime='02:00:00',                  # 2 horas
    job_timeout=7200,
)

# Analisar resultados
print("\nResultados de Convergência:")
for vac, energy in sorted(slab_conv['vacuum_results'].items()):
    print(f"  Vácuo {vac:2.0f} Å: {energy:.8f} eV/atom")

print(f"\n✓ Vácuo ótimo: {slab_conv['optimal_vacuum']:.1f} Å")
print(f"✓ Convergência alcançada: {slab_conv['converged']}")
```

---

## 🎯 Exemplo 3: Supercela (2×2) - Mais Estável

**Use case**: Quando primitiva é instável, precisa de supercela

```python
# ✅ EXEMPLO 3: Usar supercela (2×2) em vez de primitiva
slab_conv = slab_wf.run_slab_convergence(
    surface_index=(1, 1, 1),
    use_primitive_cell=False,             # ← Usa supercela 2×2
    nlayers_for_vacuum=6,
    vacuum_test=[15, 18, 20, 25],        # Vácuos maiores para supercela
    
    protocol='moderate',
    convergence_tol=0.001,
    
    label_prefix='au111_supercell',
    machine='medusa',
    walltime='03:00:00',                  # Mais tempo (supercela é maior)
)

# A supercela é ~4× maior (2×2 em xy), então:
# - SCF é ~4-6× mais lento
# - Mas pode ser mais estável (menos efeitos de tamanho finito)
```

---

## 🎯 Exemplo 4: Teste com Skip (Mock Mode)

**Use case**: Validar setup sem rodar SCF (teste rápido)

```python
# ✅ EXEMPLO 4: Mock mode - teste estrutura sem SCF real
slab_conv = slab_wf.run_slab_convergence(
    surface_index=(1, 1, 1),
    use_primitive_cell=True,
    nlayers_for_vacuum=6,
    vacuum_test=[5, 10, 15, 20, 25],
    
    protocol='moderate',
    convergence_tol=0.001,
    
    skip_calculations=True,               # ← NÃO roda SCF
    label_prefix='au111_test',
)

# Resultado: Instantâneo, energias são mock values
# Útil para: verificar estrutura, paths, logging, etc
print("Setup validado, estruturas OK!")
```

---

## 🎯 Exemplo 5: Protocolo Accurato (Ciência séria)

**Use case**: Paper, publicação, resultados com alta acurácia

```python
# ✅ EXEMPLO 5: Alta acurácia (paper quality)
slab_conv = slab_wf.run_slab_convergence(
    surface_index=(1, 1, 1),
    use_primitive_cell=True,
    nlayers_for_vacuum=7,                 # Mais camadas
    
    # Vácuos bem distribuídos
    vacuum_test=[10, 12, 15, 18, 20, 25, 30],
    
    protocol='accurate',                  # ← ecutwfc=60 Ry (mais acurado)
    convergence_tol=0.0001,               # ← Muito rigoroso (0.1 meV/atom)
    
    skip_calculations=False,
    label_prefix='au111_accurate',
    
    code_version="7.4.1",
    machine='medusa',
    walltime='04:00:00',                  # Bastante tempo
    job_timeout=14400,                    # 4 horas max
)

# Resultado: Mais preciso, mas bem mais lento
print("Convergência de alta acurácia concluída")
```

---

## 🎯 Exemplo 6: Múltiplas Superfícies em Sequência

**Use case**: Estudar várias superfícies (Au bulk com (100), (110), (111))

```python
# ✅ EXEMPLO 6: Múltiplas superfícies
surfaces = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
results = {}

for surface in surfaces:
    print(f"\n{'='*70}")
    print(f"Analisando {surface}")
    print(f"{'='*70}")
    
    slab_conv = slab_wf.run_slab_convergence(
        surface_index=surface,
        use_primitive_cell=True,
        nlayers_for_vacuum=5,
        vacuum_test=[10, 15, 20, 25],
        protocol='moderate',
        convergence_tol=0.001,
        label_prefix=f'au_{surface[0]}{surface[1]}{surface[2]}_conv',
        machine='medusa',
    )
    
    results[surface] = slab_conv

# Resumo comparativo
print("\n" + "="*70)
print("RESUMO: Vácuos Ótimos por Superfície")
print("="*70)
for surface, conv in results.items():
    print(f"{surface}: {conv['optimal_vacuum']:.1f} Å")
```

---

## 🎯 Exemplo 7: Customização de Input Data (Avançado)

**Use case**: Ajustar pseudopotenciais, vdW correction, ocupação, etc.

```python
# ✅ EXEMPLO 7: Customização avançada (herança de classe)
from xespresso.workflow.slab_workflow import SlabWorkflow

class CustomAuSlabWorkflow(SlabWorkflow):
    """SlabWorkflow customizado para Au com vdW correction"""
    
    def _prepare_input_data(self):
        """Override para adicionar vdW correction"""
        base = super()._prepare_input_data()
        
        # Adicionar Grimme D2 (vdW)
        base['system']['vdw_corr'] = 'grimme-d2'
        
        # Ajustar smearing (mais conservador)
        base['system']['degauss'] = 0.01
        base['system']['occupations'] = 'smearing'
        base['system']['smearing'] = 'gaussian'
        
        # Mais rigoroso na convergência eletrônica
        base['electrons']['conv_thr'] = 1.0e-9
        base['electrons']['mixing_beta'] = 0.4  # Mais conservador
        
        return base

# Usar classe customizada
bulk_au = bulk('Au', 'fcc', a=4.08)
slab_wf = CustomAuSlabWorkflow(
    bulk_atoms=bulk_au,
    surface_indices=[(1, 1, 1)],
    protocol='moderate',
)
slab_wf.generate_slabs()

# Agora run_slab_convergence usa seus customizations
slab_conv = slab_wf.run_slab_convergence(
    surface_index=(1, 1, 1),
    use_primitive_cell=True,
    nlayers_for_vacuum=6,
    vacuum_test=[10, 15, 20, 25, 30],
    protocol='moderate',
    machine='medusa',
)

print("✓ Cálculos com vdW correction concluídos")
```

---

## 🎯 Exemplo 8: Monitoramento Detalhado

**Use case**: Entender exatamente o que está acontecendo

```python
import logging

# ✅ EXEMPLO 8: Logging detalhado
# Habilitar debug logging
logging.basicConfig(level=logging.DEBUG)

slab_conv = slab_wf.run_slab_convergence(
    surface_index=(1, 1, 1),
    use_primitive_cell=True,
    nlayers_for_vacuum=6,
    vacuum_test=[10, 15, 20],
    protocol='moderate',
    
    # Ajustes para monitoramento
    skip_calculations=False,
    label_prefix='au111_debug',
    machine='medusa',
)

# Saída será muito detalhada, mostrando:
# - Estrutura da slab gerada
# - K-mesh calculado
# - Cada job submetido
# - Convergência de cada vácuo
```

---

## 📊 Comparação Rápida: Qual Opção Escolher?

| Cenário | Parâmetros |
|---------|-----------|
| **Teste rápido** | `protocol='fast'`, `use_primitive_cell=True`, `vacuum_test=[5,10,15]`, `skip_calculations=True` |
| **Produção normal** | `protocol='moderate'`, `use_primitive_cell=True`, `vacuum_test=[10,15,20,25,30]`, `convergence_tol=0.001` |
| **Alta acurácia (paper)** | `protocol='accurate'`, `use_primitive_cell=True`, `vacuum_test=[10,12,15,18,20,25,30]`, `convergence_tol=0.0001` |
| **Supercela (estável)** | `use_primitive_cell=False`, `protocol='moderate'`, vácuos maiores |
| **Mock (validação)** | `skip_calculations=True`, qualquer outro parâmetro |

---

## 🔍 Entendendo a Saída

```python
slab_conv = slab_wf.run_slab_convergence(...)

# slab_conv é um dicionário com:
{
    'surface_index': (1, 1, 1),           # Qual superfície
    'kmesh_calc': (18, 18, 1),            # K-mesh anisotropico
    'vacuum_results': {                   # Energias por vácuo
        5.0: -5.123456,                   # eV/atom
        7.0: -5.124123,
        9.0: -5.124567,
        11.0: -5.124589,                  # ← Mínimo aqui
        13.0: -5.124590,                  # Praticamente igual
    },
    'optimal_vacuum': 11.0,               # Vácuo que minimiza energia
    'converged': True,                    # Convergência alcançada?
    'layer_results': {},                  # Vazio (test_nlayers=False default)
    'optimal_nlayers': None,
}
```

---

## 💾 Salvando Resultados

```python
import json

# Salvar resultados
with open('slab_convergence_results.json', 'w') as f:
    # Converter numpy types para tipos Python
    safe_results = {
        'surface_index': slab_conv['surface_index'],
        'optimal_vacuum': float(slab_conv['optimal_vacuum']),
        'kmesh': slab_conv['kmesh_calc'],
        'converged': bool(slab_conv['converged']),
        'vacuum_results': {
            float(k): float(v) for k, v in slab_conv['vacuum_results'].items()
        }
    }
    json.dump(safe_results, f, indent=2)

print("✓ Resultados salvos em slab_convergence_results.json")
```

---

## 🚀 Próximos Passos Após Convergência

Depois de `run_slab_convergence()` completar:

```python
# 1. Use vácuo ótimo para relaxação
optimal_vac = slab_conv['optimal_vacuum']

# 2. Rodar Phase 4: relaxação de estrutura
# slab_relax = slab_wf.run_slab_relax(
#     surface_index=(1, 1, 1),
#     optimal_vacuum=optimal_vac,
#     nlayers_test=[3, 4, 5, 6, 7],  # Testar espessura com relaxação
#     relax_type='relax',
# )

# 3. Depois: calcular energias de superfície
# surface_energy = slab_wf.calculate_surface_energy(
#     surface_index=(1, 1, 1),
#     # ... parâmetros
# )
```
