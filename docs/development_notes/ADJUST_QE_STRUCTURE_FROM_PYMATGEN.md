# Como Ajustar a Estrutura QE Gerada pelo PyMatGen na Convergência de Slab

Quando você chama `run_slab_convergence()`, a estrutura QE é gerada automaticamente a partir da slab criada pelo PyMatGen. Aqui mostramos como ajustar essa estrutura em cada etapa.

---

## 📋 Fluxo de Geração da Estrutura

```
1. PyMatGen SlabGenerator
   ├─ Cria slab com PyMatGen
   │  └─ min_slab_size, min_vacuum_size, miller_index
   │
2. Post-processamento ASE
   ├─ Ortogonalizar célula
   ├─ Centrar slab no vácuo
   ├─ Aparar para exatamente N camadas
   └─ Aplicar constraints (FixAtoms)
   │
3. Configuração QE
   ├─ K-mesh (calculado anisotropicamente)
   ├─ Pseudopotenciais
   ├─ Parâmetros SCF (ecutwfc, convergência)
   └─ Input parameters (ocupação, smearing, etc)
   │
4. Submissão Paralela (batch_utils)
   └─ Submit → Coletar resultados
```

---

## 🔧 OPÇÃO 1: Ajustar Parâmetros de Entrada na Chamada

**Mais direto** - Todos os parâmetros importantes estão expostos:

```python
slab_conv = slab_wf.run_slab_convergence(
    surface_index=(1, 1, 1),
    use_primitive_cell=True,              # ✓ 1×1 (primitiva, mais rápido)
    vacuum_test=[5, 7, 9],                # ✓ Vácuos para testar
    nlayers_for_vacuum=6,                 # ✓ N camadas fixo para testes de vácuo
    convergence_tol=0.001,                # ✓ Critério de convergência (eV/atom)
    skip_calculations=False,               # ✓ Se True: mock mode (teste rápido)
    label_prefix='au111_convergence',     # ✓ Prefixo para diretórios
    protocol='moderate',                  # ✓ fast/moderate/accurate → ecutwfc, kpts
    machine='medusa',                     # ✓ Máquina remota
    code_version="7.4.1",                 # ✓ Versão do QE
    job_timeout=3600,                     # ✓ Timeout para jobs remotos (segundos)
    walltime='01:30:00',                  # ✓ Tempo máximo SLURM/PBS
)
```

### O que cada parâmetro ajusta:

| Parâmetro | Afeta | Exemplo |
|-----------|-------|---------|
| `use_primitive_cell` | Tamanho da supercela (1×1 vs 2×2) | `True` = primitiva, rápido |
| `protocol` | ecutwfc, occupations, smearing | 'moderate' = balanceado |
| `nlayers_for_vacuum` | Espessura da slab na convergência de vácuo | 6 camadas |
| `vacuum_test` | Vácuos a testar | `[5,7,9,11,13]` |
| `convergence_tol` | Limite de energia ∆E | 0.001 eV/atom = 1 meV/atom |
| `job_timeout` | Tempo máximo para job remoto | 3600s = 1 hora |
| `walltime` | Tempo solicitado ao scheduler | SLURM format: '1:30:00' |

---

## 🔨 OPÇÃO 2: Customizar Antes de Chamar run_slab_convergence()

Se você quer controle fino sobre a slab **antes** da convergência:

```python
# 1. Gerar a slab
slab = slab_wf._regenerate_slab_with_nlayers(
    surface_index=(1, 1, 1),
    nlayers=6,
    vacuum_size=15.0,
    use_primitive_cell=True
)

# 2. AQUI você pode ajustar a slab manualmente
# ============================================

# Ajustar posições (ex: adsorção)
slab.positions[-1, 2] += 0.5  # Levantar último átomo 0.5 Å

# Expandir célula
slab.set_cell(slab.get_cell() * [1.05, 1.05, 1.0])  # 5% expansão em xy

# Modificar constraints
from ase.constraints import FixAtoms
constraint = FixAtoms(indices=[0, 1, 2, 3])  # Fixar primeiros 4 átomos
slab.set_constraint(constraint)

# 3. Agora rodar convergência com a slab personalizada
# (você precisaria criar um método customizado para isso)
```

---

## 🎯 OPÇÃO 3: Modificar Parâmetros QE Internos (Advanced)

Se você quer ajustar **parâmetros QE específicos** como pseudopotenciais, ecutrho, etc:

```python
# 1. Preparar input_data customizado
custom_input_data = {
    'control': {
        'calculation': 'scf',
        'verbosity': 'high',
    },
    'system': {
        'occupations': 'smearing',
        'smearing': 'gaussian',
        'degauss': 0.02,             # Broadening (Ry)
        'vdw_corr': 'grimme-d2',     # Dispersion correction
    },
    'electrons': {
        'mixing_beta': 0.3,          # Mais conservador que 0.7
        'conv_thr': 1.0e-9,          # Convergência mais rigorosa
    }
}

# 2. Modificar protocolo (intermediate)
slab_wf.protocol = 'custom'  # Flag que você define
slab_wf.input_data = custom_input_data

# 3. Agora rodar convergência
slab_conv = slab_wf.run_slab_convergence(
    surface_index=(1, 1, 1),
    protocol='custom',  # Usa seu input_data customizado
    # ... outros parâmetros
)
```

---

## 📊 Exemplo Completo: Ajustar Au(111) com Vácuo Variável

```python
from ase.build import bulk
from xespresso.workflow.slab_workflow import SlabWorkflow

# 1. Bulk
bulk_au = bulk('Au', 'fcc', a=4.08)

# 2. Criar workflow
slab_wf = SlabWorkflow(
    bulk_atoms=bulk_au,
    surface_indices=[(1, 1, 1)],
    min_vacuum_size=15.0,
    use_primitive_cell=True,  # Usar primitiva (1×1)
    protocol='moderate',
)

# 3. Gerar slabs
slabs = slab_wf.generate_slabs()

# 4. Convergência com parâmetros customizados
slab_conv = slab_wf.run_slab_convergence(
    surface_index=(1, 1, 1),
    
    # === AJUSTES DA SLAB ===
    use_primitive_cell=True,          # Usar primitiva
    nlayers_for_vacuum=6,             # 6 camadas para testes de vácuo
    
    # === VÁCUOS A TESTAR ===
    vacuum_test=[5, 7, 9, 11, 13],    # Vácuos variados
    
    # === K-MESH (automático baseado em bulk) ===
    # Não precisa ajustar - calculado anisotropicamente
    
    # === PROTOCOLO QE ===
    protocol='moderate',              # ecutwfc=40 Ry, kpts automático
    code_version="7.4.1",
    
    # === CONVERGÊNCIA ===
    convergence_tol=0.001,            # 1 meV/atom
    
    # === EXECUÇÃO REMOTA ===
    machine='medusa',
    job_timeout=3600,                 # 1 hora max
    walltime='01:30:00',              # Pedir 1.5h ao SLURM
)

# 5. Resultados
print(f"Vácuo ótimo: {slab_conv['optimal_vacuum']:.1f} Å")
print(f"Energias por vácuo:")
for vac, energy in slab_conv['vacuum_results'].items():
    print(f"  {vac:.0f} Å: {energy:.6f} eV/atom")
```

---

## 🎬 Ajustes Avançados: Hooks Customizados

Se você precisa de controle **ainda mais fino**, você pode estender `SlabWorkflow`:

```python
from xespresso.workflow.slab_workflow import SlabWorkflow

class CustomSlabWorkflow(SlabWorkflow):
    
    def _prepare_input_data(self):
        """Override para customizar input_data"""
        base_input = super()._prepare_input_data()
        
        # Adicionar ajustes customizados
        base_input['system']['vdw_corr'] = 'grimme-d2'
        base_input['system']['degauss'] = 0.015
        base_input['electrons']['mixing_beta'] = 0.4
        
        return base_input
    
    def _regenerate_slab_with_nlayers(self, surface_index, nlayers, **kwargs):
        """Override para customizar slab ANTES de usar"""
        slab = super()._regenerate_slab_with_nlayers(surface_index, nlayers, **kwargs)
        
        # Seu ajuste aqui
        # Ex: adicionar adsorbato
        # slab.positions[-1, 2] += 1.5
        
        return slab


# Usar sua classe customizada
slab_wf_custom = CustomSlabWorkflow(
    bulk_atoms=bulk('Au', 'fcc', a=4.08),
    surface_indices=[(1, 1, 1)],
)

# Agora usa seus customizações
slab_conv = slab_wf_custom.run_slab_convergence(
    surface_index=(1, 1, 1),
    vacuum_test=[5, 7, 9],
)
```

---

## 📋 Tabela: Onde Cada Parâmetro Afeta

| Aspecto | Parâmetro | Onde Ajustar | Impacto |
|---------|-----------|--------------|--------|
| **Slab size** | `nlayers_for_vacuum` | `run_slab_convergence()` | Espessura (SCF cost) |
| **Slab shape** | `use_primitive_cell` | `run_slab_convergence()` | 1×1 (rápido) vs 2×2 (estável) |
| **Vácuo** | `vacuum_test` | `run_slab_convergence()` | Isolamento entre camadas |
| **QE cutoff** | `protocol` | `run_slab_convergence()` | ecutwfc, ecutrho (acurácia vs speed) |
| **QE occupações** | `protocol` ou custom | Classe customizada | smearing, degauss |
| **Convergência SCF** | `convergence_tol` | `run_slab_convergence()` | conv_thr (precisão) |
| **K-mesh** | automático | Baseado em bulk | Anisotropicamente calculado |
| **Tempo de job** | `job_timeout` | `run_slab_convergence()` | Timeout remoto |
| **Walltime SLURM** | `walltime` | `run_slab_convergence()` | Tempo solicitado |

---

## ⚡ Dicas Práticas

### Para testes rápidos (mock mode):
```python
slab_conv = slab_wf.run_slab_convergence(
    surface_index=(1, 1, 1),
    skip_calculations=True,  # ← Não roda SCF, só retorna mock energies
    vacuum_test=[5, 10, 15],
)
```

### Para convergência rigorosa:
```python
slab_conv = slab_wf.run_slab_convergence(
    surface_index=(1, 1, 1),
    protocol='accurate',      # ecutwfc=60 Ry
    convergence_tol=0.0001,   # 0.1 meV/atom (muito rigoroso!)
    vacuum_test=[10, 12, 15, 18, 20, 25, 30],
)
```

### Para supercela (mais estável que primitiva):
```python
slab_conv = slab_wf.run_slab_convergence(
    surface_index=(1, 1, 1),
    use_primitive_cell=False,  # ← Usa 2×2 (supercela)
    # ... resto dos parâmetros
)
```

---

## ✅ Checklist: Antes de Rodar

- [ ] Bulk convergence concluída? (`run_bulk_convergence()`)
- [ ] Slabs geradas? (`generate_slabs()`)
- [ ] Vácuos a testar são realistas? (10-30 Å típico)
- [ ] `nlayers_for_vacuum` é razoável? (3-7 camadas)
- [ ] `protocol` é apropriado? (fast=rápido, moderate=balanceado, accurate=lento)
- [ ] Tempo de job é suficiente? (check `walltime`)
- [ ] Máquina remota está configurada? (`machine='medusa'`)
