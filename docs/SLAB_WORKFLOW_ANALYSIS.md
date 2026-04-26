# Workflow de Criação e Otimização de Slabs

## 📋 Arquitetura Existente

### Classes Atuais

1. **OER_bulk** (xespresso/workflow/oer.py)
   - Cria slabs a partir de estrutura bulk
   - Implementa SlabGenerator do pymatgen
   - Parâmetros: min_slab_size, min_vacuum_size, nlayer, fix (quais camadas fixar)

2. **OER_site** (herda de OER_bulk)
   - Adiciona adsorbatos em sites específicos
   - Calcula adsorção de O*, OH*, OOH*
   - Diagramas de Pourbaix e OER

### Fluxo Atual (OER)
```
bulk (relax) 
  ↓
surfaces (relax) 
  ↓
adsorbates (VASP+ASE)
  ↓
thermochemistry (free energy)
```

---

## 🔧 Proposta: SlabWorkflow Integrado

### 1. Fase de Criação (novo módulo: slab_workflow.py)

```python
class SlabWorkflow:
    """Workflow completo para slabs com ConvergenceWorkflow"""
    
    def __init__(
        self,
        bulk_atoms,
        surface_indices=[(1,0,0), (1,1,0), (1,1,1)],
        min_slab_size=6.0,      # Angstrom
        min_vacuum_size=15.0,   # Angstrom
        nlayer=4,               # camadas
        fix_layers=[0, 3],      # quais camadas fixar
        pseudopotentials=None,
        pseudopotentials_config="default",
        machine="medusa",
        code_version="7.4.1"
    ):
        self.bulk_atoms = bulk_atoms
        self.surface_indices = surface_indices
        self.slabs = {}  # {(1,0,0): Atoms object, ...}
```

### 2. Fase de Convergência de Slab

**Diferenças de uma SCF normal:**

```python
def _prepare_slab_for_convergence(slab, surface_index):
    """
    Preparar slab para convergência
    
    Diferenças do bulk:
    1. K-mesh anisotrópico (x,y) > z
       - Exemplo: (12,12,1) em vez de (8,8,8)
    2. Dipole correction (edir=3)
    3. Quatro camadas: 2 primeiras fixas, 2 últimas relaxam
    """
    # k-spacing convergência em x,y (não em z)
    # Para Au(111): use (12,12,1) equivalente a kspacing=0.2 em x,y
```

### 3. Arquitetura Proposta

```
Workflow Slab
├── ETAPA 1: Bulk Convergence
│   └── run_bulk_convergence()
│       ├── ecutwfc convergence
│       ├── kspacing convergence (k-mesh 3D)
│       └── Get optimal_ecutwfc, optimal_kspacing
│
├── ETAPA 2: Gerar Slabs
│   └── generate_slabs()
│       ├── Para cada (hkl): SlabGenerator(bulk, hkl, ...)
│       └── Adicionar dipole correction
│
├── ETAPA 3: Convergência Slab (POR SUPERFÍCIE)
│   └── run_slab_convergence(surface_index)
│       ├── Usemesmos ecutwfc do bulk
│       ├── Convergência kspacing 2D (x,y apenas)
│       ├── Testar diferentes número de camadas
│       │   (4, 5, 6 layers)
│       └── Testar diferentes fix_layers
│
├── ETAPA 4: Relaxação Estrutural
│   └── relax_slab(surface_index)
│       ├── vc-relax (relaxar célula + átomos)
│       └── Fixar bottom 2 layers
│
└── ETAPA 5: Análise
    └── analyze_surfaces()
        ├── Surface energy
        ├── Reconstituição
        └── Adsorption sites
```

---

## 💻 Código de Exemplo

### Versão Simples (hoje, usando OER_bulk)

```python
from xespresso.workflow.oer import OER_bulk
from ase.build import bulk

# Criar bulk
bulk_atoms = bulk('Au', 'fcc', a=4.08)

# Gerar slabs
oer_bulk = OER_bulk(
    bulk_atoms,
    indexs=[(1,0,0), (1,1,0), (1,1,1)],
    min_slab_size=6.0,
    min_vacuum_size=15.0,
    nlayer=4,
    fix=[0, 3],  # Fixar primeiras 2 e últimas 2 camadas
    calculator={
        'pseudopotentials': ...,
        'ecutwfc': 60,  # Assumir já convergido
        'kpts': (12, 12, 1),  # anisotrópico!
        ...
    }
)

# Obter slabs
slabs = oer_bulk.get_slabs()  # {(1,0,0): slab, ...}
```

### Versão Proposta (com ConvergenceWorkflow)

```python
from xespresso.workflow.slab_workflow import SlabWorkflow
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow
from ase.build import bulk

# PASSO 1: Convergência do bulk
bulk_atoms = bulk('Au', 'fcc', a=4.08)
bulk_conv = ConvergenceWorkflow(
    atoms=bulk_atoms,
    pseudopotentials_config="default",
    precision='low',
    machine="medusa"
)
bulk_results = bulk_conv.run_convergence(
    label_prefix='bulk_conv'
)
bulk_recommendations = bulk_conv.get_recommendations()

# PASSO 2: SlabWorkflow com parâmetros otimizados
slab_wf = SlabWorkflow(
    bulk_atoms=bulk_atoms,
    surface_indices=[(1,0,0), (1,1,0), (1,1,1)],
    pseudopotentials_config="default",
    machine="medusa",
    ecutwfc=bulk_recommendations['optimal_ecutwfc'],
    # kspacing será tratado diferentemente para slabs!
)

# PASSO 3: Convergência específica para slabs
slab_conv_results = slab_wf.run_slab_convergence(
    surface_index=(1,1,1),
    # Teste diferentes:
    nlayers=[3, 4, 5, 6],
    fix_layer_schemes=[[0,1], [0,2], [0,3]],
    label_prefix='slab_conv_au111'
)

# PASSO 4: Relaxação de slabs com parâmetros otimizados
relaxed_slabs = slab_wf.relax_slabs(
    surface_indices=[(1,0,0), (1,1,0), (1,1,1)],
    relax_type='vc-relax',
    label_prefix='relax_au'
)

# PASSO 5: Análise
surface_energies = slab_wf.calculate_surface_energies(
    bulk_energy=bulk_results['final_energy'],
    relaxed_slabs=relaxed_slabs
)
```

---

## 🔑 Pontos Técnicos Críticos

### 1. K-mesh para Slabs (Anisotrópico)

**Bulk**: Isotrópico
```python
kpts = (8, 8, 8)
kspacing = 0.27 Â⁻¹
```

**Slab (fcc 111)**:Anisotrópico em z
```python
# x,y: converge normalmente (0.27 Â⁻¹)
# z: sempre 1 (não há periódico)

kpts_xy = 12  # (0.27 Â⁻¹ em x,y)
kpts = (12, 12, 1)
```

### 2. Dipole Correction

```python
slab.center(vacuum=15, axis=2)

# Adicionar campo elétrico oposto para anular dipolo
# Necessário para slabs 2D
dip_correction = {
    'dipole': True,
    'edir': 3,     # perpendicular à superfície (z)
    'emaxpos': 0.9,
    'eopreg': 0.1
}
```

### 3. Convergência de Camadas

```
∆E(5 layers) vs 4 layers < 1 meV/atom?
∆E(6 layers) vs 5 layers < 0.5 meV/atom?

→ Parar quando ∆E entre camadas < threshold
```

### 4. Estratégia de Fixação

```
Au(111) com 4 camadas:
┌─────────────────┐
│ Layer 3 (relaxa)│  ← Free
│ Layer 2 (relaxa)│  ← Free  
├─────────────────┤
│ Layer 1 (fixa)  │  ← Fixed (bulk-like)
│ Layer 0 (fixa)  │  ← Fixed (bulk-like)
└─────────────────┘

Fix com: ase.constraints.FixAtoms(indices=[...])
```

---

## 📊 Estrutura de Dados Recomendada

```python
slab_results = {
    'Au111': {
        'atoms': Atoms,
        'energy': float,
        'n_layers': int,
        'fixed_layers': [0,1,2],
        'surface_energy': float,
        'k_convergence': {
            'kpts': (12,12,1),
            'kspacing_xy': 0.27,
            'converged': True
        }
    },
    'Au100': {
        # ...
    },
    'Au110': {
        # ...
    }
}
```

---

## ⚡ Fases de Implementação

### Fase 1 (Fácil): Integração com OER existente
- Usar `OER_bulk.get_slabs()` existente
- Adaptar `ConvergenceWorkflow` para slabs (k-mesh anisotrópico)
- Tempo: ~2 horas

### Fase 2 (Média): SlabWorkflow com automação
- Classe `SlabWorkflow` que encapsula:
  - Criação de slabs
  - Convergência 2D de k-points
  - Convergência de espessura (n_layers)
- Tempo: ~8 horas

### Fase 3 (Complexa): Análise de superfícies
- Energy landscape de diferentes (hkl)
- Reconstrução automática
- Surface thermochemistry
- Tempo: ~16 horas

---

## 🎯 Checklist de Design

- [ ] K-mesh anisotrópico em `ConvergenceWorkflow`
- [ ] Dipole correction automática para slabs
- [ ] Convergência de camadas (n_layers sweep)
- [ ] Armazenamento eficiente de múltiplos slabs
- [ ] Cálculo de surface energy: `γ = (E_slab - n_bulk*E_atom) / (2*Area)`
- [ ] Support para multiple terminations
- [ ] Constraints automáticas (fix bottom layers)

---

## 🔗 Integração com Código Existente

```
SlabWorkflow
├── Herda estrutura de ConvergenceWorkflow
├── Usa CalculationWorkflow.run_relax()
├── Importa de OER_bulk.get_slabs()
└── Reutiliza ConvergenceWorkflow mas com:
    ├── K-mesh anisotrópico
    └── Constraints customizados
```

---

## Exemplo: Au(111)

```python
# Bulk
Au_bulk = bulk('Au', 'fcc', a=4.0782)
# ecutwfc=60 Ry, kspacing=0.27 (converged)

# Slab Au(111) com 4 camadas
Au111 = fcc111('Au', size=(4,4,4), a=4.0782, vacuum=15.0)
# k-mesh: (12, 12, 1)  ← Anisotrópico!
# Fix: primeiras 2 camadas

# Superfície (111) energia
γ_111 ≈ (E_slab - 4*E_atom) / (2 × Area) ≈ 0.124 J/m² (experimental)
```

