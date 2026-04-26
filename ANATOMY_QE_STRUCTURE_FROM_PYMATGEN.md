# Anatomia: Como a Estrutura QE é Gerada do PyMatGen

Documento técnico mostrando exatamente **como** a estrutura QE é gerada a partir da slab do PyMatGen, para quem quer fazer modificações avançadas no código.

---

## 📍 Fluxo Completo de Geração

```
┌─────────────────────────────────────────────────────────────┐
│ ENTRADA: bulk_atoms (ASE Atoms)                             │
└──────────────┬──────────────────────────────────────────────┘
               │
               ▼
┌─────────────────────────────────────────────────────────────┐
│ STEP 1: Converter ASE → PyMatGen                            │
│ Arquivo: slab_workflow.py, método _regenerate_slab_...     │
│                                                              │
│ bulk_atoms                                                   │
│   ↓ AseAtomsAdaptor.get_structure()                         │
│ struct (PyMatGen Structure)                                 │
│   ↓ SpacegroupAnalyzer().get_conventional...()             │
│ struct_conv (Convencional - cubic para Au)                 │
└──────────────┬──────────────────────────────────────────────┘
               │
               ▼
┌─────────────────────────────────────────────────────────────┐
│ STEP 2: Calcular d-spacing Dinâmicamente                    │
│                                                              │
│ d_hkl = a / sqrt(h² + k² + l²)                             │
│                                                              │
│ Exemplo para Au(111):                                       │
│   - a = 4.08 Å (parâmetro de rede)                         │
│   - (h,k,l) = (1,1,1)                                      │
│   - d_hkl = 4.08 / sqrt(3) = 2.357 Å                       │
│   - Para 6 camadas: min_slab_size = 6 × 2.357 - 0.5 Å    │
│                     = 14.142 - 0.5 = 13.642 Å             │
└──────────────┬──────────────────────────────────────────────┘
               │
               ▼
┌─────────────────────────────────────────────────────────────┐
│ STEP 3: PyMatGen SlabGenerator                              │
│                                                              │
│ SlabGenerator(                                              │
│     struct_conv,                                            │
│     miller_index=(1,1,1),                                   │
│     min_slab_size=13.642,  ← dinâmico!                     │
│     min_vacuum_size=15.0,                                   │
│     center_slab=False,                                      │
│ )                                                            │
│                                                              │
│ Resultado: slab_list[0] (PyMatGen Slab object)             │
└──────────────┬──────────────────────────────────────────────┘
               │
               ▼
┌─────────────────────────────────────────────────────────────┐
│ STEP 4: Converter PyMatGen → ASE                            │
│                                                              │
│ slab_pymatgen                                               │
│   ↓ AseAtomsAdaptor.get_atoms()                            │
│ slab (ASE Atoms)                                            │
└──────────────┬──────────────────────────────────────────────┘
               │
               ▼
┌─────────────────────────────────────────────────────────────┐
│ STEP 5: Pós-processamento - Ortogonalizar Célula           │
│                                                              │
│ Método: _orthogonalize_cell()                              │
│                                                              │
│ PyMatGen retorna célula não-ortogonal (com xy, yz, xz)     │
│ Converter para célula ortogonal:                            │
│   - Aplicar rotação para alinhar z com vácuo               │
│   - x,y ficam paralelos à superfície                       │
│   - z aponta para o vácuo                                  │
│                                                              │
│ Antes:                 Depois:                             │
│ [[4.08, 0.00, 0.00]   [[4.08, 0.00, 0.00]                │
│  [2.04, 3.53, 0.00] →  [0.00, 4.08, 0.00]                │
│  [0.00, 0.00, 20.0]]   [0.00, 0.00, 20.0]]               │
└──────────────┬──────────────────────────────────────────────┘
               │
               ▼
┌─────────────────────────────────────────────────────────────┐
│ STEP 6: Centrar Slab no Vácuo                              │
│                                                              │
│ slab.center(vacuum=15.0, axis=2)                           │
│                                                              │
│ Posiciona os átomos centralizados na célula com vácuo      │
│ acima e abaixo igualmente distribuído                      │
│                                                              │
│ Antes center(): átomos entre z=0-10                        │
│ Depois center(): átomos entre z=5-15 (com vácuo nas pontas)│
└──────────────┬──────────────────────────────────────────────┘
               │
               ▼
┌─────────────────────────────────────────────────────────────┐
│ STEP 7: Aparar para Exatamente N Camadas                   │
│                                                              │
│ Método: _trim_to_exact_layers()                            │
│                                                              │
│ SlabGenerator pode gerar MAIS camadas do que solicitado    │
│ (min_slab_size é um MÍNIMO)                                │
│                                                              │
│ Aqui aparamos para exatamente 6 camadas:                   │
│   1. Identificar camadas (por z-position)                  │
│   2. Manter apenas as primeiras 6 camadas                  │
│   3. Remover átomos das camadas extras                     │
│                                                              │
│ Resultado: Exatamente 6 camadas atômicas                   │
└──────────────┬──────────────────────────────────────────────┘
               │
               ▼
┌─────────────────────────────────────────────────────────────┐
│ STEP 8: Aplicar Constraints (FixAtoms)                     │
│                                                              │
│ Método: _calculate_fixed_layers()                          │
│                                                              │
│ Congelar camadas inferiores (bulk-like):                   │
│   - nlayers=6: congelar primeiras 2-3 camadas              │
│   - Deixar livres: últimas 3-4 camadas (superfície)        │
│                                                              │
│ Exemplo para nlayers=6:                                    │
│   fix_indices = [0, 1] ou [0, 1, 2] (depende do método)   │
│   constraint = FixAtoms(indices=fix_indices)               │
│   slab.set_constraint(constraint)                          │
│                                                              │
│ Isso previne relaxações não-físicas (o bulk não relaxa)    │
└──────────────┬──────────────────────────────────────────────┘
               │
               ▼
┌─────────────────────────────────────────────────────────────┐
│ SAÍDA: slab (ASE Atoms, pronta para QE)                    │
│                                                              │
│ Propriedades:                                               │
│  - cell: 3×3 ortogonal, com z para vácuo                   │
│  - positions: átomos centralizados                          │
│  - constraints: FixAtoms para camadas bulk                 │
│  - pbc: [True, True, True]                                 │
└─────────────────────────────────────────────────────────────┘
```

---

## 🔍 STEP-BY-STEP: Vendo Cada Transformação

### Código Relevante no `slab_workflow.py`:

**STEP 1-2: Converter e Calcular d-spacing**
```python
# Linhas 450-480
struct = AseAtomsAdaptor.get_structure(self.bulk_atoms)
struct_conv = SpacegroupAnalyzer(struct).get_conventional_standard_structure()

h, k, l = surface_index  # ex: (1, 1, 1)
a = struct_conv.lattice.a  # ex: 4.08 Å

# Fórmula cristalográfica
d_hkl = a / np.sqrt(h**2 + k**2 + l**2)  # ex: 4.08 / sqrt(3) = 2.357 Å

# Calcular mínimo de slab
min_slab_size_dynamic = nlayers * d_hkl - 0.5  # ex: 6 * 2.357 - 0.5
```

**STEP 3: SlabGenerator**
```python
# Linhas 483-493
slabgen = SlabGenerator(
    struct_conv,
    miller_index=surface_index,
    min_slab_size=min_slab_size_dynamic,  # ← Dinâmico!
    min_vacuum_size=vac,
    center_slab=False,  # Importante: deixamos centrar depois
)
slab_pymatgen = slabgen.get_slabs(tol=0.1)[0]
slab = AseAtomsAdaptor.get_atoms(slab_pymatgen)
```

**STEP 5: Ortogonalizar**
```python
# Linhas 501
slab = self._orthogonalize_cell(slab)
```
Ver implementação em `_orthogonalize_cell()` (linhas ~600):
```python
def _orthogonalize_cell(self, slab):
    """Converte célula não-ortogonal para ortogonal"""
    # ... rotações complexas para alinhar z com vácuo
    # Resultado: célula ortogonal, z aponta "up"
```

**STEP 6: Centrar**
```python
# Linha 502
slab.center(vacuum=self.min_vacuum_size, axis=2)
```

**STEP 7: Aparar para N exato**
```python
# Linha 506
slab = self._trim_to_exact_layers(slab, nlayers, surface_index)
```
Ver implementação em `_trim_to_exact_layers()` (linhas ~670):
```python
def _trim_to_exact_layers(self, slab, target_nlayers, surface_index):
    """Remove átomos para ter exatamente target_nlayers camadas"""
    # Calcula layer_id para cada átomo
    # Mantém apenas os primeiros target_nlayers layers
    # Remove o resto
```

**STEP 8: Constraints**
```python
# Linhas 527-542
fix_indices = self._calculate_fixed_layers(nlayers)
fixed_atom_indices = [...]  # Identifica átomos a congelar
constraint = FixAtoms(indices=fixed_atom_indices)
slab.set_constraint(constraint)
```

---

## 🎯 Valores Típicos Durante Transformação

Exemplo: Au(111) com 6 camadas, vácuo 15 Å

```
ENTRADA:
  bulk_atoms: 4 átomos Au (célula FCC)

STEP 1-2:
  struct_conv: convencional, a=4.08 Å
  d_hkl (111): 4.08 / sqrt(3) = 2.357 Å
  min_slab_size: 6 × 2.357 - 0.5 = 13.642 Å

STEP 3:
  PyMatGen SlabGenerator gera slab com ~7 camadas
  (mais do que solicitado, porque é mínimo!)

STEP 4:
  slab (ASE): ~30 átomos

STEP 5:
  Célula antes:  [[4.08, 0.00, 0.00],
                  [2.04, 3.53, 0.00],
                  [0.00, 0.00, 25.0]]
  
  Célula depois: [[4.08, 0.00, 0.00],
                  [0.00, 4.08, 0.00],
                  [0.00, 0.00, 25.0]]  ← Ortogonal!

STEP 6:
  z_positions antes: [0.0 até 10.0]
  z_positions depois: [5.0 até 15.0] (centrado, com vácuo nas pontas)

STEP 7:
  Aparar de ~7 camadas → exatamente 6 camadas
  Átomos: ~30 → ~24

STEP 8:
  FixAtoms indices: [0, 1, 2] (primeiros 3 átomos quando sorteados por z)
  Átomos livres: ~21 (para relaxar na Phase 4)
```

---

## 🔧 Onde Você Pode Fazer Ajustes Avançados

### Opção A: Modificar `_orthogonalize_cell()`
Se a célula ortogonal não está ficando boa:
```python
def _orthogonalize_cell(self, slab):
    # Você pode ajustar a rotação aqui
    # Ex: forçar um eixo específico para z
    pass
```

### Opção B: Modificar `_trim_to_exact_layers()`
Se o aparamento de camadas não está certo:
```python
def _trim_to_exact_layers(self, slab, target_nlayers, surface_index):
    # Você pode mudar como as camadas são identificadas
    # Ex: usar método diferente que não seja baseado em z
    pass
```

### Opção C: Modificar `_calculate_fixed_layers()`
Se o padrão de congelamento não é ideal:
```python
def _calculate_fixed_layers(self, nlayers):
    # Padrão: congelar ~1/3 das camadas
    # Você pode mudar para: congelar ~1/2, ou usar outra fórmula
    return [0, 1]  # Ex: sempre congelar só 2 camadas
```

### Opção D: Modificar SlabGenerator parameters
Se quer usar min_slab_size diferente:
```python
# Na linha ~483, você pode trocar:
min_slab_size_dynamic = nlayers * d_hkl - 0.5
# Para:
min_slab_size_dynamic = nlayers * d_hkl * 1.2  # 20% mais espesso
# Ou:
min_slab_size_dynamic = nlayers * d_hkl  # Exatamente nlayers
```

### Opção E: Adicionar pós-processamento customizado
```python
def _regenerate_slab_with_nlayers(self, ...):
    # ... código original até slab.set_constraint(constraint) ...
    
    # AQUI você pode adicionar customizações:
    
    # Ex 1: Adicionar adsorbato
    # from ase.build import add_vacuum
    # adsorbato = Atoms('O')
    # adsorbato.position = [slab.cell[0,0]/2, slab.cell[1,1]/2, slab.cell[2,2] + 2.0]
    # slab += adsorbato
    
    # Ex 2: Distorcer superfície deliberadamente
    # slab.positions[-1] += [0.1, 0.1, 0.5]  # Mover último átomo
    
    # Ex 3: Aplicar strain
    # slab.set_cell(slab.get_cell() * [1.05, 1.05, 1.0])
    
    return slab
```

---

## 📊 Comparação: Primitiva vs Supercela

Ambas passam pelo mesmo fluxo, mas com tamanhos diferentes:

```
PRIMITIVA (use_primitive_cell=True)
├─ SlabGenerator → 1×1 estrutura de repetição
├─ ~6-8 átomos por camada
├─ Rápido (~5-10 min SCF)
└─ Menos estável (efeitos de tamanho finito)

SUPERCELA (use_primitive_cell=False)
├─ SlabGenerator → 2×2 estrutura de repetição
├─ ~24-32 átomos por camada
├─ Lento (~20-50 min SCF)
└─ Mais estável (menos efeitos de tamanho finito)
```

---

## 🔄 Pipeline de Convergência

Toda vez que `run_slab_convergence()` chama `_test_vacuum_convergence()`:

```
Para cada vácuo em vacuum_test=[5, 10, 15, 20, 25]:
  │
  ├─ Chamar _regenerate_slab_with_nlayers(nlayers=6)
  │  └─ Passa por todo fluxo acima
  │
  ├─ Copiar slab
  │
  ├─ Aplicar novo vácuo: slab.center(vacuum=vácuo, axis=2)
  │
  ├─ Submeter para QE (batch_utils.submit_structures_parallel)
  │  ├─ Gerar input QE
  │  ├─ Submeter job remoto
  │  └─ Esperar convergência
  │
  └─ Coletar energia E(vácuo)

Final: energias[5] = ..., energies[10] = ..., etc.
```

---

## 📝 Valores de Saída Finais

Após todo o processamento, a slab está pronta para QE com:

```python
slab.get_cell()
# Array de 3×3, ortogonal, z=vácuo

slab.get_positions()
# Posições de átomos, todas posicionadas em z entre [vácuo/2 - espessura/2, vácuo/2 + espessura/2]

slab.get_constraint()
# FixAtoms constraint aplicado

slab.pbc
# [True, True, True] - condições periódicas em 3D
```

---

## 🚨 Problemas Comuns e Diagnóstico

| Problema | Causa Provável | Solução |
|----------|----------------|---------|
| Slab muito espessa | SlabGenerator gerou >nlayers | `_trim_to_exact_layers` pode ter bug |
| Slab muito fina | d_hkl calculado errado | Verificar `a` do bulk |
| Célula não ortogonal | `_orthogonalize_cell()` falhou | Checar rotação na função |
| Átomos não centrados | `slab.center()` não funcionou | Checar PBC antes de center |
| Constraint errado | `_calculate_fixed_layers()` errado | Ajustar índices de congelamento |

---

## 📚 Referências no Código

| Função | Arquivo | Linhas | Responsabilidade |
|--------|---------|--------|-----------------|
| `run_slab_convergence()` | slab_workflow.py | 1186-1385 | Orquestrador principal |
| `_regenerate_slab_with_nlayers()` | slab_workflow.py | 448-556 | Gera slab (STEP 1-8) |
| `_test_vacuum_convergence()` | slab_workflow.py | 1388-1568 | Loop de vácuos |
| `_orthogonalize_cell()` | slab_workflow.py | ~600 | Ortogonaliza célula |
| `_trim_to_exact_layers()` | slab_workflow.py | ~670 | Aparar para N camadas |
| `_calculate_fixed_layers()` | slab_workflow.py | ~730 | Quais camadas congelar |
| `submit_structures_parallel()` | batch_utils.py | ? | Submissão paralela QE |
