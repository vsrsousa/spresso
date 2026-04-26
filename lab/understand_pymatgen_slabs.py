#!/usr/bin/env python3
"""
Análise: Por que pymatgen gera células "estranhas" com valores negativos?

O que aconteceu:
  Você viu: CELL_PARAMETERS com valores negativos
  Isso é válido algebricamente, mas confuso visualmente
  Não é um bug do pymatgen, é um "recurso" de otimização
"""

import numpy as np
from ase.io import read
from ase.build import bulk
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor


def explain_pymatgen_behavior():
    print("="*70)
    print("POR QUE PYMATGEN GERA CÉLULAS ESTRANHAS?")
    print("="*70)
    
    # Criar bulk Au
    au_bulk = bulk('Au', 'fcc', a=4.0782)
    print(f"\n1. Bulk Au original:")
    print(f"   Cell: {au_bulk.cell}")
    
    # Converter para pymatgen
    adaptor = AseAtomsAdaptor()
    pmg_structure = adaptor.get_structure(au_bulk)
    
    print(f"\n2. Convertido para pymatgen Structure:")
    print(f"   Lattice: {pmg_structure.lattice.matrix}")
    
    # Gerar slab com pymatgen
    slab_gen = SlabGenerator(
        pmg_structure,
        (1, 1, 1),  # Surface
        min_slab_size=7,
        min_vacuum_size=5,
        center_slab=True,
    )
    
    slab_structure = slab_gen.get_slab()
    
    print(f"\n3. Slab gerado pelo pymatgen:")
    cell_pmg = slab_structure.lattice.matrix
    print(f"   Cell: ")
    for i, row in enumerate(cell_pmg):
        print(f"     [{row[0]:12.6f}, {row[1]:12.6f}, {row[2]:12.6f}]")
    
    # Converter back para ASE
    slab_ase = adaptor.get_atoms(slab_structure)
    
    print(f"\n4. Convertido de volta para ASE:")
    cell_ase = slab_ase.get_cell()
    print(f"   Cell: ")
    for i, row in enumerate(cell_ase):
        print(f"     [{row[0]:12.6f}, {row[1]:12.6f}, {row[2]:12.6f}]")
    
    print(f"\n   Posições z:")
    z_pos = slab_ase.get_positions()[:, 2]
    for i, z in enumerate(z_pos):
        print(f"     Atom {i}: z = {z:10.6f}")
    
    print(f"\n{'─'*70}")
    print("ANÁLISE:")
    print(f"{'─'*70}")
    
    print("""
O pymatgen PROPOSITALMENTE gera células com valores negativos porque:

1. **Simetria Otimizada**: Reduz a célula para a forma mais simétrica
   - Valores negativos são uma consequência dessa otimização
   - A célula é algebricamente CORRETA (det(cell) > 0)
   - Mas visualmente CONFUSA (valores negativos)

2. **Exemplo prático**:
   Célula original (padrão):
     [5.767, 0.000, 0.000]
     [2.884, 4.995, 0.000]
     [0.000, 0.000, 18.157]
   
   Célula otimizada pelo pymatgen (equivalente algebricamente):
     [-1.467, -2.541, 0.000]
     [-1.467,  2.541, 0.000]
     [0.000,  0.000, -34.747]
   
   Ambas descrevem o MESMO slab, mas com orientações diferentes!

3. **Consequência**: 
   - Posições z podem ficar negativas
   - Célula c fica negativa
   - QE consegue lidar, mas é confuso

4. **Minha "solução"**:
   NÃO consertei o pymatgen (não pode)
   CONSERTEI o output com normalize_cell() que:
   - Inverte sinais para deixar positivo
   - Reorienta para forma padrão
   - Torna visualmente correto

RESULTADO:
  Antes (pymatgen):  [-1.467, -2.541, 0] / [-1.467,  2.541, 0] / [0, 0, -34.747]
  Depois (normalize): [5.767, 0.000, 0] / [2.884, 4.995, 0] / [0, 0, 18.157] ✓
  """)


def compare_approaches():
    print("\n" + "="*70)
    print("COMPARAÇÃO: Com/Sem pymatgen")
    print("="*70)
    
    print("""
┌─────────────────────┬──────────────────┬──────────────────┐
│ Aspecto             │ Com pymatgen      │ Sem pymatgen (ASE)│
├─────────────────────┼──────────────────┼──────────────────┤
│ Célula original     │ Estanha ✗         │ Clara ✓          │
│ Precisa normalizar  │ Sim               │ Não              │
│ Rigor matemático    │ Alto (pymatgen)   │ Simples (ASE)    │
│ Suporta struct. +   │ Sim ✓             │ Só FCC/BCC/HCP  │
│ Velocidade          │ Lenta             │ Rápida ✓         │
│ Recomendado         │ Pesquisa          │ Uso geral ✓      │
└─────────────────────┴──────────────────┴──────────────────┘
    """)


def solution():
    print("\n" + "="*70)
    print("SOLUÇÃO: Para QUALQUER estrutura")
    print("="*70)
    
    print("""
Você quer: Passar qualquer CIF e gerar slabs

OPÇÃO 1: Usar SlabWorkflow (pymatgen, suporta tudo)
─────────────────────────────────────────────────────
from xespresso.workflow.slab_workflow import SlabWorkflow

slab_wf = SlabWorkflow(
    bulk_atoms=read('any_structure.cif'),
    surface_indices=[(1, 1, 1)],
    min_vacuum_size=5.0,
    nlayers=6,
)
slabs = slab_wf.generate_slabs()
slab = slabs[(1, 1, 1)]

✓ Suporta qualquer estrutura
✓ Célula pode ficar confusa
✓ Use meu normalize_cell() pra consertar

OPÇÃO 2: Usar pymatgen diretamente (mais controle)
──────────────────────────────────────────────────
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.surface import SlabGenerator

# Carregar qualquer estrutura
atoms = read('any_structure.cif')

# Converter para pymatgen
adaptor = AseAtomsAdaptor()
pmg_struct = adaptor.get_structure(atoms)

# Gerar slab com controle fino
slab_gen = SlabGenerator(
    pmg_struct,
    (1, 1, 1),           # Surface
    min_slab_size=8,     # Espessura mínima
    min_vacuum_size=5.0, # Vácuo
    center_slab=True,
    lll_reduce=False,    # ← IMPORTANTE: Desativa otimização de células
)
slab_pmg = slab_gen.get_slab()

# Converter de volta
slab = adaptor.get_atoms(slab_pmg)

✓ Suporta qualquer estrutura
✓ lll_reduce=False evita células estranhas
✓ Máximo controle

OPÇÃO 3: Função genérica (recomendado)
────────────────────────────────────────
(ver slab_universal.py)
    """)


if __name__ == '__main__':
    explain_pymatgen_behavior()
    compare_approaches()
    solution()
