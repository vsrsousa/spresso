#!/usr/bin/env python3
"""
Demonstração: Opções para gerar slabs com xespresso agora.

Você tem 3 caminhos:
1. SlabWorkflow (pymatgen) - completo mas complexo
2. ASE nativo - simples, rápido, sem dependências
3. Híbrido - ASE para slab + xespresso para cálculos
"""

from ase.build import bulk
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def example_1_slabworkflow_original():
    """Forma original com SlabWorkflow."""
    print("\n" + "="*70)
    print("OPÇÃO 1: SlabWorkflow (pymatgen)")
    print("="*70)
    print("""
from xespresso.workflow.slab_workflow import SlabWorkflow

slab_wf = SlabWorkflow(
    bulk_atoms=relaxed_atoms,
    surface_indices=[(1, 1, 1)],
    min_vacuum_size=5.0,
    nlayers=6,
    pseudopotentials_config='default',
    protocol='standard',
    precision='low',
    verbose=True,
    use_primitive_cell=True
)

# Gera slabs com pymatgen
slabs = slab_wf.generate_slabs()
slab_111 = slabs[(1, 1, 1)]

print(f"Slab: {len(slab_111)} atoms")
    """)
    print("✓ Completo (converge bulk, gera slab)")
    print("✗ Depende de pymatgen (maiores que ASE)")
    print("✗ Pode gerar células rotacionadas/confusas")


def example_2_ase_native():
    """Forma simplificada com ASE puro."""
    print("\n" + "="*70)
    print("OPÇÃO 2: ASE Nativo (Recomendado - NOVO)")
    print("="*70)
    print("""
from ase.build import fcc111, add_vacuum
from ase.io import read

# Opção A: Carregar CIF
bulk = read('au_bulk.cif')
a = bulk.cell[0, 0]

# Opção B: Usar bulk do ASE direto
from ase.build import bulk
bulk = bulk('Au', 'fcc', a=4.0782)

# Gerar slab
slab = fcc111('Au', size=(1, 1, 6), a=a, vacuum=0)
add_vacuum(slab, 5.0)

print(f"Slab: {len(slab)} atoms")
    """)
    print("✓ Simples e rápido")
    print("✓ Sem dependências extras (só ASE)")
    print("✓ Células sempre bem-orientadas")
    print("✗ Não converge bulk automaticamente")


def example_3_hybrid():
    """Forma híbrida: ASE para slab + xespresso para cálculos."""
    print("\n" + "="*70)
    print("OPÇÃO 3: Híbrida (ASE + xespresso)")
    print("="*70)
    print("""
from ase.build import fcc111, add_vacuum
from ase.io import read
from ase.constraints import FixAtoms
from xespresso.workflow.calculation_workflow import CalculationWorkflow

# 1. Gerar slab com ASE (rápido, sem pymatgen)
slab = fcc111('Au', size=(1, 1, 6), a=4.0782, vacuum=0)
add_vacuum(slab, 5.0)

# 2. Fixar camadas de baixo
z_pos = slab.get_positions()[:, 2]
z_min = z_pos.min()
z_range = z_pos.max() - z_min
fix_threshold = z_min + z_range * 0.4
indices = [i for i, z in enumerate(z_pos) if z < fix_threshold]
slab.set_constraint(FixAtoms(indices=indices))

# 3. Usar xespresso para relaxação
calc_wf = CalculationWorkflow(
    slab,
    pseudopotentials_config='default',
    machine='medusa',
)

# Executar relaxação
results = calc_wf.run_calculation(
    label='au111_relax',
    fmax=0.05,  # Force convergence
)

print(f"Final energy: {results['energy']} eV")
print(f"Final forces: {results['forces']} eV/Å")
    """)
    print("✓ Melhor dos dois mundos")
    print("✓ ASE simples para geometria")
    print("✓ xespresso para cálculos QE")
    print("✓ Sem dependências desnecessárias")


def example_4_comparison_table():
    """Tabela comparativa."""
    print("\n" + "="*70)
    print("COMPARAÇÃO")
    print("="*70)
    
    comparison = """
┌─────────────────┬──────────────────┬──────────────┬─────────────┐
│ Aspecto         │ SlabWorkflow     │ ASE Nativo   │ Híbrida     │
├─────────────────┼──────────────────┼──────────────┼─────────────┤
│ Complexidade    │ Alta             │ Baixa        │ Média       │
│ Dependências    │ pymatgen         │ Só ASE       │ ASE + xespresso│
│ Célula correta  │ Às vezes         │ Sempre ✓     │ Sempre ✓    │
│ Velocidade      │ Lenta            │ Rápida ✓     │ Rápida      │
│ Cálc. bulk      │ Sim ✓            │ Não          │ Opcional    │
│ Cálc. slab      │ Não              │ Não          │ Sim ✓       │
│ Recomendado     │ Pesquisa+debug   │ Testes/dev   │ Produção    │
└─────────────────┴──────────────────┴──────────────┴─────────────┘
    """
    print(comparison)


def main():
    print("="*70)
    print("OPÇÕES DE GERAÇÃO DE SLABS COM XESPRESSO AGORA")
    print("="*70)
    
    example_1_slabworkflow_original()
    example_2_ase_native()
    example_3_hybrid()
    example_4_comparison_table()
    
    # Recomendação final
    print("\n" + "="*70)
    print("RECOMENDAÇÃO")
    print("="*70)
    print("""
Para seu caso (Au(111), relaxação):

  MELHOR: Opção 3 (Híbrida)
  
  slab = fcc111('Au', size=(1,1,6), a=4.0782, vacuum=0)
  add_vacuum(slab, 5.0)
  
  # Fixar camadas
  z_pos = slab.get_positions()[:, 2]
  fix_idx = [i for i,z in enumerate(z_pos) if z < z_pos.min() + (z_pos.max()-z_pos.min())*0.4]
  slab.set_constraint(FixAtoms(indices=fix_idx))
  
  # Relaxar com xespresso
  calc_wf = CalculationWorkflow(slab, machine='medusa')
  calc_wf.run_calculation(label='au111_relax', fmax=0.05)

Vantagens:
  ✓ Simples e direto
  ✓ Sem pymatgen (célula garantida correta)
  ✓ Integração xespresso completa
  ✓ Rápido (ASE para slab)
  ✓ Preciso (xespresso para cálculos)
    """)
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
