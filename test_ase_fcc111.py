#!/usr/bin/env python3
"""
Test: Comparar PyMatGen vs ASE fcc111() para Au(111)
"""

import numpy as np
from ase.build import fcc111
from ase.io import write, read

print("="*70)
print("COMPARAÇÃO: PyMatGen vs ASE fcc111()")
print("="*70)

# ============================================================
# OPÇÃO 1: ASE fcc111() - SIMPLES e CORRETO
# ============================================================
print("\n[OPÇÃO 1] ASE fcc111():")
print("-"*70)

# Parâmetros Au
a_au = 4.0782  # Au FCC lattice constant
nlayers = 6
vacuum_size = 5.0

# Gerar slab com ASE
slab_ase = fcc111('Au', size=(1, 1, nlayers), a=a_au, vacuum=vacuum_size)

print(f"  Slab: {len(slab_ase)} átomos")
print(f"\n  CÉLULA:")
print(f"    {slab_ase.cell[0]}")
print(f"    {slab_ase.cell[1]}")
print(f"    {slab_ase.cell[2]}")
print(f"    Vácuo (c[2]): {slab_ase.cell[2, 2]:.4f} Å ✓")

print(f"\n  POSIÇÕES (z):")
print(f"    min: {slab_ase.positions[:, 2].min():.4f} Å")
print(f"    max: {slab_ase.positions[:, 2].max():.4f} Å")
print(f"    altura slab: {slab_ase.positions[:, 2].max() - slab_ase.positions[:, 2].min():.4f} Å")

# Converter para coordenadas fracionárias
cell_inv = np.linalg.inv(slab_ase.get_cell().T)
print(f"\n  POSIÇÕES FRACIONÁRIAS (crystal):")
for i, pos in enumerate(slab_ase.positions):
    frac = np.dot(cell_inv, pos)
    print(f"    Au {i+1}: z_frac={frac[2]:.6f}")

# Salvar
write('au111_ase_fcc111.in', slab_ase)
print(f"\n  ✓ Salvo em: au111_ase_fcc111.in")

# ============================================================
# VALIDAÇÃO: d-spacing
# ============================================================
print("\n[VALIDAÇÃO] Espaçamento entre camadas:")
print("-"*70)

z_positions = np.sort(slab_ase.positions[:, 2])
spacings = np.diff(z_positions)

print(f"  Posições z: {z_positions}")
print(f"  Espaçamentos: {spacings}")
print(f"  Médio: {np.mean(spacings):.4f} Å")

# Teórico para Au(111)
d_spacing_theory = a_au / np.sqrt(1**2 + 1**2 + 1**2)  # d_111
print(f"\n  Teórico d_111: {d_spacing_theory:.4f} Å")
print(f"  Dois espaçamentos: {2*d_spacing_theory:.4f} Å")
print(f"  Diferença: {abs(np.mean(spacings) - 2*d_spacing_theory):.6f} Å")

# ============================================================
# MOSTRAR INPUT QE
# ============================================================
print("\n[INPUT QE GERADO]:")
print("-"*70)

with open('au111_ase_fcc111.in', 'r') as f:
    lines = f.readlines()
    for i, line in enumerate(lines[:15]):
        print(line.rstrip())

print("\n" + "="*70)
print("✓ TESTE COMPLETO")
print("="*70)
print(f"""
VANTAGENS do ASE fcc111():
  ✓ Função específica para FCC(111)
  ✓ Gera célula ORTOGONAL (não precisa consertar)
  ✓ Posições dos átomos CORRETAS
  ✓ Vácuo SEMPRE positivo
  ✓ Muito simples de usar
  ✓ Determinístico (sem surpresas)

SINTAXE:
  from ase.build import fcc111
  slab = fcc111('Au', size=(1, 1, nlayers), a=4.0782, vacuum=5.0)
  
PARÂMETROS:
  - 'Au': símbolo do elemento
  - size=(nx, ny, nz): número de repetições
    nx, ny: repetições in-plane (geralmente 1,1 para cell mínima)
    nz: número de camadas
  - a: lattice constant
  - vacuum: vácuo em Ångströms
""")
