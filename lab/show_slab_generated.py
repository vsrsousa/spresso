#!/usr/bin/env python3
"""
Gera slab usando SlabWorkflow e mostra o resultado
"""
import sys
sys.path.insert(0, '/home/vinicius/projects/spresso')

from ase.io import read
from xespresso.workflow.slab_workflow import SlabWorkflow
import numpy as np

print("="*70)
print("GERANDO SLAB COM SLABWORKFLOW")
print("="*70)

# Carregar bulk Au
bulk_au = read('au_bulk.cif')
print(f"\n[BULK]")
print(f"  Átomos: {len(bulk_au)}")
print(f"  Cell: {bulk_au.cell.cellpar()}")

# Criar SlabWorkflow
slab_wf = SlabWorkflow(
    bulk_atoms=bulk_au,
    surface_indices=[(1, 1, 1)],
    nlayers=6,
    min_vacuum_size=5.0,
    pseudopotentials_config='default',
)

# Gerar slabs
slab_wf.generate_slabs()

slab = slab_wf.slabs[(1, 1, 1)]

print(f"\n[SLAB GERADO]")
print(f"  Átomos: {len(slab)}")
print(f"  Cell a,b,c: {slab.cell.cellpar()[:3]}")
print(f"  Cell angles: {slab.cell.cellpar()[3:]}")
print(f"  Volume: {slab.get_volume():.2f} Ų")
print(f"  Vácuo (c): {slab.cell[2,2]:.4f} Å")

print(f"\n[POSIÇÕES ATÔMICAS] (em Å)")
print(f"{'Átomo':<6} {'X':<12} {'Y':<12} {'Z':<12}")
print("-" * 50)
for i, pos in enumerate(slab.positions):
    print(f"Au{i+1:<4} {pos[0]:>11.6f} {pos[1]:>11.6f} {pos[2]:>11.6f}")

# Analisar camadas
z_coords = slab.positions[:, 2]
z_sorted = np.sort(z_coords)

print(f"\n[ANÁLISE DE CAMADAS]")
print(f"  Z mínimo: {z_coords.min():.4f} Å")
print(f"  Z máximo: {z_coords.max():.4f} Å")
print(f"  Altura do slab: {z_coords.max() - z_coords.min():.4f} Å")

print(f"\n  Coordenadas Z únicas:")
for i, z in enumerate(z_sorted):
    print(f"    Camada {i+1}: z = {z:.4f} Å")

print(f"\n  Espaçamento entre camadas:")
for i in range(len(z_sorted)-1):
    delta = z_sorted[i+1] - z_sorted[i]
    print(f"    Camada {i+1} → {i+2}: {delta:.4f} Å")

# Constraints
if slab.constraints:
    print(f"\n[CONSTRAINTS]")
    for constraint in slab.constraints:
        print(f"  {constraint}")
else:
    print(f"\n[CONSTRAINTS] Nenhum aplicado")

EOF
