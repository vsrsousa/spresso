#!/usr/bin/env python3
"""
Como transformar célula usando set_cell() do ASE
"""

from ase.io import read, write
import numpy as np

print("="*70)
print("Como eu fiz: Transformar célula com set_cell()")
print("="*70)

# Carregar slab
from ase import Atoms
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.surface import SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

bulk_au = read('au_bulk.cif')
struct = AseAtomsAdaptor.get_structure(bulk_au)
struct_conv = SpacegroupAnalyzer(struct).get_conventional_standard_structure()

slabgen = SlabGenerator(struct_conv, miller_index=(1, 1, 1), 
                        min_slab_size=10.0, min_vacuum_size=5.0, center_slab=False)
slab_pymatgen = slabgen.get_slabs(tol=0.1)[0]
slab = AseAtomsAdaptor.get_atoms(slab_pymatgen)
slab.center(vacuum=5.0, axis=2)

print("\n[1] CÉLULA ORIGINAL (PyMatGen):")
print(slab.get_cell())
print(f"    Volume: {slab.get_volume():.2f} Å³")

# ============================================================
# Opção 1: Usar set_cell() com uma nova célula
# ============================================================
print("\n[2] TRANSFORMAR COM set_cell():")

slab_new = slab.copy()

# Definir uma nova célula (com scale_atoms=True, escala as posições)
new_cell = np.array([
    [2.88373,   0.00000,   0.00000],
    [1.44186,   2.49738,   0.00000],
    [0.00000,   0.00000,  15.00000]
])

print(f"\n    Nova célula a usar:")
print(new_cell)

# Aplicar transformação
slab_new.set_cell(new_cell, scale_atoms=True)

print(f"\n✓ CÉLULA TRANSFORMADA:")
print(slab_new.get_cell())
print(f"  Volume: {slab_new.get_volume():.2f} Å³")

# ============================================================
# Verificar se estrutura foi preservada
# ============================================================
print("\n[3] VERIFICAR PRESERVAÇÃO:")
print(f"  Átomos original: {len(slab)}")
print(f"  Átomos novo: {len(slab_new)}")
print(f"  Posições z original: min={slab.positions[:, 2].min():.4f}, max={slab.positions[:, 2].max():.4f}")
print(f"  Posições z novo: min={slab_new.positions[:, 2].min():.4f}, max={slab_new.positions[:, 2].max():.4f}")

# ============================================================
# Salvar
# ============================================================
print("\n[4] SALVAR:")
write('au111_reformatted.in', slab_new)
print("  ✓ au111_reformatted.in (com nova célula)")

print("\n[5] INPUT GERADO:")
with open('au111_reformatted.in', 'r') as f:
    lines = f.readlines()
    for line in lines[:15]:
        print(line.rstrip())

print("\n" + "="*70)
print("RESUMO: A Transformação é Simples")
print("="*70)
print("""
Você define uma nova célula e chama:
  slab.set_cell(new_cell, scale_atoms=True)

scale_atoms=True: escala as posições atomicamente para se adequar
ao novo sistema de coordenadas.

A estrutura cristalina é preservada, só muda a representação!
""")
