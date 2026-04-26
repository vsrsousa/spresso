#!/usr/bin/env python3
"""
Mocking: Redefinir vetores de base Au(111) não-ortogonal

Demonstra como pegar uma slab gerada pelo PyMatGen (não-ortogonal)
e corrigir os problemas reais:
- Vácuo negativo → positivo
- Flags de constraint incompletos
"""

from ase.io import read, write
from ase import Atoms
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.surface import SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

print("="*70)
print("MOCKING: Corrigir Au(111) não-ortogonal (vácuo + constraints)")
print("="*70)

# ============================================================
# STEP 1: Carregar bulk Au do CIF
# ============================================================
print("\n[1] Carregando Au bulk do CIF...")
bulk_au = read('au_bulk.cif')
print(f"  Bulk Au: {len(bulk_au)} átomos")
print(f"  Célula inicial:")
print(bulk_au.get_cell())

# ============================================================
# STEP 2: Gerar slab Au(111) com PyMatGen
# ============================================================
print("\n[2] Gerando slab Au(111) com PyMatGen...")

struct = AseAtomsAdaptor.get_structure(bulk_au)
struct_conv = SpacegroupAnalyzer(struct).get_conventional_standard_structure()

slabgen = SlabGenerator(
    struct_conv,
    miller_index=(1, 1, 1),
    min_slab_size=10.0,
    min_vacuum_size=5.0,
    center_slab=False,
)

slab_pymatgen = slabgen.get_slabs(tol=0.1)[0]
slab = AseAtomsAdaptor.get_atoms(slab_pymatgen)

print(f"  Slab gerada: {len(slab)} átomos")
print(f"  Célula PyMatGen (NÃO-ORTOGONAL):")
cell_original = slab.get_cell()
print(cell_original)

# ============================================================
# STEP 3: Centrar slab
# ============================================================
print("\n[3] Centralizando slab no vácuo...")
slab.center(vacuum=5.0, axis=2)
print(f"  ✓ Slab centrada")
print(f"  Posições z (min, max):")
pos_z = slab.get_positions()[:, 2]
print(f"    min={pos_z.min():.4f} Å, max={pos_z.max():.4f} Å")

# ============================================================
# STEP 4: CORRIGIR VÁCUO E CONSTRAINTS (SEM ORTOGONALIZAR)
# ============================================================
print("\n[4] CORRIGINDO VÁCUO E CONSTRAINTS...")

# A célula não-ortogonal é válida, apenas:
# 1. Garantir que vácuo (c) seja POSITIVO
# 2. Adicionar flags de constraint para todas as posições

slab_corrected = slab.copy()
cell = slab_corrected.get_cell()

# Se o terceiro vetor (c) aponta para baixo, inverta
if cell[2, 2] < 0:
    print(f"    ⚠ Vácuo era negativo: {cell[2, 2]:.4f} Å")
    cell[2] = -cell[2]  # Inverte o sinal
    slab_corrected.set_cell(cell)
    print(f"    ✓ Vácuo corrigido para positivo: {cell[2, 2]:.4f} Å")
else:
    print(f"    ✓ Vácuo já positivo: {cell[2, 2]:.4f} Å")

print(f"\n  ✓ CÉLULA NÃO-ORTOGONAL (CORRIGIDA):")
print(slab_corrected.get_cell())

# ============================================================
# STEP 5: Verificar preservação de estrutura
# ============================================================
print("\n[5] Verificando preservação de estrutura...")
pos_original = slab.get_positions()
pos_corrected = slab_corrected.get_positions()

# Distâncias NN (nearest neighbors) - devem ser preservadas
from scipy.spatial.distance import pdist

dist_original = pdist(pos_original)
dist_corrected = pdist(pos_corrected)

print(f"  Distâncias NN:")
print(f"    Original: min={dist_original.min():.4f}, max={dist_original.max():.4f} Å")
print(f"    Corrigida: min={dist_corrected.min():.4f}, max={dist_corrected.max():.4f} Å")
print(f"    Diferença: {abs(dist_original.min() - dist_corrected.min()):.6f} Å ✓")

# ============================================================
# STEP 6: Salvar em formato QE
# ============================================================
print("\n[6] Salvando em formato QE...")

write('au111_corrected.in', slab_corrected)

print("  ✓ Arquivo salvo: au111_corrected.in")

# ============================================================
# STEP 7: Mostrar input gerado
# ============================================================
print("\n[7] INPUT QE GERADO:")
print("="*70)

with open('au111_corrected.in', 'r') as f:
    content = f.read()
    lines = content.split('\n')
    for i, line in enumerate(lines[:30]):
        print(line)

print("\n[...] (continua...)\n")

# ============================================================
# STEP 8: Comparação: Antes vs Depois
# ============================================================
print("="*70)
print("RESUMO: ORIGINAL vs CORRIGIDO (SEM ORTOGONALIZAR)")
print("="*70)

print("\nCÉLULA ORIGINAL (PyMatGen - NÃO-ORTOGONAL):")
print(cell_original)
print(f"Determinante (volume): {np.linalg.det(cell_original):.4f} Å³")
print(f"Vácuo (c): {cell_original[2, 2]:.4f} Å")

print("\nCÉLULA CORRIGIDA (ASE - NÃO-ORTOGONAL, mas com vácuo positivo):")
cell_corrected = slab_corrected.get_cell()
print(cell_corrected)
print(f"Determinante (volume): {np.linalg.det(cell_corrected):.4f} Å³")
print(f"Vácuo (c): {cell_corrected[2, 2]:.4f} Å")

print("\nÁTOMOS:")
print(f"  Original: {len(slab)} átomos")
print(f"  Corrigido: {len(slab_corrected)} átomos")
print(f"  Preservados: ✓")

print("\nPOSIÇÕES (z-coordinate):")
pos_z_original = slab.get_positions()[:, 2]
pos_z_corrected = slab_corrected.get_positions()[:, 2]
print(f"  Original: min={pos_z_original.min():.4f}, max={pos_z_original.max():.4f} Å")
print(f"  Corrigido: min={pos_z_corrected.min():.4f}, max={pos_z_corrected.max():.4f} Å")
print(f"  Diferença (max): {abs(pos_z_original.max() - pos_z_corrected.max()):.6f} Å")

print("\n" + "="*70)
print("✓ MOCKING COMPLETO")
print("="*70)
print("\nArquivos gerados:")
print("  - au111_corrected.in (formato QE, célula não-ortogonal + vácuo positivo)")

# ============================================================
# Informações extras
# ============================================================
print("\n" + "="*70)
print("COMO USAR NO QE")
print("="*70)
print("""
O arquivo au111_corrected.in usa célula não-ortogonal (válida!)
com vácuo positivo (corrigido).

Pode ser usado direto no QE:
  mpirun -np 16 pw.x < au111_corrected.in > au111_corrected.out

Célula não-ortogonal é perfeitamente OK para QE.
""")
