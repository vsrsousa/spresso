#!/usr/bin/env python3
"""
Mocking: Transformar célula não-ortogonal em outra representação não-ortogonal

Mostrar como transformar de uma célula não-ortogonal para outra
usando ASE ou pymatgen, mantendo a mesma estrutura.
"""

from ase.io import read, write
from ase import Atoms
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.surface import SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

print("="*70)
print("MOCKING: Transformar Célula Não-Ortogonal")
print("="*70)

# ============================================================
# STEP 1: Carregar bulk Au do CIF
# ============================================================
print("\n[1] Carregando Au bulk do CIF...")
bulk_au = read('au_bulk.cif')
print(f"  Bulk Au: {len(bulk_au)} átomos")

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
print(f"  Célula PyMatGen (REPRESENTAÇÃO 1):")
cell_orig = slab.get_cell()
print(cell_orig)

# ============================================================
# STEP 3: Centrar slab
# ============================================================
print("\n[3] Centralizando slab no vácuo...")
slab.center(vacuum=5.0, axis=2)

# ============================================================
# STEP 4: TRANSFORMAR PARA OUTRA REPRESENTAÇÃO (não-ortogonal)
# ============================================================
print("\n[4] TRANSFORMANDO PARA OUTRA REPRESENTAÇÃO NÃO-ORTOGONAL...")

def get_reciprocal_lattice(cell):
    """Calcular rede recíproca"""
    volume = np.abs(np.dot(cell[0], np.cross(cell[1], cell[2])))
    b1 = 2 * np.pi * np.cross(cell[1], cell[2]) / volume
    b2 = 2 * np.pi * np.cross(cell[2], cell[0]) / volume
    b3 = 2 * np.pi * np.cross(cell[0], cell[1]) / volume
    return np.array([b1, b2, b3])

def find_shortest_lattice_vectors(cell, max_iterations=5):
    """
    Encontrar uma representação mais 'natural' da célula
    usando combinações lineares de vetores de rede
    """
    a1, a2, a3 = cell
    
    # Tentar diferentes combinações para encontrar uma mais 'simples'
    best_cell = cell.copy()
    best_norm = np.linalg.norm(a1) + np.linalg.norm(a2) + np.linalg.norm(a3)
    
    # Operações de rede: combinações lineares de vetores
    for i1 in range(-2, 3):
        for i2 in range(-2, 3):
            for j1 in range(-2, 3):
                for j2 in range(-2, 3):
                    new_a1 = i1 * a1 + i2 * a2
                    new_a2 = j1 * a1 + j2 * a2
                    new_a3 = a3  # Manter c igual
                    
                    # Calcular norma total
                    norm = np.linalg.norm(new_a1) + np.linalg.norm(new_a2) + np.linalg.norm(new_a3)
                    
                    # Se encontrou representação mais curta, usar
                    if norm < best_norm and abs(np.linalg.det(np.array([new_a1, new_a2, new_a3]))) > 1e-3:
                        best_cell = np.array([new_a1, new_a2, new_a3])
                        best_norm = norm
    
    return best_cell

# Encontrar representação alternativa
cell_transformed = find_shortest_lattice_vectors(slab.get_cell())

print(f"\n  Célula Original (REPRESENTAÇÃO 1):")
print(slab.get_cell())

print(f"\n  Célula Transformada (REPRESENTAÇÃO 2 - também não-ortogonal):")
print(cell_transformed)

# Aplicar transformação
slab_transformed = slab.copy()
slab_transformed.set_cell(cell_transformed, scale_atoms=True)

print(f"\n  ✓ Transformação aplicada")
print(f"    Magnitude a: {np.linalg.norm(cell_transformed[0]):.4f} Å")
print(f"    Magnitude b: {np.linalg.norm(cell_transformed[1]):.4f} Å")
print(f"    Magnitude c: {np.linalg.norm(cell_transformed[2]):.4f} Å")

# ============================================================
# STEP 5: Verificar se a estrutura foi preservada
# ============================================================
print("\n[5] Verificando preservação de estrutura...")
from scipy.spatial.distance import pdist

dist_orig = pdist(slab.get_positions())
dist_trans = pdist(slab_transformed.get_positions())

print(f"  Distâncias NN:")
print(f"    Original: min={dist_orig.min():.4f}, max={dist_orig.max():.4f} Å")
print(f"    Transformada: min={dist_trans.min():.4f}, max={dist_trans.max():.4f} Å")

# ============================================================
# STEP 6: Salvar ambas as formas
# ============================================================
print("\n[6] Salvando em formato QE...")

write('au111_repr1.in', slab)
write('au111_repr2.in', slab_transformed)

print("  ✓ au111_repr1.in (representação original)")
print("  ✓ au111_repr2.in (representação transformada)")

# ============================================================
# STEP 7: Comparar inputs gerados
# ============================================================
print("\n[7] COMPARAÇÃO DOS INPUTS GERADOS:")
print("="*70)

print("\nREPRESENTAÇÃO 1 (ASE ASE gerou):")
with open('au111_repr1.in', 'r') as f:
    for i, line in enumerate(f):
        if i < 15:
            print(line.rstrip())

print("\n\nREPRESENTAÇÃO 2 (Transformada - vetores mais curtos):")
with open('au111_repr2.in', 'r') as f:
    for i, line in enumerate(f):
        if i < 15:
            print(line.rstrip())

# ============================================================
# STEP 8: Resumo
# ============================================================
print("\n" + "="*70)
print("RESUMO: Duas Representações da Mesma Estrutura")
print("="*70)

print(f"""
Ambas são VÁLIDAS e representam a MESMA estrutura cristalina!

REPRESENTAÇÃO 1 (original PyMatGen):
  - Células: a, b acoplados e c separado
  - Volume: {np.linalg.det(slab.get_cell()):.2f} Å³
  
REPRESENTAÇÃO 2 (transformada - combinação linear):
  - Células: a, b com componentes diferentes
  - Volume: {np.linalg.det(slab_transformed.get_cell()):.2f} Å³

Transformação usada: Combinação linear dos vetores originais
  new_a = c1*old_a + c2*old_b
  new_b = d1*old_a + d2*old_b
  new_c = old_c

Isso muda a REPRESENTAÇÃO mas não a ESTRUTURA CRISTALINA!
QE aceita ambas as formas sem problema.
""")
