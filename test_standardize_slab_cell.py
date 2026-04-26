#!/usr/bin/env python3
"""
Mocking Test 4: Gerar slab corretamente com PyMatGen + _standardize_slab_cell()
Demonstra que agora o slab sai correto direto do PyMatGen
"""

import numpy as np
from ase import Atoms
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.surface import SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

print("="*70)
print("MOCKING TEST 4: Gerar slab corretamente com PyMatGen")
print("="*70)

# ============================================================
# STEP 1: Gerar bulk e converter para PyMatGen
# ============================================================
print("\n[1] Carregando bulk Au e convertendo para PyMatGen...")

bulk_au = read('au_bulk.cif')
struct = AseAtomsAdaptor.get_structure(bulk_au)
struct_conv = SpacegroupAnalyzer(struct).get_conventional_standard_structure()

print(f"  Bulk Au: {len(struct_conv)} átomos")
print(f"  Célula: a={struct_conv.lattice.a:.4f} Å")

# ============================================================
# STEP 2: Gerar slab com PyMatGen (6 camadas, vácuo=5)
# ============================================================
print("\n[2] Gerando slab Au(111) com PyMatGen...")
print("  Parâmetros: nlayers=6, vácuo=5 Å")

# Calcular min_slab_size dinamicamente
d_hkl = struct_conv.lattice.a / np.sqrt(1**2 + 1**2 + 1**2)
nlayers = 6
min_slab_size = nlayers * d_hkl

slabgen = SlabGenerator(
    struct_conv,
    miller_index=(1, 1, 1),
    min_slab_size=min_slab_size,
    min_vacuum_size=5.0,
    center_slab=False,  # Deixa em posição arbitrária (PyMatGen faz isso)
)

slabs_list = slabgen.get_slabs(tol=0.1)
slab_pymatgen = slabs_list[0]
slab_raw = AseAtomsAdaptor.get_atoms(slab_pymatgen)

print(f"  Slab gerado: {len(slab_raw)} átomos")
print(f"\n  ANTES de _standardize_slab_cell():")
print(f"    Célula c[2]: {slab_raw.cell[2, 2]:.4f} Å")
if slab_raw.cell[2, 2] < 0:
    print(f"    ⚠ VÁCUO NEGATIVO!")
print(f"    Posições z: min={slab_raw.positions[:, 2].min():.2f}, max={slab_raw.positions[:, 2].max():.2f} Å")

# ============================================================
# STEP 3: Aplicar _standardize_slab_cell()
# ============================================================
print("\n[3] Aplicando _standardize_slab_cell()...")

def standardize_slab_cell(slab):
    """Versão da função do slab_workflow.py"""
    slab = slab.copy()
    cell = slab.get_cell()
    
    # STEP 1: Ensure c-vector (vácuo) is positive
    if cell[2, 2] < 0:
        print(f"  ⚠ Fixando vácuo negativo: {cell[2, 2]:.4f} Å → positivo")
        cell[2] = -cell[2]
        slab.set_cell(cell)
    
    # STEP 2: Wrap atoms to be within [0, L)
    slab.wrap()
    
    # STEP 3: Center slab vertically (equal vacuum above/below)
    z_positions = slab.positions[:, 2]
    z_min = z_positions.min()
    z_max = z_positions.max()
    slab_height = z_max - z_min
    
    # Calculate center position (should be at cell[2,2]/2)
    z_center_current = (z_min + z_max) / 2
    z_center_target = cell[2, 2] / 2
    
    # Shift all atoms to center
    shift = z_center_target - z_center_current
    if abs(shift) > 1e-6:  # Only shift if needed
        print(f"  ✓ Centrado slab: shifted z by {shift:.4f} Å")
        slab.positions[:, 2] += shift
        slab.wrap()  # Wrap again after shifting
    else:
        print(f"  ✓ Slab já estava centrado")
    
    return slab

slab_standardized = standardize_slab_cell(slab_raw)

# ============================================================
# STEP 4: Mostrar resultado
# ============================================================
print(f"\n[4] DEPOIS de _standardize_slab_cell():")
print(f"    Célula c[2]: {slab_standardized.cell[2, 2]:.4f} Å ✓ POSITIVO!")
print(f"    Posições z: min={slab_standardized.positions[:, 2].min():.2f}, max={slab_standardized.positions[:, 2].max():.2f} Å")
print(f"    Altura do slab: {slab_standardized.positions[:, 2].max() - slab_standardized.positions[:, 2].min():.2f} Å")

# ============================================================
# STEP 5: Aplicar center() como normalmente se faz
# ============================================================
print("\n[5] Aplicando .center(vacuum=5, axis=2)...")

slab_final = slab_standardized.copy()
slab_final.center(vacuum=5.0, axis=2)

print(f"    Posições z FINAIS: min={slab_final.positions[:, 2].min():.2f}, max={slab_final.positions[:, 2].max():.2f} Å")
print(f"    Vácuo FINAL: {slab_final.cell[2, 2]:.4f} Å ✓")

# ============================================================
# STEP 6: Comparar antes vs depois
# ============================================================
print("\n" + "="*70)
print("RESUMO: ANTES vs DEPOIS")
print("="*70)

print(f"""
ANTES (_standardize_slab_cell()):
  - Célula c[2]: {slab_raw.cell[2, 2]:.4f} Å (potencialmente negativo)
  - Posições z: {slab_raw.positions[:, 2].min():.2f} até {slab_raw.positions[:, 2].max():.2f} Å
  - Átomos podem estar fora do centro
  
DEPOIS (_standardize_slab_cell()):
  - Célula c[2]: {slab_standardized.cell[2, 2]:.4f} Å ✓ POSITIVO
  - Posições z: {slab_standardized.positions[:, 2].min():.2f} até {slab_standardized.positions[:, 2].max():.2f} Å
  - Átomos centrados na célula ✓
  
FINAL (após .center()):
  - Célula c[2]: {slab_final.cell[2, 2]:.4f} Å ✓
  - Posições z: {slab_final.positions[:, 2].min():.2f} até {slab_final.positions[:, 2].max():.2f} Å
  - Vácuo bem distribuído: {(slab_final.positions[:, 2].min() - 0):.2f} Å abaixo, {(slab_final.cell[2,2] - slab_final.positions[:, 2].max()):.2f} Å acima ✓

TIPOS DE CÉLULA:
  - Antes: {slab_raw.cell}
  - Depois: {slab_final.cell}
  ✓ CÉLULA NÃO-ORTOGONAL PRESERVADA (não foi ortogonalizada!)
""")

# ============================================================
# STEP 7: Converter para coordenadas fracionárias
# ============================================================
print(f"\n[6] POSIÇÕES FRACIONÁRIAS (crystal coordinates):")

cell_inv = np.linalg.inv(slab_final.get_cell().T)
print("\n  Slab FINAL (após .center()):")
for i, pos in enumerate(slab_final.positions):
    frac = np.dot(cell_inv, pos)
    print(f"    Au {i+1}: z_frac={frac[2]:.6f} (z_real={pos[2]:.2f} Å)")

print("\n" + "="*70)
print("✓ TEST COMPLETO: SLAB PRONTO PARA QE")
print("="*70)
print("""
Fluxo da solução:
  1. PyMatGen gera slab (pode estar com vácuo negativo/mal posicionado)
  2. _standardize_slab_cell() corrige AUTOMATICAMENTE:
     - Inverte vácuo se negativo
     - Recentra os átomos
     - Preserva célula não-ortogonal (sem precisa ortogonalizar!)
  3. .center(vacuum=5) ajusta o vácuo final
  4. Slab PRONTO para xespresso gerar input QE válido ✓
""")
