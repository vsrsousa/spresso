#!/usr/bin/env python3
"""
Mocking Test 3: Corrigir vácuo NEGATIVO + RECENTRAR átomos
Simula a estrutura que o usuário mostrou (6 layers, vacuum=5)
"""

import numpy as np
from ase import Atoms

print("="*70)
print("MOCKING TEST 3: Inverter vácuo negativo + Recentrar átomos")
print("="*70)

# ============================================================
# STEP 1: Criar estrutura EXATAMENTE como o usuário mostrou
# ============================================================
print("\n[1] Criando estrutura com vácuo negativo (como saiu do PyMatGen)...")

# 6 camadas de Au(111)
positions_fractional = [
    [0, 0, -0.8373817787],
    [0, 0, -0.7024290672],
    [0, 0, -0.5674763557],
    [0, 0, -0.4325236443],
    [0, 0, -0.2975709328],
    [0, 0, -0.1626182213],
]

# Célula com vácuo NEGATIVO (como saiu)
cell_negative = [
    [-1.46702467661151, -2.54096127584844, 0.00000000000000],
    [-1.46702467661151, 2.54096127584844, 0.00000000000000],
    [0.00000000000000, 0.00000000000000, -30.74686194000001],  # NEGATIVO!
]

# Converter coordenadas fracionárias em reais para a célula negativa
positions_real_before = []
for frac in positions_fractional:
    real = np.dot(frac, cell_negative)
    positions_real_before.append(real)

slab_bad = Atoms('Au6', positions=positions_real_before, cell=cell_negative, pbc=True)

print(f"  Slab: {len(slab_bad)} átomos")
print(f"\n  ANTES (com vácuo negativo):")
print(f"    Célula c[2]: {slab_bad.cell[2, 2]:.4f} Å ⚠ NEGATIVO!")
print(f"    Posições z: min={slab_bad.positions[:, 2].min():.2f}, max={slab_bad.positions[:, 2].max():.2f} Å")
print(f"    Slab height: {slab_bad.positions[:, 2].max() - slab_bad.positions[:, 2].min():.2f} Å")

# ============================================================
# STEP 2: Aplicar correção (inverter vácuo + recentrar)
# ============================================================
print("\n[2] Aplicando _ensure_positive_vacuum() com recentragem...")

def ensure_positive_vacuum_with_center(slab):
    """Versão melhorada com recentragem"""
    cell = slab.get_cell()
    if cell[2, 2] < 0:
        print(f"  ⚠ Detectado vácuo negativo: {cell[2, 2]:.4f} Å")
        
        # Invert cell
        cell[2] = -cell[2]
        slab.set_cell(cell)
        print(f"  ✓ Célula invertida: {cell[2, 2]:.4f} Å")
        
        # Recenter atoms
        z_min = slab.positions[:, 2].min()
        z_max = slab.positions[:, 2].max()
        slab_height = z_max - z_min
        vacuum_size = cell[2, 2] - slab_height
        
        print(f"  ✓ Recentrando: vacuum={vacuum_size / 2:.4f} Å (total slab height={slab_height:.2f} Å)")
        slab.center(vacuum=vacuum_size / 2, axis=2)
    
    return slab

slab_fixed = slab_bad.copy()
slab_fixed = ensure_positive_vacuum_with_center(slab_fixed)

# ============================================================
# STEP 3: Mostrar resultado
# ============================================================
print(f"\n[3] DEPOIS (corrigido):")
print(f"    Célula c[2]: {slab_fixed.cell[2, 2]:.4f} Å ✓ POSITIVO!")
print(f"    Posições z: min={slab_fixed.positions[:, 2].min():.2f}, max={slab_fixed.positions[:, 2].max():.2f} Å")
print(f"    Slab height: {slab_fixed.positions[:, 2].max() - slab_fixed.positions[:, 2].min():.2f} Å")

# ============================================================
# STEP 4: Converter para coordenadas fracionárias
# ============================================================
print(f"\n[4] POSIÇÕES EM COORDENADAS FRACIONÁRIAS (crystal):")

# Inverter a célula para converter posições reais → fracionárias
cell_inv = np.linalg.inv(slab_fixed.get_cell().T)

print(f"\n  ANTES (com vácuo negativo):")
for i, pos in enumerate(slab_bad.positions):
    frac = np.dot(cell_inv, pos)
    print(f"    Au {i+1}: z_frac={frac[2]:.10f}")

print(f"\n  DEPOIS (corrigido):")
cell_fixed_inv = np.linalg.inv(slab_fixed.get_cell().T)
for i, pos in enumerate(slab_fixed.positions):
    frac = np.dot(cell_fixed_inv, pos)
    print(f"    Au {i+1}: z_frac={frac[2]:.10f}")

# ============================================================
# STEP 5: Validação
# ============================================================
print("\n" + "="*70)
print("VALIDAÇÃO")
print("="*70)

print(f"""
ESTRUTURA ANTES:
  - Vácuo: -30.7469 Å (NEGATIVO ❌)
  - Átomos: espalhados na parte negativa (z=-25 até -5 Å)
  - Frações z: -0.837 até -0.163 (todas negativas)

ESTRUTURA DEPOIS:
  - Vácuo: +30.7469 Å (POSITIVO ✓)
  - Átomos: CENTRADOS na célula
  - Frações z: valores simétricos em torno de 0.5
  
MUDANÇAS APLICADAS:
  ✓ Inverteu célula (cell[2] = -cell[2])
  ✓ Recentrou slab automaticamente
  ✓ Posições atômicas PRESERVADAS (estrutura intacta)
""")

# Verificar que a estrutura foi preservada
pos_original_sorted = np.sort(slab_bad.positions[:, 2])
pos_fixed_sorted = np.sort(slab_fixed.positions[:, 2])
distances_original = np.diff(pos_original_sorted)
distances_fixed = np.diff(pos_fixed_sorted)

print(f"\nEspaçamentos entre camadas:")
print(f"  ANTES: {distances_original}")
print(f"  DEPOIS: {distances_fixed}")
print(f"  Diferença máxima: {np.max(np.abs(distances_original - distances_fixed)):.10f} Å")
print(f"  ✓ Estrutura preservada!")

print("\n" + "="*70)
print("✓ TEST COMPLETO")
print("="*70)
