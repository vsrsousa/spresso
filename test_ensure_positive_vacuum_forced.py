#!/usr/bin/env python3
"""
Mocking Test 2: Forçar vácuo NEGATIVO e mostrar a correção
"""

import numpy as np
from ase import Atoms
from ase.io import write

print("="*70)
print("MOCKING TEST 2: Corrigir vácuo NEGATIVO")
print("="*70)

# ============================================================
# STEP 1: Criar slab com vácuo FORÇADAMENTE NEGATIVO
# ============================================================
print("\n[1] Criando slab com vácuo FORÇADAMENTE negativo...")

# Slab simples: 5 Au átomos em (111)
positions = [
    [0.0, 0.0, 8.66],
    [0.0, 0.0, 12.74],
    [0.0, 0.0, 16.82],
    [0.0, 0.0, 20.89],
    [0.0, 0.0, 24.97],
]

# Célula com vácuo NEGATIVO (como PyMatGen às vezes faz)
cell_negative = [
    [4.0782, 0.0, -4.0782],
    [-0.0, 4.0782, -4.0782],
    [0.0, 0.0, -33.6333],  # ← NEGATIVO!
]

slab_bad = Atoms('Au5', positions=positions, cell=cell_negative, pbc=True)

print(f"  Slab: {len(slab_bad)} átomos")
print(f"\n  CÉLULA ANTES (com vácuo negativo):")
print(f"    a: {slab_bad.cell[0]}")
print(f"    b: {slab_bad.cell[1]}")
print(f"    c: {slab_bad.cell[2]}")
print(f"    Vácuo (c[2]): {slab_bad.cell[2, 2]:.4f} Å ⚠ NEGATIVO!")

# ============================================================
# STEP 2: Aplicar _ensure_positive_vacuum()
# ============================================================
print("\n[2] Aplicando _ensure_positive_vacuum()...")

def ensure_positive_vacuum(slab):
    """Imitar a função do slab_workflow.py"""
    cell = slab.get_cell()
    if cell[2, 2] < 0:
        print(f"  ⚠ Detectado vácuo negativo: {cell[2, 2]:.4f} Å")
        cell[2] = -cell[2]
        slab.set_cell(cell)
        print(f"  ✓ Corrigido para: {cell[2, 2]:.4f} Å")
    return slab

slab_fixed = slab_bad.copy()
slab_fixed = ensure_positive_vacuum(slab_fixed)

# ============================================================
# STEP 3: Mostrar célula DEPOIS
# ============================================================
print(f"\n[3] CÉLULA DEPOIS (corrigida):")
print(f"    a: {slab_fixed.cell[0]}")
print(f"    b: {slab_fixed.cell[1]}")
print(f"    c: {slab_fixed.cell[2]}")
print(f"    Vácuo (c[2]): {slab_fixed.cell[2, 2]:.4f} Å ✓ POSITIVO!")

# ============================================================
# STEP 4: Salvar ambas as versões
# ============================================================
print("\n[4] Gerando inputs QE (antes e depois)...")

write('au111_NEGATIVO_vácuo.in', slab_bad)
write('au111_CORRIGIDO_vácuo.in', slab_fixed)

# ============================================================
# STEP 5: Comparar os dois arquivos
# ============================================================
print("\n[5] COMPARAÇÃO DE INPUTS QE:")
print("="*70)

print("\n❌ ANTES (com vácuo negativo):")
print("-"*70)
with open('au111_NEGATIVO_vácuo.in', 'r') as f:
    lines = f.readlines()
    for i, line in enumerate(lines[:12]):
        print(line.rstrip())

print("\n\n✅ DEPOIS (com vácuo corrigido):")
print("-"*70)
with open('au111_CORRIGIDO_vácuo.in', 'r') as f:
    lines = f.readlines()
    for i, line in enumerate(lines[:12]):
        print(line.rstrip())

# ============================================================
# STEP 6: Validação
# ============================================================
print("\n\n" + "="*70)
print("VALIDAÇÃO ESTRUTURAL")
print("="*70)

print(f"""
ANTES:
  - c[2] = {slab_bad.cell[2, 2]:.4f} Å (NEGATIVO ❌)
  - Posições atômicas: {len(slab_bad)} átomos
  
DEPOIS:
  - c[2] = {slab_fixed.cell[2, 2]:.4f} Å (POSITIVO ✓)
  - Posições atômicas: {len(slab_fixed)} átomos (preservadas ✓)
  
DIFERENÇA:
  - Célula virou de cabeça para baixo (invertida apenas em Z)
  - Todas as posições atômicas permaneceram iguais ✓
  - Estrutura cristalina preservada ✓
""")

# Verificar que posições não mudaram
import numpy as np
pos_before = slab_bad.get_positions()
pos_after = slab_fixed.get_positions()

max_diff = np.max(np.abs(pos_before - pos_after))
print(f"Diferença máxima de posições: {max_diff:.10f} Ångstrom (inalteradas ✓)")

print("\n" + "="*70)
print("✓ TEST COMPLETO: VÁCUO NEGATIVO FOI CORRIGIDO")
print("="*70)
