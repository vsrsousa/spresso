#!/usr/bin/env python
"""
Quick EOS test - local execution only (no remote/SLURM)
"""

from ase.build import bulk
from xespresso.workflow import EOSWorkflow
import os

# Create test structure
print("=" * 70)
print("TESTE RÁPIDO: EOSWorkflow com scale_factors exatos")
print("=" * 70)

atoms = bulk('Fe', 'bcc', a=2.87)
print(f"✓ Estrutura: {atoms.get_chemical_formula()}")
print(f"  Volume: {atoms.get_volume():.4f} Ų")

# Pseudopotentials (usar config local)
pseudopotentials_config = 'default'

# Create EOS workflow (LOCAL only - no machine parameter)
eos = EOSWorkflow(
    atoms=atoms,
    pseudopotentials_config=pseudopotentials_config,
    protocol='moderate',
    machine=None  # LOCAL execution
)

print(f"✓ EOSWorkflow inicializado")
print(f"  Protocol: moderate")
print(f"  Machine: LOCAL (sem remoto)")

# Test 1: DRY RUN - apenas gera inputs
print("\n" + "=" * 70)
print("TEST 1: DRY RUN (apenas gera arquivos de input, sem rodar)")
print("=" * 70)

try:
    df = eos.run_eos_study(
        scale_factors=[0.98, 1.00, 1.02],  # ← 3 pontos exatos
        dry_run=True,  # ← DRY RUN
        parallel=False,
        label='test_eos_dry'
    )
    
    print(f"✓ DRY RUN concluído")
    print(f"  Diretórios criados em: test_eos_dry/")
    
    # List the generated files
    for root, dirs, files in os.walk('test_eos_dry'):
        for file in files:
            if file.endswith('.pwi'):
                filepath = os.path.join(root, file)
                print(f"    - {filepath}")
    
except Exception as e:
    print(f"✗ Erro no dry run: {e}")

# Test 2: Synthetic test - sem rodar Espresso
print("\n" + "=" * 70)
print("TEST 2: Synthetic data (testa EOS fitting sem rodar Espresso)")
print("=" * 70)

import numpy as np

# Create synthetic E-V data
factors = np.array([0.98, 1.00, 1.02])
volumes = atoms.get_volume() * factors**(3)  # isotropic scaling
energies = -50 + 0.5*(volumes - volumes[1])**2 / (volumes[1]**2)  # parabola

print(f"Dados sintéticos:")
print(f"  Factors:  {factors}")
print(f"  Volumes:  {volumes}")
print(f"  Energies: {energies}")

# Manually add data to EOS object
eos.results_df = __import__('pandas').DataFrame({
    'factor': factors,
    'volume': volumes,
    'energy': energies
})

# Fit EOS
try:
    eos.fit_eos()
    print(f"\n✓ EOS fitting concluído!")
    
    # Get properties
    props = eos.get_eos_properties()
    print(f"\nPropriedades EOS:")
    print(f"  V₀ (Ų):              {props['v0']:.4f}")
    print(f"  E₀ (eV):             {props['e0']:.6f}")
    print(f"  B₀ (GPa):            {props['bulk_modulus']:.2f}")
    print(f"  B₀' (adimensional):  {props['bulk_modulus_prime']:.2f}")
    print(f"  R²:                  {props['r_squared']:.6f}")
    print(f"  Convergido:          {props['converged']}")
    
except Exception as e:
    print(f"✗ Erro no fitting: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 70)
print("TESTES CONCLUÍDOS")
print("=" * 70)
