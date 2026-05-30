#!/usr/bin/env python
"""
Test EOS fitting with real data from user
"""

import numpy as np
import pandas as pd
from xespresso.workflow import EOSWorkflow
from ase.build import bulk

# User's real data
data = {
    'factor': [0.90, 0.92, 0.94, 0.96, 0.98, 1.00, 1.02, 1.04, 1.06, 1.08, 1.10],
    'energy': [-310.554809, -310.606456, -310.661272, -310.702053, -310.730165, -310.746915, -310.753498, -310.750984, -310.740326, -310.722380, -310.697792],
    'volume': [36.023177, 36.823692, 37.624207, 38.424722, 39.225237, 40.025752, 40.826267, 41.626782, 42.427297, 43.227812, 44.028327],
}

print("=" * 70)
print("TEST: EOS Fitting with Real Data (ASE algorithm)")
print("=" * 70)

# Create a dummy atoms object (not used for data but needed for workflow)
atoms = bulk('Fe', 'bcc', a=2.87)
eos = EOSWorkflow(atoms=atoms, pseudopotentials_config='default', protocol='moderate', machine=None)

# Manually inject the real data
eos.results_df = pd.DataFrame(data)

print(f"\nData points: {len(data['energy'])}")
print(f"Volume range: {np.min(data['volume']):.4f} - {np.max(data['volume']):.4f} Ų")
print(f"Energy range: {np.min(data['energy']):.6f} - {np.max(data['energy']):.6f} eV")

# Fit EOS
try:
    eos.fit_eos()
    print(f"\n✓ EOS fitting concluído!")
    
    props = eos.get_eos_properties()
    print(f"\n{'='*70}")
    print(f"EOS Fit Results (using ASE algorithm)")
    print(f"{'='*70}")
    print(f"V₀ (Ų):              {props['v0']:.6f}")
    print(f"E₀ (eV):             {props['e0']:.6f}")
    print(f"B₀ (GPa):            {props['bulk_modulus']:.2f}")
    print(f"B₀' (adimensional):  {props['bulk_modulus_prime']:.2f}")
    print(f"R²:                  {props['r_squared']:.6f}")
    print(f"Converged:           {props['converged']}")
    
    # Compare with ASE directly
    print(f"\n{'='*70}")
    print(f"Comparison with ASE EquationOfState")
    print(f"{'='*70}")
    
    from ase.eos import EquationOfState
    eos_ase = EquationOfState(data['volume'], data['energy'], eos='birchmurnaghan')
    v0_ase, e0_ase, b0_ase = eos_ase.fit()
    
    print(f"V₀ (Ų):              {v0_ase:.6f}")
    print(f"E₀ (eV):             {e0_ase:.6f}")
    print(f"B₀ (eV/Ų):           {b0_ase:.6f}")
    print(f"B₀ (GPa):            {b0_ase * 160.217662:.2f}")
    print(f"B₀' (adimensional):  {eos_ase.eos_parameters[2]:.2f}")
    
    # Check if results match
    print(f"\n{'='*70}")
    print(f"Verification (differences)")
    print(f"{'='*70}")
    print(f"ΔV₀:   {abs(props['v0'] - v0_ase):.8f} Ų ({abs(props['v0'] - v0_ase)/v0_ase*100:.4f}%)")
    print(f"ΔE₀:   {abs(props['e0'] - e0_ase):.8f} eV ({abs(props['e0'] - e0_ase)/abs(e0_ase)*100:.4f}%)")
    print(f"ΔB₀:   {abs(props['bulk_modulus'] - b0_ase*160.217662):.4f} GPa ({abs(props['bulk_modulus'] - b0_ase*160.217662)/(b0_ase*160.217662)*100:.4f}%)")
    print(f"ΔB₀':  {abs(props['bulk_modulus_prime'] - eos_ase.eos_parameters[2]):.4f}")
    
except Exception as e:
    print(f"✗ Erro: {e}")
    import traceback
    traceback.print_exc()
