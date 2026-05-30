#!/usr/bin/env python
"""Debug script for R² calculation."""

import numpy as np
from xespresso.workflow.eos_workflow import fit_birch_murnaghan, birch_murnaghan_eos

# Synthetic test data (Fe BCC, small range)
factors = np.array([0.98, 1.0, 1.02])
v0_ref = 11.8199515
# Isotropic scaling: V(factor) = V0 * factor^3
volumes = v0_ref * factors**3
# Parabolic energy (minimum at factor=1.0)
energies = -50 + 0.5*(volumes - volumes[1])**2 / (volumes[1]**2)

print("=" * 70)
print("DEBUG: R² Calculation")
print("=" * 70)
print(f"Volumes (Ų):  {volumes}")
print(f"Energies (eV): {energies}")
print(f"Mean energy:   {np.mean(energies):.10f}")
print(f"SS_tot = sum((E - mean(E))²)")

# Manual calculation
mean_e = np.mean(energies)
ss_tot_manual = np.sum((energies - mean_e)**2)
print(f"  SS_tot = {ss_tot_manual:.15f}")
print()

# Fit EOS
print("Fitting Birch-Murnaghan EOS...")
result = fit_birch_murnaghan(volumes, energies)

print(f"\nFit Results:")
print(f"  E₀:        {result['e0']:.10f} eV")
print(f"  V₀:        {result['v0']:.10f} Ų")
print(f"  B₀ (GPa):  {result['b0'] * 160.2:.10f}")  # Convert back to GPa to see
print(f"  B₀ (eV/Ų): {result['b0']:.10f}")
print(f"  B':        {result['b0_prime']:.10f}")
print(f"  R²:        {result['r_squared']:.10f}")
print(f"  Converged: {result['converged']}")
print()

# Check fitted values
print("Validation:")
E_fit = birch_murnaghan_eos(
    volumes, 
    result['e0'],
    result['v0'],
    result['b0'] * 160.2,  # Convert back to GPa
    result['b0_prime']
)
print(f"  Fitted energies: {E_fit}")
print(f"  Residuals:       {energies - E_fit}")
ss_res = np.sum((energies - E_fit)**2)
print(f"  SS_res:          {ss_res:.15f}")
print(f"  SS_tot:          {ss_tot_manual:.15f}")
r2_manual = 1 - (ss_res / ss_tot_manual) if ss_tot_manual > 0 else 0
print(f"  R² (manual):     {r2_manual:.10f}")
