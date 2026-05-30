#!/usr/bin/env python
"""Debug script for initial B0 estimation."""

import numpy as np

# Synthetic test data (Fe BCC, small range)
factors = np.array([0.98, 1.0, 1.02])
v0_ref = 11.8199515
# Isotropic scaling: V(factor) = V0 * factor^3
volumes = v0_ref * factors**3
# Parabolic energy (minimum at factor=1.0)
energies = -50 + 0.5*(volumes - volumes[1])**2 / (volumes[1]**2)

print("=" * 70)
print("DEBUG: Initial B0 Estimation")
print("=" * 70)
print(f"Volumes (Ų):  {volumes}")
print(f"Energies (eV): {energies}")
print()

# Initial parameter guesses (from fit_birch_murnaghan)
v0_initial = volumes[np.argmin(energies)]
e0_initial = energies[np.argmin(energies)]

print(f"V₀ initial:   {v0_initial}")
print(f"E₀ initial:   {e0_initial}")
print()

# Estimate bulk modulus from data spread
dE = np.max(energies) - np.min(energies)
dV = np.max(volumes) - np.min(volumes)
print(f"dE (max-min): {dE}")
print(f"dV (max-min): {dV}")
print()

b0_ev_ang3 = dE * v0_initial / (dV**2)  # In eV/Ų
b0_initial_gpa = max(50, b0_ev_ang3 * 160.2)  # Convert to GPa, minimum 50 GPa

print(f"b0_ev_ang3 (estimated):   {b0_ev_ang3:.10f} eV/Ų")
print(f"b0_initial (after ×160.2): {b0_ev_ang3 * 160.2:.10f} GPa")
print(f"b0_initial (after max(50)): {b0_initial_gpa:.10f} GPa")
print()

# The problem: for small data ranges, the estimated B0 can be too small
# and get clamped to 50 GPa, but the optimization still struggles
print("ISSUE: For synthetic data with tiny energy variations (~0.002 eV),")
print("the estimated B0 is very small (1.6 GPa), but clamped to 50 GPa.")
print("The optimizer then has difficulty finding a good fit because")
print("the initial guess is far from the actual optimal value.")
print()
print("SOLUTION: For poor fits, consider:")
print("  1. Initial B0 guess improvement (use better heuristics)")
print("  2. Multiple optimization attempts with different initial values")
print("  3. Data quality checks (ensure sufficient energy variation)")
