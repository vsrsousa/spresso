#!/usr/bin/env python
"""
Convergence study with exponential fit analysis using ConvergenceWorkflow.

This example demonstrates:
1. Running Phase 1 (ecutwfc convergence with fixed kspacing)
2. Using exponential fit to extrapolate asymptotic energy E_inf
3. Determining optimal ecutwfc based on tolerance vs E_inf (not vs max_ecutwfc)
4. Running Phase 2 (kspacing convergence with optimal ecutwfc)
5. Comparing tested values vs exponential fit extrapolations

The new method:
- Fits: E(ecutwfc) = E_inf + A * exp(-B * ecutwfc)
- Uses E_inf (asymptotic) as reference instead of ecutwfc=200 Ry
- Finds minimum ecutwfc where |E - E_inf| < tolerance
- Returns R² goodness-of-fit metric

Benefits:
- More robust: uses all data via fit (not just highest tested point)
- Theoretically sound: E_inf is the true asymptotic limit
- Extrapolation: can predict ecutwfc needed for any tolerance
- Efficiency: may find that tested ecutwfc exceeds what's needed
"""

from ase.build import bulk
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow
import os

# ============================================================================
# Setup
# ============================================================================

print("="*80)
print("CONVERGENCE STUDY WITH EXPONENTIAL FIT")
print("Au Bulk (FCC) - Phase 1 (ecutwfc) + Phase 2 (kspacing)")
print("="*80)

# Create Au bulk structure
atoms = bulk('Au', 'fcc', a=4.0782)
print(f"\nStructure: {atoms.get_chemical_formula()}")
print(f"Cell: {atoms.get_cell().cellpar()}\n")

# Create workflow with custom parameters
wf = ConvergenceWorkflow(
    atoms=atoms,
    pseudopotentials={'Au': 'Au.pbe-n-rrkjus_psl.1.0.0.UPF'},
    protocol='moderate',
    precision='low',
    min_ecutwfc=30.0,
    max_ecutwfc=100.0,  # Reduced for faster testing
    initial_kspacing=0.30,  # Å⁻¹
)

# ============================================================================
# Run convergence study (Independent Mode: Phase 1 + Phase 2)
# ============================================================================

print("\nRunning convergence study (independent mode)...")
print("-" * 80)

results = wf.run_convergence_study(
    label_prefix='au_convergence/exponential_fit',
    max_ecutwfc=100.0,
    ecut_step=10.0,
    min_kspacing_allowed=0.10,  # Å⁻¹
    kspacing_step=0.03,
    verbose=True,
    batch_timeout=3600,
)

# ============================================================================
# Get recommendations (includes exponential fit analysis)
# ============================================================================

print("\n" + "="*80)
recommendations = wf.get_recommendations(verbose=True)

# ============================================================================
# Print detailed fit analysis
# ============================================================================

if 'exponential_fit' in recommendations:
    fit = recommendations['exponential_fit']
    print("\n" + "="*80)
    print("DETAILED EXPONENTIAL FIT ANALYSIS")
    print("="*80)
    print(f"\nFitted equation: E(ecutwfc) = E_inf + A * exp(-B * ecutwfc)")
    print(f"\nParameters:")
    print(f"  E_inf = {fit['E_inf']:.8f} eV     (asymptotic energy as ecutwfc → ∞)")
    print(f"  A     = {fit['A']:.8f} eV     (exponential amplitude)")
    print(f"  B     = {fit['B']:.6f} Ry⁻¹  (decay constant)")
    print(f"  R²    = {fit['R_squared']:.6f}  (goodness of fit: 1.0 = perfect)")
    
    print(f"\nOptimal ecutwfc selection:")
    print(f"  Target tolerance: {recommendations['energy_tolerance_meV_atom']:.2f} meV/atom")
    print(f"  E_inf (extrapolated) = {fit['E_inf']:.8f} eV")
    print(f"  Min ecutwfc for tolerance = {fit['min_ecutwfc_for_tolerance']:.1f} Ry (extrapolated)")
    print(f"  Tested ecutwfc = {recommendations['optimal_ecutwfc']:.1f} Ry (recommended)")
    
    # Show margin
    if fit['R_squared'] > 0.99:
        print(f"\n  ✓ Excellent fit (R² = {fit['R_squared']:.6f})")
        print(f"    Extrapolations are highly reliable")
    elif fit['R_squared'] > 0.95:
        print(f"\n  ✓ Good fit (R² = {fit['R_squared']:.6f})")
        print(f"    Extrapolations are reasonably reliable")
    else:
        print(f"\n  ⚠ Fair fit (R² = {fit['R_squared']:.6f})")
        print(f"    Use extrapolations with caution")

print("\n" + "="*80)
print("CONVERGENCE STUDY COMPLETE")
print("="*80)
print(f"\nResults saved in: au_convergence/exponential_fit/")
print(f"  - SCF calculations at different ecutwfc and kspacing values")
print(f"  - Convergence analysis with exponential fit")
print(f"  - Optimal parameters: ecutwfc={recommendations['optimal_ecutwfc']} Ry, kspacing={recommendations['optimal_kspacing']} Å⁻¹")
