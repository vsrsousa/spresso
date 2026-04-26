#!/usr/bin/env python
"""
Exponential Fit Reuse: Estimate ecutwfc for MULTIPLE precisions from single Phase 1.

This example demonstrates the most powerful feature of exponential fit:
ONCE you have the fit (E_inf, A, B), you can estimate ecutwfc for ANY tolerance
WITHOUT refitting or recalculating!

Workflow:
1. Run Phase 1 with precision='low' (cheapest)
   └─ Gets exponential fit with E_inf, A, B
   
2. Use fit to extrapolate ecutwfc for:
   - precision='low'    (1.0 meV)
   - precision='medium' (0.5 meV)
   - precision='high'   (0.1 meV)
   - precision='ultra'  (0.01 meV)
   
3. ZERO additional calculations needed!

This is the ultimate efficiency: Get recommendations for all precision levels
from a single Phase 1 convergence study.
"""

from ase.build import bulk
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow
import numpy as np

print("="*80)
print("EXPONENTIAL FIT REUSE: Multi-Precision Recommendations from Single Phase 1")
print("="*80)

# ============================================================================
# Step 1: Run Phase 1 (single time, cheapest precision)
# ============================================================================

print("\n🔧 STEP 1: Run Phase 1 with precision='low'")
print("-" * 80)

atoms = bulk('Au', 'fcc', a=4.0782)
print(f"Structure: {atoms.get_chemical_formula()}")

# Create with low precision (cheapest Phase 1)
wf = ConvergenceWorkflow(
    atoms=atoms,
    pseudopotentials={'Au': 'Au.pbe-n-rrkjus_psl.1.0.0.UPF'},
    protocol='moderate',
    precision='low',  # ← Cheap Phase 1
    min_ecutwfc=30.0,
    max_ecutwfc=80.0,  # Reduced for demo
    initial_kspacing=0.30,
)

print("\nRunning Phase 1 (ecutwfc convergence with fixed kspacing)...")
results = wf.run_convergence_study(
    label_prefix='au_multi_precision/phase1',
    max_ecutwfc=80.0,
    ecut_step=10.0,
    min_kspacing_allowed=0.10,
    kspacing_step=0.03,
    verbose=False,  # Quiet for clarity
)

print("✅ Phase 1 complete with exponential fit!")

# ============================================================================
# Step 2: Get baseline recommendation for precision='low'
# ============================================================================

print("\n📊 STEP 2: Get baseline recommendation")
print("-" * 80)

rec_low = wf.get_recommendations(verbose=False)
print(f"\nPrecision='low' (tested):")
print(f"  optimal_ecutwfc: {rec_low['optimal_ecutwfc']:.1f} Ry")
print(f"  tolerance: {rec_low['energy_tolerance_meV_atom']:.2f} meV/atom")

if 'exponential_fit' in rec_low:
    fit = rec_low['exponential_fit']
    print(f"\n✅ Exponential fit available:")
    print(f"  E_inf = {fit['E_inf']:.8f} eV")
    print(f"  A = {fit['A']:.8f} eV")
    print(f"  B = {fit['B']:.6f} Ry⁻¹")
    print(f"  R² = {fit['R_squared']:.6f}")

# ============================================================================
# Step 3: Extrapolate for ALL precision levels (ZERO new calculations!)
# ============================================================================

print("\n" + "="*80)
print("⭐ STEP 3: EXTRAPOLATE for all precision levels (NO new calculations!)")
print("="*80)

multi_rec = wf.recommend_for_multiple_precisions(verbose=True)

# ============================================================================
# Step 4: Detailed analysis
# ============================================================================

print("\n📈 DETAILED ANALYSIS")
print("-" * 80)

# Show the exponential decay at different ecutwfc values
if 'exponential_fit' in rec_low:
    fit = rec_low['exponential_fit']
    E_inf = fit['E_inf']
    A = fit['A']
    B = fit['B']
    
    print(f"\nEnergy convergence at different ecutwfc (using fit):")
    print(f"{'ecutwfc (Ry)':<15} {'E (eV)':<20} {'ΔE (meV)':<15} {'Precision':<15}")
    print("-" * 65)
    
    test_ecutwfc = [30, 40, 50, 60, 70, 80, 100, 120, 150, 200]
    for ecut in test_ecutwfc:
        E = E_inf + A * np.exp(-B * ecut)
        dE_meV = abs(E - E_inf) * 1000
        
        # Assign precision based on tolerance
        if dE_meV < 0.01:
            precision_label = "ultra (0.01)"
        elif dE_meV < 0.1:
            precision_label = "high (0.1)"
        elif dE_meV < 0.5:
            precision_label = "medium (0.5)"
        elif dE_meV < 1.0:
            precision_label = "low (1.0)"
        else:
            precision_label = "very coarse"
        
        print(f"{ecut:<15.0f} {E:<20.8f} {dE_meV:<15.4f} {precision_label:<15}")

# ============================================================================
# Step 5: Comparison with single high-precision Phase 1
# ============================================================================

print("\n💡 COMPUTATIONAL EFFICIENCY")
print("-" * 80)

print("\nScenario A: Individual Phase 1 for each precision")
print("  precision='low'    → Run Phase 1 (5 jobs × ~5 min = 25 min)")
print("  precision='medium' → Run Phase 1 (5 jobs × ~5 min = 25 min)")
print("  precision='high'   → Run Phase 1 (5 jobs × ~5 min = 25 min)")
print("  precision='ultra'  → Run Phase 1 (5 jobs × ~5 min = 25 min)")
print("  ────────────────────────────────────────────────────────")
print("  TOTAL TIME: ~100 min (20×5 jobs)")

print("\nScenario B: Single Phase 1 + exponential fit extrapolation")
print("  precision='low'    → Run Phase 1 (5 jobs × ~5 min = 25 min)")
print("  precision='medium' → EXTRAPOLATE from fit (instant)")
print("  precision='high'   → EXTRAPOLATE from fit (instant)")
print("  precision='ultra'  → EXTRAPOLATE from fit (instant)")
print("  ────────────────────────────────────────────────────────")
print("  TOTAL TIME: ~25 min (5 jobs only!)")
print(f"\n  🚀 SPEEDUP: 4× faster!")

# ============================================================================
# Step 6: Practical recommendations
# ============================================================================

print("\n" + "="*80)
print("✅ PRACTICAL RECOMMENDATIONS (from fit)")
print("="*80)

precision_names = {
    'low': 'Quick scans, initial screening',
    'medium': 'Production calculations, good balance',
    'high': 'Publication-quality results',
    'ultra': 'High-precision, demanding applications',
}

print()
for precision, ecut in multi_rec.items():
    if ecut is not None:
        print(f"{precision.upper():8s}: ecutwfc = {ecut:6.1f} Ry")
        print(f"          {precision_names.get(precision, '')}")
        print()

print("="*80)
print("✨ KEY INSIGHT: Once you have exponential fit, you can answer:")
print("  'What ecutwfc do I need for X meV tolerance?'")
print("   WITHOUT any additional calculations!")
print("="*80)

# ============================================================================
# Advanced: Estimate for custom tolerance
# ============================================================================

print("\n🎯 ADVANCED: Custom tolerance estimates")
print("-" * 80)

custom_tolerances = [0.5, 0.2, 0.05, 0.02]
print(f"\nEstimated ecutwfc for custom tolerances:")
print(f"{'Tolerance (meV)':<20} {'Ecutwfc (Ry)':<15}")
print("-" * 35)

for tol in custom_tolerances:
    ecut = wf.estimate_ecutwfc_for_tolerance(tol)
    if ecut is not None:
        print(f"{tol:<20.2f} {ecut:<15.1f}")

print("\n" + "="*80)
print("EXPONENTIAL FIT REUSE: Complete!")
print("="*80)
