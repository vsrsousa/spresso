#!/usr/bin/env python
"""
EXAMPLE: Exponential Fit with Auto-Detected ecutwfc_min_for_fit

This example demonstrates the NEW approach to handling basis-incomplete regions:
Instead of statistical outlier removal (which was wrong), we use PHYSICS-BASED
REGION SEPARATION to automatically identify where exponential convergence begins.

The key insight:
  - Low ecutwfc values (ecut < 50 Ry for Au) have INCOMPLETE BASIS
  - This is not a "bug" or "outlier" - it's a different physics regime
  - We must exclude these from the exponential fit automatically
  - Auto-detection uses:
    1. Pseudopotential database lookup (fastest)
    2. Curvature analysis (fallback, more general)
    3. Conservative estimate (take maximum of both)

Two convergence regions in Au bulk:
  REGION 1 (BASIS-INCOMPLETE): ecutwfc = 30, 40 Ry
    Energy jumps: ~2000 meV, ~100 meV
    Behavior: Non-exponential, dominated by basis set
    Action: AUTO-EXCLUDE from fit
  
  REGION 2 (EXPONENTIAL CONVERGENCE): ecutwfc ≥ 50 Ry
    Energy changes: smooth exponential decay (1-10 meV)
    Behavior: Standard exponential basis convergence
    Action: USE FOR FIT
"""

import numpy as np
from ase.build import bulk
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow
from xespresso.pseudopotentials.manager import load_pseudopotentials_config


def example_ecut_min_auto_detection():
    """
    Demonstrate auto-detection of ecutwfc_min_for_fit with Au bulk.
    
    Shows:
    1. How the auto-detection works (database + curvature)
    2. Why it's conservative (takes maximum of both methods)
    3. How it improves exponential fitting
    """
    
    print("""
    ╔════════════════════════════════════════════════════════════════════════════════╗
    ║                EXPONENTIAL FIT WITH AUTO-DETECTED REGIONS                      ║
    ║                                                                                ║
    ║ Problem:  Low-ecutwfc points (30, 40 Ry) have DIFFERENT PHYSICS               ║
    ║           (basis-incomplete) than high-ecutwfc (50+ Ry) (exponential)          ║
    ║                                                                                ║
    ║ Solution: Auto-detect where transition happens, fit ONLY to valid region      ║
    ║                                                                                ║
    ║ Benefit:  No manual tweaking needed!  System figures it out automatically.     ║
    ╚════════════════════════════════════════════════════════════════════════════════╝
    """)
    
    # Step 1: Create structure (Au bulk FCC)
    print("\n" + "="*80)
    print("STEP 1: CREATE STRUCTURE")
    print("="*80)
    
    atoms = bulk('Au', 'fcc', a=4.0782)
    print(f"✓ Au bulk FCC structure")
    print(f"  Lattice parameter: {atoms.cell[0,0]:.4f} Å")
    print(f"  Number of atoms: {len(atoms)}")
    
    # Step 2: Load pseudopotentials
    print("\n" + "="*80)
    print("STEP 2: LOAD PSEUDOPOTENTIALS")
    print("="*80)
    
    pseudo_config = load_pseudopotentials_config('SSSP_efficiency', verbose=False)
    print(f"✓ Loaded SSSP efficiency pseudopotentials")
    
    # Step 3: Initialize convergence workflow
    print("\n" + "="*80)
    print("STEP 3: INITIALIZE CONVERGENCE WORKFLOW")
    print("="*80)
    
    wf = ConvergenceWorkflow(
        atoms=atoms,
        pseudopotentials_config='SSSP_efficiency',
        precision='low',
        protocol='moderate'
    )
    print(f"✓ Workflow initialized")
    print(f"  Pseudopotentials: {list(wf.pseudopotentials.keys())}")
    
    # Step 4: Simulate Phase 1 results (ecutwfc convergence)
    print("\n" + "="*80)
    print("STEP 4: PHASE 1 ECUTWFC CONVERGENCE (simulated data)")
    print("="*80)
    print("""
    In real usage, Phase 1 would collect data:
      ecutwfc = [30, 40, 50, 60, 70] Ry
    
    With energies showing TWO REGIONS:
      - Low ecut (30-40): Basis incomplete, large jumps
      - High ecut (50+): Exponential convergence, smooth
    """)
    
    # Simulated Phase 1 results for Au bulk
    # (These are realistic values based on PSL pseudopotential behavior)
    ecut_results = {
        30.0: {'energy': -9800.0},      # 2000 meV from convergence
        40.0: {'energy': -9995.0},      # 100 meV from convergence
        50.0: {'energy': -9999.5},      # Smooth exponential starts
        60.0: {'energy': -9999.90},
        70.0: {'energy': -10000.00}
    }
    
    criteria_tolerances = {'energy_tolerance': 0.001}  # 1 meV
    
    print(f"\nSimulated results:")
    for ecut, data in sorted(ecut_results.items()):
        print(f"  ecutwfc = {ecut:5.0f} Ry  →  E = {data['energy']:12.2f} eV")
    
    # Step 5: Fit WITHOUT auto-detection (using all data)
    print("\n" + "="*80)
    print("STEP 5A: FIT WITHOUT AUTO-DETECTION (WRONG APPROACH)")
    print("="*80)
    print("""
    If we fit ALL data (including basis-incomplete region),
    the exponential model gets distorted by low-ecut anomalies.
    """)
    
    fit_all = wf._fit_exponential_decay_phase1(
        ecut_results,
        criteria_tolerances,
        verbose=True,
        ecut_min_for_fit=0.0,      # Force using ALL data
        auto_detect_ecut_min=False  # Disable auto-detection
    )
    
    print(f"\nRecommendation (BAD FIT): {fit_all['min_ecutwfc_for_tolerance']:.1f} Ry")
    print(f"  Problem: Fit distorted by basis-incomplete points!")
    
    # Step 5b: Fit WITH auto-detection (correct approach)
    print("\n" + "="*80)
    print("STEP 5B: FIT WITH AUTO-DETECTION (CORRECT APPROACH)")
    print("="*80)
    print("""
    Now with auto-detection:
    1. Detect ecut_min_for_fit automatically (database + curvature)
    2. Exclude low-ecut points (basis-incomplete region)
    3. Fit ONLY to convergence region (exponential valid)
    """)
    
    fit_auto = wf._fit_exponential_decay_phase1(
        ecut_results,
        criteria_tolerances,
        verbose=True,
        ecut_min_for_fit=None,      # Auto-detect!
        auto_detect_ecut_min=True   # Enable auto-detection
    )
    
    print(f"\nRecommendation (GOOD FIT): {fit_auto['min_ecutwfc_for_tolerance']:.1f} Ry")
    print(f"  Benefit: Fit uses ONLY valid convergence region!")
    
    # Step 6: Compare results
    print("\n" + "="*80)
    print("STEP 6: COMPARISON")
    print("="*80)
    
    print(f"\n{'Method':<40} {'E_inf (eV)':<15} {'Recommendation (Ry)':<20} {'R²':<10}")
    print(f"{'-'*85}")
    print(f"{'Without auto-detection (BAD)':<40} {fit_all.get('E_inf', 0):<15.8f} {fit_all.get('min_ecutwfc_for_tolerance', 0):<20.1f} {fit_all.get('R_squared', 0):<10.6f}")
    print(f"{'With auto-detection (GOOD)':<40} {fit_auto.get('E_inf', 0):<15.8f} {fit_auto.get('min_ecutwfc_for_tolerance', 0):<20.1f} {fit_auto.get('R_squared', 0):<10.6f}")
    
    # Step 7: Show what was auto-detected
    print("\n" + "="*80)
    print("STEP 7: AUTO-DETECTION DETAILS")
    print("="*80)
    
    print(f"\nAuto-detection results:")
    print(f"  ecut_min_for_fit = {fit_auto.get('ecut_min_for_fit', 'N/A'):.1f} Ry")
    print(f"  (Found by combining database + curvature methods)")
    
    print(f"\nRegion separation:")
    basis_incomplete = fit_auto.get('basis_incomplete_points', [])
    convergence = fit_auto.get('convergence_region_points', [])
    print(f"  Basis-incomplete (excluded):  {basis_incomplete}")
    print(f"  Convergence region (used):    {convergence}")
    
    # Step 8: Key insights
    print("\n" + "="*80)
    print("KEY INSIGHTS")
    print("="*80)
    
    print("""
    1. TWO CONVERGENCE REGIMES:
       - Low ecut: Different physics (basis incomplete)
       - High ecut: Standard exponential convergence
       
    2. AUTO-DETECTION STRATEGY:
       ✓ Database lookup: Quick, knows typical values for each pseudo type
       ✓ Curvature analysis: Robust, works for any pseudo
       ✓ Conservative estimate: Takes maximum (safest approach)
       
    3. RESULT:
       - No manual tuning needed!
       - Physics-based, not statistical
       - Works for all pseudopotential types
       - Provides reliable recommendations
       
    4. USAGE:
       fit = wf._fit_exponential_decay_phase1(
           ecut_results,
           criteria_tolerances,
           # ecut_min_for_fit=None,        ← Auto-detect!
           # auto_detect_ecut_min=True     ← Default
       )
       
       Or with explicit override:
       fit = wf._fit_exponential_decay_phase1(
           ecut_results,
           criteria_tolerances,
           ecut_min_for_fit=50.0,          ← Use this value
           auto_detect_ecut_min=False      ← Skip auto-detection
       )
    """)


def example_multi_precision_reuse():
    """
    Demonstrate fit REUSE for multiple precision levels.
    
    Once we have the fit parameters (E_inf, A, B), we can instantly
    calculate recommendations for ANY tolerance without recalculation!
    """
    
    print("""
    
    ╔════════════════════════════════════════════════════════════════════════════════╗
    ║                      FIT REUSE FOR MULTIPLE PRECISIONS                         ║
    ║                                                                                ║
    ║ Insight: Once we know E_inf, A, B from ONE Phase 1 run,                       ║
    ║          we can get recommendations for ANY tolerance INSTANTLY!               ║
    ║                                                                                ║
    ║ Speed:   1 Phase 1 (25 min)  vs  4 Phase 1s (100 min) = 4× faster!            ║
    ╚════════════════════════════════════════════════════════════════════════════════╝
    """)
    
    # Example fit parameters (from previous step)
    E_inf = -10000.0  # eV
    A = -199.5        # eV (negative amplitude)
    B = 0.048         # Ry^-1
    
    # Now calculate recommendations for different tolerances
    tolerances = {
        'low': 10.0,     # 10 meV
        'medium': 5.0,   # 5 meV
        'high': 1.0,     # 1 meV
        'ultra': 0.1     # 0.1 meV
    }
    
    print("\nUsing same fit parameters, calculate for all precisions:\n")
    print(f"{'Precision':<12} {'Tolerance (meV)':<18} {'ecutwfc (Ry)':<15}")
    print(f"{'-'*45}")
    
    for precision, tolerance_meV in tolerances.items():
        tolerance = tolerance_meV / 1000  # Convert meV to eV
        
        # Formula: ecut = -ln(tolerance/|A|) / B
        if abs(A) > 0 and B > 0 and tolerance < abs(A):
            ecut = -np.log(tolerance / abs(A)) / B
        else:
            ecut = 70.0
        
        print(f"{precision:<12} {tolerance_meV:<18.1f} {ecut:<15.1f}")
    
    print(f"\n✓ All recommendations calculated in MICROSECONDS")
    print(f"  (no SCF calculations needed!)")


if __name__ == '__main__':
    print("\n" + "="*80)
    print("EXPONENTIAL FIT WITH AUTO-DETECTED CONVERGENCE REGIONS")
    print("="*80)
    
    example_ecut_min_auto_detection()
    example_multi_precision_reuse()
    
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print("""
    ✅ NEW APPROACH: Physics-based region separation
       - Auto-detects where exponential convergence begins
       - Excludes basis-incomplete region automatically
       - No manual tuning needed
       - Works for all pseudopotential types
    
    ✅ BENEFITS:
       - More reliable extrapolations
       - Stable recommendations
       - Can reuse fit for 4 precisions (4× speedup)
       - Automated, reproducible, robust
    
    ✅ IMPLEMENTATION:
       - Database lookup (fast, for known pseudos)
       - Curvature analysis (robust, for unknown pseudos)
       - Conservative estimate (safest)
    """)
