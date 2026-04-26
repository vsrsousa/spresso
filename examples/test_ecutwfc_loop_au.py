#!/usr/bin/env python
"""
Batch ecutwfc convergence test for Au bulk (FCC) with parallel submission.

This test demonstrates:
1. Batch submission (all jobs submitted at once, non-blocking)
2. Parallel execution via Slurm queue
3. Simultaneous monitoring of all jobs
4. Energy convergence tracking
5. Comparison of serial vs parallel execution

Run with:
    python examples/test_ecutwfc_loop_au.py
"""

from ase.build import bulk
from xespresso.workflow.calculation_workflow import CalculationWorkflow
import os
import time
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# ============================================================================
# Setup
# ============================================================================

print("="*80)
print("ECUTWFC CONVERGENCE TEST - Au Bulk (FCC) - BATCH MODE")
print("="*80)

# Create Au bulk structure
atoms = bulk('Au', 'fcc', a=4.0782)
print(f"\nStructure: {atoms.get_chemical_formula()}")
print(f"Atoms: {len(atoms)}")
print(f"Cell: {atoms.get_cell().cellpar()}\n")

# Parameters
ecut_vals = [30, 40, 50, 60, 70]
pseudopotentials = {'Au': 'Au.pbe-n-rrkjus_psl.1.0.0.UPF'}
kpts = (8, 8, 8)  # 8x8x8 k-point mesh
ecut_min = None  # Minimum ecutwfc to use for fit (e.g., 40 to exclude first point)
                 # Set to None to use all points

# ============================================================================
# Create base workflow (will be used for batch submission)
# ============================================================================

print("Creating base workflow...")
wf = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials=pseudopotentials,
    protocol='moderate',
    ecutwfc=ecut_vals[0],  # Will be overridden per job
)

# ============================================================================
# Prepare parameter sets for batch submission
# ============================================================================

parameter_sets = [
    {
        'label': f'au_test/ecut{ecutwfc}',
        'ecutwfc': ecutwfc,
    }
    for ecutwfc in ecut_vals
]

print(f"Prepared {len(parameter_sets)} parameter sets for batch submission\n")

# ============================================================================
# BATCH SUBMISSION (non-blocking)
# ============================================================================

print("="*80)
print("BATCH SUBMISSION PHASE")
print("="*80)
start_submit = time.time()

# Submit ALL jobs at once (returns immediately, jobs run in parallel via Slurm)
batch_results = wf.submit_scf_batch_multiple(parameter_sets, verbose=True)

submit_time = time.time() - start_submit
print(f"Submission completed in {submit_time:.2f} seconds\n")

# ============================================================================
# BATCH MONITORING (parallel)
# ============================================================================

print("="*80)
print("BATCH MONITORING PHASE (Parallel)")
print("="*80)
start_wait = time.time()

# Monitor ALL jobs simultaneously (Slurm queue handles parallelization)
final_results = wf.wait_for_batch_jobs(batch_results, verbose=True)

wait_time = time.time() - start_wait
total_time = time.time() - start_submit

# ============================================================================
# Extract results
# ============================================================================

print("\n" + "="*80)
print("RESULTS")
print("="*80)

results = {}
energies = []

print(f"\n{'Iteration':<12} {'ecutwfc':<12} {'Energy (eV)':<18}")
print("-" * 70)

for i, result in enumerate(final_results, 1):
    # wait_for_batch_jobs returns results without 'calc', but with 'energy' and 'success'
    if result and result.get('success') and result.get('energy') is not None:
        try:
            ecutwfc = ecut_vals[i-1]
            energy = result['energy']
            energies.append(energy)
            
            results[ecutwfc] = energy
            print(f"{i:<12} {ecutwfc:<12} {energy:<18.8f}")
            
        except Exception as e:
            print(f"{i:<12} {ecut_vals[i-1]:<12} ERROR: {str(e)[:40]}")
            results[ecut_vals[i-1]] = None
    else:
        ecutwfc = ecut_vals[i-1]
        status = "FAILED"
        if result:
            status = result.get('error', 'FAILED')
        print(f"{i:<12} {ecutwfc:<12} {status}")
        results[ecutwfc] = None

print("-" * 70)

# ============================================================================
# Convergence Analysis with Exponential Fit
# ============================================================================

print("\n" + "="*80)
print("CONVERGENCE ANALYSIS WITH EXPONENTIAL FIT")
print("="*80)

# Extract successful results
ecutwfc_vals = []
energy_vals = []
for i, ecutwfc in enumerate(ecut_vals):
    if results[ecutwfc] is not None:
        ecutwfc_vals.append(ecutwfc)
        energy_vals.append(results[ecutwfc])

ecutwfc_vals = np.array(ecutwfc_vals)
energy_vals = np.array(energy_vals)

if len(ecutwfc_vals) > 2:
    print(f"\n✓ Found {len(ecutwfc_vals)} successful calculations")
    print(f"  ecutwfc range: {ecutwfc_vals.min():.1f} - {ecutwfc_vals.max():.1f} Ry")
    print(f"  Energy range: {energy_vals.min():.8f} - {energy_vals.max():.8f} eV")
    
    # Show data summary
    print(f"\n📋 DATA SUMMARY:")
    print(f"{'-'*70}")
    print(f"{'ecutwfc (Ry)':<15} {'Energy (eV)':<20} {'ΔE (meV)':<15}")
    print(f"{'-'*70}")
    for ecut, energy in zip(ecutwfc_vals, energy_vals):
        de = (energy - energy_vals[0]) * 1000
        print(f"{ecut:<15.1f} {energy:<20.8f} {de:<15.2f}")
    
    # Filter data by ecut_min
    if ecut_min is not None:
        mask = ecutwfc_vals >= ecut_min
        ecutwfc_vals_fit = ecutwfc_vals[mask]
        energy_vals_fit = energy_vals[mask]
        excluded_mask = ecutwfc_vals < ecut_min
        print(f"\n📌 Using ecut_min = {ecut_min} Ry for fit")
        print(f"Using {len(ecutwfc_vals_fit)}/{len(ecutwfc_vals)} points for fit")
        if np.any(excluded_mask):
            excluded_ecutwfc = ecutwfc_vals[excluded_mask]
            print(f"Excluded ecutwfc values: {excluded_ecutwfc.astype(int).tolist()}")
    else:
        ecutwfc_vals_fit = ecutwfc_vals
        energy_vals_fit = energy_vals
        excluded_mask = np.zeros_like(ecutwfc_vals, dtype=bool)
        print(f"\n✓ Using all {len(ecutwfc_vals_fit)} points for fit (ecut_min = None)")
    
    if len(ecutwfc_vals_fit) > 2:
        # Define exponential decay function: E(x) = E_inf + A * exp(-B * x)
        def exponential_decay(x, E_inf, A, B):
            return E_inf + A * np.exp(-B * x)
        
        # Initial parameter guess
        E_inf_guess = energy_vals_fit[-1]  # Last (highest ecutwfc) value
        A_guess = energy_vals_fit[0] - E_inf_guess
        B_guess = 0.05
        
        try:
            # Fit exponential decay
            popt, pcov = curve_fit(
                exponential_decay, 
                ecutwfc_vals_fit, 
                energy_vals_fit,
                p0=[E_inf_guess, A_guess, B_guess],
                maxfev=10000
            )
            
            E_inf, A, B = popt
            
            # Calculate R² (goodness of fit)
            residuals = energy_vals_fit - exponential_decay(ecutwfc_vals_fit, *popt)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((energy_vals_fit - np.mean(energy_vals_fit))**2)
            r_squared = 1 - (ss_res / ss_tot)
            
            print(f"\n📊 EXPONENTIAL FIT RESULTS")
            print(f"{'-'*70}")
            print(f"Function: E(ecutwfc) = E_inf + A * exp(-B * ecutwfc)")
            print(f"\nFitted parameters:")
            print(f"  E_inf (asymptotic energy) = {E_inf:.8f} eV")
            print(f"  A (amplitude)             = {A:.8f} eV")
            print(f"  B (decay constant)        = {B:.6f} Ry⁻¹")
            print(f"  R² (goodness of fit)      = {r_squared:.6f}")
            
            # Extrapolate to higher ecutwfc
            ecutwfc_predict = np.array([100, 150, 200])
            E_predict = exponential_decay(ecutwfc_predict, *popt)
            print(f"\n🔮 Extrapolations:")
            print(f"  ecutwfc = 100 Ry  →  E = {E_predict[0]:.8f} eV  (ΔE = {(E_predict[0]-E_inf)*1000:+.2f} meV)")
            print(f"  ecutwfc = 150 Ry  →  E = {E_predict[1]:.8f} eV  (ΔE = {(E_predict[1]-E_inf)*1000:+.2f} meV)")
            print(f"  ecutwfc = 200 Ry  →  E = {E_predict[2]:.8f} eV  (ΔE = {(E_predict[2]-E_inf)*1000:+.2f} meV)")
            
            # Convergence estimate: find ecutwfc where ΔE < 1 meV
            tolerance_meV = 1.0
            # A * exp(-B * ecut_conv) = tolerance / 1000
            ecut_conv = -np.log(tolerance_meV / 1000 / abs(A)) / B if A != 0 else ecutwfc_vals_fit.max()
            print(f"\n✓ Estimated ecutwfc for ΔE < {tolerance_meV} meV: {ecut_conv:.1f} Ry")
            
            # Create convergence plot
            # Determine plot range based on data used in fit
            ecut_max_fit = ecutwfc_vals_fit.max()
            ecut_min_fit = ecutwfc_vals_fit.min()
            ecut_plot_max = max(150, ecut_max_fit * 1.5)  # Extrapolate a bit
            ecutwfc_smooth = np.linspace(ecut_min_fit, ecut_plot_max, 200)
            E_smooth = exponential_decay(ecutwfc_smooth, *popt)
            
            plt.figure(figsize=(10, 6))
            
            # Plot data points (used in fit)
            plt.scatter(ecutwfc_vals_fit, (energy_vals_fit - E_inf) * 1000, 
                       s=100, color='red', label='Data points (used in fit)', zorder=5)
            
            # Plot fit
            plt.plot(ecutwfc_smooth, (E_smooth - E_inf) * 1000, 
                    'b-', linewidth=2, label=f'Exponential fit (R²={r_squared:.4f})')
            
            # Plot tolerance line
            plt.axhline(y=tolerance_meV, color='green', linestyle='--', 
                       linewidth=1.5, label=f'Tolerance: {tolerance_meV} meV')
            plt.axhline(y=-tolerance_meV, color='green', linestyle='--', linewidth=1.5)
            
            # Auto-scale with some padding
            plt.margins(x=0.05, y=0.1)
            
            # Formatting
            plt.xlabel('ecutwfc (Ry)', fontsize=12)
            plt.ylabel('ΔE = E - E_inf (meV)', fontsize=12)
            plt.title('Au Bulk: ecutwfc Convergence with Exponential Fit', fontsize=13)
            plt.grid(True, alpha=0.3)
            plt.legend(fontsize=10)
            plt.tight_layout()
            
            # Save plot
            plt.savefig('au_test/convergence_fit.png', dpi=150)
            print(f"\n📈 Convergence plot saved: au_test/convergence_fit.png")
            plt.show()
            
        except Exception as e:
            print(f"\n✗ Fit failed: {e}")
            print(f"  Trying with fewer data points or check data quality")
    else:
        print(f"\n⚠ Not enough data points for fit after exclusion (need ≥3, got {len(ecutwfc_vals_fit)})")
else:
    print(f"\n⚠ Not enough data points for fit (need ≥3, got {len(ecutwfc_vals)})")

# ============================================================================
# Performance Analysis
# ============================================================================

print("\n" + "="*80)
print("PERFORMANCE ANALYSIS")
print("="*80)

print(f"\n{'Metric':<40} {'Time (seconds)':<15}")
print("-" * 55)
print(f"{'Submission phase (non-blocking)':<40} {submit_time:<15.2f}")
print(f"{'Monitoring phase (parallel)':<40} {wait_time:<15.2f}")
print(f"{'Total time':<40} {total_time:<15.2f}")

# Estimate serial time (assuming ~5 min per ecutwfc)
serial_estimate = len(ecut_vals) * 300
speedup = serial_estimate / total_time
print(f"\nEstimated serial time: {serial_estimate:.0f} seconds ({serial_estimate/60:.1f} min)")
print(f"Batch time: {total_time:.0f} seconds ({total_time/60:.1f} min)")
print(f"Speedup: {speedup:.1f}×")

print("\n" + "="*80)
print("Results saved in: au_test/ directory")
print("="*80)
