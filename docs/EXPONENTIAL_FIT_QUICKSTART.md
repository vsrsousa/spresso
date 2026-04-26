# Exponential Fit Integration: Quick Start Guide

## What Changed?

**ConvergenceWorkflow Phase 1 now uses exponential fit analysis instead of hard reference at ecutwfc=200 Ry**

## Why Should I Care?

### ❌ Old Method (Legacy)
```
Phase 1 selection:
  Compare all energies vs ecutwfc=200 Ry
  ↓
  "Is this point within tolerance of ecutwfc=200?"
  ↓
  Select minimum ecutwfc that passes this test
  
Problems:
- Uses only ONE test point as reference
- Ignores all other measurements
- No way to know if reference is truly "converged"
```

### ✅ New Method (Exponential Fit)
```
Phase 1 selection:
  Fit E(ecutwfc) = E_inf + A·exp(-B·ecutwfc) to ALL data
  ↓
  Extract E_inf (true asymptotic energy)
  ↓
  Compare vs E_inf instead of vs ecutwfc=200
  ↓
  Select minimum ecutwfc that meets tolerance
  
Benefits:
- Uses ALL test data simultaneously
- E_inf is the true asymptotic limit
- Provides R² quality metric
- Can extrapolate beyond tested range
```

## Quick Example

### Before
```python
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow
from ase.build import bulk

atoms = bulk('Au', 'fcc', a=4.0782)
wf = ConvergenceWorkflow(atoms, 
    pseudopotentials={'Au': 'Au.pbe-n-rrkjus_psl.1.0.0.UPF'},
    precision='low')
wf.run_convergence_study()
rec = wf.get_recommendations()

print(f"Optimal ecutwfc: {rec['optimal_ecutwfc']} Ry")
print(f"Optimal kspacing: {rec['optimal_kspacing']} Å⁻¹")
# That's it!
```

### After (with fit analysis)
```python
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow
from ase.build import bulk

atoms = bulk('Au', 'fcc', a=4.0782)
wf = ConvergenceWorkflow(atoms, 
    pseudopotentials={'Au': 'Au.pbe-n-rrkjus_psl.1.0.0.UPF'},
    precision='low')
wf.run_convergence_study()
rec = wf.get_recommendations()

print(f"Optimal ecutwfc: {rec['optimal_ecutwfc']} Ry")
print(f"Optimal kspacing: {rec['optimal_kspacing']} Å⁻¹")

# NEW: Exponential fit analysis!
if 'exponential_fit' in rec:
    fit = rec['exponential_fit']
    print(f"\n📊 Fit Analysis:")
    print(f"  E_inf = {fit['E_inf']:.8f} eV")
    print(f"  R² = {fit['R_squared']:.6f}")
    print(f"  Min ecutwfc for tolerance = {fit['min_ecutwfc_for_tolerance']:.1f} Ry")
```

## Understanding the Output

### New Recommendation Structure

```python
rec = wf.get_recommendations(verbose=False)

# Always available
rec['optimal_ecutwfc']           # ← Min ecutwfc that meets tolerance
rec['optimal_kspacing']          # ← Min kspacing that meets tolerance
rec['precision']                 # ← Level used (low/medium/high/ultra)
rec['energy_tolerance_meV_atom'] # ← Target tolerance

# NEW: Only if fit succeeded
if 'exponential_fit' in rec:
    fit = rec['exponential_fit']
    
    fit['E_inf']                        # Extrapolated asymptotic energy
    fit['A']                            # Exponential amplitude
    fit['B']                            # Decay constant
    fit['R_squared']                    # Goodness of fit (0 = bad, 1 = perfect)
    fit['min_ecutwfc_for_tolerance']    # Estimated min ecutwfc for tolerance
    fit['tolerance_meV']                # Target tolerance in meV
    fit['method']                       # = 'exponential_decay'
```

## Interpreting Fit Results

### R² (Goodness of Fit)
```
R² = 1.0  ✅ Perfect fit - extrapolations are highly reliable
R² > 0.99 ✅ Excellent fit - extrapolations are very reliable  
R² > 0.95 ✓ Good fit - extrapolations are reasonably reliable
R² < 0.95 ⚠ Poor fit - use extrapolations with caution
```

### E_inf (Asymptotic Energy)
```
E_inf = -19.253456 eV

This is the energy value as ecutwfc → ∞
Used as reference instead of ecutwfc=200 Ry
More physically meaningful than any single test point
```

### Interpretation Examples

```python
rec = wf.get_recommendations(verbose=False)
fit = rec['exponential_fit']

# Example 1: Excellent fit with safety margin
if fit['R_squared'] > 0.99 and rec['optimal_ecutwfc'] > fit['min_ecutwfc_for_tolerance']:
    margin = (rec['optimal_ecutwfc'] - fit['min_ecutwfc_for_tolerance']) / fit['min_ecutwfc_for_tolerance'] * 100
    print(f"✅ Tested ecutwfc is {margin:.1f}% above minimum")
    print(f"   Provides good safety margin")

# Example 2: Borderline fit
elif fit['R_squared'] < 0.95:
    print(f"⚠ Fit quality is marginal (R²={fit['R_squared']:.3f})")
    print(f"   Recommendation is still valid (tested, not extrapolated)")

# Example 3: Could use lower ecutwfc
elif rec['optimal_ecutwfc'] ≈ fit['min_ecutwfc_for_tolerance']:
    print(f"💡 Tested ecutwfc is very close to minimum required")
    print(f"   Could potentially use lower value")
```

## When Does Exponential Fit Apply?

### ✅ Fits Well
- Ecutwfc convergence (typical R² = 0.999+)
- Smooth exponential-like decay behavior
- 5+ test points
- No spurious measurements

### ⚠ Fits Okay
- 3-4 test points
- Some noise in measurements
- Mostly exponential behavior

### ❌ Falls Back to Legacy
- <3 test points
- Non-exponential behavior
- If fit fails for any reason

## Backward Compatibility

**Don't worry - everything still works!**

```python
# Old code without fit analysis
rec = wf.get_recommendations(verbose=False)
ecut = rec['optimal_ecutwfc']
ksp = rec['optimal_kspacing']
# Works exactly as before!

# New code with fit analysis
if 'exponential_fit' in rec:
    # Use fit information if available
    fit_info = rec['exponential_fit']
else:
    # Fit not available (fallback used)
    pass
```

## Advanced: Extrapolate Your Own Values

```python
import numpy as np

rec = wf.get_recommendations(verbose=False)
fit = rec['exponential_fit']

if fit:
    E_inf = fit['E_inf']
    A = fit['A']
    B = fit['B']
    
    # Predict energy for ecutwfc = 150 Ry (without testing!)
    ecut = 150
    E_150 = E_inf + A * np.exp(-B * ecut)
    delta_E = abs(E_150 - E_inf) * 1000  # meV
    
    print(f"Predicted for ecutwfc={ecut} Ry:")
    print(f"  E = {E_150:.8f} eV")
    print(f"  ΔE = {delta_E:.2f} meV from asymptotic")
    
    # For which ecutwfc do you get 0.5 meV?
    target_tolerance = 0.0005  # eV
    ecut_for_tolerance = -np.log(target_tolerance / abs(A)) / B
    print(f"\nEcutwfc needed for 0.5 meV tolerance: {ecut_for_tolerance:.1f} Ry")
```

## Documentation

For more details:
- **Visual comparison**: `docs/EXPONENTIAL_FIT_INTEGRATION.md`
- **Workflow diagrams**: `docs/PHASE1_EXPONENTIAL_FIT_WORKFLOW.md`
- **Full summary**: `docs/EXPONENTIAL_FIT_IMPLEMENTATION_SUMMARY.md`
- **Working example**: `examples/convergence_workflow_with_exponential_fit.py`
- **Code reference**: `xespresso/workflow/convergence_workflow.py` (search for `_fit_exponential_decay_phase1`)

## FAQs

### Q: Why don't I see fit information in my output?
**A:** The fit may have failed (e.g., R² threshold not met). Check:
```python
rec = wf.get_recommendations(verbose=True)  # Shows why
if 'exponential_fit' not in rec:
    print("Fit unavailable, fell back to legacy method")
```

### Q: Can I force use of legacy method?
**A:** Not directly, but if fit fails, legacy method is used automatically.

### Q: Are the `optimal_ecutwfc` values different?
**A:** Possibly! The new method uses E_inf, the old used ecutwfc=200. Results may differ slightly.

### Q: Which is more accurate?
**A:** The exponential fit is more physically meaningful because:
- E_inf is the true asymptotic limit
- Uses all data, not just one point
- Validated by R² metric

### Q: Can I use fit to avoid testing at ecutwfc=200?
**A:** Yes! If fit is good (R² > 0.99), you could extrapolate.
But testing provides validation, so continue testing range recommended.

---

**That's it! The exponential fit integration is transparent and automatic.**
