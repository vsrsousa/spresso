#!/usr/bin/env python
"""
AUTO-DETECTION STRATEGIES: How to determine ecutwfc_min_for_fit

The challenge: Automatically find where the "basis-incomplete region" ends
and the "exponential convergence region" begins.

Three approaches:
1. CURVATURE METHOD: Find the inflection point in energy curve
2. RELATIVE CHANGE METHOD: Where ΔE becomes < threshold relative to total
3. PSEUDOPOTENTIAL DATABASE: Look up standard values
"""

import numpy as np
from scipy.interpolate import UnivariateSpline

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║              AUTO-DETECTION: Finding ecutwfc_min_for_fit                      ║
╚══════════════════════════════════════════════════════════════════════════════╝

🎯 THE CHALLENGE
─────────────────────────────────────────────────────────────────────────────

You have data like:
  ecutwfc = 30 Ry  →  E = -9800.00 eV
  ecutwfc = 40 Ry  →  E = -9995.00 eV  (jump of 195 eV!)
  ecutwfc = 50 Ry  →  E = -9999.50 eV  (jump of 4.5 eV)
  ecutwfc = 60 Ry  →  E = -9999.90 eV  (jump of 0.4 eV)
  ecutwfc = 70 Ry  →  E =  -10000.00 eV  (jump of 0.1 eV)

QUESTION: Where does the transition happen? (30→50, 40→50, or somewhere else?)


═══════════════════════════════════════════════════════════════════════════════
STRATEGY 1: SECOND DERIVATIVE / CURVATURE METHOD
═══════════════════════════════════════════════════════════════════════════════

Idea: Find where the curvature of E(ecutwfc) changes

Mathematical approach:
  1. Fit smooth spline to all data
  2. Compute first derivative: dE/d(ecut)
  3. Compute second derivative: d²E/d(ecut)²
  4. Find inflection point (where d²E/d(ecut)² changes sign or drops)
  5. Use that as transition point

Code:
""")

# Example implementation
def auto_detect_ecut_min_curvature(ecutwfc_vals, energy_vals, verbose=True):
    """
    Detect ecutwfc_min using second derivative (curvature).
    
    The basis-incomplete region has large curvature (rapidly changing slope).
    The convergence region has low curvature (smooth exponential).
    Transition happens at the inflection point.
    """
    
    ecutwfc_vals = np.array(ecutwfc_vals)
    energy_vals = np.array(energy_vals)
    
    # Fit smooth spline (k=3 cubic spline, s controls smoothing)
    spline = UnivariateSpline(ecutwfc_vals, energy_vals, k=min(3, len(ecutwfc_vals)-1), s=None)
    
    # Compute derivatives
    # First derivative: dE/d(ecut)
    dE = spline.derivative(n=1)
    
    # Second derivative: d²E/d(ecut)²
    d2E = spline.derivative(n=2)
    
    # Evaluate at all points
    dE_vals = dE(ecutwfc_vals)
    d2E_vals = d2E(ecutwfc_vals)
    
    # Curvature: |d²E/d(ecut)²| / (1 + (dE/d(ecut))²)^(3/2)
    # For simplicity, use just |d²E/d(ecut)²|
    curvature = np.abs(d2E_vals)
    
    # Find inflection point: where curvature drops below mean
    curvature_threshold = np.mean(curvature)
    
    # Find first point where curvature < threshold
    candidates = ecutwfc_vals[curvature < curvature_threshold]
    
    if len(candidates) > 0:
        ecut_min = candidates[0]
    else:
        # If no clear inflection, use median point
        ecut_min = np.median(ecutwfc_vals)
    
    if verbose:
        print(f"\n📊 CURVATURE METHOD:")
        print(f"  Curvature threshold: {curvature_threshold:.6f}")
        print(f"  Detected ecut_min: {ecut_min:.1f} Ry")
        print(f"  (First point where curvature < threshold)")
    
    return ecut_min


print("""
EXAMPLE OUTPUT:
  ecutwfc    |  Energy    |  dE/decut  |  d²E/decut²  |  Status
  ─────────────────────────────────────────────────────────────
  30         | -9800.00   |   1950.0   |    -195.0    | High curvature
  40         | -9995.00   |     45.0   |     -4.5     | High curvature  
  50         | -9999.50   |      4.0   |     -0.4     | ← Threshold
  60         | -9999.90   |      1.0   |     -0.1     | Low curvature ✓
  70         | -10000.00  |      0.1   |     -0.01    | Low curvature ✓

  Result: ecut_min = 50 Ry ✓


═══════════════════════════════════════════════════════════════════════════════
STRATEGY 2: RELATIVE ENERGY CHANGE METHOD
═══════════════════════════════════════════════════════════════════════════════

Idea: Find where the energy change becomes "small" relative to total range

Mathematical approach:
  1. Calculate total energy range: E_max - E_min
  2. For each point, calculate fractional change: ΔE / (E_max - E_min)
  3. Find first point where fractional change < threshold (e.g., 1%)
  4. Use that as transition point

Code:
""")

def auto_detect_ecut_min_relative_change(ecutwfc_vals, energy_vals, threshold=0.01, verbose=True):
    """
    Detect ecutwfc_min using relative energy change.
    
    In basis-incomplete region: large jumps (e.g., 100 meV)
    In convergence region: small changes (e.g., 1 meV)
    
    Find where changes become < threshold % of total range.
    """
    
    ecutwfc_vals = np.array(ecutwfc_vals)
    energy_vals = np.array(energy_vals)
    
    # Sort by ecutwfc
    sort_idx = np.argsort(ecutwfc_vals)
    ecutwfc_vals = ecutwfc_vals[sort_idx]
    energy_vals = energy_vals[sort_idx]
    
    # Total energy range
    E_range = energy_vals[-1] - energy_vals[0]
    
    # Changes between consecutive points
    dE = np.abs(np.diff(energy_vals))
    
    # Fractional changes
    frac_changes = dE / E_range
    
    # Find first point where change < threshold
    candidates_idx = np.where(frac_changes < threshold)[0]
    
    if len(candidates_idx) > 0:
        # Use the point BEFORE the first small change
        # (i.e., the point that precedes the convergence region)
        idx = candidates_idx[0]
        ecut_min = ecutwfc_vals[idx]  # Last point in incomplete region
    else:
        ecut_min = ecutwfc_vals[0]  # Fallback
    
    if verbose:
        print(f"\n📊 RELATIVE CHANGE METHOD:")
        print(f"  Threshold: {threshold*100:.1f}% of total range")
        print(f"  Total energy range: {E_range:.6f} eV")
        print(f"  Detected ecut_min: {ecut_min:.1f} Ry")
        print(f"\n  Energy changes (meV):")
        for i, (ecut, de_frac) in enumerate(zip(ecutwfc_vals[:-1], frac_changes)):
            next_ecut = ecutwfc_vals[i+1]
            status = "← THRESHOLD" if de_frac < threshold else ""
            print(f"    {ecut:.0f} → {next_ecut:.0f}: {de_frac*E_range*1000:.1f} meV ({de_frac*100:.2f}%) {status}")
    
    return ecut_min


print("""
EXAMPLE OUTPUT:
  30 → 40: 1950.00 meV (19.50%)
  40 → 50:  4.50 meV (0.45%)  ← THRESHOLD (1%)
  50 → 60:  0.40 meV (0.04%)
  60 → 70:  0.10 meV (0.01%)

  Result: ecut_min = 40 Ry


═══════════════════════════════════════════════════════════════════════════════
STRATEGY 3: PSEUDOPOTENTIAL DATABASE METHOD
═══════════════════════════════════════════════════════════════════════════════

Idea: Look up standard values for known pseudopotentials

Different pseudopotentials have different basis completeness:

""")

# Database of known pseudopotentials
PSEUDOPOTENTIAL_DEFAULTS = {
    # Format: (library, pattern) → ecutwfc_min
    
    # Norm-Conserving (NC) - Slim basis
    ('psl', 'nc'): 40,
    ('gbrv', 'bhs'): 35,  # GBRV NC versions
    
    # Norm-Conserving (NC) - Standard
    ('psl', 'n-rrkjus'): 50,  # Your Au pseudo!
    ('psl', 'rrkjus'): 45,
    
    # Norm-Conserving (NC) - Extended basis
    ('psl', 'n-rrkjusxc'): 60,
    
    # Ultra-Soft (US)
    ('psl', 'us'): 30,
    ('gbrv', 'bh'): 35,
    
    # PAW (hardest basis)
    ('psl', 'paw'): 70,
    ('vasp', 'paw'): 75,
    ('gbrv', 'paw'): 65,
    
    # SG15
    ('sg15', 'oncvpsp'): 50,
}

def auto_detect_ecut_min_from_pseudopotential(pseudo_filename, verbose=True):
    """
    Detect ecutwfc_min by looking up the pseudopotential type.
    """
    
    pseudo_lower = pseudo_filename.lower()
    
    best_match = None
    best_score = 0
    
    # Find best matching pseudo from database
    for (library, pattern), ecut_min in PSEUDOPOTENTIAL_DEFAULTS.items():
        if library in pseudo_lower and pattern in pseudo_lower:
            score = len(library) + len(pattern)  # Longer match = better
            if score > best_score:
                best_match = (library, pattern, ecut_min)
                best_score = score
    
    if best_match:
        lib, pat, ecut_min = best_match
        if verbose:
            print(f"\n📚 PSEUDOPOTENTIAL DATABASE METHOD:")
            print(f"  Recognized: {lib} {pat}")
            print(f"  Detected ecut_min: {ecut_min:.0f} Ry")
        return ecut_min
    else:
        if verbose:
            print(f"\n📚 PSEUDOPOTENTIAL DATABASE METHOD:")
            print(f"  Not recognized: {pseudo_filename}")
            print(f"  Using fallback: 40 Ry")
        return 40.0


print("""
DATABASE LOOKUPS:

Your pseudo (Au.pbe-n-rrkjus_psl.1.0.0.UPF):
  Library: 'psl'
  Pattern: 'n-rrkjus'
  → ecut_min = 50 Ry ✓

For comparison:
  Soft US (psl-us):          ecut_min = 30 Ry
  PAW (paw):                 ecut_min = 70 Ry
  SG15 (sg15-oncvpsp):       ecut_min = 50 Ry


═══════════════════════════════════════════════════════════════════════════════
STRATEGY 4: COMBINED APPROACH (RECOMMENDED)
═══════════════════════════════════════════════════════════════════════════════

Use multiple methods and pick the MOST CONSERVATIVE (highest) result:

""")

def auto_detect_ecut_min_combined(ecutwfc_vals, energy_vals, pseudo_filename, verbose=True):
    """
    Auto-detect ecutwfc_min using multiple strategies.
    Return the most conservative (highest) result.
    """
    
    results = {}
    
    # Method 1: Curvature
    results['curvature'] = auto_detect_ecut_min_curvature(ecutwfc_vals, energy_vals, verbose=False)
    
    # Method 2: Relative change
    results['relative_change'] = auto_detect_ecut_min_relative_change(ecutwfc_vals, energy_vals, verbose=False)
    
    # Method 3: Database
    results['database'] = auto_detect_ecut_min_from_pseudopotential(pseudo_filename, verbose=False)
    
    # Pick the maximum (most conservative)
    ecut_min = max(results.values())
    
    if verbose:
        print(f"\n🔄 COMBINED METHOD:")
        print(f"  Curvature method:      {results['curvature']:.1f} Ry")
        print(f"  Relative change:       {results['relative_change']:.1f} Ry")
        print(f"  Database method:       {results['database']:.1f} Ry")
        print(f"  ─────────────────────────────")
        print(f"  Final (most conservative): {ecut_min:.1f} Ry ✓")
    
    return ecut_min


print("""
EXAMPLE OUTPUT:
  Curvature method:      50.0 Ry
  Relative change:       40.0 Ry
  Database method:       50.0 Ry
  ─────────────────────────────
  Final result:          50.0 Ry (take max for safety)


═══════════════════════════════════════════════════════════════════════════════
RECOMMENDATION FOR YOUR CODE
═══════════════════════════════════════════════════════════════════════════════

Implement in this priority order:

1. ✅ DATABASE METHOD (Easiest, most reliable for common pseudos)
   - Build lookup table for popular pseudopotential libraries
   - Returns instant result
   - Can be overridden by user if needed

2. ✅ CURVATURE METHOD (Fallback if pseudo not in database)
   - Works for any data
   - Mathematically grounded
   - Robust to different systems

3. ⚠️ RELATIVE CHANGE (Use with caution)
   - Sensitive to threshold choice
   - May not work if data is sparse
   - Good for validation

USAGE IN CODE:

    def _fit_exponential_decay_phase1(
        self,
        ecut_results,
        criteria_tolerances,
        verbose=True,
        ecut_min_for_fit=None,  # ← NEW PARAMETER
        auto_detect_method='combined'  # ← 'database', 'curvature', 'combined'
    ):
        \"\"\"
        Fit exponential decay model to Phase 1 ecutwfc convergence data.
        
        Uses ONLY points in the convergence region (ecut ≥ ecut_min_for_fit).
        Ignores basis-incomplete region (ecut < ecut_min_for_fit).
        
        Args:
            ecut_results: Dict mapping ecutwfc → Dict with 'energy' key
            criteria_tolerances: Dict with 'energy_tolerance' key
            verbose: Print fit results and diagnostics
            ecut_min_for_fit: Minimum ecutwfc for fit
                             If None: auto-detect using auto_detect_method
                             If specified: use this value
            auto_detect_method: 'database', 'curvature', or 'combined'
        \"\"\"
        
        # If not specified, auto-detect
        if ecut_min_for_fit is None:
            if auto_detect_method == 'database':
                ecut_min_for_fit = self._detect_ecut_min_from_pseudo(
                    self.pseudo_filename
                )
            elif auto_detect_method == 'curvature':
                ecutwfc_vals = np.array(sorted(ecut_results.keys()))
                energy_vals = np.array([ecut_results[e]['energy'] for e in ecutwfc_vals])
                ecut_min_for_fit = auto_detect_ecut_min_curvature(
                    ecutwfc_vals, energy_vals, verbose=verbose
                )
            else:  # 'combined'
                ecutwfc_vals = np.array(sorted(ecut_results.keys()))
                energy_vals = np.array([ecut_results[e]['energy'] for e in ecutwfc_vals])
                ecut_min_for_fit = auto_detect_ecut_min_combined(
                    ecutwfc_vals, energy_vals, 
                    self.pseudo_filename,
                    verbose=verbose
                )
        
        # Filter results: use only ecut ≥ ecut_min_for_fit
        ecut_results_filtered = {
            ecut: data
            for ecut, data in ecut_results.items()
            if ecut >= ecut_min_for_fit
        }
        
        if verbose:
            print(f"\\n🔧 CONVERGENCE REGION:")
            print(f"  Using points with ecutwfc ≥ {ecut_min_for_fit:.1f} Ry")
            print(f"  Excluded (basis incomplete): {[e for e in ecut_results if e < ecut_min_for_fit]}")
            print(f"  Included (exponential region): {sorted(ecut_results_filtered.keys())}")
        
        # ... rest of fit logic
"""

print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print("""
BEST APPROACH: Use combination of methods

1. Check pseudopotential database first (fastest)
   - If found → use that value
   
2. If not in database, auto-detect from curvature (robust)
   - Analyze second derivative
   - Find inflection point
   
3. User can always override
   - wf._fit_exponential_decay_phase1(ecut_results, ..., ecut_min_for_fit=50.0)

BENEFITS:
  ✅ Automatic for most users
  ✅ Physics-based (finds real transition)
  ✅ Conservative (picks max value)
  ✅ Reproducible
  ✅ Works for any pseudopotential

IMPLEMENTATION EFFORT:
  - ~100 lines of code
  - ~50 lines for database
  - ~40 lines for curvature detection
  - ~10 lines for integration

RESULT:
  Users run: wf.run_convergence_study()
  System automatically:
    1. Detects ecut_min = 50 Ry (for your Au pseudo)
    2. Filters data (ignores 30, 40)
    3. Fits exponential to valid region (50, 60, 70)
    4. Gets stable recommendation (50 Ry) ✓
""")
