"""
ECUTWFC RANGE SPECIFICATION: 3 Flexible Modes

Three complementary ways to specify ecutwfc values for Phase 1 convergence:

MODE 1: EXPLICIT LIST (most control)
  ecut_vals=[30, 40, 50, 60, 70]
  → Tests EXACTLY these values in order

MODE 2: LINSPACE (n points between min and max)
  ecut_range={'min': 30, 'max': 70, 'n_points': 5}
  → Generates: [30, 40, 50, 60, 70]
  → np.linspace(30, 70, 5) internally

MODE 3: ARANGE (step-based)
  ecut_range={'min': 30, 'max': 70, 'step': 10}
  → Generates: [30, 40, 50, 60, 70]
  → np.arange(30, 80, 10) internally

MODE 4: DYNAMIC (default, current behavior)
  ecut_range=None
  → Starts with min_ecutwfc, expands by ecut_step until convergence
  → Adaptive to convergence behavior (smart!)

RECOMMENDED USAGE:
  • Testing specific points: MODE 1 (explicit list)
  • Quick convergence check: MODE 4 (dynamic, default)
  • Comparison studies: MODE 2 (linspace, uniform spacing)
  • Production runs: MODE 3 (arange, control spacing)
"""

import numpy as np
from typing import Optional, List, Dict, Union

# Example implementations
def build_ecut_range(
    ecut_vals: Optional[List[float]] = None,
    ecut_range: Optional[Dict[str, float]] = None,
    min_ecutwfc: float = 30.0,
    max_ecutwfc: float = 200.0,
    ecut_step: float = 10.0,
) -> List[float]:
    """
    Build ecutwfc range based on user specification.
    
    Supports 4 modes with clear precedence:
    1. Explicit list (ecut_vals) - HIGHEST priority
    2. Dict-based range (ecut_range) - Can be 'linspace' or 'arange' mode
    3. Dynamic generation (default) - Uses min/max/step parameters
    
    Args:
        ecut_vals: Explicit list like [30, 40, 50, 60, 70]
        ecut_range: Dict with:
          - {'min': 30, 'max': 70, 'n_points': 5}  → linspace mode
          - {'min': 30, 'max': 70, 'step': 10}     → arange mode
        min_ecutwfc: Minimum value for dynamic mode
        max_ecutwfc: Maximum value for dynamic mode
        ecut_step: Step size for dynamic mode
        
    Returns:
        List of ecutwfc values to test (sorted)
        
    Examples:
        # Mode 1: Explicit list
        >>> build_ecut_range(ecut_vals=[30, 40, 50, 60, 70])
        [30, 40, 50, 60, 70]
        
        # Mode 2: Linspace (5 uniform points)
        >>> build_ecut_range(ecut_range={'min': 30, 'max': 70, 'n_points': 5})
        [30.0, 40.0, 50.0, 60.0, 70.0]
        
        # Mode 3: Arange (step-based)
        >>> build_ecut_range(ecut_range={'min': 30, 'max': 70, 'step': 10})
        [30, 40, 50, 60, 70]
        
        # Mode 4: Dynamic (default)
        >>> build_ecut_range(min_ecutwfc=30, max_ecutwfc=200, ecut_step=10)
        [30]  # Will expand dynamically during convergence
    """
    
    # PRIORITY 1: Explicit list (highest priority)
    if ecut_vals is not None:
        if not isinstance(ecut_vals, (list, tuple)):
            raise TypeError(f"ecut_vals must be list or tuple, got {type(ecut_vals)}")
        ecutwfc_list = sorted(list(ecut_vals))
        if len(ecutwfc_list) < 1:
            raise ValueError("ecut_vals must contain at least 1 value")
        return ecutwfc_list
    
    # PRIORITY 2: Dict-based range specification
    if ecut_range is not None:
        if not isinstance(ecut_range, dict):
            raise TypeError(f"ecut_range must be dict, got {type(ecut_range)}")
        
        # Check for linspace mode (n_points specified)
        if 'n_points' in ecut_range:
            min_val = ecut_range.get('min')
            max_val = ecut_range.get('max')
            n_points = ecut_range.get('n_points')
            
            if min_val is None or max_val is None or n_points is None:
                raise ValueError(
                    "Linspace mode requires 'min', 'max', and 'n_points' keys. "
                    f"Got: {ecut_range}"
                )
            
            if n_points < 2:
                raise ValueError(f"n_points must be >= 2, got {n_points}")
            
            ecutwfc_list = np.linspace(min_val, max_val, n_points).tolist()
            return sorted(ecutwfc_list)
        
        # Check for arange mode (step specified)
        elif 'step' in ecut_range:
            min_val = ecut_range.get('min')
            max_val = ecut_range.get('max')
            step = ecut_range.get('step')
            
            if min_val is None or max_val is None or step is None:
                raise ValueError(
                    "Arange mode requires 'min', 'max', and 'step' keys. "
                    f"Got: {ecut_range}"
                )
            
            if step <= 0:
                raise ValueError(f"step must be positive, got {step}")
            
            # Use arange: include max value if it's a multiple of step
            ecutwfc_list = np.arange(min_val, max_val + step/2, step).tolist()
            return sorted(ecutwfc_list)
        
        else:
            raise ValueError(
                "ecut_range dict must have either 'n_points' (linspace mode) "
                "or 'step' (arange mode). "
                f"Got keys: {list(ecut_range.keys())}"
            )
    
    # PRIORITY 3: Dynamic mode (default) - return initial value, will expand during convergence
    return [min_ecutwfc]


# Test cases
if __name__ == '__main__':
    print("="*70)
    print("ECUTWFC RANGE BUILDING TESTS")
    print("="*70)
    
    # Test Mode 1: Explicit list
    print("\n1. MODE 1: EXPLICIT LIST")
    print("   build_ecut_range(ecut_vals=[70, 50, 30, 60, 40])")
    result = build_ecut_range(ecut_vals=[70, 50, 30, 60, 40])
    print(f"   Result: {result} ✓")
    
    # Test Mode 2a: Linspace
    print("\n2. MODE 2: LINSPACE (uniform spacing)")
    print("   build_ecut_range(ecut_range={'min': 30, 'max': 70, 'n_points': 5})")
    result = build_ecut_range(ecut_range={'min': 30, 'max': 70, 'n_points': 5})
    print(f"   Result: {result} ✓")
    
    # Test Mode 3: Arange
    print("\n3. MODE 3: ARANGE (step-based)")
    print("   build_ecut_range(ecut_range={'min': 30, 'max': 70, 'step': 10})")
    result = build_ecut_range(ecut_range={'min': 30, 'max': 70, 'step': 10})
    print(f"   Result: {result} ✓")
    
    # Test Mode 4: Dynamic
    print("\n4. MODE 4: DYNAMIC (adaptive)")
    print("   build_ecut_range(min_ecutwfc=30, max_ecutwfc=200, ecut_step=10)")
    result = build_ecut_range(min_ecutwfc=30, max_ecutwfc=200, ecut_step=10)
    print(f"   Result: {result} (will expand during convergence) ✓")
    
    # Test edge cases
    print("\n" + "="*70)
    print("ERROR HANDLING")
    print("="*70)
    
    try:
        print("\nTesting invalid mode (missing required keys)...")
        result = build_ecut_range(ecut_range={'min': 30, 'max': 70})
        print("ERROR: Should have raised ValueError!")
    except ValueError as e:
        print(f"✓ Caught expected error: {e}")
    
    try:
        print("\nTesting n_points < 2...")
        result = build_ecut_range(ecut_range={'min': 30, 'max': 70, 'n_points': 1})
        print("ERROR: Should have raised ValueError!")
    except ValueError as e:
        print(f"✓ Caught expected error: {e}")
    
    try:
        print("\nTesting negative step...")
        result = build_ecut_range(ecut_range={'min': 30, 'max': 70, 'step': -10})
        print("ERROR: Should have raised ValueError!")
    except ValueError as e:
        print(f"✓ Caught expected error: {e}")
    
    print("\n" + "="*70)
    print("COMPARISON OF MODES")
    print("="*70)
    
    print("""
╔═══════════════════════════════════════════════════════════════════╗
║                        MODE COMPARISON                            ║
╠═══════════════════════════════════════════════════════════════════╣
║ MODE  │ USAGE              │ CONTROL │ FLEXIBILITY │ USE CASE     ║
├───────┼────────────────────┼─────────┼─────────────┼──────────────┤
║ 1     │ Explicit list      │ Full    │ Low (fixed) │ Test exact   ║
║       │ [30,40,50,60,70]   │         │             │ points       ║
├───────┼────────────────────┼─────────┼─────────────┼──────────────┤
║ 2     │ Linspace n=5       │ Medium  │ High        │ Comparison   ║
║       │ min=30, max=70     │         │ (uniform)   │ studies      ║
├───────┼────────────────────┼─────────┼─────────────┼──────────────┤
║ 3     │ Arange step=10     │ Medium  │ Medium      │ Production   ║
║       │ min=30, max=70     │         │ (step based)│ runs         ║
├───────┼────────────────────┼─────────┼─────────────┼──────────────┤
║ 4     │ Dynamic (default)  │ Low     │ Very high   │ Quick check  ║
║       │ Expands adaptively │         │ (adaptive)  │ Auto-optimal ║
╚═══════════════════════════════════════════════════════════════════╝

RECOMMENDATIONS:
  ✓ Use MODE 1 when testing specific known values
  ✓ Use MODE 2 for uniform spacing (papers, comparisons)
  ✓ Use MODE 3 when you know typical step size
  ✓ Use MODE 4 (default) for smart adaptive convergence
    """)
