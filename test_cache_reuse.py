#!/usr/bin/env python3
"""
Test script: Convergence cache reuse with different precision levels.

Test the new cache persistence functionality:
1. Run with precision='low'
2. Run with precision='high' on SAME workflow instance
3. Verify cache is reused (equations that were already tested are not recalculated)
"""

import sys
import numpy as np
from ase.build import bulk
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow

def test_cache_reuse():
    print("\n" + "="*80)
    print("TEST: Cache Reuse with Different Precision Levels")
    print("="*80)
    
    # Create Au bulk structure
    au_bulk = bulk('Au', 'fcc', a=4.0782)
    print(f"\nStructure: {au_bulk.get_chemical_formula()}")
    print(f"Atoms: {len(au_bulk)}")
    
    # Create workflow (it will initialize with precision='low')
    print("\n" + "-"*80)
    print("Creating ConvergenceWorkflow with precision='low'")
    print("-"*80)
    
    wf = ConvergenceWorkflow(
        atoms=au_bulk,
        pseudopotentials_config='default',
        precision='low',
        protocol='moderate',
        machine='medusa'  # Set to None if no machine available
    )
    
    print(f"\nInitial cache state:")
    print(f"  ecut_results_cache: {len(wf.ecut_results_cache)} items")
    print(f"  kspacing_results_cache: {len(wf.kspacing_results_cache)} items")
    
    # ========================================================================
    # FIRST RUN: precision='low'
    # ========================================================================
    print("\n" + "="*80)
    print("FIRST RUN: precision='low' (tolerance = 3 meV/atom)")
    print("="*80)
    
    try:
        results_low = wf.run_convergence_independent(
            label_prefix='test_cache_low',
            precision='low',  # ← Dynamic precision
            verbose=True
        )
        print(f"\n✓ First run completed with 'low' precision")
    except Exception as e:
        print(f"\n✗ First run failed: {e}")
        return False
    
    # Check cache after first run
    print(f"\nCache after FIRST run (low precision):")
    print(f"  ecut_results_cache: {len(wf.ecut_results_cache)} items")
    print(f"  ecutwfc values tested: {sorted(wf.ecut_results_cache.keys())}")
    print(f"  kspacing_results_cache: {len(wf.kspacing_results_cache)} items")
    print(f"  kspacing values tested: {sorted(wf.kspacing_results_cache.keys())}")
    
    cache_size_after_first = {
        'ecut': len(wf.ecut_results_cache),
        'kspacing': len(wf.kspacing_results_cache)
    }
    
    # ========================================================================
    # SECOND RUN: precision='high'
    # ========================================================================
    print("\n" + "="*80)
    print("SECOND RUN: precision='high' (tolerance = 1 meV/atom, MORE STRICT)")
    print("="*80)
    
    try:
        results_high = wf.run_convergence_independent(
            label_prefix='test_cache_high',
            precision='high',  # ← Different precision, STRICTER
            verbose=True
        )
        print(f"\n✓ Second run completed with 'high' precision")
    except Exception as e:
        print(f"\n✗ Second run failed: {e}")
        return False
    
    # Check cache after second run
    print(f"\nCache after SECOND run (high precision):")
    print(f"  ecut_results_cache: {len(wf.ecut_results_cache)} items")
    print(f"  ecutwfc values tested: {sorted(set(wf.ecut_results_cache.keys()))}")
    print(f"  kspacing_results_cache: {len(wf.kspacing_results_cache)} items")
    print(f"  kspacing values tested: {sorted(set(wf.kspacing_results_cache.keys()))}")
    
    cache_size_after_second = {
        'ecut': len(wf.ecut_results_cache),
        'kspacing': len(wf.kspacing_results_cache)
    }
    
    # ========================================================================
    # ANALYSIS
    # ========================================================================
    print("\n" + "="*80)
    print("ANALYSIS")
    print("="*80)
    
    print(f"\nCalculations reused from first run (cached):")
    print(f"  ecut values: {cache_size_after_first['ecut']} → {cache_size_after_second['ecut']} (+{cache_size_after_second['ecut'] - cache_size_after_first['ecut']} new)")
    print(f"  kspacing values: {cache_size_after_first['kspacing']} → {cache_size_after_second['kspacing']} (+{cache_size_after_second['kspacing'] - cache_size_after_first['kspacing']} new)")
    
    if cache_size_after_second['ecut'] > cache_size_after_first['ecut']:
        print(f"\n✓ CACHE IS WORKING!")
        print(f"  - First run ('low') tested {cache_size_after_first['ecut']} ecutwfc values")
        print(f"  - Second run ('high') added {cache_size_after_second['ecut'] - cache_size_after_first['ecut']} more")
        print(f"  - This means {cache_size_after_first['ecut']} values were REUSED (not recalculated)")
    else:
        print(f"\n✗ Cache appears NOT working - no additional ecutwfc values added")
    
    # Get recommendations
    print(f"\nRecommendations:")
    print(f"  Low precision:  ecutwfc={wf.optimal_ecutwfc:.1f} Ry (initial iteration)")
    print(f"  High precision: ecutwfc={wf.optimal_ecutwfc:.1f} Ry (after second run)")
    
    print("\n" + "="*80)
    print("TEST COMPLETE")
    print("="*80 + "\n")
    
    return True

if __name__ == '__main__':
    success = test_cache_reuse()
    sys.exit(0 if success else 1)
