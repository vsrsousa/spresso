#!/usr/bin/env python3
"""
Test: Verify nlayers with both supercell and primitive options
"""

import numpy as np
from ase.build import bulk
from xespresso.workflow.slab_workflow import SlabWorkflow

au_bulk = bulk('Au', 'fcc', a=4.08)

print("="*80)
print("FULL TEST: nlayers CORRESPONDENCE (SUPERCELL vs PRIMITIVE)")
print("="*80)

for use_primitive in [False, True]:
    mode = "PRIMITIVE (1×1)" if use_primitive else "SUPERCELL (2×2)"
    
    slab_wf = SlabWorkflow(
        bulk_atoms=au_bulk,
        surface_indices=[(1, 1, 1)],
        min_slab_size=6.0,
        min_vacuum_size=5.0,
        use_primitive_cell=use_primitive,
    )
    
    print(f"\n{mode}")
    print("-" * 80)
    print(f"{'nlayers':<10} {'Total atoms':<15} {'Atoms/layer':<15}")
    print("-" * 80)
    
    for nlayers in [3, 4, 5, 6]:
        slab = slab_wf._regenerate_slab_with_nlayers((1, 1, 1), nlayers)
        natoms = len(slab)
        atoms_per_layer = natoms / nlayers
        
        print(f"{nlayers:<10} {natoms:<15} {atoms_per_layer:<15.1f}")

print("\n" + "="*80)
print("EXPECTED")
print("="*80)
print(f"\nSupercell (2×2): atoms_per_layer = 4.0 (constant)")
print(f"Primitive (1×1): atoms_per_layer = 1.0 (constant)")
print(f"\nReduction factor: 4.0×")
