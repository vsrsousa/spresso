#!/usr/bin/env python3
"""
Detailed analysis of vertical layer structure
"""

import numpy as np
from ase.build import bulk
from xespresso.workflow.slab_workflow import SlabWorkflow

au_bulk = bulk('Au', 'fcc', a=4.08)

slab_wf = SlabWorkflow(
    bulk_atoms=au_bulk,
    surface_indices=[(1, 1, 1)],
    min_vacuum_size=15.0,
)

print("="*80)
print("DETAILED LAYER ANALYSIS")
print("="*80)

for nlayers in [3, 4]:
    slab = slab_wf._regenerate_slab_with_nlayers((1, 1, 1), nlayers)
    
    positions = slab.get_positions()
    z_positions = positions[:, 2]
    
    print(f"\nnlayers={nlayers} ({len(slab)} total atoms)")
    print("-" * 80)
    
    # Get unique z-positions
    z_unique_tolerance = 0.1  # Angstrom
    z_sorted = np.sort(z_positions)
    
    layers = []
    current_layer = [z_sorted[0]]
    
    for z in z_sorted[1:]:
        if z - current_layer[-1] < z_unique_tolerance:
            current_layer.append(z)
        else:
            layers.append(np.mean(current_layer))
            current_layer = [z]
    if current_layer:
        layers.append(np.mean(current_layer))
    
    print(f"Unique z-layers found: {len(layers)}")
    print(f"Expected z-layers: {nlayers}")
    
    for i, z_mean in enumerate(layers):
        atoms_in_layer = np.sum(np.abs(z_positions - z_mean) < z_unique_tolerance)
        print(f"  Layer {i}: z={z_mean:.4f} Å, atoms={atoms_in_layer}")
    
    # Check individual positions
    print(f"\nAll z-positions (sorted):")
    for i, z in enumerate(z_sorted):
        print(f"  Atom {i}: z={z:.4f} Å")
