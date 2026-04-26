#!/usr/bin/env python3
"""
Test: Verify that nlayers corresponds to actual geometric layers
"""

import numpy as np
from ase.build import bulk
from xespresso.workflow.slab_workflow import SlabWorkflow

# Generate Au bulk
au_bulk = bulk('Au', 'fcc', a=4.08)

# Create workflow
slab_wf = SlabWorkflow(
    bulk_atoms=au_bulk,
    surface_indices=[(1, 1, 1)],
    min_slab_size=6.0,
    min_vacuum_size=5.0,
    use_primitive_cell=False,  # Test with supercell first
)

print("="*80)
print("TEST: nlayers CORRESPONDENCE")
print("="*80)

# Expected d-spacing for Au(111): 4.08 / sqrt(3) = 2.3556 Å
# So for nlayers=N, we expect N distinct z-positions with spacing ~2.3556 Å

for nlayers in [3, 4, 5, 6, 7, 8]:
    slab = slab_wf._regenerate_slab_with_nlayers((1, 1, 1), nlayers)
    
    # Get z-positions and find unique layers
    z_positions = slab.get_positions()[:, 2]
    z_min = z_positions.min()
    z_max = z_positions.max()
    
    # Expected d-spacing
    d_hkl = 4.08 / np.sqrt(1**2 + 1**2 + 1**2)  # ~2.3556 Å
    
    # Count layers by grouping atoms with similar z (within 0.5 Å)
    z_sorted = np.sort(z_positions)
    layers = []
    current_layer = [z_sorted[0]]
    
    for z in z_sorted[1:]:
        if z - current_layer[-1] < 0.5:  # Same layer
            current_layer.append(z)
        else:  # New layer
            layers.append(np.mean(current_layer))
            current_layer = [z]
    if current_layer:
        layers.append(np.mean(current_layer))
    
    actual_layers = len(layers)
    
    # Calculate spacing between layers
    if len(layers) > 1:
        spacings = np.diff(layers)
        mean_spacing = np.mean(spacings)
    else:
        mean_spacing = 0
    
    match = "✓" if actual_layers == nlayers else "✗"
    
    print(f"\nnlayers={nlayers} {match}")
    print(f"  Total atoms: {len(slab)}")
    print(f"  Actual layers found: {actual_layers}")
    print(f"  Z range: {z_min:.4f} to {z_max:.4f} Å (height: {z_max-z_min:.4f} Å)")
    print(f"  Expected layer spacing: {d_hkl:.4f} Å")
    if len(layers) > 1:
        print(f"  Actual spacing: {mean_spacing:.4f} Å")
    print(f"  Layer z-positions: {[f'{z:.4f}' for z in layers]}")
    
print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print("\nIf all have ✓, nlayers correctly corresponds to geometric layers!")
