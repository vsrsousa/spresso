#!/usr/bin/env python3
"""
Generate and validate a slab CORRECTLY before relaxation.
Uses SlabWorkflow but with mocked calculation backend.
"""

import os
import sys
import tempfile
import shutil
from ase.build import bulk
from ase.io import write
from ase.constraints import FixAtoms
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def analyze_slab(slab, nlayers=None):
    """Analyze slab geometry."""
    
    print(f"\n{'='*70}")
    print("SLAB GEOMETRY ANALYSIS")
    print(f"{'='*70}")
    
    # Basic info
    print(f"\nBASIC STRUCTURE:")
    print(f"  Total atoms: {len(slab)}")
    print(f"  Element: {slab.get_chemical_formula()}")
    
    # Cell
    cell = slab.cell
    print(f"\nCELL PARAMETERS:")
    print(f"  a1: {np.linalg.norm(cell[0]):.6f} Å")
    print(f"  a2: {np.linalg.norm(cell[1]):.6f} Å") 
    print(f"  a3 (z): {np.linalg.norm(cell[2]):.6f} Å")
    print(f"  α: {np.degrees(np.arccos(np.dot(cell[1], cell[2]) / (np.linalg.norm(cell[1]) * np.linalg.norm(cell[2])))):.2f}°")
    print(f"  β: {np.degrees(np.arccos(np.dot(cell[0], cell[2]) / (np.linalg.norm(cell[0]) * np.linalg.norm(cell[2])))):.2f}°")
    print(f"  γ: {np.degrees(np.arccos(np.dot(cell[0], cell[1]) / (np.linalg.norm(cell[0]) * np.linalg.norm(cell[1])))):.2f}°")
    
    # Atomic positions
    positions = slab.get_positions()
    z_pos = positions[:, 2]
    z_min, z_max = z_pos.min(), z_pos.max()
    
    print(f"\nATOMIC Z-POSITIONS:")
    print(f"  Min z: {z_min:.6f} Å")
    print(f"  Max z: {z_max:.6f} Å")
    print(f"  Slab thickness: {z_max - z_min:.6f} Å")
    
    # Vacuum
    cell_z = np.linalg.norm(cell[2])
    slab_thickness = z_max - z_min
    vacuum = cell_z - slab_thickness
    print(f"  Cell z-length: {cell_z:.6f} Å")
    print(f"  Vacuum: {vacuum:.6f} Å")
    
    # Layers
    unique_z = np.unique(np.round(z_pos, 4))
    print(f"\nLAYERS:")
    print(f"  Number of layers: {len(unique_z)}")
    
    if nlayers and len(unique_z) != nlayers:
        print(f"  ⚠ Expected {nlayers} layers, got {len(unique_z)}")
    
    for i, z in enumerate(unique_z, 1):
        n_atoms = np.sum(np.abs(z_pos - z) < 0.01)
        print(f"    Layer {i}: z={z:.6f} Å ({n_atoms} atoms)")
    
    # Inter-atomic distances
    from ase.geometry import get_distances
    dists, _ = get_distances(slab, mic=True)
    dists_flat = dists[np.triu_indices_from(dists, k=1)]
    dists_nonzero = dists_flat[dists_flat > 0.01]
    
    print(f"\nATOMIC DISTANCES:")
    print(f"  Min: {dists_nonzero.min():.6f} Å (Au-Au in bulk ≈ 2.878 Å)")
    print(f"  Max: {dists_nonzero.max():.6f} Å")
    print(f"  Mean: {dists_nonzero.mean():.6f} Å")
    
    # Constraints
    print(f"\nCONSTRAINTS:")
    if hasattr(slab, 'constraints') and slab.constraints:
        print(f"  ✓ Active constraints:")
        for constraint in slab.constraints:
            print(f"    - {constraint}")
    else:
        print(f"  ✗ NO CONSTRAINTS - Slab will collapse!")
        return False
    
    # Validation checks
    print(f"\nVALIDATION:")
    issues = []
    
    if vacuum < 3.0:
        issues.append(f"  ✗ Vacuum {vacuum:.2f} Å < 3.0 Å (too small)")
    elif vacuum < 5.0:
        issues.append(f"  ⚠ Vacuum {vacuum:.2f} Å < 5.0 Å (marginal)")
    else:
        print(f"  ✓ Vacuum {vacuum:.2f} Å OK")
    
    if dists_nonzero.min() < 2.0:
        issues.append(f"  ✗ Min distance {dists_nonzero.min():.2f} Å < 2.0 Å (atoms too close)")
    else:
        print(f"  ✓ Min distance {dists_nonzero.min():.2f} Å OK")
    
    if not (hasattr(slab, 'constraints') and slab.constraints):
        issues.append(f"  ✗ No constraints - atoms free to move unrealistically")
    else:
        print(f"  ✓ Constraints present")
    
    # Layer spacing
    if len(unique_z) > 1:
        spacings = np.diff(unique_z)
        spacing_std = np.std(spacings)
        if spacing_std > 0.1:
            issues.append(f"  ⚠ Layer spacing irregular (σ={spacing_std:.4f})")
        else:
            print(f"  ✓ Layer spacing regular (σ={spacing_std:.4f})")
    
    if issues:
        print("\nISSUES:")
        for issue in issues:
            print(issue)
        return False
    else:
        print("\n✓ ALL CHECKS PASSED")
        return True


def main():
    print("="*70)
    print("SLAB VALIDATION - Without machine dependency")
    print("="*70)
    
    # Step 1: Generate bulk
    print("\n[1/3] Generating bulk Au...")
    bulk_au = bulk('Au', 'fcc', a=4.0782)
    print(f"  ✓ {len(bulk_au)} atom(s)")
    
    # Step 2: Import and use SlabWorkflow
    print("\n[2/3] Creating SlabWorkflow...")
    try:
        from xespresso.workflow.slab_workflow import SlabWorkflow
        
        # Create workflow with mocked machine
        wf = SlabWorkflow(
            bulk_atoms=bulk_au,
            surface_indices=[(1, 1, 1)],
            nlayers=3,
            vacuum=5.0,
            pseudopotentials_config='default',
            protocol='standard',
            machine='medusa',  # Real machine from config
        )
        print(f"  ✓ SlabWorkflow created")
        
    except Exception as e:
        print(f"  ✗ Failed to create SlabWorkflow: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # Step 3: Generate slab (no calculation, just geometry)
    print("\n[3/3] Generating Au(111) slab...")
    try:
        # This doesn't run calculations, just creates geometry
        slabs = wf.generate_slabs()
        
        if (1, 1, 1) not in slabs:
            print(f"  ✗ Slab (1,1,1) not created")
            return False
        
        slab = slabs[(1, 1, 1)]
        print(f"  ✓ Slab: {len(slab)} atoms")
        
    except Exception as e:
        print(f"  ✗ Failed to generate slab: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # Step 4: Analyze
    valid = analyze_slab(slab, nlayers=3)
    
    # Step 5: Save for visualization
    print(f"\n[4/4] Saving structure...")
    try:
        write('slab_validation.xyz', slab)
        print(f"  ✓ Saved to slab_validation.xyz")
        print(f"    Open with: VESTA, Ovito, or jmol")
    except Exception as e:
        print(f"  ✗ Failed to save: {e}")
    
    print(f"\n{'='*70}")
    if valid:
        print("✓ Slab is ready for relaxation!")
        print("="*70)
        return True
    else:
        print("✗ Fix geometry issues before relaxation")
        print("="*70)
        return False


if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
