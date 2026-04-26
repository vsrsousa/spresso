#!/usr/bin/env python3
"""
Quick script to generate and inspect a single slab geometry before relaxation.

No external machine dependency - uses pure ASE + pymatgen.

Usage:
    python validate_slab_before_relax.py
"""

import os
import sys
import numpy as np
from ase.build import bulk
from ase.io import write
from ase.constraints import FixAtoms

# Add xespresso to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from pymatgen.core import Structure
    from pymatgen.analysis.surface import SlabGenerator as PMGSlabGenerator
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    PYMATGEN_AVAILABLE = True
except ImportError:
    PYMATGEN_AVAILABLE = False


def generate_slab_simple(bulk_atoms, surface=(1, 1, 1), nlayers=3, vacuum=5.0):
    """
    Simple slab generator using ASE + pymatgen.
    
    Does NOT require machine configuration.
    """
    from ase.build import fcc111, surface
    
    if not PYMATGEN_AVAILABLE:
        # Fall back to simple manual construction
        print("  ⚠ pymatgen not available, using simple construction")
        # Just duplicate along z with vacuum
        slab = bulk_atoms.repeat((1, 1, nlayers))
        # Add vacuum
        cell = slab.cell.copy()
        cell[2, 2] += vacuum
        slab.set_cell(cell, scale_atoms=False)
        slab.center(vacuum=vacuum, axis=2)
        return slab
    
    # Use pymatgen for proper slab generation
    try:
        # Convert ASE to pymatgen Structure
        from pymatgen.core import Structure
        structure = Structure(
            lattice=bulk_atoms.cell,
            species=[bulk_atoms.get_chemical_symbols()[0]] * len(bulk_atoms),
            coords=bulk_atoms.get_scaled_positions(),
        )
        
        # Generate slab
        slab_gen = PMGSlabGenerator(
            structure,
            surface,
            min_slab_size=nlayers * 2.4,  # Approximate
            min_vacuum_size=vacuum,
            center_slab=True,
        )
        
        slab_structure = slab_gen.get_slab()
        
        # Convert back to ASE
        from ase import Atoms
        slab = Atoms(
            symbols=[site.species_string for site in slab_structure],
            positions=slab_structure.cart_coords,
            cell=slab_structure.lattice.matrix,
        )
        slab.center(vacuum=vacuum, axis=2)
        
        return slab
    except Exception as e:
        print(f"  ⚠ pymatgen slab generation failed: {e}")
        print("    Using simple construction")
        slab = bulk_atoms.repeat((1, 1, nlayers))
        cell = slab.cell.copy()
        cell[2, 2] += vacuum
        slab.set_cell(cell, scale_atoms=False)
        slab.center(vacuum=vacuum, axis=2)
        return slab


def main():
    print("="*70)
    print("SLAB GEOMETRY VALIDATION BEFORE RELAXATION")
    print("="*70)
    print("(No machine dependency - pure Python/ASE)\n")
    
    # Step 1: Generate bulk (quick)
    print("[1/3] Generating bulk Au...")
    bulk_au = bulk('Au', 'fcc', a=4.0782)
    print(f"  ✓ Bulk: {len(bulk_au)} atoms, a = 4.0782 Å")
    
    # Step 2: Mock bulk convergence results
    print("\n[2/3] Bulk convergence (mocked)...")
    conv_results = {
        'optimal_ecutwfc': 60.0,
        'optimal_kspacing': 0.24,
        'energy_per_atom': -3.73657,  # Au bulk energy
    }
    print(f"  ✓ Optimal ecutwfc: {conv_results['optimal_ecutwfc']} Ry")
    print(f"  ✓ Optimal k-spacing: {conv_results['optimal_kspacing']} Å⁻¹")
    
    # Step 3: Generate slab (Phase 2)
    print("\n[3/3] Generating Au(111) slab...")
    try:
        slab = generate_slab_simple(bulk_au, surface=(1, 1, 1), nlayers=3, vacuum=5.0)
        print(f"  ✓ Slab generated: {len(slab)} atoms")
    except Exception as e:
        print(f"  ✗ Slab generation failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # ========================================================================
    # GEOMETRY INSPECTION
    # ========================================================================
    print("\n" + "="*70)
    print("GEOMETRY ANALYSIS")
    print("="*70)
    
    cell = slab.cell
    positions = slab.get_positions()
    
    # Cell info
    print(f"\nCELL PARAMETERS:")
    print(f"  a1: {np.linalg.norm(cell[0]):.6f} Å")
    print(f"  a2: {np.linalg.norm(cell[1]):.6f} Å")
    print(f"  a3 (z-direction): {np.linalg.norm(cell[2]):.6f} Å")
    
    # Atomic structure
    print(f"\nATOMIC STRUCTURE:")
    print(f"  Total atoms: {len(slab)}")
    
    z_positions = positions[:, 2]
    z_min, z_max = z_positions.min(), z_positions.max()
    slab_thickness = z_max - z_min
    vacuum = np.linalg.norm(cell[2]) - slab_thickness
    
    print(f"  Z-coordinate range: {z_min:.6f} to {z_max:.6f} Å")
    print(f"  Slab thickness: {slab_thickness:.6f} Å")
    print(f"  Vacuum size: {vacuum:.6f} Å")
    
    # Layer structure
    z_tol = 0.5
    unique_z = []
    for z in z_positions:
        if not unique_z or min(abs(z - uz) for uz in unique_z) > z_tol:
            unique_z.append(z)
    
    unique_z.sort()
    print(f"\nLAYER STRUCTURE:")
    print(f"  Number of layers: {len(unique_z)}")
    print(f"  Layer positions:")
    for i, z in enumerate(unique_z, 1):
        n_atoms = sum(1 for pos_z in z_positions if abs(pos_z - z) < z_tol)
        print(f"    Layer {i}: z = {z:.6f} Å ({n_atoms} atoms)")
    
    if len(unique_z) > 1:
        spacings = [unique_z[i] - unique_z[i-1] for i in range(1, len(unique_z))]
        print(f"  Layer spacings: {[f'{s:.6f}' for s in spacings]} Å")
        print(f"  Average spacing: {np.mean(spacings):.6f} Å")
        print(f"  Spacing variation: {np.std(spacings):.6f} Å")
    
    # Distances
    print(f"\nATOMIC DISTANCES:")
    min_dist = float('inf')
    max_dist = 0
    
    for i in range(len(slab)):
        for j in range(i+1, len(slab)):
            dist = slab.get_distance(i, j, mic=True)
            min_dist = min(min_dist, dist)
            max_dist = max(max_dist, dist)
    
    print(f"  Minimum: {min_dist:.6f} Å")
    print(f"  Maximum: {max_dist:.6f} Å")
    print(f"  (Au-Au in bulk ≈ 2.878 Å)")
    
    # Constraints
    print(f"\nCONSTRAINTS:")
    if hasattr(slab, 'constraints') and slab.constraints:
        for c in slab.constraints:
            print(f"  {c}")
    else:
        print(f"  None (ISSUE: atoms will move unexpectedly)")
    
    # Validation
    print(f"\n" + "="*70)
    print(f"VALIDATION")
    print(f"="*70)
    
    issues = []
    warnings = []
    
    # Check vacuum
    if vacuum < 3.0:
        issues.append(f"Vacuum too small: {vacuum:.2f} Å < 3.0 Å")
    elif vacuum < 5.0:
        warnings.append(f"Vacuum is small: {vacuum:.2f} Å (ideal: > 5 Å)")
    
    # Check min distance
    if min_dist < 2.0:
        issues.append(f"Atoms too close: {min_dist:.2f} Å < 2.0 Å")
    elif min_dist < 2.5:
        warnings.append(f"Some atoms close: {min_dist:.2f} Å")
    
    # Check layer spacing
    if len(unique_z) > 1:
        avg_spacing = np.mean(spacings)
        if avg_spacing < 1.5:
            warnings.append(f"Layer spacing very small: {avg_spacing:.2f} Å")
        spacing_var = np.std(spacings)
        if spacing_var > 0.2:
            warnings.append(f"Irregular spacing: σ = {spacing_var:.3f} Å")
    
    # Check constraints
    if not (hasattr(slab, 'constraints') and slab.constraints):
        issues.append("No constraints! Slab will collapse during relaxation")
    
    # Print results
    if warnings:
        print("\n⚠  WARNINGS:")
        for w in warnings:
            print(f"   {w}")
    
    if issues:
        print("\n❌ ISSUES - FIX BEFORE RELAXATION:")
        for i in issues:
            print(f"   {i}")
    else:
        print("\n✓ GEOMETRY OK - Ready for relaxation")
    
    # Save for inspection
    print(f"\n" + "="*70)
    print("OUTPUT")
    print(f"="*70)
    
    output_file = 'slab_inspection.xyz'
    write(output_file, slab)
    print(f"✓ Saved slab to: {output_file}")
    print(f"  (Open with VESTA or similar to visualize)")
    
    if issues:
        print(f"\n⚠ Fix issues above before running relaxation!")
        sys.exit(1)
    else:
        print(f"\n✓ Geometry validated - safe to relax!")
        sys.exit(0)


if __name__ == '__main__':
    main()
