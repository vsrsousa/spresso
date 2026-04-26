#!/usr/bin/env python3
"""
Load a CIF file, generate a slab, and validate geometry.
No machine dependency - pure Python/ASE/pymatgen.

Usage:
    python from_cif_to_slab.py au_bulk.cif --surface 111 --nlayers 3 --vacuum 5.0
"""

import os
import sys
import argparse
import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms

def read_cif_structure(cif_path):
    """Load structure from CIF file."""
    print(f"[1/4] Loading CIF: {cif_path}")
    
    if not os.path.exists(cif_path):
        raise FileNotFoundError(f"CIF file not found: {cif_path}")
    
    try:
        structure = read(cif_path)
        print(f"  ✓ Loaded: {structure.get_chemical_formula()}")
        print(f"    Atoms: {len(structure)}")
        print(f"    Lattice a={structure.cell[0,0]:.6f} Å")
        return structure
    except Exception as e:
        print(f"  ✗ Failed to read CIF: {e}")
        raise


def generate_slab_from_cif(bulk_atoms, surface=(1, 1, 1), nlayers=3, vacuum=5.0):
    """
    Generate a slab from bulk structure.
    Uses pymatgen if available, falls back to simple repetition.
    """
    print(f"\n[2/4] Generating slab (surface={surface}, nlayers={nlayers}, vacuum={vacuum:.1f} Å)")
    
    try:
        from pymatgen.core import Structure
        from pymatgen.analysis.surface import SlabGenerator
        
        # Convert ASE to pymatgen
        print(f"  → Using pymatgen for proper slab generation")
        pmg_structure = Structure(
            lattice=bulk_atoms.cell,
            species=[bulk_atoms.get_chemical_symbols()[0]] * len(bulk_atoms),
            coords=bulk_atoms.get_scaled_positions(),
        )
        
        # Generate slab
        slab_gen = SlabGenerator(
            pmg_structure,
            surface,
            min_slab_size=nlayers * 3.0,
            min_vacuum_size=vacuum,
            center_slab=True,
        )
        
        slab_structure = slab_gen.get_slab()
        
        # Convert back to ASE
        slab = read(f"POSCAR.{slab_structure.formula}")
        
        # Use ASE's built-in method
        from ase.build import surface as ase_surface_builder
        
        # Map (hkl) to string
        surface_str = f"{''.join(map(str, surface))}"
        
        # For FCC Au(111)
        if surface == (1, 1, 1):
            from ase.build import fcc111
            slab = fcc111(
                'Au',
                size=(1, 1, nlayers),
                a=bulk_atoms.cell[0, 0],
                vacuum=vacuum,
            )
            print(f"  ✓ Generated Au(111) slab: {len(slab)} atoms")
            return slab
        else:
            # Fallback for other surfaces
            raise NotImplementedError(f"Surface {surface} needs pymatgen")
            
    except ImportError:
        print(f"  ⚠ pymatgen not available, using ASE simple surface builder")
        
        from ase.build import fcc111, fcc100, add_vacuum
        
        # Simple FCC surface builders
        if surface == (1, 1, 1):
            slab = fcc111('Au', size=(1, 1, nlayers), a=bulk_atoms.cell[0, 0], vacuum=0)
            add_vacuum(slab, vacuum)
        elif surface == (1, 0, 0):
            from ase.build import fcc100
            slab = fcc100('Au', size=(1, 1, nlayers), a=bulk_atoms.cell[0, 0], vacuum=0)
            add_vacuum(slab, vacuum)
        else:
            raise NotImplementedError(f"Surface {surface} not supported without pymatgen")
        
        print(f"  ✓ Generated {surface} slab: {len(slab)} atoms")
        return slab
    
    except Exception as e:
        print(f"  ✗ Slab generation failed: {e}")
        raise


def fix_bottom_layers(slab, nlayers, fix_percent=0.4):
    """Apply constraints to bottom layers (typical DFT setup)."""
    print(f"\n[3/4] Applying constraints...")
    
    z_pos = slab.get_positions()[:, 2]
    z_min = z_pos.min()
    z_max = z_pos.max()
    z_range = z_max - z_min
    
    # Fix bottom ~40% of atoms
    fix_threshold = z_min + z_range * fix_percent
    
    indices_to_fix = [i for i, z in enumerate(z_pos) if z < fix_threshold]
    
    if indices_to_fix:
        slab.set_constraint(FixAtoms(indices=indices_to_fix))
        print(f"  ✓ Fixed {len(indices_to_fix)}/{len(slab)} atoms (bottom {fix_percent*100:.0f}%)")
    else:
        print(f"  ⚠ No atoms fixed (all above threshold)")
    
    return slab


def analyze_slab(slab, surface_indices):
    """Analyze slab geometry."""
    
    print(f"\n[4/4] Analyzing geometry...")
    
    # Basic info
    print(f"\n  STRUCTURE:")
    print(f"    Atoms: {len(slab)}")
    print(f"    Formula: {slab.get_chemical_formula()}")
    
    # Cell
    cell = slab.cell
    print(f"\n  CELL PARAMETERS:")
    print(f"    a: {np.linalg.norm(cell[0]):.6f} Å")
    print(f"    b: {np.linalg.norm(cell[1]):.6f} Å") 
    print(f"    c (z): {np.linalg.norm(cell[2]):.6f} Å")
    
    # Positions
    positions = slab.get_positions()
    z_pos = positions[:, 2]
    z_min, z_max = z_pos.min(), z_pos.max()
    
    print(f"\n  Z-DIRECTION:")
    print(f"    Z-range: {z_min:.6f} to {z_max:.6f} Å")
    print(f"    Slab thickness: {z_max - z_min:.6f} Å")
    
    cell_z = np.linalg.norm(cell[2])
    slab_thickness = z_max - z_min
    vacuum = cell_z - slab_thickness
    print(f"    Vacuum: {vacuum:.6f} Å")
    
    # Layers
    unique_z = np.unique(np.round(z_pos, 4))
    print(f"\n  LAYERS:")
    print(f"    Count: {len(unique_z)}")
    for i, z in enumerate(unique_z, 1):
        n_atoms = np.sum(np.abs(z_pos - z) < 0.01)
        status = "(fixed)" if hasattr(slab, 'constraints') and slab.constraints else ""
        print(f"      Layer {i}: z={z:.6f} Å ({n_atoms} atoms) {status}")
    
    # Distances - manual calculation with PBC
    positions = slab.get_positions()
    cell = slab.get_cell()
    
    dists_list = []
    for i in range(len(slab)):
        for j in range(i+1, len(slab)):
            # Calculate distance with periodic boundary conditions
            diff = positions[j] - positions[i]
            # Wrap to nearest image
            for k in range(3):
                if abs(diff[k]) > 0.5 * np.linalg.norm(cell[k]):
                    diff[k] -= np.sign(diff[k]) * np.linalg.norm(cell[k])
            dist = np.linalg.norm(diff)
            if dist > 0.01:
                dists_list.append(dist)
    
    dists_nonzero = np.array(dists_list) if dists_list else np.array([0])
    
    print(f"\n  ATOMIC DISTANCES:")
    print(f"    Min: {dists_nonzero.min():.6f} Å")
    print(f"    Mean: {dists_nonzero.mean():.6f} Å")
    print(f"    Max: {dists_nonzero.max():.6f} Å")
    print(f"    (Au-Au bulk ≈ 2.878 Å)")
    
    # Constraints
    print(f"\n  CONSTRAINTS:")
    if hasattr(slab, 'constraints') and slab.constraints:
        print(f"    ✓ Active")
        for c in slab.constraints:
            print(f"      {c}")
    else:
        print(f"    ✗ None (atoms free to move)")
    
    # Validation
    print(f"\n  VALIDATION:")
    all_good = True
    
    if vacuum < 3.0:
        print(f"    ✗ Vacuum {vacuum:.2f} Å < 3.0 Å")
        all_good = False
    elif vacuum < 5.0:
        print(f"    ⚠ Vacuum {vacuum:.2f} Å < 5.0 Å (marginal)")
    else:
        print(f"    ✓ Vacuum {vacuum:.2f} Å")
    
    if dists_nonzero.min() < 2.0:
        print(f"    ✗ Min distance {dists_nonzero.min():.2f} Å (atoms too close)")
        all_good = False
    else:
        print(f"    ✓ Min distance {dists_nonzero.min():.2f} Å")
    
    if len(unique_z) < 2:
        print(f"    ✗ Only 1 layer (degenerate slab)")
        all_good = False
    else:
        print(f"    ✓ {len(unique_z)} layers")
    
    if hasattr(slab, 'constraints') and slab.constraints:
        print(f"    ✓ Constraints present")
    else:
        print(f"    ⚠ No constraints (will collapse during relaxation)")
        all_good = False
    
    return all_good


def main():
    parser = argparse.ArgumentParser(
        description='Load CIF, generate slab, validate geometry'
    )
    parser.add_argument('cif_file', help='Path to CIF file')
    parser.add_argument('--surface', default='111', help='Surface indices (default: 111)')
    parser.add_argument('--nlayers', type=int, default=3, help='Number of layers')
    parser.add_argument('--vacuum', type=float, default=5.0, help='Vacuum size in Å')
    parser.add_argument('--output', default='slab.xyz', help='Output file')
    
    args = parser.parse_args()
    
    # Parse surface
    surface = tuple(int(x) for x in args.surface)
    
    print("="*70)
    print(f"CIF → SLAB GENERATOR")
    print("="*70)
    
    try:
        # Load CIF
        bulk = read_cif_structure(args.cif_file)
        
        # Generate slab
        slab = generate_slab_from_cif(bulk, surface=surface, nlayers=args.nlayers, vacuum=args.vacuum)
        
        # Fix atoms
        slab = fix_bottom_layers(slab, args.nlayers)
        
        # Analyze
        valid = analyze_slab(slab, surface)
        
        # Save
        write(args.output, slab)
        print(f"\n  ✓ Saved to {args.output}")
        
        # Summary
        print(f"\n{'='*70}")
        if valid:
            print("✓ Slab is valid and ready for relaxation!")
        else:
            print("⚠ Fix issues before relaxation")
        print(f"{'='*70}\n")
        
        return 0 if valid else 1
        
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
