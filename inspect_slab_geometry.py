#!/usr/bin/env python3
"""
Inspect and validate slab geometry before relaxation.

Checks:
- Layer thickness and spacing
- Vacuum size
- Atomic positions
- Cell parameters
- Constraints
"""

import numpy as np
from ase import Atoms
from ase.io import read, write
import sys

def analyze_slab_geometry(slab: Atoms, surface_index: tuple = (1, 1, 1)) -> dict:
    """
    Comprehensive analysis of slab geometry.
    
    Returns dict with validation checks.
    """
    results = {
        'valid': True,
        'issues': [],
        'warnings': [],
        'info': [],
    }
    
    # Cell parameters
    cell = slab.cell
    a1, a2, a3 = cell[0], cell[1], cell[2]
    
    results['info'].append(f"Cell parameters:")
    results['info'].append(f"  a1 = {np.linalg.norm(a1):.4f} Å, angle = {np.degrees(np.arccos(np.dot(a1, a2)/(np.linalg.norm(a1)*np.linalg.norm(a2)))):.2f}°")
    results['info'].append(f"  a2 = {np.linalg.norm(a2):.4f} Å, angle = {np.degrees(np.arccos(np.dot(a1, a3)/(np.linalg.norm(a1)*np.linalg.norm(a3)))):.2f}°")
    results['info'].append(f"  a3 (perpendicular/vacuum) = {np.linalg.norm(a3):.4f} Å")
    
    # Atomic positions
    positions = slab.get_positions()
    z_positions = positions[:, 2]
    z_min, z_max = z_positions.min(), z_positions.max()
    slab_thickness = z_max - z_min
    vacuum = np.linalg.norm(a3) - slab_thickness
    
    results['info'].append(f"\nSlab structure:")
    results['info'].append(f"  Number of atoms: {len(slab)}")
    results['info'].append(f"  Z position range: {z_min:.4f} to {z_max:.4f} Å")
    results['info'].append(f"  Slab thickness: {slab_thickness:.4f} Å")
    results['info'].append(f"  Vacuum size: {vacuum:.4f} Å")
    
    # Check vacuum
    if vacuum < 3.0:
        results['issues'].append(f"⚠ VACUUM TOO SMALL: {vacuum:.2f} Å < 3.0 Å")
        results['valid'] = False
    elif vacuum < 5.0:
        results['warnings'].append(f"Vacuum is small ({vacuum:.2f} Å), consider > 5.0 Å")
    
    # Analyze layer structure
    # Group atoms by z-coordinate (within tolerance)
    z_tol = 0.5  # Å
    unique_z = []
    for z in z_positions:
        if not unique_z or min(abs(z - uz) for uz in unique_z) > z_tol:
            unique_z.append(z)
    
    unique_z.sort()
    nlayers = len(unique_z)
    
    if nlayers > 0:
        layer_spacing = []
        for i in range(1, len(unique_z)):
            spacing = unique_z[i] - unique_z[i-1]
            layer_spacing.append(spacing)
        
        avg_spacing = np.mean(layer_spacing) if layer_spacing else 0
        
        results['info'].append(f"\nLayer analysis:")
        results['info'].append(f"  Number of layers: {nlayers}")
        results['info'].append(f"  Layer positions (z): {[f'{z:.4f}' for z in unique_z]}")
        if layer_spacing:
            results['info'].append(f"  Layer spacing: {[f'{s:.4f}' for s in layer_spacing]} Å")
            results['info'].append(f"  Average spacing: {avg_spacing:.4f} Å")
            
            # Check for irregular spacing
            spacing_variation = np.std(layer_spacing)
            if spacing_variation > 0.1:
                results['warnings'].append(f"Irregular layer spacing (σ = {spacing_variation:.4f} Å)")
    
    # Check for atoms too close
    results['info'].append(f"\nAtomic distances:")
    min_dist = float('inf')
    max_dist = 0
    
    for i in range(len(slab)):
        for j in range(i+1, len(slab)):
            dist = slab.get_distance(i, j, mic=True)
            min_dist = min(min_dist, dist)
            max_dist = max(max_dist, dist)
    
    results['info'].append(f"  Minimum distance: {min_dist:.4f} Å")
    results['info'].append(f"  Maximum distance: {max_dist:.4f} Å")
    
    if min_dist < 2.0:
        results['issues'].append(f"⚠ ATOMS TOO CLOSE: min_dist = {min_dist:.2f} Å < 2.0 Å")
        results['valid'] = False
    elif min_dist < 2.5:
        results['warnings'].append(f"Some atoms quite close (min = {min_dist:.2f} Å)")
    
    # Check for structural defects
    # For Au, typical distance is ~2.88 Å
    if nlayers > 1 and layer_spacing:
        expected_spacing = np.linalg.norm(a3) / (nlayers + 1)  # Approximate
        if abs(avg_spacing - expected_spacing) > 0.5:
            results['warnings'].append(f"Layer spacing differs from expected geometry")
    
    # Check constraints
    if hasattr(slab, 'constraints') and slab.constraints:
        results['info'].append(f"\nConstraints:")
        for constraint in slab.constraints:
            results['info'].append(f"  {constraint}")
    else:
        results['warnings'].append("No constraints found (atoms may move unexpectedly)")
    
    # Summary
    results['info'].append(f"\n{'='*70}")
    if results['valid']:
        results['info'].append("✓ GEOMETRY VALID - OK to proceed with relaxation")
    else:
        results['info'].append("✗ GEOMETRY ISSUES FOUND - FIX BEFORE RELAXATION")
    
    results['info'].append(f"{'='*70}")
    
    return results


def main():
    if len(sys.argv) < 2:
        print("Usage: python inspect_slab_geometry.py <structure.in> [surface_index]")
        print("Examples:")
        print("  python inspect_slab_geometry.py relax/111/nlayers_3/nlayers_3.in")
        print("  python inspect_slab_geometry.py au111_convergence/vacuum_5/vacuum_5.in 1 1 1")
        sys.exit(1)
    
    structure_file = sys.argv[1]
    surface_index = tuple(int(x) for x in sys.argv[2:5]) if len(sys.argv) > 2 else (1, 1, 1)
    
    try:
        # Try to read ASE structure file
        if structure_file.endswith('.in'):
            # QE input - try to read as generic
            print(f"Reading QE input: {structure_file}")
            # This is tricky - we'd need a parser
            print("✗ QE input format not directly supported")
            print("  Convert to structure.in or use ASE-compatible format (e.g., .xyz)")
            sys.exit(1)
        else:
            slab = read(structure_file)
    except Exception as e:
        print(f"✗ Error reading {structure_file}: {e}")
        sys.exit(1)
    
    print(f"\n{'='*70}")
    print(f"SLAB GEOMETRY INSPECTION")
    print(f"{'='*70}\n")
    
    results = analyze_slab_geometry(slab, surface_index)
    
    # Print info
    for line in results['info']:
        print(line)
    
    # Print warnings
    if results['warnings']:
        print(f"\n{'⚠ WARNINGS':")
        for warning in results['warnings']:
            print(f"  {warning}")
    
    # Print issues
    if results['issues']:
        print(f"\n{'❌ ISSUES':")
        for issue in results['issues']:
            print(f"  {issue}")
    
    # Exit code
    if not results['valid']:
        sys.exit(1)


if __name__ == '__main__':
    main()
