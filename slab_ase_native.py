#!/usr/bin/env python3
"""
Generate slabs usando APENAS ASE - sem pymatgen.
Com posições em CRYSTAL (fractional).

Usage:
    python slab_ase_native.py au_bulk.cif --surface 111 --nlayers 3 --vacuum 5.0
"""

import os
import sys
import argparse
import numpy as np
from ase.io import read
from ase.build import surface, bulk
from ase.constraints import FixAtoms

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def slab_from_cif_ase(cif_file, surface_indices, nlayers, vacuum):
    """
    Generate slab from CIF using ONLY ASE (no pymatgen).
    
    Surface indices map to:
    (1,0,0) → fcc100
    (1,1,0) → fcc110  
    (1,1,1) → fcc111
    """
    print(f"[1/3] Loading bulk from {cif_file}...")
    bulk_atoms = read(cif_file)
    print(f"  ✓ Loaded: {bulk_atoms.get_chemical_formula()}")
    
    # Map (hkl) to ASE surface function
    hkl = surface_indices[0]  # e.g., (1,1,1)
    element = bulk_atoms.get_chemical_symbols()[0]
    a = bulk_atoms.get_cell()[0, 0]  # Lattice parameter
    
    print(f"\n[2/3] Generating {element}{hkl} slab with ASE...")
    print(f"  Lattice parameter: {a:.6f} Å")
    print(f"  Layers: {nlayers}")
    print(f"  Vacuum: {vacuum} Å")
    
    # Use ASE native surface builder
    try:
        from ase.build import fcc111, fcc100, fcc110
        
        if hkl == (1, 1, 1):
            slab = fcc111(element, size=(1, 1, nlayers), a=a, vacuum=0)
        elif hkl == (1, 0, 0):
            slab = fcc100(element, size=(1, 1, nlayers), a=a, vacuum=0)
        elif hkl == (1, 1, 0):
            slab = fcc110(element, size=(1, 1, nlayers), a=a, vacuum=0)
        else:
            raise ValueError(f"Surface {hkl} not supported by ASE (use 100, 110, 111)")
        
        # Add vacuum
        from ase.build import add_vacuum
        add_vacuum(slab, vacuum)
        
        print(f"  ✓ Generated: {len(slab)} atoms")
        
    except ImportError:
        print(f"  ✗ ASE old version, trying generic surface()...")
        slab = surface(element, hkl, nlayers, vacuum=vacuum, a=a)
        print(f"  ✓ Generated: {len(slab)} atoms")
    
    return slab


def normalize_cell(atoms):
    """Normalize cell to standard orientation."""
    positions = atoms.get_positions()
    cell = atoms.get_cell()
    
    # Ensure c (z) is positive
    if cell[2, 2] < 0:
        cell[2, 2] = abs(cell[2, 2])
        positions[:, 2] = abs(positions[:, 2])
    
    # Ensure a and b have positive x
    if cell[0, 0] < 0:
        cell[0] = -cell[0]
    if cell[1, 0] < 0:
        cell[1] = -cell[1]
    
    atoms.set_cell(cell, scale_atoms=False)
    atoms.set_positions(positions)
    
    return atoms


def apply_constraints(slab, nlayers, fix_fraction=0.4):
    """Fix bottom layers."""
    z_pos = slab.get_positions()[:, 2]
    z_min, z_max = z_pos.min(), z_pos.max()
    z_range = z_max - z_min
    fix_threshold = z_min + z_range * fix_fraction
    
    indices_to_fix = [i for i, z in enumerate(z_pos) if z < fix_threshold]
    
    if indices_to_fix:
        slab.set_constraint(FixAtoms(indices=indices_to_fix))
    
    return slab, indices_to_fix


def main():
    parser = argparse.ArgumentParser(description='CIF → Slab (ASE native, no pymatgen)')
    parser.add_argument('cif_file', help='CIF file')
    parser.add_argument('--surface', default='111', help='Surface (100/110/111)')
    parser.add_argument('--nlayers', type=int, default=3, help='Layers')
    parser.add_argument('--vacuum', type=float, default=5.0, help='Vacuum (Å)')
    parser.add_argument('--output', default='slab_ase.in', help='QE input file')
    parser.add_argument('--crystal', action='store_true', default=True, help='Use CRYSTAL coordinates')
    
    args = parser.parse_args()
    
    surface = tuple(int(x) for x in args.surface)
    
    print("="*70)
    print("ASE NATIVE SLAB GENERATOR (no pymatgen)")
    print("="*70 + "\n")
    
    try:
        # Generate slab
        slab = slab_from_cif_ase(args.cif_file, [surface], args.nlayers, args.vacuum)
        slab = normalize_cell(slab)
        
        # Apply constraints
        slab, indices_to_fix = apply_constraints(slab, args.nlayers)
        print(f"  ✓ Fixed {len(indices_to_fix)} atoms\n")
        
        # Prepare positions and cell
        cell = slab.get_cell()
        positions_angstrom = slab.get_positions()
        symbols = slab.get_chemical_symbols()
        
        # Convert to fractional (CRYSTAL) coordinates
        inv_cell = np.linalg.inv(cell)
        positions_crystal = positions_angstrom @ inv_cell.T
        
        print(f"[3/3] Generating QE input...\n")
        
        # Show both formats
        print(f"  ANGSTROM coordinates:")
        for i, (sym, pos) in enumerate(zip(symbols, positions_angstrom)):
            print(f"    {sym:2s}  {pos[0]:12.6f}  {pos[1]:12.6f}  {pos[2]:12.6f}")
        
        print(f"\n  CRYSTAL (fractional) coordinates:")
        for i, (sym, pos) in enumerate(zip(symbols, positions_crystal)):
            print(f"    {sym:2s}  {pos[0]:12.6f}  {pos[1]:12.6f}  {pos[2]:12.6f}")
        
        # Build QE input with CRYSTAL coordinates
        qe_input = f"""&CONTROL
    calculation = 'relax'
    restart_mode = 'from_scratch'
    prefix = 'slab'
    outdir = './'
    wfcdir = './'
    pseudo_dir = './pseudos'
    nstep = 150
    upscale = 150.0
    etot_conv_thr = 1.0d-5
    forc_conv_thr = 1.0d-5
    verbosity = 'high'
/
&SYSTEM
    ibrav = 0
    nat = {len(slab)}
    ntyp = 1
    ecutwfc = 60.0
    ecutrho = 240.0
    occupations = 'smearing'
    smearing = 'mv'
    degauss = 0.01
    nspin = 1
/
&ELECTRONS
    diagonalization = 'david'
    mixing_mode = 'plain'
    mixing_beta = 0.7
    conv_thr = 1.0d-9
/
&IONS
    ion_dynamics = 'bfgs'
    upscale = 150.0
/
ATOMIC_SPECIES
Au  196.966  Au_ONCV_PBE-1.0.oncvpsp.upf

CELL_PARAMETERS (angstrom)
  {cell[0,0]:12.6f}  {cell[0,1]:12.6f}  {cell[0,2]:12.6f}
  {cell[1,0]:12.6f}  {cell[1,1]:12.6f}  {cell[1,2]:12.6f}
  {cell[2,0]:12.6f}  {cell[2,1]:12.6f}  {cell[2,2]:12.6f}

ATOMIC_POSITIONS (crystal)
"""
        
        # Add atoms with constraints in CRYSTAL coordinates
        for i, (symbol, pos_frac) in enumerate(zip(symbols, positions_crystal)):
            is_fixed = False
            if hasattr(slab, 'constraints') and slab.constraints:
                for constraint in slab.constraints:
                    if hasattr(constraint, 'indices') and i in constraint.indices:
                        is_fixed = True
            
            fix_flag = "0 0 0" if is_fixed else "1 1 1"
            qe_input += f"{symbol:2s}  {pos_frac[0]:12.6f}  {pos_frac[1]:12.6f}  {pos_frac[2]:12.6f}  {fix_flag}\n"
        
        qe_input += f"\nK_POINTS (automatic)\n 7 7 1  0 0 0\n"
        
        with open(args.output, 'w') as f:
            f.write(qe_input)
        
        print(f"\n  ✓ Saved to {args.output}\n")
        
        # Summary
        print(f"{'='*70}")
        print("✓ SLAB READY FOR QE")
        print(f"{'='*70}")
        print(f"Atoms: {len(slab)}")
        a_norm = np.linalg.norm(cell[0])
        b_norm = np.linalg.norm(cell[1])
        c_norm = np.linalg.norm(cell[2])
        print(f"Cell (Å): {a_norm:.2f} × {b_norm:.2f} × {c_norm:.2f}")
        print(f"Fixed: {len(indices_to_fix)}/{len(slab)} atoms")
        print(f"Coordinates: CRYSTAL (fractional) ✓")
        print(f"\nRun: pw.x < {args.output}")
        print(f"{'='*70}\n")
        
        return 0
        
    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
