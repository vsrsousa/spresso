#!/usr/bin/env python3
"""
Generate slab and QE input using xespresso - with internal mocking.
No machine dependency - just generates the geometry and input file.

Usage:
    python slab_from_cif_xespresso.py au_bulk.cif --surface 111 --nlayers 3
"""

import os
import sys
import argparse
import tempfile
from ase.io import read
from ase.constraints import FixAtoms

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from xespresso.workflow.slab_workflow import SlabWorkflow


def normalize_cell(atoms, surface=(1, 1, 1)):
    """
    Normalize slab cell to standard orientation:
    - a, b parallel to surface (xy-plane)
    - c perpendicular to surface (z-direction)
    - All values positive
    - No negative z coordinates
    """
    print(f"  Normalizing cell orientation...")
    
    slab = atoms.copy()
    
    # Get positions
    positions = slab.get_positions()
    cell = slab.get_cell()
    
    # Find z-range
    z_coords = positions[:, 2]
    z_min = z_coords.min()
    z_max = z_coords.max()
    z_range = z_max - z_min
    
    # Get cell z-length
    cell_z = abs(cell[2, 2])  # Use absolute value
    
    # If cell_z is negative, flip it
    if cell[2, 2] < 0:
        cell[2, 2] = cell_z
        positions[:, 2] = cell_z - (positions[:, 2] - z_min)  # Flip z-coords
    
    # Shift z so atoms start from reasonable position
    # Put first layer around z = vacuum_size/2 for symmetry
    z_offset = cell[2, 2] / 2 - z_range / 2
    positions[:, 2] = positions[:, 2] - z_min + z_offset
    
    # Ensure a and b vectors have positive x-component
    if cell[0, 0] < 0:
        cell[0] = -cell[0]
    if cell[1, 0] < 0:
        cell[1] = -cell[1]
    
    slab.set_cell(cell, scale_atoms=False)
    slab.set_positions(positions)
    
    # Check if all positions are positive
    if slab.get_positions().min() < -0.1:
        slab.center()
    
    print(f"    Cell after normalization:")
    cell = slab.get_cell()
    print(f"      a = [{cell[0,0]:10.6f}, {cell[0,1]:10.6f}, {cell[0,2]:10.6f}]")
    print(f"      b = [{cell[1,0]:10.6f}, {cell[1,1]:10.6f}, {cell[1,2]:10.6f}]")
    print(f"      c = [{cell[2,0]:10.6f}, {cell[2,1]:10.6f}, {cell[2,2]:10.6f}]")
    
    return slab


def main():
    parser = argparse.ArgumentParser(description='CIF → Slab with xespresso')
    parser.add_argument('cif_file', help='CIF file with bulk structure')
    parser.add_argument('--surface', default='111', help='Surface (default: 111)')
    parser.add_argument('--nlayers', type=int, default=3, help='Layers')
    parser.add_argument('--vacuum', type=float, default=5.0, help='Vacuum (Å)')
    parser.add_argument('--output', default='slab_qe.in', help='QE input file')
    
    args = parser.parse_args()
    
    # Parse surface
    surface = tuple(int(x) for x in args.surface)
    
    print("="*70)
    print("XESPRESSO: CIF → SLAB → QE INPUT")
    print("="*70)
    
    # Step 1: Load bulk from CIF
    print(f"\n[1/3] Loading bulk from {args.cif_file}...")
    try:
        bulk = read(args.cif_file)
        print(f"  ✓ {bulk.get_chemical_formula()}: {len(bulk)} atom(s)")
    except Exception as e:
        print(f"  ✗ Error: {e}")
        return 1
    
    # Step 2: Create SlabWorkflow (no calculation, just geometry)
    print(f"\n[2/3] Creating xespresso SlabWorkflow...")
    try:
        wf = SlabWorkflow(
            bulk_atoms=bulk,
            surface_indices=[surface],
            nlayers=args.nlayers,
            min_vacuum_size=args.vacuum,
            pseudopotentials_config='default',
            machine='medusa',
        )
        print(f"  ✓ Workflow initialized")
    except Exception as e:
        print(f"  ✗ Error: {e}")
        return 1
    
    # Step 3: Generate slab (geometry only, no calculation)
    print(f"\n[3/3] Generating slab and QE input...")
    try:
        slabs = wf.generate_slabs()
        
        if surface not in slabs:
            print(f"  ✗ Slab {surface} not generated")
            return 1
        
        slab = slabs[surface]
        print(f"  ✓ Slab: {len(slab)} atoms, {slab.get_chemical_formula()}")
        
        # Normalize cell orientation
        slab = normalize_cell(slab, surface=surface)
        
        # Apply constraints (bottom 40% fixed)
        z_pos = slab.get_positions()[:, 2]
        z_min, z_max = z_pos.min(), z_pos.max()
        z_range = z_max - z_min
        fix_threshold = z_min + z_range * 0.4
        
        indices_to_fix = [i for i, z in enumerate(z_pos) if z < fix_threshold]
        if indices_to_fix:
            slab.set_constraint(FixAtoms(indices=indices_to_fix))
            print(f"  ✓ Fixed {len(indices_to_fix)} atoms")
        
        # Get QE input from CalculationWorkflow
        from xespresso.workflow.calculation_workflow import CalculationWorkflow
        
        calc_wf = CalculationWorkflow(
            slab,
            pseudopotentials_config='default',
            machine='medusa',
        )
        
        # Get the input data (this is what would be sent to QE)
        print(f"\n  Input data keys: {list(calc_wf.input_data.keys())}")
        
        # Save to file (use ASE + manual QE format)
        from ase.io import write
        import numpy as np
        
        # Create QE input manually
        cell = slab.get_cell()
        positions = slab.get_positions()
        symbols = slab.get_chemical_symbols()
        
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

ATOMIC_POSITIONS (angstrom)
"""
        
        # Add atoms with constraints
        for i, (symbol, pos) in enumerate(zip(symbols, positions)):
            is_fixed = False
            if hasattr(slab, 'constraints') and slab.constraints:
                for constraint in slab.constraints:
                    if hasattr(constraint, 'indices') and i in constraint.indices:
                        is_fixed = True
            
            fix_flag = "0 0 0" if is_fixed else "1 1 1"
            qe_input += f"{symbol:2s}  {pos[0]:12.6f}  {pos[1]:12.6f}  {pos[2]:12.6f}  {fix_flag}\n"
        
        qe_input += f"\nK_POINTS (automatic)\n 7 7 1  0 0 0\n"
        
        with open(args.output, 'w') as f:
            f.write(qe_input)
        
        print(f"\n  ✓ Saved to {args.output}")
        
        # Summary
        print(f"\n{'='*70}")
        print("OUTPUT GEOMETRY")
        print(f"{'='*70}")
        print(f"\nCell z-direction: {cell[2,2]:.6f} Å")
        print(f"Slab thickness: {z_max - z_min:.6f} Å")
        print(f"Vacuum: {cell[2,2] - (z_max - z_min):.6f} Å")
        print(f"Fixed atoms: {len(indices_to_fix) if indices_to_fix else 0}/{len(slab)}")
        print(f"\n✓ Ready to run: pw.x < {args.output}")
        print(f"{'='*70}\n")
        
        return 0
        
    except Exception as e:
        print(f"  ✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
