#!/usr/bin/env python3
"""
Convert slab geometry (XYZ) to Quantum ESPRESSO input format.
Shows exactly how the mocking becomes a real QE calculation.
"""

import sys
import os
from ase.io import read, write

def xyz_to_qe_input(xyz_file, output_file, ecutwfc=60, ecutrho=240, kspacing=0.30):
    """Convert XYZ slab to QE input file."""
    
    print(f"{'='*70}")
    print("XYZ → QUANTUM ESPRESSO INPUT CONVERSION")
    print(f"{'='*70}\n")
    
    # Read structure
    print(f"[1/4] Reading slab from {xyz_file}...")
    slab = read(xyz_file)
    print(f"  ✓ {len(slab)} atoms, formula: {slab.get_chemical_formula()}")
    
    # Get cell
    cell = slab.get_cell()
    positions = slab.get_positions()
    symbols = slab.get_chemical_symbols()
    
    print(f"\n[2/4] Cell parameters:")
    print(f"  a = {cell[0]}")
    print(f"  b = {cell[1]}")
    print(f"  c = {cell[2]}")
    
    # Get constraints info
    print(f"\n[3/4] Atomic positions & constraints:")
    for i, (symbol, pos) in enumerate(zip(symbols, positions)):
        # Check if atom is fixed
        is_fixed = False
        if hasattr(slab, 'constraints') and slab.constraints:
            for constraint in slab.constraints:
                if hasattr(constraint, 'index') and i in constraint.index:
                    is_fixed = True
                    break
                elif hasattr(constraint, 'indices') and i in constraint.indices:
                    is_fixed = True
                    break
        
        fix_str = "0 0 0 (FIXED)" if is_fixed else "1 1 1 (FREE)"
        print(f"  {i+1}. {symbol:2s}  {pos[0]:10.6f}  {pos[1]:10.6f}  {pos[2]:10.6f}  {fix_str}")
    
    # Generate QE input
    print(f"\n[4/4] Generating QE input...")
    
    # Calculate k-points from spacing
    import numpy as np
    a = np.linalg.norm(cell[0])
    b = np.linalg.norm(cell[1])
    nk_a = max(1, int(a / kspacing))
    nk_b = max(1, int(b / kspacing))
    
    # Build input
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
    ecutwfc = {ecutwfc:.1f}
    ecutrho = {ecutrho:.1f}
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
{symbols[0]:2s}  {196.966:8.3f}  Au_ONCV_PBE-1.0.oncvpsp.upf

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
                if hasattr(constraint, 'index') and i in constraint.index:
                    is_fixed = True
                elif hasattr(constraint, 'indices') and i in constraint.indices:
                    is_fixed = True
        
        fix_flag = "0 0 0" if is_fixed else "1 1 1"
        qe_input += f"{symbol:2s}  {pos[0]:12.6f}  {pos[1]:12.6f}  {pos[2]:12.6f}  {fix_flag}\n"
    
    # Add k-points
    qe_input += f"\nK_POINTS (automatic)\n {nk_a} {nk_b} 1  0 0 0\n"
    
    # Write file
    with open(output_file, 'w') as f:
        f.write(qe_input)
    
    print(f"  ✓ Wrote {output_file}\n")
    
    # Show content
    print(f"{'='*70}")
    print("GENERATED QE INPUT FILE")
    print(f"{'='*70}\n")
    print(qe_input)
    
    return qe_input


if __name__ == '__main__':
    # Use slab.xyz from previous step
    xyz_to_qe_input('slab.xyz', 'slab.in', ecutwfc=60, ecutrho=240)
