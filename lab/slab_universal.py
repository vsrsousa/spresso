#!/usr/bin/env python3
"""
Gerador Universal de Slabs - Funciona com QUALQUER estrutura.

Usa pymatgen para máxima compatibilidade, mas com lll_reduce=False
para evitar células estranhas.

Usage:
    python slab_universal.py au_bulk.cif --surface 111 --nlayers 3
    python slab_universal.py cu_bulk.cif --surface 100 --nlayers 4
"""

import os
import sys
import argparse
import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def generate_slab_universal(
    structure_file,
    surface_indices,
    nlayers=4,
    vacuum=5.0,
    fix_fraction=0.4,
):
    """
    Gera slab para QUALQUER estrutura usando pymatgen.
    
    Evita o problema de células estranhas com lll_reduce=False.
    """
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.core.surface import SlabGenerator
    
    print(f"[1/4] Loading structure: {structure_file}")
    
    # Carregar estrutura (ASE - suporta CIF, POSCAR, etc)
    atoms = read(structure_file)
    formula = atoms.get_chemical_formula()
    print(f"  ✓ Loaded: {formula} ({len(atoms)} atoms)")
    
    # Converter para pymatgen
    print(f"\n[2/4] Preparing for slab generation...")
    adaptor = AseAtomsAdaptor()
    pmg_structure = adaptor.get_structure(atoms)
    print(f"  ✓ Structure prepared")
    print(f"    Lattice: {pmg_structure.lattice}")
    print(f"    Formula: {pmg_structure.composition}")
    
    # Gerar slab
    print(f"\n[3/4] Generating slab {surface_indices}...")
    try:
        slab_gen = SlabGenerator(
            pmg_structure,
            surface_indices,
            min_slab_size=nlayers * 2.5,
            min_vacuum_size=vacuum,
            center_slab=True,
            lll_reduce=False,  # ← IMPORTANTE: Evita células rotacionadas/estranhas
        )
        
        slab_pmg = slab_gen.get_slab()
        print(f"  ✓ Slab generated: {len(slab_pmg)} atoms")
        
        # Converter de volta para ASE
        slab = adaptor.get_atoms(slab_pmg)
        
    except Exception as e:
        print(f"  ✗ Error: {e}")
        raise
    
    # Normalizar célula (remove negativos se houver)
    print(f"\n[4/4] Normalizing cell...")
    slab = normalize_cell(slab)
    
    # Aplicar constraints
    slab, indices_fixed = apply_constraints(slab, fix_fraction)
    print(f"  ✓ Fixed {len(indices_fixed)}/{len(slab)} atoms")
    
    return slab, indices_fixed


def normalize_cell(atoms):
    """Remove valores negativos e orienta corretamente."""
    cell = atoms.get_cell()
    positions = atoms.get_positions()
    
    # Se z-cell for negativo, inverte
    if cell[2, 2] < 0:
        cell[2, 2] = abs(cell[2, 2])
        positions[:, 2] = abs(positions[:, 2])
    
    # Se a ou b tiver x negativo, inverte
    for i in range(2):
        if cell[i, 0] < 0:
            cell[i] = -cell[i]
    
    atoms.set_cell(cell, scale_atoms=False)
    atoms.set_positions(positions)
    
    return atoms


def apply_constraints(slab, fix_fraction=0.4):
    """Fixa camadas de baixo."""
    z_pos = slab.get_positions()[:, 2]
    z_min, z_max = z_pos.min(), z_pos.max()
    z_range = z_max - z_min
    
    fix_threshold = z_min + z_range * fix_fraction
    indices = [i for i, z in enumerate(z_pos) if z < fix_threshold]
    
    if indices:
        slab.set_constraint(FixAtoms(indices=indices))
    
    return slab, indices


def main():
    parser = argparse.ArgumentParser(
        description='Universal Slab Generator - Works with ANY structure'
    )
    parser.add_argument('structure_file', help='Structure file (CIF, POSCAR, etc)')
    parser.add_argument('--surface', default='111', help='Surface (h,k,l)')
    parser.add_argument('--nlayers', type=int, default=4, help='Number of layers')
    parser.add_argument('--vacuum', type=float, default=5.0, help='Vacuum (Å)')
    parser.add_argument('--output', default='slab_universal.in', help='Output file')
    parser.add_argument('--format', default='crystal', choices=['crystal', 'angstrom'],
                       help='Coordinate format')
    
    args = parser.parse_args()
    
    # Parse surface
    surface = tuple(int(x) for x in args.surface)
    
    print("="*70)
    print("UNIVERSAL SLAB GENERATOR")
    print("="*70 + "\n")
    
    try:
        # Generate slab
        slab, indices_fixed = generate_slab_universal(
            args.structure_file,
            surface,
            nlayers=args.nlayers,
            vacuum=args.vacuum,
        )
        
        # Get positions and cell
        cell = slab.get_cell()
        positions_angstrom = slab.get_positions()
        symbols = slab.get_chemical_symbols()
        
        # Convert to crystal coordinates if requested
        if args.format == 'crystal':
            inv_cell = np.linalg.inv(cell)
            positions = positions_angstrom @ inv_cell.T
            coord_type = 'crystal'
        else:
            positions = positions_angstrom
            coord_type = 'angstrom'
        
        print(f"\n✓ Generated {len(slab)}-atom {surface} slab")
        print(f"  Cell: {np.linalg.norm(cell[0]):.3f} × {np.linalg.norm(cell[1]):.3f} × {np.linalg.norm(cell[2]):.3f} Ų")
        print(f"  Vacuum: {slab.get_cell()[2,2] - (slab.get_positions()[:,2].max() - slab.get_positions()[:,2].min()):.3f} Å")
        print(f"  Fixed atoms: {len(indices_fixed)}/{len(slab)}")
        
        # Build QE input
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
    ntyp = {len(set(symbols))}
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
"""
        
        # Add atomic species
        unique_symbols = sorted(set(symbols))
        for symbol in unique_symbols:
            qe_input += f"{symbol:2s}  None  {symbol}.psp8\n"
        
        qe_input += f"""
CELL_PARAMETERS (angstrom)
  {cell[0,0]:12.6f}  {cell[0,1]:12.6f}  {cell[0,2]:12.6f}
  {cell[1,0]:12.6f}  {cell[1,1]:12.6f}  {cell[1,2]:12.6f}
  {cell[2,0]:12.6f}  {cell[2,1]:12.6f}  {cell[2,2]:12.6f}

ATOMIC_POSITIONS ({coord_type})
"""
        
        # Add atoms
        for i, (symbol, pos) in enumerate(zip(symbols, positions)):
            is_fixed = False
            if hasattr(slab, 'constraints') and slab.constraints:
                for constraint in slab.constraints:
                    if hasattr(constraint, 'indices') and i in constraint.indices:
                        is_fixed = True
            
            fix_flag = "0 0 0" if is_fixed else "1 1 1"
            qe_input += f"{symbol:2s}  {pos[0]:12.6f}  {pos[1]:12.6f}  {pos[2]:12.6f}  {fix_flag}\n"
        
        qe_input += f"\nK_POINTS (automatic)\n 7 7 1  0 0 0\n"
        
        # Write file
        with open(args.output, 'w') as f:
            f.write(qe_input)
        
        # Save structure for visualization
        write('slab_universal.xyz', slab)
        
        print(f"\n✓ Output files:")
        print(f"  QE input: {args.output}")
        print(f"  Visualization: slab_universal.xyz")
        print(f"\nRun: pw.x < {args.output}\n")
        
        return 0
        
    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
