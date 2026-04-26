#!/usr/bin/env python3
"""
Add SpgLib symmetry analysis to SlabWorkflow.

This shows how to analyze and optionally reduce slab size using SpgLib.
"""

import sys
sys.path.insert(0, '/home/vinicius/projects/spresso')

import numpy as np
from ase.atoms import Atoms
from typing import Optional, Dict, Tuple

def analyze_slab_symmetry(
    slab: Atoms,
    symprec: float = 1e-4,
    verbose: bool = True
) -> Dict:
    """
    Analyze slab symmetry using SpgLib.
    
    Parameters
    ----------
    slab : ase.Atoms
        Slab structure to analyze
    symprec : float, default=1e-4
        Symmetry tolerance in Ångströms
        - 1e-5: Very tight (few reductions)
        - 1e-4: Recommended for slabs
        - 1e-3: Liberal (maximum reduction)
    verbose : bool
        Print symmetry information
        
    Returns
    -------
    Dict with keys:
        - 'space_group': Space group number (e.g., 164)
        - 'pointgroup': Point group symbol (e.g., '-3m')
        - 'n_operations': Number of symmetry operations
        - 'primitive_atoms': Number of atoms in primitive cell
        - 'reduction_factor': natoms_slab / natoms_primitive
        - 'can_reduce': Boolean, True if reduction is possible
        - 'primitive_cell': Primitive cell array (if reduction possible)
        - 'primitive_positions': Primitive positions (if reduction possible)
        - 'primitive_numbers': Primitive atomic numbers (if reduction possible)
    """
    try:
        import spglib
    except ImportError:
        raise ImportError("SpgLib not installed. Install with: pip install spglib")
    
    # Prepare data for SpgLib
    cell = slab.cell.array
    positions = slab.get_scaled_positions()
    numbers = slab.get_atomic_numbers()
    
    result = {
        'space_group': None,
        'pointgroup': None,
        'n_operations': None,
        'primitive_atoms': None,
        'reduction_factor': 1.0,
        'can_reduce': False,
        'primitive_cell': None,
        'primitive_positions': None,
        'primitive_numbers': None,
    }
    
    # Get space group info
    dataset = spglib.get_symmetry_dataset((cell, positions, numbers), symprec=symprec)
    
    if dataset:
        result['space_group'] = int(dataset.number)
        result['pointgroup'] = dataset.pointgroup
        result['n_operations'] = len(dataset.rotations)
        
        # Get primitive cell
        primitive_cell = spglib.find_primitive((cell, positions, numbers), symprec=symprec)
        if primitive_cell:
            prim_cell, prim_pos, prim_nums = primitive_cell
            prim_natoms = len(prim_nums)
            reduction = len(slab) / prim_natoms
            
            result['primitive_atoms'] = prim_natoms
            result['reduction_factor'] = reduction
            result['can_reduce'] = reduction > 1.01  # Only if significant reduction
            result['primitive_cell'] = prim_cell
            result['primitive_positions'] = prim_pos
            result['primitive_numbers'] = prim_nums
    
    if verbose:
        print(f"\n[SPGLIB SYMMETRY ANALYSIS]")
        print(f"  Original slab: {len(slab)} atoms")
        print(f"  Space group: #{result['space_group']} ({result['pointgroup']})")
        print(f"  Symmetry operations: {result['n_operations']}")
        if result['can_reduce']:
            print(f"  ✓ Primitive cell: {result['primitive_atoms']} atoms ({result['reduction_factor']:.1f}× reduction)")
            print(f"  → Can use smaller supercell for faster calculations")
        else:
            print(f"  ✗ No significant reduction possible (already primitive)")
    
    return result


def create_primitive_slab_from_analysis(analysis_result: Dict) -> Optional[Atoms]:
    """Create ASE Atoms object from SpgLib analysis result."""
    if not analysis_result['can_reduce']:
        return None
    
    from ase import Atoms
    
    cell = analysis_result['primitive_cell']
    positions = analysis_result['primitive_positions']
    numbers = analysis_result['primitive_numbers']
    
    # Map atomic numbers to symbols
    from ase.data import atomic_numbers
    symbols = {v: k for k, v in atomic_numbers.items()}
    
    symbols_list = [symbols.get(int(n), f'X{n}') for n in numbers]
    
    prim_slab = Atoms(
        symbols=symbols_list,
        positions=positions,
        cell=cell,
        pbc=(True, True, True)  # Periodic in all directions
    )
    
    return prim_slab


# Test
if __name__ == "__main__":
    from ase.io import read
    from ase.build import surface
    
    print("="*70)
    print("INTEGRATING SPGLIB INTO SLABWORKFLOW")
    print("="*70)
    
    # Test with different slab sizes
    bulk_au = read('au_bulk.cif')
    
    for nlayers in [3, 4, 6, 8, 12]:
        slab = surface(bulk_au, (1, 1, 1), layers=nlayers, vacuum=5)
        print(f"\n{'='*70}")
        print(f"Au(111) with {nlayers} layers:")
        
        analysis = analyze_slab_symmetry(slab, symprec=1e-4)
        
        if analysis['can_reduce']:
            prim = create_primitive_slab_from_analysis(analysis)
            print(f"  Primitive structure: {len(prim)} atoms")
            print(f"  Cell a,b,c: {prim.cell.cellpar()[:3]}")
        
    print("\n" + "="*70)
    print("INTEGRATION TIPS")
    print("="*70)
    print("""
✓ Add to SlabWorkflow as optional method:
    slab_wf.analyze_symmetry(surface_index, symprec=1e-4)
    slab_wf.reduce_to_primitive(surface_index)

✓ Use cases:
  1. Initial geometry optimization (no surface relaxation)
  2. Bulk properties (EOS, bulk modulus)
  3. Fast test runs before full convergence
  
✗ Don't use for:
  1. Surface energy calculations (need full slab)
  2. Surface relaxation (simetria quebra)
  3. Adsorption studies (superficial atoms não simétricos)

✓ Recommended workflow:
  Phase 1: Bulk convergence (with primitive or supercell, doesn't matter)
  Phase 3: Slab convergence (use FULL slab, not primitive!)
  Phase 4: Relaxation (use FULL slab for accurate surface)
  
  → Primitive cell mostly useful for Phase 1 bulk and initial tests
""")
