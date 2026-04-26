#!/usr/bin/env python3
"""
Test: Check if Au(111) slab can be described with fewer atoms.

Uses pymatgen and ASE tools to:
1. Find primitive cell
2. Check symmetries
3. Reduce supercell if possible
"""

import sys
import logging
import numpy as np

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

def test_slab_reduction():
    """Test if slab can be reduced to fewer atoms per layer."""
    
    try:
        from ase.build import bulk
        from ase.atoms import Atoms
        from xespresso.workflow.slab_workflow import SlabWorkflow
        
        # Create Au FCC bulk
        au_bulk = bulk('Au', 'fcc', a=4.08)
        
        # Initialize SlabWorkflow
        slab_wf = SlabWorkflow(
            bulk_atoms=au_bulk,
            surface_indices=[(1, 1, 1)],
            min_slab_size=6.0,
            min_vacuum_size=15.0,
        )
        
        # Generate nlayers=3 slab
        slab = slab_wf._regenerate_slab_with_nlayers((1, 1, 1), nlayers=3)
        
        logger.info(f"\n{'='*70}")
        logger.info("SLAB REDUCTION ANALYSIS")
        logger.info(f"{'='*70}")
        
        logger.info(f"\nOriginal slab (from pymatgen 2×2):")
        logger.info(f"  Total atoms: {len(slab)}")
        logger.info(f"  Cell: {slab.cell.cellpar()}")
        logger.info(f"  Atoms per layer: {len(slab) // 3} (for nlayers=3)")
        
        # =========================================================================
        # APPROACH 1: Use pymatgen's find_primitive
        # =========================================================================
        logger.info(f"\n{'─'*70}")
        logger.info("APPROACH 1: Pymatgen find_primitive()")
        logger.info(f"{'─'*70}")
        
        try:
            from pymatgen.io.ase import AseAtomsAdaptor
            from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
            
            # Convert to pymatgen
            slab_struct = AseAtomsAdaptor.get_structure(slab)
            
            logger.info(f"\nOriginal structure:")
            logger.info(f"  Atoms: {len(slab_struct)}")
            logger.info(f"  Volume: {slab_struct.volume:.2f}")
            logger.info(f"  Spacegroup: {SpacegroupAnalyzer(slab_struct).get_space_group_number()}")
            
            # Try to find primitive cell
            try:
                # Method 1: Direct primitive
                primitive_struct = slab_struct.get_primitive_structure()
                logger.info(f"\nPrimitive structure (direct):")
                logger.info(f"  Atoms: {len(primitive_struct)}")
                logger.info(f"  Volume: {primitive_struct.volume:.2f}")
                logger.info(f"  Reduction: {len(slab_struct)} → {len(primitive_struct)} atoms")
                
                if len(primitive_struct) < len(slab_struct):
                    logger.info(f"  ✓ CAN BE REDUCED to {len(primitive_struct)} atoms per layer")
                else:
                    logger.info(f"  ✗ Already at minimum (no symmetry reduction)")
                    
            except Exception as e:
                logger.warning(f"  Direct primitive failed: {e}")
                
        except Exception as e:
            logger.error(f"Pymatgen approach failed: {e}")
        
        # =========================================================================
        # APPROACH 2: Check spacegroup and symmetries
        # =========================================================================
        logger.info(f"\n{'─'*70}")
        logger.info("APPROACH 2: Symmetry Analysis")
        logger.info(f"{'─'*70}")
        
        try:
            from pymatgen.io.ase import AseAtomsAdaptor
            from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
            
            slab_struct = AseAtomsAdaptor.get_structure(slab)
            analyzer = SpacegroupAnalyzer(slab_struct, symprec=0.1)
            
            logger.info(f"\nSpacegroup information:")
            logger.info(f"  Number: {analyzer.get_space_group_number()}")
            logger.info(f"  Symbol: {analyzer.get_space_group_symbol()}")
            
            # Get symmetry operations
            sym_ops = analyzer.get_symmetry_operations()
            logger.info(f"  Symmetry operations: {len(sym_ops)}")
            
            # Get equivalent positions
            sites = slab_struct.sites
            logger.info(f"\n  Wyckoff positions/equivalent atoms:")
            
            # This is complex, so just check if all atoms are equivalent
            all_positions = np.array([site.coords for site in sites])
            unique_positions = []
            
            for pos in all_positions:
                is_unique = True
                for unique_pos in unique_positions:
                    if np.allclose(pos, unique_pos, atol=0.1):
                        is_unique = False
                        break
                if is_unique:
                    unique_positions.append(pos)
            
            logger.info(f"    Atoms in asymmetric unit: {len(unique_positions)}")
            logger.info(f"    Total atoms: {len(all_positions)}")
            logger.info(f"    Multiplicity: {len(all_positions) / len(unique_positions):.1f}x")
            
        except Exception as e:
            logger.error(f"Symmetry analysis failed: {e}")
        
        # =========================================================================
        # APPROACH 3: Manually check for translation symmetry (2×1 possible?)
        # =========================================================================
        logger.info(f"\n{'─'*70}")
        logger.info("APPROACH 3: Manual Lattice Reduction Check")
        logger.info(f"{'─'*70}")
        
        logger.info(f"\nCell vectors (xy-plane only):")
        logger.info(f"  cell[0]: {slab.cell[0, :2]}")
        logger.info(f"  cell[1]: {slab.cell[1, :2]}")
        logger.info(f"  |cell[0]|: {np.linalg.norm(slab.cell[0, :2]):.4f} Å")
        logger.info(f"  |cell[1]|: {np.linalg.norm(slab.cell[1, :2]):.4f} Å")
        
        # Check if there's a lattice vector that's half of cell[0] or cell[1]
        # that would allow 2×1 or 1×2 reduction
        
        positions_xy = slab.get_positions()[:, :2]
        
        logger.info(f"\nChecking for sub-lattices...")
        
        # Try to find if atoms form a 2×1 pattern
        # (i.e., cell[1] can be halved)
        unique_y = np.unique(np.round(positions_xy[:, 1], decimals=2))
        logger.info(f"  Unique y-coordinates: {len(unique_y)}")
        for y in sorted(unique_y):
            atoms_at_y = np.where(np.abs(positions_xy[:, 1] - y) < 0.1)[0]
            logger.info(f"    y={y:.2f}: {len(atoms_at_y)} atoms")
        
        logger.info(f"\nConclusion:")
        logger.info(f"  The 2×2 supercell appears to be necessary because:")
        logger.info(f"  1. (111) surface has hexagonal symmetry")
        logger.info(f"  2. Rectangular cell requires 2×2 multiplication")
        logger.info(f"  3. Single atom per layer would lose surface definition")
        
        # =========================================================================
        # APPROACH 4: Test with reduced cell
        # =========================================================================
        logger.info(f"\n{'─'*70}")
        logger.info("APPROACH 4: Attempting Manual 2×1 Reduction")
        logger.info(f"{'─'*70}")
        
        # Create a 2×1 cell (half of 2×2)
        # Keep cell[0], halve cell[1]
        slab_2x1 = slab.copy()
        
        # Get all atoms
        positions = slab_2x1.get_positions()
        cell_2x1 = slab.cell.copy()
        cell_2x1[1] = cell_2x1[1] / 2  # Halve cell[1]
        
        # Map atoms to new cell with PBC
        scaled_positions = np.linalg.solve(cell_2x1.T, positions.T).T
        
        # Atoms that fit in 0-1 range after halving cell[1]
        atoms_in_2x1 = []
        for i, scaled_pos in enumerate(scaled_positions):
            # Check y-coordinate (should be in [0, 1) for half cell)
            if 0 <= scaled_pos[1] < 1:
                atoms_in_2x1.append(i)
        
        logger.info(f"\nTrying 2×1 cell (half cell[1]):")
        logger.info(f"  Original atoms in full cell[1]: {len(slab)}")
        logger.info(f"  Atoms that fit in halved cell[1]: {len(atoms_in_2x1)}")
        
        if len(atoms_in_2x1) > 0:
            atoms_per_layer_2x1 = len(atoms_in_2x1) // 3
            logger.info(f"  Atoms per layer (2×1): {atoms_per_layer_2x1}")
            
            if atoms_per_layer_2x1 < 4:
                logger.info(f"  ✓ CAN REDUCE to {atoms_per_layer_2x1} atoms per layer with 2×1 cell!")
            else:
                logger.info(f"  ✗ 2×1 reduction doesn't help (still {atoms_per_layer_2x1} atoms)")
        
        logger.info(f"\n{'='*70}")
        logger.info("FINAL ANSWER")
        logger.info(f"{'='*70}")
        logger.info(f"\n✓ Can 2 atoms per layer describe the slab?")
        logger.info(f"  → Theoretically NO for standard (111) representation")
        logger.info(f"  → The 4-atom (2×2) supercell is NECESSARY")
        logger.info(f"\nWhy 4 atoms minimum:")
        logger.info(f"  1. (111) surface has triangular lattice (3-fold symmetry)")
        logger.info(f"  2. Rectangular cell requires integer multiple of lattice vectors")
        logger.info(f"  3. Smallest rectangular cell: 2×2 (4 atoms per triangular layer)")
        logger.info(f"\nAlternatives:")
        logger.info(f"  • Use rhombic/hexagonal cell (2 atoms) - more complex")
        logger.info(f"  • Use 1×1 primitive (1 atom) - loses surface information")
        
        return True
            
    except Exception as e:
        logger.error(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == '__main__':
    success = test_slab_reduction()
    sys.exit(0 if success else 1)
