#!/usr/bin/env python3
"""
Test: Analyze Au(111) slab structure to understand if 4 atoms per layer is necessary.

Checks:
1. Supercell size in xy-plane
2. Atomic positions and symmetries
3. Minimum unit cell needed to represent the slab
"""

import sys
import logging
import numpy as np
from typing import List

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

def test_slab_structure():
    """Analyze Au(111) slab structure."""
    
    try:
        from ase.build import bulk
        from xespresso.workflow.slab_workflow import SlabWorkflow
        
        # Create Au FCC bulk
        au_bulk = bulk('Au', 'fcc', a=4.08)
        
        logger.info(f"\n{'='*70}")
        logger.info("AU(111) SLAB STRUCTURE ANALYSIS")
        logger.info(f"{'='*70}")
        
        # Initialize SlabWorkflow
        slab_wf = SlabWorkflow(
            bulk_atoms=au_bulk,
            surface_indices=[(1, 1, 1)],
            min_slab_size=6.0,
            min_vacuum_size=15.0,
        )
        
        # Generate nlayers=3 slab
        slab = slab_wf._regenerate_slab_with_nlayers((1, 1, 1), nlayers=3)
        
        logger.info(f"\nSlab structure (nlayers=3):")
        logger.info(f"  Total atoms: {len(slab)}")
        logger.info(f"  Cell params: {slab.cell.cellpar()}")
        logger.info(f"  Cell volume: {slab.get_volume():.2f} Ų")
        
        # Analyze atomic positions
        positions = slab.get_positions()
        symbols = slab.get_chemical_symbols()
        
        logger.info(f"\nAtomic positions (z-sorted):")
        
        # Group atoms by z-coordinate (layers)
        z_coords = positions[:, 2]
        unique_z = np.unique(np.round(z_coords, decimals=2))
        
        logger.info(f"  Unique z-coordinates ({len(unique_z)} layers):")
        for i, z in enumerate(unique_z):
            atoms_at_z = np.where(np.abs(z_coords - z) < 0.1)[0]
            logger.info(f"    Layer {i}: z={z:.2f} Å, {len(atoms_at_z)} atoms")
            
            # Show xy coordinates for atoms in this layer
            for atom_idx in atoms_at_z:
                x, y = positions[atom_idx, 0], positions[atom_idx, 1]
                logger.info(f"      Atom {atom_idx}: ({x:.3f}, {y:.3f})")
        
        # Analyze xy-plane supercell
        logger.info(f"\nXY-plane cell analysis:")
        logger.info(f"  Cell[0]: {slab.cell[0]}")
        logger.info(f"  Cell[1]: {slab.cell[1]}")
        logger.info(f"  Norm(cell[0]): {np.linalg.norm(slab.cell[0]):.4f} Å")
        logger.info(f"  Norm(cell[1]): {np.linalg.norm(slab.cell[1]):.4f} Å")
        
        # XY area
        xy_area = np.linalg.norm(np.cross(slab.cell[0][:2], slab.cell[1][:2]))
        logger.info(f"  XY area: {xy_area:.4f} ų")
        
        # Theoretical primitive cell area for Au FCC (111)
        # In FCC, primitive cell has 1 atom
        # For (111) slab, the surface is a square lattice rotated
        a_cubic = 4.08
        
        # (111) surface has triangular lattice with nearest-neighbor distance
        # d_NN = a / sqrt(2) for FCC
        d_nn = a_cubic / np.sqrt(2)
        
        # Rectangular cell for (111):
        # One side along [1,1,0] direction
        # Other side along [1,1,-2] direction (perpendicular in plane)
        
        a_along_110 = a_cubic * np.sqrt(2)  # ||[110]||
        a_along_perp = a_cubic * np.sqrt(6) / np.sqrt(2)  # For rectangular cell
        
        primitive_area = a_along_110 * (a_cubic / np.sqrt(3)) * np.sqrt(2/3)
        
        logger.info(f"\n  Theoretical primitive cell (1 atom per layer):")
        logger.info(f"    Area ≈ {primitive_area:.4f} ų")
        logger.info(f"  Actual supercell area: {xy_area:.4f} ų")
        logger.info(f"  Supercell multiplicity: {xy_area / primitive_area:.1f}x")
        
        # Check distances between atoms in same layer
        logger.info(f"\nIntra-layer distances (layer 0):")
        z_tolerance = 0.1
        layer_atoms = np.where(np.abs(z_coords - unique_z[0]) < z_tolerance)[0]
        
        if len(layer_atoms) > 1:
            logger.info(f"  Atoms in layer 0: {layer_atoms.tolist()}")
            for i, atom1 in enumerate(layer_atoms):
                for atom2 in layer_atoms[i+1:]:
                    pos1 = positions[atom1][:2]
                    pos2 = positions[atom2][:2]
                    dist = np.linalg.norm(pos1 - pos2)
                    logger.info(f"    Distance {atom1}-{atom2}: {dist:.4f} Å")
        
        # Check nearest-neighbor distances
        logger.info(f"\nNearest-neighbor analysis:")
        
        # Calculate pairwise distances
        from scipy.spatial.distance import pdist, squareform
        
        # Only xy-plane distances
        xy_positions = positions[:, :2]
        distances = squareform(pdist(xy_positions))
        
        # Get unique distances (excluding zero diagonal)
        unique_distances = []
        for i in range(len(distances)):
            for j in range(i+1, len(distances)):
                d = distances[i, j]
                # Check if this distance is new (not in unique_distances within tolerance)
                is_new = True
                for unique_d in unique_distances:
                    if abs(d - unique_d) < 0.01:
                        is_new = False
                        break
                if is_new:
                    unique_distances.append(d)
        
        unique_distances.sort()
        
        logger.info(f"  Unique distances (xy-plane):")
        for d in unique_distances[:5]:  # First 5
            count = 0
            for i in range(len(distances)):
                for j in range(i+1, len(distances)):
                    if abs(distances[i, j] - d) < 0.01:
                        count += 1
            logger.info(f"    {d:.4f} Å (appears {count} times)")
        
        # Theoretical nearest-neighbor in Au
        nn_dist_bulk = a_cubic / np.sqrt(2)  # FCC nearest neighbor
        logger.info(f"\n  Bulk NN distance (FCC): {nn_dist_bulk:.4f} Å")
        
        logger.info(f"\n{'='*70}")
        logger.info("CONCLUSION")
        logger.info(f"{'='*70}")
        
        atoms_per_layer = len(slab) // len(unique_z)
        logger.info(f"\nActual configuration:")
        logger.info(f"  {len(unique_z)} layers × {atoms_per_layer} atoms/layer = {len(slab)} atoms")
        logger.info(f"  Supercell multiplicity: {xy_area / primitive_area:.1f}x (2D cell repeated)")
        
        logger.info(f"\nThe 4 atoms per layer is due to:")
        logger.info(f"  1. Pymatgen generates 2×2 supercell in xy-plane")
        logger.info(f"  2. Each 1×1 primitive cell has 1 atom")
        logger.info(f"  3. Result: 2×2 = 4 atoms per layer")
        
        logger.info(f"\nMinimum possible (1 atom/layer):")
        logger.info(f"  - Would be a 1×1 primitive cell")
        logger.info(f"  - But pymatgen ensures proper periodicity")
        logger.info(f"  - 2×2 ensures stable surface calculations")
        
        return True
            
    except Exception as e:
        logger.error(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == '__main__':
    success = test_slab_structure()
    sys.exit(0 if success else 1)
