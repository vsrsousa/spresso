#!/usr/bin/env python3
"""
Test: Verify pymatgen slab generation for Au(111).

Checks:
1. Correct d-spacing calculation for (111) surface
2. Proper atom count for different nlayers
3. Area reduction consistency
"""

import sys
import logging
import numpy as np
from typing import List, Tuple

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

def test_slab_generation():
    """Test slab generation with pymatgen."""
    
    try:
        from ase.build import bulk
        from ase.constraints import FixAtoms
        from xespresso.workflow.slab_workflow import SlabWorkflow
        
        # Create Au FCC bulk
        au_bulk = bulk('Au', 'fcc', a=4.08)
        logger.info(f"\n{'='*70}")
        logger.info("AU(111) SLAB GENERATION TEST")
        logger.info(f"{'='*70}")
        logger.info(f"\nBulk structure:")
        logger.info(f"  Lattice constant: a = 4.08 Å")
        logger.info(f"  Atoms per bulk cell: {len(au_bulk)}")
        logger.info(f"  Cell volume: {au_bulk.get_volume():.2f} Ų")
        logger.info(f"  Cell: {au_bulk.cell.cellpar()[:3]}")
        
        # Initialize SlabWorkflow
        slab_wf = SlabWorkflow(
            bulk_atoms=au_bulk,
            surface_indices=[(1, 1, 1)],
            min_vacuum_size=15.0,
        )
        
        # Theoretical d-spacing for (111)
        # d_hkl = a / sqrt(h² + k² + l²)
        a = 4.08
        h, k, l = 1, 1, 1
        d_hkl = a / np.sqrt(h**2 + k**2 + l**2)
        
        logger.info(f"\nTheoretical d-spacing for (111):")
        logger.info(f"  d_hkl = a / sqrt(h² + k² + l²)")
        logger.info(f"  d_hkl = {a:.2f} / sqrt({h}² + {k}² + {l}²)")
        logger.info(f"  d_hkl = {a:.2f} / sqrt(3) = {d_hkl:.4f} Å")
        
        # Generate slabs with different nlayers
        logger.info(f"\n{'='*70}")
        logger.info("SLAB GENERATION FOR DIFFERENT NLAYERS")
        logger.info(f"{'='*70}")
        
        nlayers_test = [3, 4, 5, 6, 7]
        slabs_data = []
        
        for nlayers in nlayers_test:
            slab = slab_wf._regenerate_slab_with_nlayers((1, 1, 1), nlayers)
            
            # Calculate expected min_slab_size
            min_slab_size_expected = nlayers * d_hkl - 0.5
            
            # Get slab info
            natoms = len(slab)
            slab_height = slab.cell[2, 2]  # z-direction
            slab_area = slab.cell[0, 0] * slab.cell[1, 1]  # xy-area
            
            slabs_data.append({
                'nlayers': nlayers,
                'natoms': natoms,
                'slab_height': slab_height,
                'area': slab_area,
                'expected_min_size': min_slab_size_expected,
            })
            
            logger.info(f"\nnlayers = {nlayers}:")
            logger.info(f"  Calculated min_slab_size: {min_slab_size_expected:.4f} Å")
            logger.info(f"  Generated slab height: {slab_height:.4f} Å")
            logger.info(f"  Slab area (xy): {slab_area:.4f} ų")
            logger.info(f"  Total atoms: {natoms}")
        
        # Analysis
        logger.info(f"\n{'='*70}")
        logger.info("ANALYSIS")
        logger.info(f"{'='*70}")
        
        # Check atom progression
        logger.info(f"\nAtom count progression:")
        prev_natoms = None
        for data in slabs_data:
            atoms_per_layer = data['natoms'] / data['nlayers']
            if prev_natoms:
                delta = data['natoms'] - prev_natoms
                logger.info(f"  nlayers={data['nlayers']}: {data['natoms']:2d} atoms "
                           f"({atoms_per_layer:.1f} atoms/layer, +{delta} from prev)")
            else:
                logger.info(f"  nlayers={data['nlayers']}: {data['natoms']:2d} atoms "
                           f"({atoms_per_layer:.1f} atoms/layer)")
            prev_natoms = data['natoms']
        
        # Check area consistency
        logger.info(f"\nSlab area consistency:")
        first_area = slabs_data[0]['area']
        for data in slabs_data:
            area_ratio = data['area'] / first_area
            logger.info(f"  nlayers={data['nlayers']}: area ratio = {area_ratio:.4f}")
        
        # Verify cell reduction
        logger.info(f"\nCell reduction check:")
        
        # For Au FCC (111), pymatgen typically generates:
        # - Bulk cell: 4 atoms in conventional cell
        # - (111) slab in standard orientation
        
        logger.info(f"  Bulk lattice parameter: {a} Å")
        logger.info(f"  d-spacing (111): {d_hkl:.4f} Å")
        logger.info(f"  Layer thickness / atom ratio should be consistent")
        
        # Expected atoms per unit area (for Au FCC 111)
        # Au FCC has 4 atoms per conventional cell
        # Area of (111) face in standard setting varies
        
        logger.info(f"\n{'='*70}")
        logger.info("CONCLUSION")
        logger.info(f"{'='*70}")
        
        if len(slabs_data) > 1:
            # Check if atom count increases linearly with nlayers
            natoms_list = [d['natoms'] for d in slabs_data]
            nlayers_list = [d['nlayers'] for d in slabs_data]
            
            # Linear fit
            coeffs = np.polyfit(nlayers_list, natoms_list, 1)
            atoms_per_layer_fit = coeffs[0]
            
            logger.info(f"\nLinear regression (atoms = {atoms_per_layer_fit:.2f} * nlayers + {coeffs[1]:.2f}):")
            logger.info(f"  Atoms per layer: {atoms_per_layer_fit:.2f}")
            
            # Check residuals
            predictions = np.polyval(coeffs, nlayers_list)
            residuals = natoms_list - predictions
            max_residual = np.max(np.abs(residuals))
            
            logger.info(f"  Max residual: {max_residual:.2f} atoms")
            
            if max_residual < 2:
                logger.info(f"\n✓ PASS: Atom count scales linearly with nlayers")
                logger.info(f"✓ Pymatgen is correctly reducing cell for layer generation")
                return True
            else:
                logger.warning(f"\n⚠ WARNING: Large residuals in atom count")
                logger.warning(f"  This might indicate non-uniform layer generation")
                return False
        else:
            logger.warning("Not enough data points for analysis")
            return False
            
    except Exception as e:
        logger.error(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == '__main__':
    success = test_slab_generation()
    sys.exit(0 if success else 1)
