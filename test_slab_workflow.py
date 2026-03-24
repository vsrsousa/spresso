#!/usr/bin/env python3
"""
Simple test of SlabWorkflow Phase 2 (slab generation).

Tests:
1. SlabWorkflow initialization
2. Slab generation for Au(100), Au(110), Au(111)
3. Validation of generated slabs
4. Saving slabs to disk
"""

import logging
from pathlib import Path
from ase.build import bulk
from xespresso.workflow.slab_workflow import SlabWorkflow

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s - %(name)s - %(message)s'
)
logger = logging.getLogger(__name__)


def test_slab_generation():
    """Test Phase 2: slab generation."""
    
    logger.info("="*70)
    logger.info("TEST: SlabWorkflow Phase 2 (Slab Generation)")
    logger.info("="*70)
    
    # Create bulk structure: Au FCC
    logger.info("\n1. Creating Au FCC bulk structure...")
    bulk_au = bulk('Au', 'fcc', a=4.078)
    logger.info(f"   ✓ Created: {bulk_au.get_chemical_formula()} {len(bulk_au)} atoms")
    
    # Initialize SlabWorkflow
    logger.info("\n2. Initializing SlabWorkflow...")
    try:
        wf = SlabWorkflow(
            bulk_atoms=bulk_au,
            surface_indices=[(1, 0, 0), (1, 1, 0), (1, 1, 1)],
            min_slab_size=6.0,
            min_vacuum_size=15.0,
            nlayers=4,
            pseudopotentials_config='default',
            protocol='moderate',
            precision='low',
            verbose=True,
        )
        logger.info(f"   ✓ Initialized: {repr(wf)}")
    except Exception as e:
        logger.error(f"   ✗ Failed to initialize: {e}")
        return False
    
    # Generate slabs
    logger.info("\n3. Generating slabs...")
    try:
        slabs = wf.generate_slabs(save_slabs=True, save_dir='./test_slabs/')
        logger.info(f"   ✓ Generated {len(slabs)} slabs")
        
        # Show slab info
        for hkl, slab in slabs.items():
            cell_lengths = slab.cell.lengths()
            logger.info(
                f"   - {hkl}: {len(slab)} atoms, "
                f"cell: {cell_lengths[0]:.3f} × {cell_lengths[1]:.3f} × {cell_lengths[2]:.3f} Å"
            )
    except Exception as e:
        logger.error(f"   ✗ Failed to generate slabs: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    # Validate slabs
    logger.info("\n4. Validating slabs...")
    try:
        assert len(slabs) > 0, "No slabs generated"
        for hkl, slab in slabs.items():
            assert len(slab) > 0, f"Slab {hkl} is empty"
            # Check cell dimensions
            cell = slab.cell.lengths()
            assert cell[0] > 0 and cell[1] > 0 and cell[2] > 0, f"Invalid cell for {hkl}"
        logger.info("   ✓ All slabs validated")
    except AssertionError as e:
        logger.error(f"   ✗ Validation failed: {e}")
        return False
    
    # Save results
    logger.info("\n5. Saving results...")
    try:
        wf.save_results(output_dir='./test_slabs/')
        logger.info("   ✓ Results saved to ./test_slabs/")
    except Exception as e:
        logger.error(f"   ✗ Failed to save results: {e}")
        return False
    
    logger.info("\n" + "="*70)
    logger.info("✓ TEST PASSED: Phase 2 (Slab Generation) works correctly")
    logger.info("="*70)
    
    return True


if __name__ == '__main__':
    success = test_slab_generation()
    exit(0 if success else 1)
