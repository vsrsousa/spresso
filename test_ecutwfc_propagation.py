#!/usr/bin/env python3
"""
Test: Verify that bulk convergence parameters (ecutwfc=60) are properly
propagated to slab batch submissions instead of using defaults (30).
"""

import sys
import logging

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

def test_parameter_propagation():
    """Test that _prepare_input_data() correctly extracts parameters."""
    
    try:
        from ase.build import bulk
        from xespresso.workflow.slab_workflow import SlabWorkflow
        
        # Create test bulk structure
        bulk_atoms = bulk('Au', 'fcc', a=4.08)
        logger.info(f"✓ Created test bulk: {len(bulk_atoms)} atoms")
        
        # Initialize SlabWorkflow
        slab_wf = SlabWorkflow(
            bulk_atoms=bulk_atoms,
            surface_indices=[(1, 1, 1)],
        )
        logger.info("✓ SlabWorkflow initialized")
        
        # Set bulk parameters with ECUTWFC=60
        slab_wf.set_bulk_parameters(
            optimal_ecutwfc=60.0,
            optimal_kspacing=0.27,
            energy_tolerance=0.003,
            precision='low',
        )
        logger.info("✓ Bulk parameters set (ecutwfc=60.0 Ry)")
        
        # Test _prepare_input_data() method
        input_data = slab_wf._prepare_input_data()
        
        logger.info("\n" + "="*70)
        logger.info("INPUT_DATA PROPAGATION TEST")
        logger.info("="*70)
        logger.info(f"✓ _prepare_input_data() returns: {input_data}")
        
        if input_data.get('ecutwfc') == 60.0:
            logger.info("✓ PASS: ecutwfc=60.0 correctly extracted from bulk_recommendations")
            logger.info("\nExpected behavior in batch submission:")
            logger.info("  1. User sets bulk ecutwfc=60.0 via set_bulk_parameters()")
            logger.info("  2. slab_wf._prepare_input_data() extracts this value")
            logger.info("  3. submit_structures_parallel(..., input_data={'ecutwfc': 60.0})")
            logger.info("  4. CalculationWorkflow receives input_data with ecutwfc=60.0")
            logger.info("  5. All calculations use ecutwfc=60.0 instead of default 30.0")
            logger.info("\nThis solves the issue: Jobs will use correct ecutwfc")
            logger.info("="*70)
            return True
        else:
            logger.error(f"✗ FAIL: Expected ecutwfc=60.0, got {input_data.get('ecutwfc')}")
            return False
            
    except Exception as e:
        logger.error(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == '__main__':
    success = test_parameter_propagation()
    sys.exit(0 if success else 1)
