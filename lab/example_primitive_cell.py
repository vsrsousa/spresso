#!/usr/bin/env python3
"""
Example: Using primitive cell option in SlabWorkflow

Shows how to generate Au(111) slabs with reduced (primitive) cells
instead of supercells.
"""

import sys
import logging
from ase.build import bulk
from xespresso.workflow.slab_workflow import SlabWorkflow

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

def example_primitive_cell():
    """Example: Generate slab with primitive cell."""
    
    logger.info("\n" + "="*70)
    logger.info("EXAMPLE: USING PRIMITIVE CELL OPTION")
    logger.info("="*70)
    
    au_bulk = bulk('Au', 'fcc', a=4.08)
    
    # =========================================================================
    # OPTION 1: SUPERCELL (2×2) - DEFAULT
    # =========================================================================
    logger.info("\n1️⃣ SUPERCELL (2×2) - Standard option:")
    
    slab_wf_supercell = SlabWorkflow(
        bulk_atoms=au_bulk,
        surface_indices=[(1, 1, 1)],
        use_primitive_cell=False,  # ← DEFAULT
    )
    
    logger.info("\n   Generating Au(111) slab with nlayers=[3,4,6,8]:")
    
    for nlayers in [3, 4, 6, 8]:
        slab = slab_wf_supercell._regenerate_slab_with_nlayers((1, 1, 1), nlayers)
        logger.info(f"   nlayers={nlayers}: {len(slab):2d} atoms (supercell 2×2)")
    
    # =========================================================================
    # OPTION 2: PRIMITIVE CELL (1×1) - REDUCED
    # =========================================================================
    logger.info("\n2️⃣ PRIMITIVE CELL (1×1) - Reduced (4× faster!):")
    
    slab_wf_primitive = SlabWorkflow(
        bulk_atoms=au_bulk,
        surface_indices=[(1, 1, 1)],
        use_primitive_cell=True,  # ← ACTIVATE THIS!
    )
    
    logger.info("\n   Generating Au(111) slab with nlayers=[3,4,6,8]:")
    
    for nlayers in [3, 4, 6, 8]:
        slab = slab_wf_primitive._regenerate_slab_with_nlayers((1, 1, 1), nlayers)
        logger.info(f"   nlayers={nlayers}: {len(slab):2d} atoms (primitive 1×1)")
    
    # =========================================================================
    # COMPARISON
    # =========================================================================
    logger.info("\n" + "="*70)
    logger.info("COMPARISON")
    logger.info("="*70)
    
    logger.info(f"\nSupercell (2×2):  atoms = 4 × nlayers")
    logger.info(f"Primitive (1×1):  atoms = 1 × nlayers")
    logger.info(f"\nComputational cost reduction: ~4× faster with primitive")
    
    logger.info("\n" + "="*70)
    logger.info("USAGE IN YOUR CODE")
    logger.info("="*70)
    
    logger.info("""
# For supercell (default, stable):
slab_wf = SlabWorkflow(
    bulk_atoms=au_bulk,
    surface_indices=[(1, 1, 1)],
    use_primitive_cell=False,  # ← explicit (not necessary)
)

# For primitive cell (faster, 75% fewer atoms):
slab_wf = SlabWorkflow(
    bulk_atoms=au_bulk,
    surface_indices=[(1, 1, 1)],
    use_primitive_cell=True,   # ← SET THIS!
)

# Then use normally:
slab = slab_wf._regenerate_slab_with_nlayers((1, 1, 1), nlayers=4)
# With primitive: 4 atoms instead of 16!
""")


if __name__ == '__main__':
    try:
        example_primitive_cell()
        logger.info("\n✓ Example completed successfully!")
    except Exception as e:
        logger.error(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
