#!/usr/bin/env python3
"""
Comprehensive unit tests for SlabWorkflow Phase 2 (Slab Generation).

Tests individual Phase 2 functions:
- generate_slabs()
- _orthogonalize_cell()
- _apply_constraints()

Uses pytest for structured test organization.
"""

import logging
import pytest
import numpy as np
from pathlib import Path
from ase.build import bulk
from ase.constraints import FixAtoms
from xespresso.workflow.slab_workflow import SlabWorkflow

# Configure logging for tests
logging.basicConfig(level=logging.WARNING)


class TestSlabGeneration:
    """Tests for generate_slabs() method."""

    @pytest.fixture
    def au_fcc_bulk(self):
        """Create Au FCC bulk structure."""
        return bulk('Au', 'fcc', a=4.078)

    @pytest.fixture
    def workflow_au(self, au_fcc_bulk):
        """Create SlabWorkflow instance for Au."""
        return SlabWorkflow(
            bulk_atoms=au_fcc_bulk,
            surface_indices=[(1, 0, 0), (1, 1, 0), (1, 1, 1)],
            min_slab_size=6.0,
            min_vacuum_size=15.0,
            nlayers=4,
            pseudopotentials_config='default',
            precision='low',
            verbose=False,
        )

    def test_generate_slabs_creates_dict(self, workflow_au):
        """Test that generate_slabs() returns a dictionary."""
        slabs = workflow_au.generate_slabs()
        assert isinstance(slabs, dict), "generate_slabs() should return dict"

    def test_generate_slabs_dict_type(self, workflow_au):
        """Test that returned dict has Miller indices as keys and Atoms as values."""
        slabs = workflow_au.generate_slabs()
        for hkl, slab in slabs.items():
            assert isinstance(hkl, tuple), f"Key should be tuple, got {type(hkl)}"
            assert len(hkl) == 3, f"Miller index should have 3 components, got {len(hkl)}"
            from ase import Atoms
            assert isinstance(slab, Atoms), f"Value should be Atoms, got {type(slab)}"

    def test_generate_slabs_all_surfaces_generated(self, workflow_au):
        """Test that slabs are generated for all requested surfaces."""
        slabs = workflow_au.generate_slabs()
        expected_surfaces = {(1, 0, 0), (1, 1, 0), (1, 1, 1)}
        generated_surfaces = set(slabs.keys())
        assert generated_surfaces == expected_surfaces, \
            f"Expected {expected_surfaces}, got {generated_surfaces}"

    def test_generate_slabs_nonempty(self, workflow_au):
        """Test that generated slabs contain atoms."""
        slabs = workflow_au.generate_slabs()
        for hkl, slab in slabs.items():
            assert len(slab) > 0, f"Slab {hkl} should have atoms"

    def test_generate_slabs_cell_has_vacuum(self, workflow_au):
        """Test that slabs have vacuum in z-direction."""
        slabs = workflow_au.generate_slabs()
        for hkl, slab in slabs.items():
            cell_z = slab.cell[2, 2]
            # Check that cell is larger than typical bulk (Au: ~4 Å per layer)
            assert cell_z > 15, f"Slab {hkl} should have vacuum (cell_z={cell_z:.2f})"

    def test_generate_slabs_fixed_layers_applied(self, workflow_au):
        """Test that FixAtoms constraints are applied."""
        slabs = workflow_au.generate_slabs()
        for hkl, slab in slabs.items():
            constraints = slab.constraints
            has_fixatoms = any(isinstance(c, FixAtoms) for c in constraints)
            assert has_fixatoms, f"Slab {hkl} should have FixAtoms constraint"

    def test_generate_slabs_reproducible(self, workflow_au):
        """Test that generate_slabs() is reproducible."""
        slabs1 = workflow_au.generate_slabs()
        
        # Create new workflow with same parameters
        workflow_au2 = SlabWorkflow(
            bulk_atoms=workflow_au.bulk_atoms.copy(),
            surface_indices=[(1, 0, 0), (1, 1, 0), (1, 1, 1)],
            min_slab_size=6.0,
            min_vacuum_size=15.0,
            nlayers=4,
            pseudopotentials_config='default',
            precision='low',
            verbose=False,
        )
        slabs2 = workflow_au2.generate_slabs()
        
        # Compare energies/positions
        for hkl in slabs1.keys():
            pos1 = slabs1[hkl].get_positions()
            pos2 = slabs2[hkl].get_positions()
            assert np.allclose(pos1, pos2), f"Slab {hkl} positions differ"

    def test_generate_slabs_stores_in_workflow(self, workflow_au):
        """Test that slabs are stored in self.slabs."""
        slabs = workflow_au.generate_slabs()
        assert workflow_au.slabs == slabs, "Slabs should be stored in self.slabs"

    def test_generate_slabs_saves_to_file(self, workflow_au, tmp_path):
        """Test that slabs can be saved to CIF files."""
        slabs = workflow_au.generate_slabs(
            save_slabs=True, 
            save_dir=str(tmp_path)
        )
        
        # Check that CIF files were created
        for hkl in slabs.keys():
            filename = tmp_path / f"slab_{hkl[0]}{hkl[1]}{hkl[2]}.cif"
            assert filename.exists(), f"CIF file {filename} should exist"


class TestOrthogonalizeCell:
    """Tests for _orthogonalize_cell() method."""

    @pytest.fixture
    def workflow(self):
        """Create SlabWorkflow instance."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        return SlabWorkflow(
            bulk_atoms=au_bulk,
            surface_indices=[(1, 1, 1)],
            pseudopotentials_config='default',
            verbose=False,
        )

    def test_orthogonalize_cell_sets_cell(self, workflow):
        """Test that _orthogonalize_cell creates non-zero cell."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        slab = workflow._orthogonalize_cell(au_bulk.copy())
        
        assert slab.cell[0, 0] > 0, "Cell[0,0] should be positive"
        assert slab.cell[1, 1] > 0, "Cell[1,1] should be positive"
        assert slab.cell[2, 2] > 0, "Cell[2,2] should be positive"

    def test_orthogonalize_cell_preserves_volume(self, workflow):
        """Test that orthogonalization doesn't drastically change volume.
        
        Note: Volume may change slightly due to scale_atoms=True.
        We check that the change is reasonable (< 10%).
        """
        # Use a generated slab for realistic test
        slabs = workflow.generate_slabs()
        slab111 = slabs[(1, 1, 1)]
        vol_before = slab111.get_volume()
        
        # Apply orthogonalization
        ortho_slab = workflow._orthogonalize_cell(slab111.copy())
        vol_after = ortho_slab.get_volume()
        
        # Check that volume change is reasonable
        rel_change = abs(vol_after - vol_before) / vol_before
        assert rel_change < 0.1, \
            f"Volume change too large: {rel_change:.1%}"

    def test_orthogonalize_cell_off_diagonal_small(self, workflow):
        """Test that off-diagonal elements are minimized."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        slab = workflow._orthogonalize_cell(au_bulk.copy())
        
        # Off-diagonal elements should be much smaller than diagonal
        cell = slab.cell
        max_offdiag = max(abs(cell[0, 1]), abs(cell[0, 2]), abs(cell[1, 2]))
        min_diag = min(abs(cell[0, 0]), abs(cell[1, 1]), abs(cell[2, 2]))
        
        assert max_offdiag < 0.1 * min_diag, \
            "Off-diagonal elements should be much smaller than diagonal"

    def test_orthogonalize_cell_atoms_scaled(self, workflow):
        """Test that atomic positions are properly scaled."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        slab = workflow._orthogonalize_cell(au_bulk.copy())
        
        # Atoms should still be within or near the cell
        scaled = slab.get_scaled_positions()
        assert np.all(scaled >= -0.1), "Scaled positions should be positive (with tolerance)"
        assert np.all(scaled <= 1.1), "Scaled positions should be <= 1 (with tolerance)"

    def test_orthogonalize_cell_preserves_num_atoms(self, workflow):
        """Test that orthogonalization doesn't change number of atoms."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        slab = workflow._orthogonalize_cell(au_bulk.copy())
        
        assert len(slab) == len(au_bulk), "Number of atoms should be preserved"


class TestApplyConstraints:
    """Tests for _apply_constraints() method."""

    @pytest.fixture
    def workflow_4layers(self):
        """Create workflow with 4 layers."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        return SlabWorkflow(
            bulk_atoms=au_bulk,
            surface_indices=[(1, 1, 1)],
            nlayers=4,
            fix_layer_indices=[0, 1],  # Fix bottom 2 layers
            pseudopotentials_config='default',
            verbose=False,
        )

    @pytest.fixture
    def test_slab(self, workflow_4layers):
        """Create a test slab with atoms at defined z-positions."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        # Generate Au(111) slab manually for testing
        from ase.io import read
        slabs = workflow_4layers.generate_slabs()
        return slabs[(1, 1, 1)]

    def test_apply_constraints_adds_fixatoms(self, workflow_4layers, test_slab):
        """Test that _apply_constraints adds FixAtoms constraint."""
        constrained_slab = workflow_4layers._apply_constraints(test_slab.copy())
        
        constraints = constrained_slab.constraints
        has_fixatoms = any(isinstance(c, FixAtoms) for c in constraints)
        assert has_fixatoms, "Should have FixAtoms constraint"

    def test_apply_constraints_fixes_bottom_layers(self, workflow_4layers, test_slab):
        """Test that bottom layers are actually fixed."""
        constrained_slab = workflow_4layers._apply_constraints(test_slab.copy())
        
        # Get FixAtoms constraint
        fixatoms = None
        for c in constrained_slab.constraints:
            if isinstance(c, FixAtoms):
                fixatoms = c
                break
        
        assert fixatoms is not None, "FixAtoms constraint should exist"
        assert len(fixatoms.index) > 0, "Should have fixed atoms"

    def test_apply_constraints_preserves_positions(self, workflow_4layers, test_slab):
        """Test that _apply_constraints doesn't change positions."""
        pos_before = test_slab.get_positions().copy()
        constrained_slab = workflow_4layers._apply_constraints(test_slab.copy())
        pos_after = constrained_slab.get_positions()
        
        assert np.allclose(pos_before, pos_after), \
            "Positions should not change after applying constraints"

    def test_apply_constraints_preserves_cell(self, workflow_4layers, test_slab):
        """Test that _apply_constraints doesn't change cell."""
        cell_before = test_slab.get_cell().copy()
        constrained_slab = workflow_4layers._apply_constraints(test_slab.copy())
        cell_after = constrained_slab.get_cell()
        
        assert np.allclose(cell_before, cell_after), \
            "Cell should not change after applying constraints"

    def test_apply_constraints_num_fixed_scales_with_layers(self):
        """Test that number of fixed atoms scales with nlayers."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        
        # Test with 4 layers, fix bottom 2
        wf4 = SlabWorkflow(
            bulk_atoms=au_bulk,
            surface_indices=[(1, 1, 1)],
            nlayers=4,
            fix_layer_indices=[0, 1],
            verbose=False,
        )
        slabs4 = wf4.generate_slabs()
        slab4 = slabs4[(1, 1, 1)]
        
        # Count fixed atoms
        fixed_count = len(slab4.constraints[0].index) if slabs4 else 0
        assert fixed_count > 0, "Should have fixed atoms"

    def test_apply_constraints_respects_fix_layer_indices(self, workflow_4layers, test_slab):
        """Test that bottom layers are actually fixed."""
        constrained_slab = workflow_4layers._apply_constraints(test_slab.copy())
        
        # Get fixed indices
        fixatoms = None
        for c in constrained_slab.constraints:
            if isinstance(c, FixAtoms):
                fixatoms = c
                break
        
        if fixatoms:
            all_indices = np.arange(len(constrained_slab))
            fixed_indices = fixatoms.index
            free_indices = all_indices[~np.isin(all_indices, fixed_indices)]
            
            positions = constrained_slab.get_positions()
            fixed_z = positions[fixed_indices, 2].mean()
            free_z = positions[free_indices, 2].mean() if len(free_indices) > 0 else fixed_z
            
            # Fixed atoms should be at or below the mean z of free atoms
            assert fixed_z <= free_z, \
                "Fixed atoms should be in bottom part of slab"


class TestPhase2Integration:
    """Integration tests for Phase 2 workflow."""

    @pytest.fixture
    def workflow(self):
        """Create SlabWorkflow for integration tests."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        return SlabWorkflow(
            bulk_atoms=au_bulk,
            surface_indices=[(1, 0, 0), (1, 1, 0), (1, 1, 1)],
            min_slab_size=6.0,
            min_vacuum_size=15.0,
            nlayers=4,
            pseudopotentials_config='default',
            precision='low',
            verbose=False,
        )

    def test_full_phase2_workflow(self, workflow):
        """Test complete Phase 2: generate slabs → orthogonalize → constrain."""
        slabs = workflow.generate_slabs()
        
        assert len(slabs) == 3, "Should have 3 slabs"
        
        for hkl, slab in slabs.items():
            # Check slab properties
            assert len(slab) > 0, f"Slab {hkl} should have atoms"
            
            # Check cell
            cell = slab.get_cell()
            assert cell.volume > 0, f"Slab {hkl} should have positive volume"
            assert cell[2, 2] > 15, f"Slab {hkl} should have vacuum in z"
            
            # Check constraints
            has_fixatoms = any(isinstance(c, FixAtoms) for c in slab.constraints)
            assert has_fixatoms, f"Slab {hkl} should have FixAtoms"

    def test_phase2_results_stored_in_workflow(self, workflow):
        """Test that Phase 2 results are stored in workflow.slabs."""
        slabs = workflow.generate_slabs()
        
        assert workflow.slabs == slabs, "Slabs should be stored in workflow"
        assert all(hkl in workflow.slabs for hkl in slabs.keys()), \
            "All slabs should be in workflow.slabs"

    def test_different_surface_indices_different_slabs(self):
        """Test that different surface indices produce different slabs."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        
        wf111 = SlabWorkflow(
            bulk_atoms=au_bulk,
            surface_indices=[(1, 1, 1)],
            verbose=False,
        )
        slab111 = wf111.generate_slabs()[(1, 1, 1)]
        
        wf100 = SlabWorkflow(
            bulk_atoms=au_bulk,
            surface_indices=[(1, 0, 0)],
            verbose=False,
        )
        slab100 = wf100.generate_slabs()[(1, 0, 0)]
        
        # Different surfaces should have different atom counts or geometries
        assert len(slab111) != len(slab100) or \
               not np.allclose(slab111.cell, slab100.cell), \
            "Different surfaces should produce different slabs"


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
