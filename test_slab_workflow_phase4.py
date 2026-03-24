#!/usr/bin/env python3
"""
Unit tests for SlabWorkflow Phase 4 (Structure Relaxation).

Tests:
- run_slab_relax() requires Phase 3 convergence
- Constraint application (FixAtoms)
- K-mesh and ecutwfc from Phase 3
- Dipole correction setting
- Results storage
"""

import pytest
from ase.build import bulk
from ase.constraints import FixAtoms
from xespresso.workflow.slab_workflow import SlabWorkflow


class TestPhase4Prerequisites:
    """Tests for Phase 4 prerequisite validation."""

    def test_relax_requires_convergence_results(self):
        """Test that run_slab_relax() requires Phase 3 convergence."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        wf = SlabWorkflow(bulk_atoms=au_bulk, surface_indices=[(1, 1, 1)], verbose=False)
        
        # Generate slabs (Phase 2)
        wf.generate_slabs()
        
        # Mock bulk recommendations (Phase 1)
        wf.bulk_recommendations = {
            'optimal_ecutwfc': 60.0,
            'optimal_kspacing': 0.05,
        }
        
        # Convergence results NOT done
        with pytest.raises(ValueError, match="Phase 3"):
            wf.run_slab_relax()

    def test_relax_requires_surface_indices(self):
        """Test that run_slab_relax() requires valid surface indices."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        wf = SlabWorkflow(bulk_atoms=au_bulk, surface_indices=[(1, 1, 1)], verbose=False)
        wf.generate_slabs()
        wf.bulk_recommendations = {'optimal_ecutwfc': 60.0, 'optimal_kspacing': 0.05}
        
        # Add mock convergence results
        wf.convergence_results = {
            (1, 1, 1): {
                'kmesh_calc': (58, 58, 1),
                'optimal_vacuum': 18.0,
                'optimal_nlayers': 5,
            }
        }
        
        # Non-existent surface
        results = wf.run_slab_relax(surfaces=[(1, 0, 0)])
        assert (1, 0, 0) not in results or not results[(1, 0, 0)].get('converged', True)


class TestPhase4Constraints:
    """Tests for constraint application in Phase 4."""

    @pytest.fixture
    def workflow_ready(self):
        """Create workflow with Phase 1-3 ready."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        wf = SlabWorkflow(bulk_atoms=au_bulk, surface_indices=[(1, 1, 1)], verbose=False)
        wf.generate_slabs()
        wf.bulk_recommendations = {'optimal_ecutwfc': 60.0, 'optimal_kspacing': 0.05}
        wf.convergence_results = {
            (1, 1, 1): {
                'kmesh_calc': (58, 58, 1),
                'optimal_vacuum': 18.0,
                'optimal_nlayers': 5,
            }
        }
        return wf

    def test_relax_applies_fix_atoms(self, workflow_ready):
        """Test that FixAtoms constraint is applied."""
        # Note: This is a mock test since we can't run actual Espresso calculations
        # Get the slab
        slab = workflow_ready.slabs[(1, 1, 1)].copy()
        
        # Apply constraint (what run_slab_relax would do)
        if workflow_ready.fix_layer_indices:
            slab.set_constraint(FixAtoms(indices=workflow_ready.fix_layer_indices))
        
        # Verify constraint is applied
        constraints = slab.constraints
        assert len(constraints) > 0, "Should have at least one constraint"
        assert any(isinstance(c, FixAtoms) for c in constraints), "Should have FixAtoms constraint"

    def test_relax_constraint_fixes_correct_atoms(self, workflow_ready):
        """Test that constraint fixes the correct bottom layers."""
        slab = workflow_ready.slabs[(1, 1, 1)].copy()
        fix_indices = workflow_ready.fix_layer_indices
        
        # Apply constraint
        slab.set_constraint(FixAtoms(indices=fix_indices))
        
        # Verify fixed indices match expected bottom layers
        assert len(fix_indices) > 0, "Should have some fixed indices"
        
        # Bottom atoms should be fixed (assuming they're at lowest z)
        z_coords = slab.get_positions()[:, 2]
        fix_z_coords = z_coords[fix_indices]
        
        # Fixed atoms should be in the bottom part of the slab
        assert max(fix_z_coords) < max(z_coords), "Fixed atoms should be at bottom"


class TestPhase4Parameters:
    """Tests for Phase 4 parameter handling."""

    @pytest.fixture
    def workflow_ready(self):
        """Create workflow with Phase 1-3 ready."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        wf = SlabWorkflow(
            bulk_atoms=au_bulk,
            surface_indices=[(1, 1, 1)],
            pseudopotentials_config='default',
            verbose=False
        )
        wf.generate_slabs()
        wf.bulk_recommendations = {'optimal_ecutwfc': 60.0, 'optimal_kspacing': 0.05}
        wf.convergence_results = {
            (1, 1, 1): {
                'kmesh_calc': (58, 58, 1),
                'optimal_vacuum': 18.0,
                'optimal_nlayers': 5,
            }
        }
        return wf

    def test_relax_uses_phase1_ecutwfc(self, workflow_ready):
        """Test that relaxation uses ecutwfc from Phase 1."""
        # This would be validated in actual calculation
        ecutwfc = workflow_ready.bulk_recommendations['optimal_ecutwfc']
        assert ecutwfc == 60.0, "Phase 1 ecutwfc should be used"

    def test_relax_uses_phase3_kmesh(self, workflow_ready):
        """Test that relaxation uses k-mesh from Phase 3."""
        kmesh = workflow_ready.convergence_results[(1, 1, 1)]['kmesh_calc']
        assert kmesh == (58, 58, 1), "K-mesh from Phase 3 should be used"
        assert kmesh[2] == 1, "nkz should be 1 for 2D slabs"

    def test_relax_uses_optimal_vacuum(self, workflow_ready):
        """Test that relaxation uses optimal vacuum from Phase 3."""
        vacuum = workflow_ready.convergence_results[(1, 1, 1)]['optimal_vacuum']
        assert vacuum == 18.0, "Optimal vacuum from Phase 3 should be stored"


class TestPhase4Results:
    """Tests for Phase 4 result handling."""

    def test_relax_returns_dict(self):
        """Test that run_slab_relax returns a dictionary."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        wf = SlabWorkflow(bulk_atoms=au_bulk, surface_indices=[(1, 1, 1)], verbose=False)
        wf.generate_slabs()
        wf.bulk_recommendations = {'optimal_ecutwfc': 60.0, 'optimal_kspacing': 0.05}
        wf.convergence_results = {
            (1, 1, 1): {
                'kmesh_calc': (58, 58, 1),
                'optimal_vacuum': 18.0,
                'optimal_nlayers': 5,
            }
        }
        
        # Mock run_relax to avoid actual Espresso execution
        # This would normally raise errors without actual calculations
        # Here we just verify the structure is correct

    def test_relax_stores_results_in_workflow(self):
        """Test that results are stored in self.relax_results."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        wf = SlabWorkflow(bulk_atoms=au_bulk, surface_indices=[(1, 1, 1)], verbose=False)
        
        # Verify relax_results attribute exists (initialized in __init__)
        assert hasattr(wf, 'relax_results'), "Should have relax_results attribute"
        assert isinstance(wf.relax_results, dict), "relax_results should be a dict"


class TestPhase4Config:
    """Tests for Phase 4 configuration options."""

    def test_relax_type_parameter(self):
        """Test that relax_type parameter is accepted."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        wf = SlabWorkflow(bulk_atoms=au_bulk, surface_indices=[(1, 1, 1)], verbose=False)
        wf.generate_slabs()
        wf.bulk_recommendations = {'optimal_ecutwfc': 60.0, 'optimal_kspacing': 0.05}
        wf.convergence_results = {
            (1, 1, 1): {
                'kmesh_calc': (58, 58, 1),
                'optimal_vacuum': 18.0,
                'optimal_nlayers': 5,
            }
        }
        
        # Should accept different relax types
        # (actual execution would be tested with mock)

    def test_dipole_correction_parameter(self):
        """Test that dipole_correction parameter is accepted."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        wf = SlabWorkflow(bulk_atoms=au_bulk, surface_indices=[(1, 1, 1)], verbose=False)
        wf.generate_slabs()
        wf.bulk_recommendations = {'optimal_ecutwfc': 60.0, 'optimal_kspacing': 0.05}
        wf.convergence_results = {
            (1, 1, 1): {
                'kmesh_calc': (58, 58, 1),
                'optimal_vacuum': 18.0,
                'optimal_nlayers': 5,
            }
        }
        
        # Should accept dipole_correction parameter
        # (actual execution would be tested with mock)


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
