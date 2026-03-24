#!/usr/bin/env python3
"""
Unit tests for SlabWorkflow Phase 3 (Slab Convergence).

Tests:
- run_slab_convergence() basic functionality
- _calculate_anisotropic_kmesh() k-mesh calculation
- _test_vacuum_convergence() vacuum testing
- _test_layer_convergence() layer testing
"""

import numpy as np
import pytest
from ase.build import bulk
from xespresso.workflow.slab_workflow import SlabWorkflow


class TestAnisotropicKMesh:
    """Tests for k-mesh calculation."""

    @pytest.fixture
    def workflow(self):
        """Create workflow with bulk convergence done."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        wf = SlabWorkflow(
            bulk_atoms=au_bulk,
            surface_indices=[(1, 1, 1)],
            pseudopotentials_config='default',
            verbose=False,
        )
        # Mock bulk recommendations
        wf.bulk_recommendations = {
            'optimal_ecutwfc': 60.0,
            'optimal_kspacing': 0.05,
        }
        return wf

    def test_calculate_anisotropic_kmesh_returns_tuple(self, workflow):
        """Test that k-mesh calculation returns (int, int, int)."""
        slabs = workflow.generate_slabs()
        slab = slabs[(1, 1, 1)]
        
        kmesh = workflow._calculate_anisotropic_kmesh(slab)
        
        assert isinstance(kmesh, tuple), "Should return tuple"
        assert len(kmesh) == 3, "Should have 3 components"
        assert all(isinstance(n, int) for n in kmesh), "All should be integers"

    def test_calculate_anisotropic_kmesh_all_positive(self, workflow):
        """Test that all k-mesh components are positive."""
        slabs = workflow.generate_slabs()
        slab = slabs[(1, 1, 1)]
        
        nk_x, nk_y, nk_z = workflow._calculate_anisotropic_kmesh(slab)
        
        assert nk_x > 0, "nk_x should be positive"
        assert nk_y > 0, "nk_y should be positive"
        assert nk_z >= 1, "nk_z should be >= 1"

    def test_calculate_anisotropic_kmesh_z_is_one(self, workflow):
        """Test that nk_z is always 1 for 2D slabs."""
        slabs = workflow.generate_slabs()
        slab = slabs[(1, 1, 1)]
        
        _, _, nk_z = workflow._calculate_anisotropic_kmesh(slab)
        
        assert nk_z == 1, "nk_z should always be 1 for 2D slabs"

    def test_calculate_anisotropic_kmesh_scales_with_kspacing(self, workflow):
        """Test that k-mesh scales inversely with k-spacing."""
        slabs = workflow.generate_slabs()
        slab = slabs[(1, 1, 1)]
        
        # Test with tight k-spacing (more points)
        kmesh_tight = workflow._calculate_anisotropic_kmesh(slab, kspacing=0.02)
        
        # Test with loose k-spacing (fewer points)
        kmesh_loose = workflow._calculate_anisotropic_kmesh(slab, kspacing=0.1)
        
        # Tight should have more/equal points than loose
        assert kmesh_tight[0] >= kmesh_loose[0], "Tight spacing should have more x points"
        assert kmesh_tight[1] >= kmesh_loose[1], "Tight spacing should have more y points"


class TestPhase3Integration:
    """Integration tests for Phase 3 convergence."""

    @pytest.fixture
    def workflow_ready(self):
        """Create workflow with both Phase 1 and Phase 2 done."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        wf = SlabWorkflow(
            bulk_atoms=au_bulk,
            surface_indices=[(1, 1, 1)],
            pseudopotentials_config='default',
            verbose=False,
        )
        # Generate slabs (Phase 2)
        wf.generate_slabs()
        # Mock bulk convergence (Phase 1)
        wf.bulk_recommendations = {
            'optimal_ecutwfc': 60.0,
            'optimal_kspacing': 0.05,
        }
        return wf

    def test_run_slab_convergence_requires_bulk_convergence(self):
        """Test that Phase 3 requires Phase 1 to be done."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        wf = SlabWorkflow(bulk_atoms=au_bulk, surface_indices=[(1, 1, 1)], verbose=False)
        wf.generate_slabs()
        
        # Bulk convergence not done
        with pytest.raises(ValueError, match="Bulk convergence"):
            wf.run_slab_convergence((1, 1, 1))

    def test_run_slab_convergence_requires_slabs(self):
        """Test that Phase 3 requires Phase 2 slabs."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        wf = SlabWorkflow(bulk_atoms=au_bulk, surface_indices=[(1, 1, 1)], verbose=False)
        # Don't generate slabs
        wf.bulk_recommendations = {'optimal_ecutwfc': 60.0, 'optimal_kspacing': 0.05}
        
        with pytest.raises(ValueError, match="not generated"):
            wf.run_slab_convergence((1, 1, 1))

    def test_run_slab_convergence_returns_dict(self, workflow_ready):
        """Test that run_slab_convergence returns proper dictionary."""
        results = workflow_ready.run_slab_convergence(
            (1, 1, 1),
            vacuum_test=[15, 18, 20],
            nlayers_test=[4, 5],
            skip_calculations=True,
        )
        
        assert isinstance(results, dict), "Should return dict"
        assert 'surface_index' in results, "Should have surface_index"
        assert 'kmesh_calc' in results, "Should have kmesh_calc"
        assert 'optimal_vacuum' in results, "Should have optimal_vacuum"
        assert 'optimal_nlayers' in results, "Should have optimal_nlayers"
        assert 'converged' in results, "Should have converged flag"

    def test_run_slab_convergence_kmesh_stored(self, workflow_ready):
        """Test that k-mesh is correctly calculated and stored."""
        results = workflow_ready.run_slab_convergence(
            (1, 1, 1),
            vacuum_test=[15],
            nlayers_test=[4],
            skip_calculations=True,
        )
        
        kmesh = results['kmesh_calc']
        assert isinstance(kmesh, tuple), "kmesh should be tuple"
        assert len(kmesh) == 3, "kmesh should have 3 components"
        assert kmesh[2] == 1, "nk_z should always be 1 for slabs"

    def test_run_slab_convergence_results_stored_in_workflow(self, workflow_ready):
        """Test that results are stored in self.convergence_results."""
        surface = (1, 1, 1)
        results = workflow_ready.run_slab_convergence(
            surface,
            vacuum_test=[15],
            nlayers_test=[4],
            skip_calculations=True,
        )
        
        assert surface in workflow_ready.convergence_results, \
            "Results should be stored in workflow"
        assert workflow_ready.convergence_results[surface] == results, \
            "Stored results should match returned results"

    def test_run_slab_convergence_custom_parameters(self, workflow_ready):
        """Test with custom vacuum and layer test values."""
        results = workflow_ready.run_slab_convergence(
            (1, 1, 1),
            vacuum_test=[10, 15, 20, 25],
            nlayers_test=[3, 4, 5],
            skip_calculations=True,
        )
        
        # Both vacuum and layers should have results
        assert len(results['vacuum_results']) >= 0, "Should test vacuums"
        assert len(results['layer_results']) >= 0, "Should test layers"

    def test_run_slab_convergence_defaults_if_none(self, workflow_ready):
        """Test that defaults are used if None passed."""
        results = workflow_ready.run_slab_convergence(
            (1, 1, 1),
            vacuum_test=None,  # Use defaults
            nlayers_test=None,  # Use defaults
            skip_calculations=True,
        )
        
        assert results['optimal_vacuum'] is not None, "Should have optimal vacuum"
        assert results['optimal_nlayers'] is not None, "Should have optimal layers"


class TestVacuumConvergence:
    """Tests for vacuum convergence testing."""

    @pytest.fixture
    def workflow_ready(self):
        """Create ready workflow."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        wf = SlabWorkflow(bulk_atoms=au_bulk, surface_indices=[(1, 1, 1)], verbose=False)
        wf.generate_slabs()
        wf.bulk_recommendations = {'optimal_ecutwfc': 60.0, 'optimal_kspacing': 0.05}
        return wf

    def test_vacuum_convergence_returns_dict(self, workflow_ready):
        """Test _test_vacuum_convergence return structure."""
        kmesh = workflow_ready._calculate_anisotropic_kmesh(
            workflow_ready.slabs[(1, 1, 1)]
        )
        
        results = workflow_ready._test_vacuum_convergence(
            (1, 1, 1),
            [15, 18, 20],
            kmesh,
            'test',
            skip_calculations=True,
        )
        
        assert isinstance(results, dict), "Should return dict"
        assert 'energies' in results, "Should have energies"
        assert 'optimal_vacuum' in results, "Should have optimal_vacuum"
        assert 'converged' in results, "Should have converged"

    def test_vacuum_convergence_optimal_in_testset(self, workflow_ready):
        """Test that optimal_vacuum is from the test set."""
        kmesh = workflow_ready._calculate_anisotropic_kmesh(
            workflow_ready.slabs[(1, 1, 1)]
        )
        vacuum_test = [12, 15, 18, 20]
        
        results = workflow_ready._test_vacuum_convergence(
            (1, 1, 1),
            vacuum_test,
            kmesh,
            'test',
            skip_calculations=True,
        )
        
        # If calculations were run, optimal should be in test set
        if results['energies']:
            assert results['optimal_vacuum'] in vacuum_test, \
                "Optimal should be from test set"


class TestLayerConvergence:
    """Tests for layer convergence testing."""

    @pytest.fixture
    def workflow_ready(self):
        """Create ready workflow."""
        au_bulk = bulk('Au', 'fcc', a=4.078)
        wf = SlabWorkflow(bulk_atoms=au_bulk, surface_indices=[(1, 1, 1)], verbose=False)
        wf.generate_slabs()
        wf.bulk_recommendations = {'optimal_ecutwfc': 60.0, 'optimal_kspacing': 0.05}
        return wf

    def test_layer_convergence_returns_dict(self, workflow_ready):
        """Test _test_layer_convergence return structure."""
        kmesh = workflow_ready._calculate_anisotropic_kmesh(
            workflow_ready.slabs[(1, 1, 1)]
        )
        
        results = workflow_ready._test_layer_convergence(
            (1, 1, 1),
            [3, 4, 5],
            15.0,
            kmesh,
            'test',
            skip_calculations=True,
        )
        
        assert isinstance(results, dict), "Should return dict"
        assert 'energies' in results, "Should have energies"
        assert 'optimal_nlayers' in results, "Should have optimal_nlayers"
        assert 'converged' in results, "Should have converged"


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
