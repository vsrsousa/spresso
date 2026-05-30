"""
Unit tests for EOS workflow module.

Tests volume scaling, EOS fitting, and utility functions.
"""

import pytest
import numpy as np
from pathlib import Path
from ase.build import bulk
import tempfile
import os

from xespresso.workflow.eos_workflow import (
    EOSWorkflow,
    birch_murnaghan_eos,
    fit_birch_murnaghan,
    DEFAULT_VOLUME_MIN_FACTOR,
    DEFAULT_VOLUME_MAX_FACTOR,
    DEFAULT_EOS_POINTS,
)


# ═════════════════════════════════════════════════════════════════════════════
# FIXTURES
# ═════════════════════════════════════════════════════════════════════════════

@pytest.fixture
def atoms_si():
    """Diamond structure Si (small unit cell)."""
    return bulk('Si', 'diamond', a=5.43)


@pytest.fixture
def atoms_fe():
    """BCC structure Fe."""
    return bulk('Fe', 'bcc', a=2.87)


@pytest.fixture
def atoms_au():
    """FCC structure Au."""
    return bulk('Au', 'fcc', a=4.08)


@pytest.fixture
def pseudo_dict_si():
    """Pseudopotential dict for Si."""
    return {'Si': 'Si.pbe.UPF'}


@pytest.fixture
def pseudo_dict_fe():
    """Pseudopotential dict for Fe."""
    return {'Fe': 'Fe.pbe-spn.UPF'}


# ═════════════════════════════════════════════════════════════════════════════
# TEST: Birch-Murnaghan EOS Function
# ═════════════════════════════════════════════════════════════════════════════

class TestBirchMurnaghanEOS:
    """Tests for birch_murnaghan_eos() function."""
    
    def test_eos_at_equilibrium_volume(self):
        """EOS value should equal E₀ at V = V₀."""
        E0, V0, B0, BP = -10.5, 20.0, 150.0, 4.0
        
        E_at_v0 = birch_murnaghan_eos(np.array([V0]), E0, V0, B0, BP)
        
        assert np.isclose(E_at_v0[0], E0, rtol=1e-6)
    
    def test_eos_symmetry(self):
        """EOS should be approximately symmetric around V₀."""
        E0, V0, B0, BP = -10.5, 20.0, 150.0, 4.0
        
        dV = 0.5
        E_below = birch_murnaghan_eos(np.array([V0 - dV]), E0, V0, B0, BP)
        E_above = birch_murnaghan_eos(np.array([V0 + dV]), E0, V0, B0, BP)
        
        # Should be approximately symmetric (within 5%)
        assert np.isclose(E_below[0], E_above[0], rtol=0.05)
    
    def test_eos_with_negative_b0(self):
        """EOS should handle negative B₀ gracefully (return large value)."""
        E0, V0, B0, BP = -10.5, 20.0, -150.0, 4.0
        
        E = birch_murnaghan_eos(np.array([V0]), E0, V0, B0, BP)
        
        assert E[0] > 1e9  # Large unphysical value
    
    def test_eos_vector_input(self):
        """EOS should handle array of volumes."""
        E0, V0, B0, BP = -10.5, 20.0, 150.0, 4.0
        volumes = np.array([19.0, 19.5, 20.0, 20.5, 21.0])
        
        energies = birch_murnaghan_eos(volumes, E0, V0, B0, BP)
        
        assert len(energies) == len(volumes)
        assert energies.shape == volumes.shape
        # Minimum should be at V0
        assert np.isclose(volumes[np.argmin(energies)], V0, rtol=0.01)


# ═════════════════════════════════════════════════════════════════════════════
# TEST: EOS Fitting
# ═════════════════════════════════════════════════════════════════════════════

class TestBirchMurnaghanFitting:
    """Tests for fit_birch_murnaghan() function."""
    
    def test_fit_synthetic_data(self):
        """Fit should recover parameters from synthetic data."""
        # Generate synthetic data
        E0_true, V0_true, B0_true, BP_true = -10.5, 20.0, 150.0, 4.0
        volumes = np.array([19.0, 19.5, 20.0, 20.5, 21.0])
        energies = birch_murnaghan_eos(volumes, E0_true, V0_true, B0_true, BP_true)
        
        # Fit
        params = fit_birch_murnaghan(volumes, energies)
        
        # Check parameters (should be close to original)
        assert np.isclose(params['v0'], V0_true, rtol=0.02)
        assert np.isclose(params['e0'], E0_true, rtol=0.02)
        assert np.isclose(params['b0'], B0_true, rtol=0.10)
        assert params['r_squared'] > 0.99
        assert params['converged']
    
    def test_fit_with_noise(self):
        """Fit should work with noisy data."""
        E0_true, V0_true, B0_true, BP_true = -10.5, 20.0, 150.0, 4.0
        volumes = np.array([19.0, 19.5, 20.0, 20.5, 21.0])
        energies = birch_murnaghan_eos(volumes, E0_true, V0_true, B0_true, BP_true)
        
        # Add 1% noise
        noise = np.random.normal(0, 0.01 * np.max(np.abs(energies)), len(energies))
        energies_noisy = energies + noise
        
        # Fit
        params = fit_birch_murnaghan(volumes, energies_noisy)
        
        # Should still be reasonable
        assert params['r_squared'] > 0.95
        assert np.isclose(params['v0'], V0_true, rtol=0.10)
    
    def test_fit_insufficient_points(self):
        """Fit should raise error with < 3 points."""
        volumes = np.array([19.0, 20.0])
        energies = np.array([-10.0, -10.5])
        
        with pytest.raises(ValueError, match="Need at least"):
            fit_birch_murnaghan(volumes, energies)


# ═════════════════════════════════════════════════════════════════════════════
# TEST: EOSWorkflow - Volume Scaling
# ═════════════════════════════════════════════════════════════════════════════

class TestEOSWorkflowVolumeScaling:
    """Tests for volume scaling methods."""
    
    def test_scale_volume_uniformly(self, atoms_si):
        """Volume should scale exactly as specified."""
        eos = EOSWorkflow(atoms_si, pseudo_dict_si, protocol='fast')
        
        scale_factor = 1.05
        atoms_scaled = eos.scale_volume_uniformly(atoms_si, scale_factor)
        
        v_original = atoms_si.get_volume()
        v_scaled = atoms_scaled.get_volume()
        ratio = v_scaled / v_original
        
        assert np.isclose(ratio, scale_factor, rtol=1e-10)
    
    def test_scale_volume_multiple_factors(self, atoms_fe):
        """Multiple scaling factors should work correctly."""
        eos = EOSWorkflow(atoms_fe, pseudo_dict_fe, protocol='fast')
        
        factors = [0.95, 1.0, 1.05]
        v_original = atoms_fe.get_volume()
        
        for factor in factors:
            atoms_scaled = eos.scale_volume_uniformly(atoms_fe, factor)
            ratio = atoms_scaled.get_volume() / v_original
            assert np.isclose(ratio, factor, rtol=1e-10)
    
    def test_scale_volume_invalid_factor(self, atoms_si):
        """Should raise error for invalid scale factor."""
        eos = EOSWorkflow(atoms_si, pseudo_dict_si, protocol='fast')
        
        with pytest.raises(ValueError):
            eos.scale_volume_uniformly(atoms_si, 0.0)
        
        with pytest.raises(ValueError):
            eos.scale_volume_uniformly(atoms_si, -1.0)
    
    def test_generate_volume_range(self, atoms_si):
        """Volume range generation should be correct."""
        eos = EOSWorkflow(atoms_si, pseudo_dict_si, protocol='fast')
        
        factors = eos.generate_volume_range((0.95, 1.05), n_points=7)
        
        assert len(factors) == 7
        assert np.isclose(factors[0], 0.95)
        assert np.isclose(factors[-1], 1.05)
        assert np.isclose(factors[3], 1.0)  # Middle should be 1.0
        # Should be monotonically increasing
        assert np.all(np.diff(factors) > 0)
    
    def test_generate_volume_range_invalid_points(self, atoms_si):
        """Should raise error with too few points."""
        eos = EOSWorkflow(atoms_si, pseudo_dict_si, protocol='fast')
        
        with pytest.raises(ValueError, match="Need at least"):
            eos.generate_volume_range((0.95, 1.05), n_points=2)
    
    def test_create_scaled_structures(self, atoms_au):
        """Should create correct number of structures."""
        eos = EOSWorkflow(atoms_au, {'Au': 'Au.pbe.UPF'}, protocol='fast')
        
        factors = np.array([0.95, 1.0, 1.05])
        structures = eos.create_scaled_structures(factors)
        
        assert len(structures) == 3
        assert all(f in structures for f in factors)
        
        # Check volumes scale correctly
        v_original = atoms_au.get_volume()
        for factor, atoms in structures.items():
            v_scaled = atoms.get_volume()
            assert np.isclose(v_scaled / v_original, factor, rtol=1e-10)


# ═════════════════════════════════════════════════════════════════════════════
# TEST: EOSWorkflow - Initialization
# ═════════════════════════════════════════════════════════════════════════════

class TestEOSWorkflowInitialization:
    """Tests for EOSWorkflow initialization."""
    
    def test_init_with_atoms_object(self, atoms_si, pseudo_dict_si):
        """Should initialize with Atoms object."""
        eos = EOSWorkflow(atoms_si, pseudo_dict_si)
        
        assert eos.atoms is not None
        assert eos.protocol == 'moderate'
        assert eos.results_df is None
        assert eos.eos_params is None
    
    def test_init_with_custom_protocol(self, atoms_fe, pseudo_dict_fe):
        """Should accept custom protocol."""
        eos = EOSWorkflow(atoms_fe, pseudo_dict_fe, protocol='accurate')
        
        assert eos.protocol == 'accurate'


# ═════════════════════════════════════════════════════════════════════════════
# TEST: EOSWorkflow - Properties
# ═════════════════════════════════════════════════════════════════════════════

class TestEOSWorkflowProperties:
    """Tests for property extraction."""
    
    def test_predict_energy_without_fit(self, atoms_si, pseudo_dict_si):
        """Should raise error if EOS not fitted."""
        eos = EOSWorkflow(atoms_si, pseudo_dict_si)
        
        with pytest.raises(ValueError, match="not fitted"):
            eos.predict_energy(20.0)
    
    def test_get_eos_properties_without_fit(self, atoms_fe, pseudo_dict_fe):
        """Should raise error if EOS not fitted."""
        eos = EOSWorkflow(atoms_fe, pseudo_dict_fe)
        
        with pytest.raises(ValueError, match="not fitted"):
            eos.get_eos_properties()


# ═════════════════════════════════════════════════════════════════════════════
# TEST: Backward Compatibility
# ═════════════════════════════════════════════════════════════════════════════

def test_eos_workflow_import():
    """EOSWorkflow should be importable from xespresso.workflow."""
    from xespresso.workflow import EOSWorkflow as EOSWorkflowImported
    assert EOSWorkflowImported is EOSWorkflow


def test_quick_eos_import():
    """quick_eos should be importable from xespresso.workflow."""
    from xespresso.workflow import quick_eos as quick_eos_imported
    from xespresso.workflow.eos_workflow import quick_eos as quick_eos_direct
    assert quick_eos_imported is quick_eos_direct


# ═════════════════════════════════════════════════════════════════════════════
# TEST: Error Handling
# ═════════════════════════════════════════════════════════════════════════════

class TestErrorHandling:
    """Tests for error handling and edge cases."""
    
    def test_eos_curve_without_data(self, atoms_si, pseudo_dict_si):
        """plot_eos_curve should raise error without data."""
        eos = EOSWorkflow(atoms_si, pseudo_dict_si)
        
        with pytest.raises(ValueError, match="No EOS data"):
            eos.plot_eos_curve()
    
    def test_summary_without_data(self, atoms_fe, pseudo_dict_fe):
        """summary() should handle incomplete workflow."""
        eos = EOSWorkflow(atoms_fe, pseudo_dict_fe)
        summary = eos.summary()
        
        assert "not completed" in summary


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
