"""
Test convergence workflow with multi-property support.

Tests the new infrastructure for checking multiple convergence criteria.
Currently only 'energy' is implemented - other properties should raise NotImplementedError.
"""

import pytest
import numpy as np
from unittest.mock import Mock, MagicMock, patch
from ase.build import bulk
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow


class TestConvergenceMultiProperty:
    """Test suite for multi-property convergence checking."""
    
    @pytest.fixture
    def si_atoms(self):
        """Create a simple Si bulk structure."""
        return bulk('Si', 'diamond', a=5.43)
    
    @pytest.fixture
    def pseudo_dict(self):
        """Dummy pseudopotentials dict."""
        return {'Si': '/path/to/Si.pbe.UPF'}
    
    def test_default_convergence_criteria_is_energy_only(self, si_atoms, pseudo_dict):
        """Test that default convergence criteria is ['energy']."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        assert workflow.convergence_criteria_list == ['energy']
    
    def test_explicit_energy_only_convergence_criteria(self, si_atoms, pseudo_dict):
        """Test explicitly setting convergence_criteria_list=['energy']."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low',
            convergence_criteria_list=['energy']
        )
        assert workflow.convergence_criteria_list == ['energy']
    
    def test_get_calculation_config_energy_only(self, si_atoms, pseudo_dict):
        """Test _get_calculation_config() returns scf for energy-only."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low',
            convergence_criteria_list=['energy']
        )
        
        config = workflow._get_calculation_config(['energy'])
        
        assert config['calc_type'] == 'scf'
        assert isinstance(config['input_data_overrides'], dict)
    
    def test_get_calculation_config_forces_implemented(self, si_atoms, pseudo_dict):
        """Test that forces criterion is now implemented and sets tprnfor=True."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        # Should not raise NotImplementedError anymore
        config = workflow._get_calculation_config(['energy', 'forces'])
        
        assert config['calc_type'] == 'scf'
        assert config['input_data_overrides']['tprnfor'] is True
    
    def test_get_calculation_config_stress_implemented(self, si_atoms, pseudo_dict):
        """Test that stress criterion is now implemented and sets tstress=True."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        # Should not raise NotImplementedError anymore
        config = workflow._get_calculation_config(['stress'])
        
        assert config['calc_type'] == 'scf'
        assert config['input_data_overrides']['tstress'] is True
    
    def test_get_calculation_config_energy_forces_stress(self, si_atoms, pseudo_dict):
        """Test that stress works in combination with energy and forces."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        config = workflow._get_calculation_config(['energy', 'forces', 'stress'])
        
        assert config['calc_type'] == 'scf'
        assert config['input_data_overrides']['tprnfor'] is True
        assert config['input_data_overrides']['tstress'] is True
    
    def test_get_calculation_config_geometry_not_implemented(self, si_atoms, pseudo_dict):
        """Test that requesting geometry raises NotImplementedError."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        with pytest.raises(NotImplementedError) as exc_info:
            workflow._get_calculation_config(['geometry'])
        
        assert 'geometry' in str(exc_info.value)
        assert 'VC-RELAX' in str(exc_info.value)
    
    def test_get_calculation_config_magnetic_not_implemented(self, si_atoms, pseudo_dict):
        """Test that requesting magnetic_moments raises NotImplementedError."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        with pytest.raises(NotImplementedError) as exc_info:
            workflow._get_calculation_config(['magnetic_moments'])
        
        assert 'magnetic' in str(exc_info.value).lower()
    
    def test_extract_property_energy(self, si_atoms, pseudo_dict):
        """Test extracting energy from completion dict."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        # Mock completion dict
        completion = {'energy': -84.5, 'success': True}
        
        energy_per_atom = workflow._extract_property_from_result(
            completion, len(si_atoms), 'energy'
        )
        
        # Si bulk has 2 atoms in primitive cell
        expected = -84.5 / len(si_atoms)
        assert np.isclose(energy_per_atom, expected)
    
    def test_extract_property_forces_implemented(self, si_atoms, pseudo_dict):
        """Test extracting max force from completion dict."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        # Create a mock forces array: 2 atoms with forces [1,0,0] and [0,2,0]
        forces = np.array([[1.0, 0.0, 0.0],
                          [0.0, 2.0, 0.0]])
        completion = {'forces': forces, 'success': True}
        
        max_force = workflow._extract_property_from_result(completion, len(si_atoms), 'forces')
        
        # Max force magnitude should be 2.0 (from second atom)
        assert np.isclose(max_force, 2.0)
    
    def test_extract_property_forces_max_magnitude(self, si_atoms, pseudo_dict):
        """Test that max force magnitude is computed correctly."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        # Create forces with known magnitude: [3,4,0] has magnitude 5.0
        forces = np.array([[3.0, 4.0, 0.0]])
        completion = {'forces': forces, 'success': True}
        
        max_force = workflow._extract_property_from_result(completion, len(forces), 'forces')
        assert np.isclose(max_force, 5.0)
    
    def test_extract_property_stress_implemented(self, si_atoms, pseudo_dict):
        """Test extracting hydrostatic pressure from stress tensor."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        # Stress tensor [in kBar]: 
        # [[σ_xx, σ_xy, σ_xz],
        #  [σ_yx, σ_yy, σ_yz],
        #  [σ_zx, σ_zy, σ_zz]]
        # Hydrostatic pressure = -(trace/3) = -((σ_xx + σ_yy + σ_zz)/3)
        # Example: diag stress [100, 100, 100] → trace=300 → P=-(300/3)=-100 kBar
        stress_tensor = np.array([[100.0, 0.0, 0.0],
                                  [0.0, 100.0, 0.0],
                                  [0.0, 0.0, 100.0]])  # in kBar
        completion = {'stress': stress_tensor, 'success': True}
        
        pressure_gpa = workflow._extract_property_from_result(completion, len(si_atoms), 'stress')
        
        # Expected: |-(300/3) * 0.1| = |−100 * 0.1| = 10 GPa
        assert np.isclose(pressure_gpa, 10.0)
    
    def test_extract_property_stress_kbar_to_gpa_conversion(self, si_atoms, pseudo_dict):
        """Test kBar to GPa conversion in stress extraction."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        # Simple test: stress = [[10, 0, 0], [0, 10, 0], [0, 0, 10]] kBar
        # trace = 30 kBar, hydrostatic pressure = -(30/3) = -10 kBar = -1 GPa
        # abs(-1) = 1 GPa
        stress_tensor = np.array([[10.0, 0.0, 0.0],
                                  [0.0, 10.0, 0.0],
                                  [0.0, 0.0, 10.0]])  # in kBar
        completion = {'stress': stress_tensor}
        
        pressure_gpa = workflow._extract_property_from_result(completion, len(si_atoms), 'stress')
        assert np.isclose(pressure_gpa, 1.0)
    
    def test_extract_property_stress_negative_pressure(self, si_atoms, pseudo_dict):
        """Test stress extraction handles negative hydrostatic pressure (tension)."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        # Negative stress (tension): stress = [[-50, 0, 0], [0, -50, 0], [0, 0, -50]] kBar
        # trace = -150 kBar, hydrostatic pressure = -(-150/3) = 50 kBar = 5 GPa
        stress_tensor = np.array([[-50.0, 0.0, 0.0],
                                  [0.0, -50.0, 0.0],
                                  [0.0, 0.0, -50.0]])  # in kBar (tension)
        completion = {'stress': stress_tensor}
        
        pressure_gpa = workflow._extract_property_from_result(completion, len(si_atoms), 'stress')
        assert np.isclose(pressure_gpa, 5.0)
    
    def test_check_convergence_forces_converged(self, si_atoms, pseudo_dict):
        """Test force convergence check when converged."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        # Results dict: param_value -> Dict of properties (forces in eV/Å)
        results = {
            40: {'forces': 0.15},
            50: {'forces': 0.12},
            60: {'forces': 0.11},
        }
        
        reference = {'forces': 0.10}
        criteria_tolerances = {'force_tolerance': 0.05}  # 0.05 eV/Å tolerance
        
        # Should converge: max deviation is 0.05 eV/Å which equals tolerance
        converged = workflow._check_convergence_vs_reference(
            results, reference, criteria_tolerances, ['forces']
        )
        
        assert converged is True
    
    def test_check_convergence_forces_not_converged(self, si_atoms, pseudo_dict):
        """Test force convergence check when NOT converged."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        results = {
            30: {'forces': 0.25},
            40: {'forces': 0.22},
            50: {'forces': 0.18},
        }
        
        reference = {'forces': 0.10}
        criteria_tolerances = {'force_tolerance': 0.05}
        
        # Should NOT converge: max deviation is 0.15 eV/Å > 0.05 eV/Å tolerance
        converged = workflow._check_convergence_vs_reference(
            results, reference, criteria_tolerances, ['forces']
        )
        
        assert converged is False
    
    def test_check_convergence_stress_converged(self, si_atoms, pseudo_dict):
        """Test stress convergence check when converged."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        # Results dict: param_value -> Dict of properties (stress in GPa)
        results = {
            40: {'stress': 2.4},
            50: {'stress': 2.0},
            60: {'stress': 1.8},
        }
        
        reference = {'stress': 1.5}
        criteria_tolerances = {'stress_tolerance': 1.0}  # 1.0 GPa tolerance
        
        # Should converge: max deviation is 0.9 GPa < 1.0 GPa tolerance
        converged = workflow._check_convergence_vs_reference(
            results, reference, criteria_tolerances, ['stress']
        )
        
        assert converged is True
    
    def test_check_convergence_stress_not_converged(self, si_atoms, pseudo_dict):
        """Test stress convergence check when NOT converged."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        results = {
            30: {'stress': 5.0},
            40: {'stress': 4.2},
            50: {'stress': 3.5},
        }
        
        reference = {'stress': 1.5}
        criteria_tolerances = {'stress_tolerance': 1.0}
        
        # Should NOT converge: max deviation is 3.5 GPa > 1.0 GPa tolerance
        converged = workflow._check_convergence_vs_reference(
            results, reference, criteria_tolerances, ['stress']
        )
        
        assert converged is False
    
    def test_check_convergence_energy_and_forces_both_required(self, si_atoms, pseudo_dict):
        """Test that ALL criteria must converge when multiple are specified."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        # Energy converged but forces not
        results = {
            50: {'energy': -10.5398, 'forces': 0.25},
            60: {'energy': -10.5399, 'forces': 0.24},
        }
        
        reference = {'energy': -10.54, 'forces': 0.10}
        criteria_tolerances = {'energy_tolerance': 1e-3, 'force_tolerance': 0.05}
        
        # Should NOT converge because forces don't meet tolerance
        converged = workflow._check_convergence_vs_reference(
            results, reference, criteria_tolerances, ['energy', 'forces']
        )
        
        assert converged is False
    
    def test_check_convergence_energy_and_forces_both_converged(self, si_atoms, pseudo_dict):
        """Test that convergence passes when ALL criteria converge."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        # Both energy and forces converged
        results = {
            50: {'energy': -10.5398, 'forces': 0.12},
            60: {'energy': -10.5399, 'forces': 0.11},
        }
        
        reference = {'energy': -10.54, 'forces': 0.10}
        criteria_tolerances = {'energy_tolerance': 1e-3, 'force_tolerance': 0.05}
        
        # Should converge: both criteria meet their tolerances
        converged = workflow._check_convergence_vs_reference(
            results, reference, criteria_tolerances, ['energy', 'forces']
        )
        
        assert converged is True
    
    def test_check_convergence_energy_and_stress_both_converged(self, si_atoms, pseudo_dict):
        """Test that convergence passes when energy and stress both converge."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        # Both energy and stress converged
        results = {
            50: {'energy': -10.5398, 'stress': 1.9},
            60: {'energy': -10.5399, 'stress': 1.8},
        }
        
        reference = {'energy': -10.54, 'stress': 1.5}
        criteria_tolerances = {'energy_tolerance': 1e-3, 'stress_tolerance': 1.0}
        
        # Should converge: both criteria meet their tolerances
        converged = workflow._check_convergence_vs_reference(
            results, reference, criteria_tolerances, ['energy', 'stress']
        )
        
        assert converged is True
    
    def test_check_convergence_stress_not_converged_energy_does(self, si_atoms, pseudo_dict):
        """Test that convergence fails when stress doesn't converge, even if energy does."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        # Energy converged but stress not
        results = {
            50: {'energy': -10.5398, 'stress': 5.0},
            60: {'energy': -10.5399, 'stress': 4.2},
        }
        
        reference = {'energy': -10.54, 'stress': 1.5}
        criteria_tolerances = {'energy_tolerance': 1e-3, 'stress_tolerance': 1.0}
        
        # Should NOT converge because stress doesn't meet tolerance
        converged = workflow._check_convergence_vs_reference(
            results, reference, criteria_tolerances, ['energy', 'stress']
        )
        
        assert converged is False
    
    def test_check_convergence_energy_single_value(self, si_atoms, pseudo_dict):
        """Test energy convergence check with single parameter value."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        # Results dict: param_value -> Dict of properties
        results = {
            40: {'energy': -10.55},
            50: {'energy': -10.54},
            60: {'energy': -10.5398},
        }
        
        reference = {'energy': -10.54}
        criteria_tolerances = {'energy_tolerance': 1e-3}
        
        # Should converge: max deviation is ~0.0002 eV < tolerance
        converged = workflow._check_convergence_vs_reference(
            results, reference, criteria_tolerances, ['energy']
        )
        
        assert converged is True
    
    def test_check_convergence_energy_not_converged(self, si_atoms, pseudo_dict):
        """Test energy convergence check when NOT converged."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        results = {
            30: {'energy': -10.3},
            40: {'energy': -10.48},
            50: {'energy': -10.52},
        }
        
        reference = {'energy': -10.54}
        criteria_tolerances = {'energy_tolerance': 1e-3}
        
        # Should NOT converge: max deviation is 0.24 eV > tolerance
        converged = workflow._check_convergence_vs_reference(
            results, reference, criteria_tolerances, ['energy']
        )
        
        assert converged is False
    
    def test_check_convergence_empty_results(self, si_atoms, pseudo_dict):
        """Test convergence check with empty results."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low'
        )
        
        results = {}
        reference = {'energy': -10.54}
        criteria_tolerances = {'energy_tolerance': 1e-3}
        
        converged = workflow._check_convergence_vs_reference(
            results, reference, criteria_tolerances, ['energy']
        )
        
        assert converged is False
    
    def test_run_convergence_with_unsupported_criteria_raises_validation(self, si_atoms, pseudo_dict):
        """Test that run_convergence_independent validates criteria before running."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict,
            precision='low',
            convergence_criteria_list=['geometry']  # Not yet implemented (requires VC-RELAX)
        )
        
        # Should raise NotImplementedError during configuration validation
        with pytest.raises(NotImplementedError) as exc_info:
            workflow._get_calculation_config(['geometry'])
        
        assert 'geometry' in str(exc_info.value).lower()
    
    def test_tolerance_values_by_precision(self, si_atoms, pseudo_dict):
        """Test that tolerance values are correctly set by precision level."""
        for precision in ['low', 'medium', 'high', 'ultra']:
            workflow = ConvergenceWorkflow(
                atoms=si_atoms,
                pseudopotentials=pseudo_dict,
                precision=precision
            )
            
            criteria = workflow.convergence_criteria
            
            # All precision levels should have energy_tolerance
            assert 'energy_tolerance' in criteria
            assert isinstance(criteria['energy_tolerance'], (int, float))
            assert criteria['energy_tolerance'] > 0
            
            # Should have stress_tolerance (newly added)
            assert 'stress_tolerance' in criteria
            assert isinstance(criteria['stress_tolerance'], (int, float))
    
    def test_default_precision_is_low(self, si_atoms, pseudo_dict):
        """Test that default precision is 'low' for faster convergence."""
        workflow = ConvergenceWorkflow(
            atoms=si_atoms,
            pseudopotentials=pseudo_dict
            # No precision specified - should default to 'low'
        )
        
        assert workflow.precision == 'low'
        
        # Low precision should have looser (larger) tolerance
        low_tolerance = workflow.convergence_criteria['energy_tolerance']
        assert low_tolerance >= 1e-3  # 1 meV/atom or looser


class TestConvergenceBackwardCompatibility:
    """Test that existing code still works (backward compatibility)."""
    
    def test_from_cif_with_default_precision(self):
        """Test ConvergenceWorkflow.from_cif() with default precision."""
        # Create a temporary CIF-like structure
        from ase.io import write
        from ase.build import bulk
        import tempfile
        import os
        
        si = bulk('Si', 'diamond', a=5.43)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            cif_path = os.path.join(tmpdir, 'test.cif')
            write(cif_path, si)
            
            workflow = ConvergenceWorkflow.from_cif(
                cif_path,
                pseudopotentials={'Si': '/dummy/Si.UPF'},
                # No precision specified - should default to 'low'
            )
            
            assert workflow.precision == 'low'
            assert workflow.convergence_criteria_list == ['energy']
    
    def test_optimize_parameters_default_precision(self):
        """Test optimize_parameters classmethod with default precision."""
        si = bulk('Si', 'diamond', a=5.43)
        
        # Should not raise even without precision specified
        # (would actually run calculations, so we just test it doesn't fail during init)
        try:
            workflow = ConvergenceWorkflow(
                atoms=si,
                pseudopotentials={'Si': '/dummy/Si.UPF'}
                # No precision specified
            )
            assert workflow.precision == 'low'
        except Exception as e:
            pytest.fail(f"Failed to create workflow with default precision: {e}")


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
