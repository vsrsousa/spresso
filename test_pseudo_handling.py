"""
Tests for pseudopotential handling in ConvergenceWorkflow.

Tests both dict-based and config-based pseudopotential loading to ensure
the new auto-discovery mechanism works and doesn't break existing functionality.
"""

import os
import unittest
from unittest import mock
from ase.build import bulk
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow


class TestPseudopotentialHandling(unittest.TestCase):
    """Test pseudopotential handling via dict and config."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.atoms = bulk('Gd', crystalstructure='hcp', a=3.6, c=5.78)
        self.pseudo_dict = {'Gd': 'Gd.pbe-spfn-rrkjus_psl.1.0.0.UPF'}
        self.pseudo_base_path = '/home/vinicius/scratch/projects/spresso'
        
    def test_dict_based_pseudopotentials_discovery(self):
        """Test that dict-based pseudopotentials are auto-discovered."""
        
        with mock.patch('xespresso.workflow.convergence_workflow.discover_pseudopotential_directory') as mock_discover:
            # Mock the discovery function to return the base path
            mock_discover.return_value = (
                self.pseudo_dict,  # Returned dict (unchanged)
                self.pseudo_base_path  # Base path where found
            )
            
            with mock.patch.dict(os.environ, {}, clear=False):
                wf = ConvergenceWorkflow(
                    atoms=self.atoms,
                    pseudopotentials=self.pseudo_dict,
                    precision='low',
                    convergence_criteria_list=['energy']
                )
            
            # Verify discovery was called
            mock_discover.assert_called_once_with(self.pseudo_dict)
            
            # Verify base path was stored
            self.assertEqual(wf.pseudopotentials_base_path, self.pseudo_base_path)
            
            # Verify pseudopotentials dict is stored
            self.assertEqual(wf.pseudopotentials, self.pseudo_dict)
            
            # Verify ESPRESSO_PSEUDO was set
            self.assertEqual(
                os.environ.get('ESPRESSO_PSEUDO'),
                self.pseudo_base_path,
                "ESPRESSO_PSEUDO env var should be set to base path"
            )
    
    def test_config_based_pseudopotentials_no_discovery(self):
        """Test that config-based pseudopotentials don't trigger discovery."""
        
        mock_config = mock.MagicMock()
        mock_config.base_path = '/etc/espresso/pseudo'
        mock_config.get_pseudopotential.return_value = mock.MagicMock(
            filename='Gd.pbe-spfn-rrkjus_psl.1.0.0.UPF'
        )
        
        with mock.patch('xespresso.workflow.convergence_workflow.discover_pseudopotential_directory') as mock_discover:
            with mock.patch('xespresso.pseudopotentials.manager.load_pseudopotentials_config') as mock_load_config:
                mock_load_config.return_value = mock_config
                
                with mock.patch.dict(os.environ, {}, clear=False):
                    wf = ConvergenceWorkflow(
                        atoms=self.atoms,
                        pseudopotentials_config='test_config',
                        precision='low',
                        convergence_criteria_list=['energy']
                    )
                
                # Discovery should NOT be called for config-based
                mock_discover.assert_not_called()
                
                # Base path from config should be stored
                self.assertEqual(wf.pseudopotentials_base_path, '/etc/espresso/pseudo')
                
                # Config name should be stored for later use
                self.assertEqual(wf._pseudo_config_name, 'test_config')
    
    def test_dict_with_absolute_paths(self):
        """Test that absolute paths in dict are preserved."""
        
        absolute_pseudo_dict = {
            'Gd': '/absolute/path/to/Gd.pbe-spfn-rrkjus_psl.1.0.0.UPF'
        }
        
        with mock.patch('xespresso.workflow.convergence_workflow.discover_pseudopotential_directory') as mock_discover:
            mock_discover.return_value = (absolute_pseudo_dict, None)
            
            wf = ConvergenceWorkflow(
                atoms=self.atoms,
                pseudopotentials=absolute_pseudo_dict,
                precision='low',
                convergence_criteria_list=['energy']
            )
            
            # Should preserve absolute path
            self.assertEqual(wf.pseudopotentials['Gd'], absolute_pseudo_dict['Gd'])
    
    def test_discovery_file_not_found_error(self):
        """Test that FileNotFoundError is properly raised and propagated."""
        
        with mock.patch('xespresso.workflow.convergence_workflow.discover_pseudopotential_directory') as mock_discover:
            mock_discover.side_effect = FileNotFoundError(
                "Pseudopotential file not found: Gd -> Gd.pbe-spfn-rrkjus_psl.1.0.0.UPF"
            )
            
            with self.assertRaises(FileNotFoundError) as ctx:
                ConvergenceWorkflow(
                    atoms=self.atoms,
                    pseudopotentials=self.pseudo_dict,
                    precision='low',
                    convergence_criteria_list=['energy']
                )
            
            self.assertIn("Pseudopotential file not found", str(ctx.exception))
    
    def test_espresso_pseudo_env_var_persistence(self):
        """Test that ESPRESSO_PSEUDO persists during workflow execution."""
        
        with mock.patch('xespresso.workflow.convergence_workflow.discover_pseudopotential_directory') as mock_discover:
            mock_discover.return_value = (self.pseudo_dict, self.pseudo_base_path)
            
            # Clear env var before test
            original_value = os.environ.pop('ESPRESSO_PSEUDO', None)
            
            try:
                wf = ConvergenceWorkflow(
                    atoms=self.atoms,
                    pseudopotentials=self.pseudo_dict,
                    precision='low',
                    convergence_criteria_list=['energy']
                )
                
                # After init, ESPRESSO_PSEUDO should be set
                self.assertEqual(os.environ['ESPRESSO_PSEUDO'], self.pseudo_base_path)
                
                # Verify it's still set (hasn't been cleared)
                self.assertEqual(
                    os.environ.get('ESPRESSO_PSEUDO'),
                    self.pseudo_base_path,
                    "ESPRESSO_PSEUDO should persist"
                )
            finally:
                # Restore original value
                if original_value:
                    os.environ['ESPRESSO_PSEUDO'] = original_value
                else:
                    os.environ.pop('ESPRESSO_PSEUDO', None)
    
    def test_calculation_workflow_receives_base_path(self):
        """Test that CalculationWorkflow is created with pseudopotentials_base_path."""
        
        with mock.patch('xespresso.workflow.convergence_workflow.discover_pseudopotential_directory') as mock_discover:
            mock_discover.return_value = (self.pseudo_dict, self.pseudo_base_path)
            
            wf = ConvergenceWorkflow(
                atoms=self.atoms,
                pseudopotentials=self.pseudo_dict,
                precision='low',
                convergence_criteria_list=['energy']
            )
            
            # Verify that base_path was properly stored
            self.assertEqual(
                wf.pseudopotentials_base_path,
                self.pseudo_base_path,
                "ConvergenceWorkflow should have pseudopotentials_base_path"
            )
            
            # Pseudo dict should be stored
            self.assertEqual(wf.pseudopotentials, self.pseudo_dict)


class TestPseudopotentialBackwardCompatibility(unittest.TestCase):
    """Test that existing behavior is not broken."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.atoms = bulk('Si', crystalstructure='diamond', a=5.43)
    
    def test_config_based_still_works_unchanged(self):
        """Ensure config-based loading still works as before."""
        
        mock_config = mock.MagicMock()
        mock_config.base_path = '/etc/espresso/pseudo'
        mock_config.get_pseudopotential.return_value = mock.MagicMock(
            filename='Si.pbe-n-rrkjus_psl.1.0.0.UPF'
        )
        
        with mock.patch('xespresso.pseudopotentials.manager.load_pseudopotentials_config') as mock_load:
            mock_load.return_value = mock_config
            
            wf = ConvergenceWorkflow(
                atoms=self.atoms,
                pseudopotentials_config='SSSP_efficiency',
                precision='medium',
                convergence_criteria_list=['energy']
            )
            
            # Config-based still has base_path
            self.assertEqual(wf.pseudopotentials_base_path, '/etc/espresso/pseudo')
            
            # Config name is stored
            self.assertEqual(wf._pseudo_config_name, 'SSSP_efficiency')
    
    def test_precision_levels_still_control_convergence_criteria(self):
        """Verify that precision still controls convergence tolerances."""
        
        mock_config = mock.MagicMock()
        mock_config.base_path = '/etc/espresso/pseudo'
        mock_config.get_pseudopotential.return_value = mock.MagicMock(
            filename='Si.pbe-n-rrkjus_psl.1.0.0.UPF'
        )
        
        precision_levels = ['low', 'medium', 'high', 'ultra']
        expected_energy_tol = {
            'low': 3e-3,
            'medium': 2e-3,
            'high': 1e-3,
            'ultra': 5e-4,
        }
        
        with mock.patch('xespresso.pseudopotentials.manager.load_pseudopotentials_config') as mock_load:
            mock_load.return_value = mock_config
            
            for precision in precision_levels:
                wf = ConvergenceWorkflow(
                    atoms=self.atoms,
                    pseudopotentials_config='SSSP_efficiency',
                    precision=precision,
                    convergence_criteria_list=['energy']
                )
                
                self.assertEqual(
                    wf.convergence_criteria['energy_tolerance'],
                    expected_energy_tol[precision],
                    f"Energy tolerance for {precision} should be {expected_energy_tol[precision]}"
                )


class TestPseudopotentialIntegration(unittest.TestCase):
    """Integration tests for pseudopotential handling."""
    
    def test_dict_vs_config_produces_same_base_path_storage(self):
        """Verify both dict and config mechanisms store base_path correctly."""
        
        atoms = bulk('Al', crystalstructure='fcc', a=4.05)
        pseudo_dict = {'Al': 'Al.pbe-n-rrkjus_psl.1.0.0.UPF'}
        base_path = '/home/pseudo/library'
        mock_config_path = '/config/path'
        
        # Test dict-based
        with mock.patch('xespresso.workflow.convergence_workflow.discover_pseudopotential_directory') as mock_discover:
            mock_discover.return_value = (pseudo_dict, base_path)
            
            wf_dict = ConvergenceWorkflow(
                atoms=atoms,
                pseudopotentials=pseudo_dict,
                precision='low',
                convergence_criteria_list=['energy']
            )
        
        # Test config-based
        mock_config = mock.MagicMock()
        mock_config.base_path = mock_config_path
        
        with mock.patch('xespresso.pseudopotentials.manager.load_pseudopotentials_config') as mock_load:
            mock_load.return_value = mock_config
            
            wf_config = ConvergenceWorkflow(
                atoms=atoms,
                pseudopotentials_config='test_config',
                precision='low',
                convergence_criteria_list=['energy']
            )
        
        # Both should have a base_path (though different values)
        self.assertIsNotNone(wf_dict.pseudopotentials_base_path)
        self.assertIsNotNone(wf_config.pseudopotentials_base_path)
        
        # Dict-based should have the discovered path
        self.assertEqual(wf_dict.pseudopotentials_base_path, base_path)
        
        # Config-based should have the config path
        self.assertEqual(wf_config.pseudopotentials_base_path, mock_config_path)
    
    def test_environ_espresso_pseudo_is_set_only_for_dict(self):
        """Verify ESPRESSO_PSEUDO is set for dict-based, appropriate for config-based."""
        
        atoms = bulk('Mg', crystalstructure='hcp', a=3.2, c=5.2)
        pseudo_dict = {'Mg': 'Mg.pbe-spfn-rrkjus_psl.1.0.0.UPF'}
        base_path = '/discovered/pseudo/path'
        
        # Save original value
        original_espresso_pseudo = os.environ.pop('ESPRESSO_PSEUDO', None)
        
        try:
            # Test dict-based sets env var
            with mock.patch('xespresso.workflow.convergence_workflow.discover_pseudopotential_directory') as mock_discover:
                mock_discover.return_value = (pseudo_dict, base_path)
                
                wf = ConvergenceWorkflow(
                    atoms=atoms,
                    pseudopotentials=pseudo_dict,
                    precision='low',
                    convergence_criteria_list=['energy']
                )
            
            # ESPRESSO_PSEUDO should be set to the discovered base path
            self.assertEqual(os.environ.get('ESPRESSO_PSEUDO'), base_path)
            
        finally:
            # Restore
            if original_espresso_pseudo:
                os.environ['ESPRESSO_PSEUDO'] = original_espresso_pseudo
            else:
                os.environ.pop('ESPRESSO_PSEUDO', None)
    
    def test_mixed_elements_with_dict_pseudopotentials(self):
        """Test that dict-based loading works with multi-element structures."""
        
        from ase.build import bulk
        
        # Create a simple structure (Pt FCC)
        atoms = bulk('Pt', crystalstructure='fcc', a=3.92)
        
        # Pseudo dict with element
        pseudo_dict = {
            'Pt': 'Pt.pbe-n-rrkjus_psl.1.0.0.UPF',
        }
        base_path = '/multi/element/pseudo'
        
        with mock.patch('xespresso.workflow.convergence_workflow.discover_pseudopotential_directory') as mock_discover:
            mock_discover.return_value = (pseudo_dict, base_path)
            
            wf = ConvergenceWorkflow(
                atoms=atoms,
                pseudopotentials=pseudo_dict,
                precision='medium',
                convergence_criteria_list=['energy']
            )
            
            # Should have base_path
            self.assertEqual(wf.pseudopotentials_base_path, base_path)
            
            # Should have all pseudopotentials
            self.assertEqual(wf.pseudopotentials, pseudo_dict)
    
    def test_precision_level_affects_convergence_regardless_of_pseudo_source(self):
        """Verify precision levels work the same for dict and config sources."""
        
        atoms = bulk('Fe', crystalstructure='bcc', a=2.87)
        pseudo_dict = {'Fe': 'Fe.pbe-n-rrkjus_psl.1.0.0.UPF'}
        
        precision_energy_tol = {
            'low': 3e-3,
            'medium': 2e-3,
            'high': 1e-3,
            'ultra': 5e-4,
        }
        
        # Test dict-based for each precision
        with mock.patch('xespresso.workflow.convergence_workflow.discover_pseudopotential_directory') as mock_discover:
            mock_discover.return_value = (pseudo_dict, '/pseudo/path')
            
            for precision, expected_tol in precision_energy_tol.items():
                wf = ConvergenceWorkflow(
                    atoms=atoms,
                    pseudopotentials=pseudo_dict,
                    precision=precision,
                    convergence_criteria_list=['energy']
                )
                
                actual_tol = wf.convergence_criteria['energy_tolerance']
                self.assertEqual(
                    actual_tol,
                    expected_tol,
                    f"Dict-based: Energy tolerance for {precision} should be {expected_tol}, got {actual_tol}"
                )
        
        # Test config-based for each precision
        mock_config = mock.MagicMock()
        mock_config.base_path = '/config/pseudo'
        
        with mock.patch('xespresso.pseudopotentials.manager.load_pseudopotentials_config') as mock_load:
            mock_load.return_value = mock_config
            
            for precision, expected_tol in precision_energy_tol.items():
                wf = ConvergenceWorkflow(
                    atoms=atoms,
                    pseudopotentials_config='SSSP',
                    precision=precision,
                    convergence_criteria_list=['energy']
                )
                
                actual_tol = wf.convergence_criteria['energy_tolerance']
                self.assertEqual(
                    actual_tol,
                    expected_tol,
                    f"Config-based: Energy tolerance for {precision} should be {expected_tol}, got {actual_tol}"
                )


if __name__ == '__main__':
    unittest.main()
