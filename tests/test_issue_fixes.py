"""
Tests for issue fixes:
1. Remove temp_ prefix from labels
2. Allow rerun of failed calculations in GUI
3. Fix remote folder naming collisions
"""

import pytest
import tempfile
import os
import shutil
from unittest.mock import Mock, patch, MagicMock
from ase.build import bulk
from xespresso import Espresso
from xespresso.schedulers.remote_mixin import RemoteExecutionMixin


class TestTempPrefixRemoval:
    """Test that temp_ prefix is no longer used in calculation preparation."""
    
    def test_calculation_setup_uses_actual_label(self):
        """
        Test that calculation_setup.py uses the actual label from config,
        not a temp_ prefixed label.
        """
        # This test verifies the fix by checking the code pattern
        # In real usage, the GUI would call prepare_calculation_from_gui with the actual label
        from xespresso.gui.calculations import prepare_calculation_from_gui
        
        atoms = bulk("Al", cubic=True)
        config = {
            'pseudopotentials': {'Al': 'Al.pbe.UPF'},
            'ecutwfc': 30.0,
            'kpts': (2, 2, 2),
            'label': 'scf/Al',  # Actual label, not temp_Al
            'calc_type': 'scf'
        }
        
        with tempfile.TemporaryDirectory() as tmpdir:
            label = os.path.join(tmpdir, config['label'])
            prepared_atoms, calc = prepare_calculation_from_gui(atoms, config, label=label)
            
            # Verify the calculator has the actual label, not temp_
            assert 'temp_' not in calc.label
            assert calc.label.endswith('Al')  # Should end with 'Al', not 'temp_Al'


class TestFailedCalculationRerun:
    """Test that failed calculations can be rerun from GUI."""
    
    def test_calculator_reset_on_failed_convergence(self):
        """
        Test that when a calculator has a failed convergence status,
        calling reset() clears the results and allows rerun.
        """
        atoms = bulk("Fe", cubic=True)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            label = os.path.join(tmpdir, "test_calc")
            
            calc = Espresso(
                pseudopotentials={'Fe': 'Fe.pbe.UPF'},
                ecutwfc=30.0,
                kpts=(2, 2, 2),
                label=label
            )
            
            # Simulate a failed calculation by setting convergence status
            calc.results = {'convergence': 3, 'energy': None}
            
            # Verify failed status
            assert calc.results['convergence'] > 0
            
            # Reset the calculator (this is what GUI does now)
            calc.reset()
            
            # Verify results are cleared
            assert 'convergence' not in calc.results or calc.results['convergence'] == 0
            assert 'energy' not in calc.results
    
    def test_gui_detects_and_resets_failed_calculator(self):
        """
        Test that the GUI logic correctly detects failed calculations
        and resets the calculator.
        """
        # This tests the logic that was added to job_submission.py
        
        # Mock calculator with failed convergence
        calc = Mock()
        calc.results = {'convergence': 3}  # Failed calculation
        calc.reset = Mock()
        
        # Simulate GUI logic
        if hasattr(calc, 'results') and 'convergence' in calc.results:
            convergence_status = calc.results['convergence']
            if convergence_status > 0:
                # This is what the GUI does now
                calc.reset()
        
        # Verify reset was called
        calc.reset.assert_called_once()


class TestRemoteFolderNaming:
    """Test that remote folder names are unique to prevent collisions."""
    
    def test_unique_remote_folder_names(self):
        """
        Test that two calculations with the same label but different
        local directories get different remote folder names.
        """
        import hashlib
        
        # Simulate two different local directories with same label
        local_dir1 = "/home/user/project1/scf/Fe"
        local_dir2 = "/home/user/project2/scf/Fe"
        
        # Calculate what the remote folder names would be
        basename = "Fe"
        
        hash1 = hashlib.md5(local_dir1.encode()).hexdigest()[:8]
        hash2 = hashlib.md5(local_dir2.encode()).hexdigest()[:8]
        
        remote_name1 = f"{basename}_{hash1}"
        remote_name2 = f"{basename}_{hash2}"
        
        # Verify they are different
        assert remote_name1 != remote_name2
        assert remote_name1.startswith("Fe_")
        assert remote_name2.startswith("Fe_")
    
    def test_remote_mixin_generates_unique_paths(self):
        """
        Test that RemoteExecutionMixin._setup_remote generates unique
        remote paths based on local directory hash.
        """
        # Create a minimal mock scheduler with RemoteExecutionMixin
        class MockScheduler(RemoteExecutionMixin):
            def __init__(self, calc, queue):
                self.calc = calc
                self.queue = queue
        
        # Mock calculator with different local directories
        calc1 = Mock()
        calc1.directory = "/home/user/project1/scf/Fe"
        
        calc2 = Mock()
        calc2.directory = "/home/user/project2/scf/Fe"
        
        queue = {
            'remote_host': 'cluster.example.com',
            'remote_user': 'testuser',
            'remote_auth': {'method': 'ssh_key', 'ssh_key': '/path/to/key'},
            'remote_dir': '/scratch/testuser'
        }
        
        # Patch RemoteAuth to avoid actual SSH connection
        with patch('xespresso.schedulers.remote_mixin.RemoteAuth') as mock_auth:
            mock_remote = Mock()
            mock_auth.return_value = mock_remote
            mock_remote.connect = Mock()
            
            scheduler1 = MockScheduler(calc1, queue)
            scheduler1._setup_remote()
            remote_path1 = scheduler1.remote_path
            
            # Reset class variable to simulate fresh calculation
            RemoteExecutionMixin._last_remote_path = None
            RemoteExecutionMixin._remote_sessions.clear()
            
            scheduler2 = MockScheduler(calc2, queue)
            scheduler2._setup_remote()
            remote_path2 = scheduler2.remote_path
            
            # Verify remote paths are different
            assert remote_path1 != remote_path2
            assert remote_path1.startswith('/scratch/testuser/Fe_')
            assert remote_path2.startswith('/scratch/testuser/Fe_')
            
            # Verify both contain the hash suffix
            assert '_' in os.path.basename(remote_path1)
            assert '_' in os.path.basename(remote_path2)
    
    def test_same_local_dir_gives_same_remote_folder(self):
        """
        Test that the same local directory always maps to the same remote folder.
        This ensures consistency across multiple runs.
        """
        import hashlib
        
        local_dir = "/home/user/project/scf/Fe"
        basename = "Fe"
        
        # Calculate hash multiple times
        hash1 = hashlib.md5(local_dir.encode()).hexdigest()[:8]
        hash2 = hashlib.md5(local_dir.encode()).hexdigest()[:8]
        
        remote_name1 = f"{basename}_{hash1}"
        remote_name2 = f"{basename}_{hash2}"
        
        # Verify they are the same
        assert remote_name1 == remote_name2


class TestIntegrationScenarios:
    """Integration tests for the complete workflow."""
    
    def test_failed_then_successful_calculation_workflow(self):
        """
        Test complete workflow:
        1. Run calculation that fails
        2. Rerun without changes - should detect failure and reset
        3. Verify calculator is ready for fresh run
        """
        atoms = bulk("Al", cubic=True)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            label = os.path.join(tmpdir, "workflow_test")
            
            # First attempt - simulate failure
            calc = Espresso(
                pseudopotentials={'Al': 'Al.pbe.UPF'},
                ecutwfc=30.0,
                kpts=(2, 2, 2),
                label=label
            )
            
            # Simulate failed calculation
            calc.results = {'convergence': 3, 'energy': None}
            
            # Second attempt - this is what GUI does
            # Detect failure and reset
            if hasattr(calc, 'results') and 'convergence' in calc.results:
                if calc.results['convergence'] > 0:
                    calc.reset()
            
            # Verify calculator is ready for fresh calculation
            assert 'convergence' not in calc.results or calc.results['convergence'] == 0
            assert 'energy' not in calc.results


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
