"""
Test that remote execution properly uses scheduler.run() instead of local command execution.

This test verifies the fix for the issue where remote SLURM jobs were attempted
to be executed locally, causing "sbatch: command not found" errors (exit code 127).

The issue occurs in ASE 3.22.1 where FileIOCalculator.calculate() directly calls
subprocess.Popen(command) locally, which fails when command is "sbatch job_file"
and sbatch is not available on the local machine.

The fix overrides calculate() to use scheduler.run() for remote execution,
which properly handles SSH connection, file transfer, and remote job submission.
"""

import pytest
from unittest.mock import Mock, MagicMock, patch
from ase import Atoms
from xespresso import Espresso
from xespresso.scheduler import set_queue


class TestRemoteExecutionFix:
    """Test suite for remote execution fix (ASE 3.22.1 compatibility)."""

    def test_calculate_calls_scheduler_run_for_remote(self):
        """Test that calculate() calls scheduler.run() for remote execution."""
        # Create a simple atoms object
        atoms = Atoms('H', positions=[(0, 0, 0)], cell=[10, 10, 10], pbc=True)
        
        # Create queue config for remote execution
        queue_config = {
            "execution": "remote",
            "scheduler": "slurm",
            "remote_host": "test.cluster.edu",
            "remote_user": "testuser",
            "remote_dir": "/home/testuser/calc",
            "remote_auth": {"method": "key", "ssh_key": "~/.ssh/id_rsa"},
            "resources": {
                "nodes": 1,
                "ntasks-per-node": 4,
                "time": "01:00:00"
            }
        }
        
        # Create calculator with mock parameters
        with patch.dict('os.environ', {'ESPRESSO_PSEUDO': '/tmp/pseudo'}):
            calc = Espresso(
                label='/tmp/test_remote_exec',
                queue=queue_config,
                pseudopotentials={'H': 'H.pbe-rrkjus_psl.1.0.0.UPF'},
                input_data={'ecutwfc': 30}
            )
        
        # Mock the scheduler
        calc.scheduler = Mock()
        calc.scheduler.run = Mock()
        
        # Mock read_results to avoid errors
        calc.read_results = Mock()
        
        # Call calculate which should use scheduler.run() for remote execution
        calc.calculate(atoms)
        
        # Verify that scheduler.run() was called
        calc.scheduler.run.assert_called_once()

    def test_calculate_falls_back_to_subprocess_for_local(self):
        """Test that calculate() falls back to subprocess execution for local jobs."""
        # Create a simple atoms object
        atoms = Atoms('H', positions=[(0, 0, 0)], cell=[10, 10, 10], pbc=True)
        
        # Create queue config for local execution
        queue_config = {
            "execution": "local",
            "scheduler": "direct"
        }
        
        # Create calculator
        with patch.dict('os.environ', {'ESPRESSO_PSEUDO': '/tmp/pseudo'}):
            calc = Espresso(
                label='/tmp/test_local_exec',
                queue=queue_config,
                pseudopotentials={'H': 'H.pbe-rrkjus_psl.1.0.0.UPF'},
                input_data={'ecutwfc': 30}
            )
        
        # Mock the scheduler
        calc.scheduler = Mock()
        calc.scheduler.run = Mock()
        
        # Mock read_results to avoid errors
        calc.read_results = Mock()
        
        # Mock subprocess.Popen to verify it gets called for local execution
        with patch('subprocess.Popen') as mock_popen:
            mock_proc = Mock()
            mock_proc.wait.return_value = 0  # Success
            mock_popen.return_value = mock_proc
            
            calc.calculate(atoms)
            
            # Verify that subprocess.Popen was called for local execution
            mock_popen.assert_called_once()
            
            # Verify that scheduler.run() was NOT called for local execution
            calc.scheduler.run.assert_not_called()

    def test_set_queue_doesnt_set_command_for_remote(self):
        """Test that set_queue doesn't set an executable command for remote execution."""
        # Create mock calculator
        mock_calc = Mock()
        mock_calc.prefix = "test"
        mock_calc.package = "pw"
        mock_calc.parallel = ""
        mock_calc.directory = "/tmp/test"
        
        # Create remote queue config
        queue_config = {
            "execution": "remote",
            "scheduler": "slurm",
            "remote_host": "test.cluster.edu",
            "remote_user": "testuser",
            "remote_dir": "/home/testuser/calc",
            "remote_auth": {"method": "key", "ssh_key": "~/.ssh/id_rsa"}
        }
        
        mock_calc.queue = queue_config
        
        # Mock the scheduler factory
        with patch('xespresso.scheduler.get_scheduler') as mock_get_scheduler:
            mock_scheduler = Mock()
            mock_scheduler.write_script = Mock()
            mock_scheduler.submit_command = Mock(return_value="sbatch job_file")
            mock_get_scheduler.return_value = mock_scheduler
            
            # Call set_queue
            with patch.dict('os.environ', {'ASE_ESPRESSO_COMMAND': 'pw.x -in test.pwi'}):
                set_queue(mock_calc)
            
            # Verify that command was set to placeholder for remote execution
            assert "Remote execution" in mock_calc.command or "echo" in mock_calc.command
            # The command should NOT be the actual sbatch command for remote execution
            assert mock_calc.command != "sbatch job_file"

    def test_set_queue_sets_command_for_local(self):
        """Test that set_queue sets the actual command for local execution."""
        # Create mock calculator
        mock_calc = Mock()
        mock_calc.prefix = "test"
        mock_calc.package = "pw"
        mock_calc.parallel = ""
        mock_calc.directory = "/tmp/test"
        
        # Create local queue config
        queue_config = {
            "execution": "local",
            "scheduler": "direct"
        }
        
        mock_calc.queue = queue_config
        
        # Mock the scheduler factory
        with patch('xespresso.scheduler.get_scheduler') as mock_get_scheduler:
            mock_scheduler = Mock()
            mock_scheduler.write_script = Mock()
            mock_scheduler.submit_command = Mock(return_value="bash job_file")
            mock_get_scheduler.return_value = mock_scheduler
            
            # Mock SLURM check to prevent it from being called
            with patch('xespresso.scheduler.check_slurm_available'):
                # Call set_queue
                with patch.dict('os.environ', {'ASE_ESPRESSO_COMMAND': 'pw.x -in test.pwi'}):
                    set_queue(mock_calc)
                
                # Verify that command was set to actual command for local execution
                assert mock_calc.command == "bash job_file"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
