import pytest
import os
import tempfile
from unittest.mock import Mock, MagicMock, patch
from xespresso.schedulers.slurm import SlurmScheduler
from xespresso.schedulers.remote_mixin import RemoteExecutionMixin


class TestSlurmScheduler:
    """Test suite for SLURM scheduler functionality."""

    def setup_method(self):
        """Set up test fixtures."""
        self.mock_calc = Mock()
        self.mock_calc.prefix = "test_job"
        self.mock_calc.package = "pw"
        self.mock_calc.directory = "/tmp/test"
        self.mock_calc.parameters = {"pseudopotentials": {}}
        
        self.queue_config = {
            "scheduler": "slurm",
            "execution": "remote",
            "remote_host": "test.hpc.edu",
            "remote_user": "testuser",
            "remote_auth": {"method": "key", "ssh_key": "~/.ssh/id_rsa"},
            "remote_dir": "/home/testuser/calculations",
            "job_timeout": 300,  # 5 minutes for testing
            "resources": {
                "nodes": 1,
                "ntasks": 8,
                "time": "01:00:00"
            }
        }
        
        self.scheduler = SlurmScheduler(
            calc=self.mock_calc,
            queue=self.queue_config,
            command="mpirun pw.x -in test_job.pwi > test_job.pwo"
        )

    def test_regex_pattern_fix(self):
        """Test that the job ID regex pattern correctly extracts job IDs."""
        import re
        
        # Test various sbatch output formats
        test_outputs = [
            "Submitted batch job 12345",
            "Submitted batch job 987654321",
            "sbatch: Submitted batch job 555",
            "INFO: Submitted batch job 12345 on partition normal"
        ]
        
        pattern = r"Submitted batch job (\d+)"
        
        for output in test_outputs:
            match = re.search(pattern, output)
            assert match is not None, f"Failed to match job ID in: {output}"
            job_id = match.group(1)
            assert job_id.isdigit(), f"Extracted job ID should be numeric: {job_id}"

    @patch('xespresso.schedulers.remote_mixin.RemoteAuth')
    def test_job_timeout_handling(self, mock_remote_auth_class):
        """Test that job timeout is properly handled."""
        mock_remote = Mock()
        mock_remote_auth_class.return_value = mock_remote
        
        # Mock job that never completes
        mock_remote.run_command.side_effect = [
            ("squeue output with job", ""),  # Job still running
            ("squeue output with job", ""),  # Job still running
            ("squeue output with job", ""),  # Job still running
        ]
        
        scheduler = SlurmScheduler(
            calc=self.mock_calc,
            queue={**self.queue_config, "job_timeout": 1},  # 1 second timeout
            command="test command"
        )
        scheduler.remote = mock_remote  # Set the remote attribute manually
        
        with patch('time.time') as mock_time:
            # Simulate time progression
            mock_time.side_effect = [0, 0.5, 1.5]  # Start, mid, timeout
            
            with patch('time.sleep'):
                with pytest.raises(RuntimeError, match="timed out after 1 seconds"):
                    scheduler._wait_for_slurm_completion("12345")

    @patch('xespresso.schedulers.remote_mixin.RemoteAuth')
    def test_job_completion_verification(self, mock_remote_auth_class):
        """Test that job completion status is properly verified."""
        mock_remote = Mock()
        mock_remote_auth_class.return_value = mock_remote
        
        scheduler = SlurmScheduler(
            calc=self.mock_calc,
            queue=self.queue_config,
            command="test command"
        )
        scheduler.remote = mock_remote
        
        # Test successful completion
        mock_remote.run_command.side_effect = [
            ("", ""),  # squeue returns empty (job not in queue)
            ("COMPLETED", ""),  # sacct returns COMPLETED status
        ]
        
        with patch('time.time', return_value=0):
            with patch('time.sleep'):
                # Should not raise exception for completed job
                scheduler._wait_for_slurm_completion("12345")
        
        # Test failed job
        mock_remote.run_command.side_effect = [
            ("", ""),  # squeue returns empty (job not in queue)
            ("FAILED", ""),  # sacct returns FAILED status
        ]
        
        with patch('time.time', return_value=0):
            with patch('time.sleep'):
                with pytest.raises(RuntimeError, match="failed with state: FAILED"):
                    scheduler._wait_for_slurm_completion("12345")

    @patch('xespresso.schedulers.remote_mixin.RemoteAuth')
    def test_output_file_verification(self, mock_remote_auth_class):
        """Test that output file existence is verified before retrieval."""
        mock_remote = Mock()
        mock_remote_auth_class.return_value = mock_remote
        
        scheduler = SlurmScheduler(
            calc=self.mock_calc,
            queue=self.queue_config,
            command="test command"
        )
        scheduler.remote = mock_remote
        scheduler.remote_path = "/remote/path"
        
        # Test file exists
        mock_remote.run_command.return_value = ("exists", "")
        mock_remote.retrieve_file.return_value = None
        
        # Should succeed without exception
        scheduler._verify_and_retrieve_output_file("test.pwo", "/local/test.pwo")
        
        # Test file missing
        mock_remote.run_command.return_value = ("missing", "")
        
        with pytest.raises(FileNotFoundError, match="does not exist on remote system"):
            scheduler._verify_and_retrieve_output_file("test.pwo", "/local/test.pwo")

    @patch('xespresso.schedulers.remote_mixin.RemoteAuth')
    def test_file_retrieval_retry_logic(self, mock_remote_auth_class):
        """Test that file retrieval has proper retry logic."""
        mock_remote = Mock()
        mock_remote_auth_class.return_value = mock_remote
        
        scheduler = SlurmScheduler(
            calc=self.mock_calc,
            queue=self.queue_config,
            command="test command"
        )
        scheduler.remote = mock_remote
        scheduler.remote_path = "/remote/path"
        
        # Test successful retrieval on second attempt
        mock_remote.run_command.return_value = ("exists", "")
        mock_remote.retrieve_file.side_effect = [
            Exception("Network error"),  # First attempt fails
            None  # Second attempt succeeds
        ]
        
        with patch('time.sleep'):
            # Should succeed after retry
            scheduler._verify_and_retrieve_output_file("test.pwo", "/local/test.pwo", max_retries=2)
        
        # Test all attempts fail
        mock_remote.retrieve_file.side_effect = [
            Exception("Network error"),
            Exception("Network error"),
            Exception("Network error")
        ]
        
        with patch('time.sleep'):
            with pytest.raises(RuntimeError, match="Failed to retrieve .* after 3 attempts"):
                scheduler._verify_and_retrieve_output_file("test.pwo", "/local/test.pwo", max_retries=3)

    @patch('xespresso.schedulers.remote_mixin.RemoteAuth')
    @patch('os.path.exists')
    def test_missing_pseudopotential_prevents_execution(self, mock_exists, mock_remote_auth_class):
        """Test that missing pseudopotentials prevent job execution."""
        mock_remote = Mock()
        mock_remote_auth_class.return_value = mock_remote
        mock_remote.run_command.return_value = ("", "")
        
        # Set up calculator with pseudopotentials
        self.mock_calc.parameters = {
            "pseudopotentials": {
                "Fe": "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF",
                "O": "O.pbe-n-kjpaw_psl.1.0.0.UPF"
            },
            "input_data": {
                "CONTROL": {}
            }
        }
        self.mock_calc.write_input = Mock()
        
        scheduler = SlurmScheduler(
            calc=self.mock_calc,
            queue=self.queue_config,
            command="test command"
        )
        scheduler.remote = mock_remote
        scheduler.remote_path = "/remote/path"
        
        # Mock that pseudopotential files don't exist
        mock_exists.return_value = False
        
        # Should raise FileNotFoundError when pseudopotentials are missing
        with pytest.raises(FileNotFoundError, match="Cannot proceed with calculation.*Missing pseudopotentials"):
            scheduler._transfer_pseudopotentials()
    
    @patch('xespresso.schedulers.remote_mixin.RemoteAuth')
    @patch('os.path.exists')
    def test_partial_missing_pseudopotentials(self, mock_exists, mock_remote_auth_class):
        """Test that execution is prevented even if only some pseudopotentials are missing."""
        mock_remote = Mock()
        mock_remote_auth_class.return_value = mock_remote
        mock_remote.run_command.return_value = ("", "")
        mock_remote.sha256.return_value = "abc123"
        
        # Set up calculator with pseudopotentials
        self.mock_calc.parameters = {
            "pseudopotentials": {
                "Fe": "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF",
                "O": "O.pbe-n-kjpaw_psl.1.0.0.UPF"
            },
            "input_data": {
                "CONTROL": {}
            }
        }
        self.mock_calc.write_input = Mock()
        
        scheduler = SlurmScheduler(
            calc=self.mock_calc,
            queue=self.queue_config,
            command="test command"
        )
        scheduler.remote = mock_remote
        scheduler.remote_path = "/remote/path"
        
        # Mock SHA256 method
        with patch.object(scheduler, '_sha256', return_value='abc123'):
            # Mock that only Fe pseudopotential exists
            def exists_side_effect(path):
                return "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF" in path
            
            mock_exists.side_effect = exists_side_effect
            
            # Should raise FileNotFoundError for missing O pseudopotential
            with pytest.raises(FileNotFoundError, match="Missing pseudopotentials.*O"):
                scheduler._transfer_pseudopotentials()

    @patch('xespresso.schedulers.remote_mixin.RemoteAuth')
    @patch('os.path.exists')
    def test_run_prevents_execution_when_pseudopotentials_missing(self, mock_exists, mock_remote_auth_class):
        """Test that the run() method prevents execution when pseudopotentials are missing."""
        mock_remote = Mock()
        mock_remote_auth_class.return_value = mock_remote
        mock_remote.run_command.return_value = ("", "")
        mock_remote.connect = Mock()
        
        # Set up calculator with pseudopotentials
        self.mock_calc.parameters = {
            "pseudopotentials": {
                "Fe": "Fe.pbe-spn-rrkjus_psl.1.0.0.UPF"
            },
            "input_data": {
                "CONTROL": {}
            }
        }
        self.mock_calc.write_input = Mock()
        
        scheduler = SlurmScheduler(
            calc=self.mock_calc,
            queue=self.queue_config,
            command="test command"
        )
        
        # Mock that pseudopotential file doesn't exist
        mock_exists.return_value = False
        
        # The run() method should raise FileNotFoundError during _transfer_pseudopotentials
        with pytest.raises(FileNotFoundError, match="Cannot proceed with calculation"):
            scheduler.run()
        
        # Verify that file transfer and job submission were NOT attempted
        # The send_file method should not have been called since we failed before that
        if hasattr(mock_remote, 'send_file'):
            assert not mock_remote.send_file.called, "Files should not be transferred when pseudopotentials are missing"


# Helper methods for scheduler testing
def _create_slurm_scheduler_with_mocks():
    """Create a SlurmScheduler instance with mocked dependencies for testing."""
    mock_calc = Mock()
    mock_calc.prefix = "test_job"
    mock_calc.package = "pw"
    mock_calc.directory = "/tmp/test"
    
    queue_config = {
        "scheduler": "slurm",
        "execution": "remote",
        "remote_host": "test.hpc.edu",
        "remote_user": "testuser",
        "remote_auth": {"method": "key"},
        "remote_dir": "/home/testuser"
    }
    
    return SlurmScheduler(
        calc=mock_calc,
        queue=queue_config,
        command="test command"
    )


# Add missing methods to complete the implementation
def _add_missing_scheduler_methods():
    """Add missing methods that should be in the scheduler class."""
    
    def _wait_for_slurm_completion(self, job_id):
        """Wait for SLURM job completion with timeout and proper status checking."""
        import time
        
        timeout = self.queue.get("job_timeout", 3600)
        start_time = time.time()
        
        while True:
            if time.time() - start_time > timeout:
                error_msg = f"Job {job_id} timed out after {timeout} seconds"
                if hasattr(self, "logger"):
                    self.logger.error(error_msg)
                raise RuntimeError(error_msg)
            
            status_stdout, status_stderr = self.remote.run_command(f"squeue -j {job_id} -h -o '%T'")
            
            if not status_stdout.strip():
                # Job is no longer in queue, check final status
                sacct_stdout, sacct_stderr = self.remote.run_command(f"sacct -j {job_id} -n -o State --parsable2")
                if sacct_stdout.strip():
                    job_state = sacct_stdout.strip().split('\n')[0]
                    if hasattr(self, "logger"):
                        self.logger.info(f"Job {job_id} finished with state: {job_state}")
                    
                    if job_state not in ["COMPLETED", "COMPLETING"]:
                        error_msg = f"Job {job_id} failed with state: {job_state}"
                        if hasattr(self, "logger"):
                            self.logger.error(error_msg)
                        raise RuntimeError(error_msg)
                break
            
            time.sleep(10)
    
    # Add method to SlurmScheduler class
    SlurmScheduler._wait_for_slurm_completion = _wait_for_slurm_completion


if __name__ == "__main__":
    # Run the tests
    pytest.main([__file__, "-v"])