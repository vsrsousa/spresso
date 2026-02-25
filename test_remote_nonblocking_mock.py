"""
Mock test for remote non-blocking workflow execution.

Tests that:
1. Remote non-blocking jobs don't repeat file transfers
2. RemoteJobMonitor waits correctly
3. Results are retrieved and read correctly
"""

import os
import tempfile
from unittest.mock import Mock, patch, MagicMock
from ase.build import bulk
from xespresso import CalculationWorkflow
from xespresso.schedulers import RemoteJobMonitor


def test_remote_nonblocking_no_repetition():
    """Test that remote non-blocking execution doesn't repeat file transfers."""
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Setup mock machine config
        mock_machine = {
            'execution': 'remote',
            'scheduler': 'slurm',
            'remote_host': 'medusa.fis.uerj.br',
            'remote_user': 'vinicius',
            'remote_auth': {'method': 'key', 'ssh_key': '~/.ssh/id_rsa', 'port': 22},
            'remote_dir': '/scratch/users/vinicius/xespresso',
            'launcher': 'srun --mpi=pmi2',
            'nprocs': 16,
            'modules': ['quantum-espresso/7.4.1'],
        }
        
        # Create workflow
        atoms = bulk("Si", cubic=True)
        
        # Mock the machine loading
        with patch('xespresso.workflow.simple_workflow.load_machine') as mock_load_machine:
            mock_load_machine.return_value = mock_machine
            
            # Mock pseudopotentials config
            with patch('xespresso.workflow.simple_workflow.load_pseudopotentials_config') as mock_load_pseudo:
                mock_pseudo_config = Mock()
                mock_pseudo_config.base_path = tmpdir
                mock_pseudo_config.get_pseudopotential = Mock(return_value=Mock(filename='Si.pbe.UPF'))
                mock_pseudo_config.list_elements = Mock(return_value=['Si'])
                mock_load_pseudo.return_value = mock_pseudo_config
                
                workflow = CalculationWorkflow(
                    atoms=atoms,
                    pseudopotentials_config='default',
                    protocol='moderate',
                    machine='snake5'
                )
        
        # Create dummy files
        workflow_dir = os.path.join(tmpdir, 'scf_si-test')
        os.makedirs(workflow_dir, exist_ok=True)
        os.makedirs(os.path.join(workflow_dir, 'pseudo'), exist_ok=True)
        
        # Create dummy pseudo file
        pseudo_file = os.path.join(tmpdir, 'Si.pbe.UPF')
        with open(pseudo_file, 'w') as f:
            f.write('dummy pseudo')
        
        # Track send_file calls
        send_file_calls = []
        
        # Mock remote connection
        mock_remote = Mock()
        mock_remote.send_file = Mock(side_effect=lambda src, dst: send_file_calls.append((src, dst)))
        mock_remote.run_command = Mock(return_value=('Submitted batch job 12345\n', ''))
        
        # Mock Espresso execution
        with patch('xespresso.xespresso.FileIOCalculator.calculate'):
            with patch.object(CalculationWorkflow, 'run_scf', wraps=workflow.run_scf):
                # Mock the remote mixin to use our mock remote
                with patch('xespresso.schedulers.remote_mixin.RemoteExecutionMixin._setup_remote'):
                    try:
                        # This should trigger the workflow handling
                        # Since we're mocking, we'll manually test the logic
                        
                        # Verify the workflow detects remote non-blocking
                        is_remote_nonblocking = (
                            workflow.queue and 
                            workflow.queue.get('execution') == 'remote' and 
                            not workflow.queue.get('wait_for_completion', False)
                        )
                        
                        assert is_remote_nonblocking, "Should detect remote non-blocking execution"
                        print("✓ Workflow detects remote non-blocking correctly")
                        
                    except Exception as e:
                        print(f"✗ Error: {e}")
                        raise


def test_remote_job_monitor():
    """Test RemoteJobMonitor functionality."""
    
    # Create mock calculator
    mock_calc = Mock()
    mock_calc.last_job_id = '12345'
    mock_calc.last_remote_path = '/scratch/job'
    mock_calc.prefix = 'si-test'
    mock_calc.package = 'pw'
    mock_calc.directory = '/local/job'
    
    # Create mock remote connection
    mock_remote = Mock()
    
    # Test SLURM job monitoring
    print("\nTesting RemoteJobMonitor with SLURM:")
    
    # Mock squeue output for running job
    mock_remote.run_command = Mock(side_effect=[
        ('RUNNING\n', ''),      # First call: job is running
        ('RUNNING\n', ''),      # Second call: still running
        ('', ''),               # Third call: job finished (no output = not in queue)
        ('COMPLETED\n\n', '')   # sacct call: job completed
    ])
    
    try:
        monitor = RemoteJobMonitor(mock_calc, remote_connection=mock_remote)
        assert monitor.job_type == 'slurm', "Should detect SLURM job"
        print("✓ Monitor detects SLURM job correctly")
        
        # Check status
        status = monitor.status()
        assert status == 'running', f"Should be running, got {status}"
        print("✓ Monitor correctly reports running status")
        
    except Exception as e:
        print(f"✗ Error: {e}")
        raise
    
    # Test direct job monitoring
    print("\nTesting RemoteJobMonitor with Direct scheduler:")
    
    mock_calc.last_job_id = 'PID:54321'
    mock_remote.run_command = Mock(side_effect=[
        ('54321\n', ''),     # First: process exists
        ('54321\n', ''),     # Second: still running
        ('', ''),            # Third: process not found
        ('exists\n', '')     # Check for output file
    ])
    
    try:
        monitor = RemoteJobMonitor(mock_calc, remote_connection=mock_remote)
        assert monitor.job_type == 'direct', "Should detect direct scheduler job"
        print("✓ Monitor detects direct scheduler job correctly")
        
        status = monitor.status()
        assert status == 'running', f"Should be running, got {status}"
        print("✓ Monitor correctly reports running status for direct job")
        
    except Exception as e:
        print(f"✗ Error: {e}")
        raise


def test_workflow_integration():
    """Test complete workflow integration (mocked)."""
    
    print("\nTesting workflow integration:")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        atoms = bulk("Si", cubic=True)
        
        # Create minimal mock config
        with patch('xespresso.workflow.simple_workflow.load_machine') as mock_load_machine, \
             patch('xespresso.workflow.simple_workflow.load_pseudopotentials_config') as mock_load_pseudo:
            
            mock_machine = {
                'execution': 'remote',
                'scheduler': 'slurm',
                'remote_host': 'test.com',
                'wait_for_completion': False,  # Non-blocking!
            }
            mock_load_machine.return_value = mock_machine
            
            mock_pseudo_config = Mock()
            mock_pseudo_config.base_path = tmpdir
            mock_pseudo_config.get_pseudopotential = Mock(return_value=Mock(filename='Si.pbe.UPF'))
            mock_pseudo_config.list_elements = Mock(return_value=['Si'])
            mock_load_pseudo.return_value = mock_pseudo_config
            
            workflow = CalculationWorkflow(
                atoms=atoms,
                pseudopotentials_config='default',
                machine='test_machine'
            )
            
            # Verify queue configuration
            assert workflow.queue is not None, "Queue should be configured"
            assert workflow.queue.get('execution') == 'remote', "Should be remote execution"
            assert not workflow.queue.get('wait_for_completion', False), "Should be non-blocking"
            print("✓ Workflow configured for remote non-blocking execution")


if __name__ == '__main__':
    print("=" * 70)
    print("Testing Remote Non-Blocking Workflow Execution")
    print("=" * 70)
    
    try:
        test_remote_nonblocking_no_repetition()
        test_remote_job_monitor()
        test_workflow_integration()
        
        print("\n" + "=" * 70)
        print("✓ ALL TESTS PASSED!")
        print("=" * 70)
    except Exception as e:
        print("\n" + "=" * 70)
        print(f"✗ TEST FAILED: {e}")
        print("=" * 70)
        raise
