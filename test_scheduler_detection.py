#!/usr/bin/env python
"""
Test scheduler type detection in wait_for_batch_jobs
Tests both direct (bash) and SLURM schedulers with mocked remote connections.
"""

import os
import tempfile
from unittest.mock import Mock, MagicMock, patch
from ase.build import bulk
from xespresso.workflow.calculation_workflow import CalculationWorkflow


def test_direct_scheduler_detection():
    """Test that direct scheduler uses ps instead of squeue/sacct"""
    print("\n" + "="*70)
    print("TEST 1: Direct Scheduler Detection")
    print("="*70)
    
    atoms = bulk('Si', 'diamond', a=5.43)
    
    # Create mock calculator
    calc = Mock()
    calc.prefix = 'test_calc'
    calc.package = 'pw'
    calc.directory = '/tmp/test_calc'
    calc.queue = {'scheduler': 'direct', 'execution': 'remote'}
    calc.results = {'energy': -100.5}
    calc.read_results = Mock()
    calc.last_remote_path = '/remote/path'
    
    # Create batch result (not cached, submitted)
    batch_result = {
        'calc': calc,
        'job_id': 'PID:12345',
        'label': 'test_job',
        'submitted': True,
        'completed': False,
    }
    
    # Mock remote connection
    mock_remote = Mock()
    
    # Mock ps command - job completed
    mock_remote.run_command = Mock(side_effect=[
        ('COMPLETED', ''),  # ps check returns COMPLETED
    ])
    
    # Mock file retrieval
    mock_remote.retrieve_file = Mock()
    
    # Attach mock remote to calc
    calc.remote = mock_remote
    
    # Create workflow and test wait_for_batch_jobs
    wf = CalculationWorkflow(
        atoms,
        pseudopotentials_config='SSSP_efficiency',
        machine='snake5',
    )
    
    # Test the wait function with direct scheduler
    results = wf.wait_for_batch_jobs([batch_result], timeout=10, verbose=True)
    
    # Verify ps was called (not squeue)
    calls = [str(call) for call in mock_remote.run_command.call_args_list]
    print(f"\nRemote commands called: {calls}")
    
    # Should have called ps for direct scheduler
    assert any('ps -p' in str(call) for call in calls), "Should use 'ps -p' for direct scheduler"
    assert not any('squeue' in str(call) for call in calls), "Should NOT use 'squeue' for direct scheduler"
    
    # Check result
    assert results[0]['completed'] == True, "Job should be marked as completed"
    assert results[0]['success'] == True, "Job should be successful"
    assert results[0]['job_id'] == 'PID:12345', "Job ID should match"
    
    print("✅ Direct scheduler test PASSED")
    print(f"   - Correctly used 'ps -p' for job status check")
    print(f"   - Job marked as completed: {results[0]['completed']}")
    print(f"   - Job marked as success: {results[0]['success']}")


def test_slurm_scheduler_detection():
    """Test that SLURM scheduler uses squeue/sacct"""
    print("\n" + "="*70)
    print("TEST 2: SLURM Scheduler Detection")
    print("="*70)
    
    atoms = bulk('Si', 'diamond', a=5.43)
    
    # Create mock calculator
    calc = Mock()
    calc.prefix = 'test_calc'
    calc.package = 'pw'
    calc.directory = '/tmp/test_calc'
    calc.queue = {'scheduler': 'slurm', 'execution': 'remote'}
    calc.results = {'energy': -100.5}
    calc.read_results = Mock()
    calc.last_remote_path = '/remote/path'
    
    # Create batch result (not cached, submitted)
    batch_result = {
        'calc': calc,
        'job_id': '12345',  # Regular SLURM job ID
        'label': 'test_job',
        'submitted': True,
        'completed': False,
    }
    
    # Mock remote connection
    mock_remote = Mock()
    
    # Mock squeue and sacct commands
    mock_remote.run_command = Mock(side_effect=[
        ('', ''),  # squeue returns empty (job not in queue anymore)
        ('COMPLETED\n', ''),  # sacct returns COMPLETED
    ])
    
    # Mock file retrieval
    mock_remote.retrieve_file = Mock()
    
    # Attach mock remote to calc
    calc.remote = mock_remote
    
    # Create workflow and test wait_for_batch_jobs
    wf = CalculationWorkflow(
        atoms,
        pseudopotentials_config='SSSP_efficiency',
        machine='snake5',
    )
    
    # Test the wait function with SLURM scheduler
    results = wf.wait_for_batch_jobs([batch_result], timeout=10, verbose=True)
    
    # Verify squeue and sacct were called
    calls = [str(call) for call in mock_remote.run_command.call_args_list]
    print(f"\nRemote commands called: {calls}")
    
    # Should have called squeue for SLURM scheduler
    assert any('squeue' in str(call) for call in calls), "Should use 'squeue' for SLURM scheduler"
    # Should have called sacct when job not in queue
    assert any('sacct' in str(call) for call in calls), "Should use 'sacct' for completed SLURM jobs"
    
    # Check result
    assert results[0]['completed'] == True, "Job should be marked as completed"
    assert results[0]['success'] == True, "Job should be successful"
    assert results[0]['job_id'] == '12345', "Job ID should match"
    
    print("✅ SLURM scheduler test PASSED")
    print(f"   - Correctly used 'squeue' for job status check")
    print(f"   - Correctly used 'sacct' for completed job status")
    print(f"   - Job marked as completed: {results[0]['completed']}")
    print(f"   - Job marked as success: {results[0]['success']}")


def test_energy_extraction():
    """Test that energy is extracted from results"""
    print("\n" + "="*70)
    print("TEST 3: Energy Extraction")
    print("="*70)
    
    atoms = bulk('Si', 'diamond', a=5.43)
    
    # Create mock calculator
    calc = Mock()
    calc.prefix = 'test_calc'
    calc.package = 'pw'
    calc.directory = '/tmp/test_calc'
    calc.queue = {'scheduler': 'direct', 'execution': 'remote'}
    calc.results = {'energy': -100.5}
    calc.read_results = Mock()
    calc.last_remote_path = '/remote/path'
    
    # Create batch result
    batch_result = {
        'calc': calc,
        'job_id': 'PID:12345',
        'label': 'test_job',
        'submitted': True,
        'completed': False,
    }
    
    # Mock remote connection
    mock_remote = Mock()
    mock_remote.run_command = Mock(return_value=('COMPLETED', ''))
    mock_remote.retrieve_file = Mock()
    
    calc.remote = mock_remote
    
    # Create workflow
    wf = CalculationWorkflow(
        atoms,
        pseudopotentials_config='SSSP_efficiency',
        machine='snake5',
    )
    
    # Test wait_for_batch_jobs
    results = wf.wait_for_batch_jobs([batch_result], timeout=10, verbose=True)
    
    # Check energy was extracted
    print(f"\nExtracted energy: {results[0]['energy']}")
    assert results[0]['energy'] == -100.5, "Energy should be extracted from calc.results"
    
    # Check that retrieve_file was called
    assert mock_remote.retrieve_file.called, "Should retrieve output file from remote"
    
    # Check that read_results was called
    assert calc.read_results.called, "Should call calc.read_results()"
    
    print("✅ Energy extraction test PASSED")
    print(f"   - Energy correctly extracted: {results[0]['energy']} eV")
    print(f"   - retrieve_file called: {mock_remote.retrieve_file.called}")
    print(f"   - read_results() called: {calc.read_results.called}")


if __name__ == '__main__':
    try:
        test_direct_scheduler_detection()
        test_slurm_scheduler_detection()
        test_energy_extraction()
        
        print("\n" + "="*70)
        print("✅ ALL TESTS PASSED")
        print("="*70)
    except AssertionError as e:
        print(f"\n❌ TEST FAILED: {e}")
        exit(1)
    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
