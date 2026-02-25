#!/usr/bin/env python3
"""
Test the _monitor_remote_job() method
"""

from ase.build import bulk
from xespresso import CalculationWorkflow
from unittest.mock import Mock, patch
import subprocess

print("Testing _monitor_remote_job() method...")
print("=" * 70)

# Create a workflow
atoms = bulk("Al", "fcc", a=4.0)
wf = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials={"Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF"},
    protocol="moderate",
)

# Create mock calculator
mock_calc = Mock()
mock_calc.results = {'output': 'JOB DONE\nconvergence has been achieved\n'}

# Test 1: Job already completed (not in queue)
print("\n[TEST 1] Job already completed (not in squeue)")
print("-" * 70)
with patch('subprocess.run') as mock_run:
    # squeue returns error code (job not found = already done)
    mock_run.return_value = Mock(returncode=1, stdout='', stderr='Job not found')
    
    result = wf._monitor_remote_job(mock_calc, job_id='12345', timeout=60, poll_interval=5)
    print(f"Result: {result['message']}")
    assert result['state'] == 'COMPLETED', "Should detect job as COMPLETED"
    assert result['success'] == True, "Should be successful"
    print("✓ PASS")

# Test 2: Job running successfully
print("\n[TEST 2] Job running successfully")
print("-" * 70)
call_count = 0
def mock_run_side_effect(*args, **kwargs):
    global call_count
    call_count += 1
    
    if call_count == 1:
        # First call: PENDING
        return Mock(returncode=0, stdout='PENDING,None,0-00:10', stderr='')
    elif call_count == 2:
        # Second call: RUNNING
        return Mock(returncode=0, stdout='RUNNING,,0-00:30', stderr='')
    else:
        # Third call: COMPLETED (not in queue)
        return Mock(returncode=1, stdout='', stderr='Job not found')

with patch('subprocess.run', side_effect=mock_run_side_effect):
    result = wf._monitor_remote_job(mock_calc, job_id='12346', timeout=120, poll_interval=1)
    print(f"Result: {result['message']}")
    assert result['success'] == True, "Should be successful"
    print("✓ PASS")

# Test 3: Job failed
print("\n[TEST 3] Job failed with error")
print("-" * 70)
with patch('subprocess.run') as mock_run:
    mock_run.return_value = Mock(returncode=0, stdout='FAILED,OutOfMemory,0-00:45', stderr='')
    
    result = wf._monitor_remote_job(mock_calc, job_id='12347', timeout=60, poll_interval=5)
    print(f"Result: {result['message']}")
    assert 'FAILED' in result['state'], "Should detect job as FAILED"
    assert 'OutOfMemory' in result['reason'], "Should capture error reason"
    print("✓ PASS")

# Test 4: Job cancelled
print("\n[TEST 4] Job cancelled by user")
print("-" * 70)
with patch('subprocess.run') as mock_run:
    mock_run.return_value = Mock(returncode=0, stdout='CANCELLED,Cancelled by user,0-00:20', stderr='')
    
    result = wf._monitor_remote_job(mock_calc, job_id='12348', timeout=60, poll_interval=5)
    print(f"Result: {result['message']}")
    assert 'CANCELLED' in result['state'], "Should detect job as CANCELLED"
    print("✓ PASS")

print("\n" + "=" * 70)
print("✓ ALL TESTS PASSED!")
print("=" * 70)
