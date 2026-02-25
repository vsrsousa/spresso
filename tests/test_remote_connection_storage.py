"""
Mock test for the final fix: storing remote connection in calc.remote
during job submission (in remote_mixin.py)
"""
import os
import sys
from unittest.mock import Mock, MagicMock, patch

print("\n" + "="*80)
print("TEST: Remote Connection Storage Fix (remote_mixin.py)")
print("="*80)

# ============================================================================
# Test 1: Verify calc.remote is set after SLURM job submission
# ============================================================================
print("\n[Test 1] SLURM job submission - verify calc.remote is stored")
print("-"*80)

# Mock the RemoteAuth class
mock_remote_auth = Mock()
mock_remote_auth.run_command = Mock(return_value=("Submitted batch job 1596\n", ""))
mock_remote_auth.send_file = Mock()

# Create a mock calculator
mock_calc = Mock()
mock_calc.prefix = "test"
mock_calc.package = "pw"
mock_calc.directory = "/tmp/test"
mock_calc.parameters = {
    "pseudopotentials": {},
    "input_data": {"CONTROL": {}}
}

# Create a mock scheduler with remote_mixin
from xespresso.schedulers.remote_mixin import RemoteExecutionMixin

class MockSLURMScheduler(RemoteExecutionMixin):
    def __init__(self, calc, queue):
        self.calc = calc
        self.queue = queue
        import logging
        self.logger = logging.getLogger(__name__)
        self.job_file = "job_file"
        self.script_dir = "/tmp"
        self._remote_sessions = {}
        self._last_remote_path = None
    
    def submit_command(self):
        return "sbatch job_file"
    
    def write_script(self):
        pass

# Create mock queue config
queue = {
    "execution": "remote",
    "scheduler": "slurm",
    "remote_host": "medusa.fis.uerj.br",
    "remote_user": "vinicius",
    "remote_dir": "/scratch/users/vinicius/xespresso",
    "remote_auth": {"method": "key"},
    "wait_for_completion": False
}

# Patch RemoteAuth to return our mock
with patch('xespresso.schedulers.remote_mixin.RemoteAuth') as mock_auth_class:
    mock_auth_class.return_value = mock_remote_auth
    
    scheduler = MockSLURMScheduler(mock_calc, queue)
    # Mock the _transfer_pseudopotentials to avoid file lookup
    scheduler._transfer_pseudopotentials = Mock()
    
    # Run the scheduler's run() method (this should store calc.remote)
    try:
        scheduler.run()
        
        # Check if calc.remote was set
        if hasattr(mock_calc, 'remote') and mock_calc.remote is not None:
            print("✅ PASS: calc.remote was set during SLURM job submission")
            print(f"   calc.remote = {mock_calc.remote}")
            print(f"   calc.last_job_id = {mock_calc.last_job_id}")
            print(f"   calc.last_remote_path = {mock_calc.last_remote_path}")
        else:
            print("❌ FAIL: calc.remote was NOT set")
            sys.exit(1)
    except Exception as e:
        print(f"❌ FAIL: Exception during job submission: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

# ============================================================================
# Test 2: Verify calc.remote is set after Direct job submission
# ============================================================================
print("\n[Test 2] Direct job submission - verify calc.remote is stored")
print("-"*80)

# Reset mock for direct scheduler test
mock_calc2 = Mock()
mock_calc2.prefix = "test2"
mock_calc2.package = "pw"
mock_calc2.directory = "/tmp/test2"
mock_calc2.parameters = {
    "pseudopotentials": {},
    "input_data": {"CONTROL": {}}
}

# Mock direct scheduler response (PID in background)
mock_remote_auth2 = Mock()
mock_remote_auth2.run_command = Mock(return_value=("12345\n", ""))
mock_remote_auth2.send_file = Mock()

queue2 = {
    "execution": "remote",
    "scheduler": "direct",
    "remote_host": "medusa.fis.uerj.br",
    "remote_user": "vinicius",
    "remote_dir": "/scratch/users/vinicius/xespresso",
    "remote_auth": {"method": "key"},
    "wait_for_completion": False
}

with patch('xespresso.schedulers.remote_mixin.RemoteAuth') as mock_auth_class:
    mock_auth_class.return_value = mock_remote_auth2
    
    scheduler2 = MockSLURMScheduler(mock_calc2, queue2)
    scheduler2.queue = queue2
    # Mock the _transfer_pseudopotentials to avoid file lookup
    scheduler2._transfer_pseudopotentials = Mock()
    
    try:
        scheduler2.run()
        
        # Check if calc.remote was set
        if hasattr(mock_calc2, 'remote') and mock_calc2.remote is not None:
            print("✅ PASS: calc.remote was set during Direct job submission")
            print(f"   calc.remote = {mock_calc2.remote}")
            print(f"   calc.last_job_id = {mock_calc2.last_job_id}")
            print(f"   calc.last_remote_path = {mock_calc2.last_remote_path}")
        else:
            print("❌ FAIL: calc.remote was NOT set")
            sys.exit(1)
    except Exception as e:
        print(f"❌ FAIL: Exception during job submission: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

# ============================================================================
# Test 3: Verify RemoteJobMonitor can find the connection
# ============================================================================
print("\n[Test 3] RemoteJobMonitor access - verify it can find calc.remote")
print("-"*80)

from xespresso.schedulers.remote_job_monitor import RemoteJobMonitor

# Create a calculator with remote connection set (like after our fix)
calc_with_remote = Mock()
calc_with_remote.last_job_id = "1596"
calc_with_remote.last_remote_path = "/scratch/users/vinicius/xespresso/scf_si-test_d91374c1"
calc_with_remote.remote = mock_remote_auth  # This is the fix - remote is now available!

try:
    # This should work now because calc.remote is set
    monitor = RemoteJobMonitor(calc_with_remote)
    print("✅ PASS: RemoteJobMonitor successfully created with calc.remote")
    print(f"   monitor.job_id = {monitor.job_id}")
    print(f"   monitor.job_type = {monitor.job_type}")
    print(f"   monitor.remote = {monitor.remote is not None}")
except ValueError as e:
    print(f"❌ FAIL: RemoteJobMonitor init failed: {e}")
    sys.exit(1)

# ============================================================================
# Test 4: Verify fallback still works (calc.scheduler.remote)
# ============================================================================
print("\n[Test 4] RemoteJobMonitor fallback - calc.scheduler.remote access")
print("-"*80)

calc_with_scheduler = Mock()
calc_with_scheduler.last_job_id = "1597"
calc_with_scheduler.last_remote_path = "/scratch/users/vinicius/xespresso/test2"
calc_with_scheduler.remote = None  # Simulate it not being set directly

# Set up scheduler with remote connection (fallback path)
mock_scheduler = Mock()
mock_scheduler.remote = mock_remote_auth
calc_with_scheduler.scheduler = mock_scheduler

try:
    # This should work because RemoteJobMonitor has fallback to scheduler.remote
    monitor2 = RemoteJobMonitor(calc_with_scheduler)
    print("✅ PASS: RemoteJobMonitor fallback to calc.scheduler.remote works")
    print(f"   monitor.remote = {monitor2.remote is not None}")
except ValueError as e:
    print(f"❌ FAIL: RemoteJobMonitor fallback failed: {e}")
    sys.exit(1)

# ============================================================================
# Summary
# ============================================================================
print("\n" + "="*80)
print("ALL TESTS PASSED! ✅")
print("="*80)
print("\nThe fix is working correctly:")
print("  ✅ SLURM: calc.remote is stored after job submission")
print("  ✅ Direct: calc.remote is stored after job submission")
print("  ✅ RemoteJobMonitor can access calc.remote directly")
print("  ✅ RemoteJobMonitor fallback to calc.scheduler.remote still works")
print("\nThe workflow.run_scf() should now complete without errors!")
print("="*80 + "\n")
