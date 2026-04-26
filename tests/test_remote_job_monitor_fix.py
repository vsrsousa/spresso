#!/usr/bin/env python
"""
Test to verify RemoteJobMonitor correctly detects completed jobs.

Simulates the scenario where jobs complete quickly but squeue no longer 
shows them (need to check sacct instead).
"""

import time
from unittest.mock import Mock, patch
from xespresso.workflow.remote_job_monitor import RemoteJobMonitor


def test_job_completed_not_in_squeue():
    """
    Test the scenario where job completed and doesn't appear in squeue.
    
    This was the BUG:
    - squeue returned nothing (job already done)
    - old _get_job_status() returned 'RUNNING' (infinite loop)
    - New version checks sacct for final status
    """
    
    print("=" * 80)
    print("TEST: Job completed but not in squeue (must check sacct)")
    print("=" * 80)
    
    # Create mock calculator
    mock_calc = Mock()
    mock_calc.last_job_id = '4428'
    mock_calc.last_remote_path = '/scratch/users/vinicius/xespresso/111_nlayers_3_d30b8eb1'
    mock_calc.directory = '/home/vinicius/data/spresso/au-fcc/relax/111/nlayers_3'
    
    # Create mock remote connection
    mock_remote = Mock()
    mock_calc.remote = mock_remote
    
    # Simulate the actual behavior:
    # 1. First poll: squeue returns empty (job already done)
    # 2. sacct shows COMPLETED
    call_count = [0]
    
    def mock_run_command(cmd):
        call_count[0] += 1
        
        print(f"\n📝 Call {call_count[0]}: {cmd}")
        
        # First call: squeue returns nothing (job not in queue)
        if 'squeue' in cmd and call_count[0] == 1:
            print("  → squeue: [empty, job already done]")
            return ('', '')  # Empty: job not in queue
        
        # Second call: sacct shows completed
        elif 'sacct' in cmd and call_count[0] == 2:
            print("  → sacct: [COMPLETED]")
            return ('COMPLETED\n', '')  # Job completed
        
        return ('', '')
    
    mock_remote.run_command = mock_run_command
    
    print("\n🔧 Monitor setup:")
    print(f"  Job ID: {mock_calc.last_job_id}")
    print(f"  Remote path: {mock_calc.last_remote_path}")
    
    # Create monitor
    monitor = RemoteJobMonitor(mock_calc)
    
    print("\n⏱️  Waiting for job completion (timeout=10s, poll=2s)...")
    
    # Wait for completion (should succeed quickly)
    result = monitor.wait(timeout=10, poll_interval=2)
    
    print(f"\n✅ Result: {'COMPLETED ✓' if result else 'TIMEOUT ✗'}")
    print(f"📊 Total squeue/sacct calls: {call_count[0]}")
    
    assert result == True, "Job should be detected as completed!"
    assert call_count[0] >= 2, f"Should have called squeue and sacct, got {call_count[0]} calls"
    
    print("\n✅ TEST PASSED: Monitor correctly detected completed job via sacct!")


def test_job_still_running():
    """Test that running jobs are detected correctly."""
    
    print("\n" + "=" * 80)
    print("TEST: Job still running (detected via squeue)")
    print("=" * 80)
    
    mock_calc = Mock()
    mock_calc.last_job_id = '4429'
    mock_calc.directory = '/scratch/test'
    
    mock_remote = Mock()
    mock_calc.remote = mock_remote
    
    call_count = [0]
    
    def mock_run_command(cmd):
        call_count[0] += 1
        
        if call_count[0] <= 3:
            # First 3 calls: job is RUNNING
            if 'squeue' in cmd:
                print(f"  Call {call_count[0]}: squeue → [RUNNING]")
                return ('RUNNING\n', '')
        else:
            # 4th call: job finished
            if 'squeue' in cmd:
                print(f"  Call {call_count[0]}: squeue → [empty, checking sacct]")
                return ('', '')
            elif 'sacct' in cmd:
                print(f"  Call {call_count[0]}: sacct → [COMPLETED]")
                return ('COMPLETED\n', '')
        
        return ('', '')
    
    mock_remote.run_command = mock_run_command
    
    monitor = RemoteJobMonitor(mock_calc)
    
    print("\n⏱️  Waiting for job completion (timeout=15s, poll=1s)...")
    result = monitor.wait(timeout=15, poll_interval=1)
    
    print(f"\n✅ Result: {'COMPLETED ✓' if result else 'TIMEOUT ✗'}")
    print(f"📊 Total calls: {call_count[0]}")
    
    assert result == True, "Job should complete!"
    assert call_count[0] >= 4, "Should have multiple polls"
    
    print("\n✅ TEST PASSED: Monitor correctly handled running → completed!")


def test_job_failed():
    """Test that failed jobs are detected correctly."""
    
    print("\n" + "=" * 80)
    print("TEST: Job failed (detected via sacct)")
    print("=" * 80)
    
    mock_calc = Mock()
    mock_calc.last_job_id = '4430'
    mock_calc.directory = '/scratch/test'
    
    mock_remote = Mock()
    mock_calc.remote = mock_remote
    
    def mock_run_command(cmd):
        if 'squeue' in cmd:
            print(f"  squeue → [empty, job finished]")
            return ('', '')
        elif 'sacct' in cmd:
            print(f"  sacct → [FAILED]")
            return ('FAILED\n', '')
        return ('', '')
    
    mock_remote.run_command = mock_run_command
    
    monitor = RemoteJobMonitor(mock_calc)
    
    print("\n⏱️  Waiting for job completion...")
    result = monitor.wait(timeout=5, poll_interval=1)
    
    print(f"\n✅ Result: {'FAILED (correctly detected) ✓' if result == False else 'Unexpected ✗'}")
    
    assert result == False, "Failed job should return False!"
    
    print("\n✅ TEST PASSED: Monitor correctly detected failed job!")


def test_timeout_behavior():
    """Test that timeout works correctly."""
    
    print("\n" + "=" * 80)
    print("TEST: Job timeout (never completes)")
    print("=" * 80)
    
    mock_calc = Mock()
    mock_calc.last_job_id = '4431'
    mock_calc.directory = '/scratch/test'
    
    mock_remote = Mock()
    mock_calc.remote = mock_remote
    
    # Always report RUNNING
    mock_remote.run_command = Mock(return_value=('RUNNING\n', ''))
    
    monitor = RemoteJobMonitor(mock_calc)
    
    print("\n⏱️  Waiting with 3s timeout...")
    start = time.time()
    result = monitor.wait(timeout=3, poll_interval=1)
    elapsed = time.time() - start
    
    print(f"\n✅ Result: {'TIMEOUT ✓' if result == False else 'Unexpected ✗'}")
    print(f"⏱️  Elapsed: {elapsed:.1f}s")
    
    assert result == False, "Should timeout!"
    assert elapsed >= 3, f"Should wait at least 3s, only waited {elapsed:.1f}s"
    
    print("\n✅ TEST PASSED: Monitor correctly handled timeout!")


if __name__ == "__main__":
    test_job_completed_not_in_squeue()
    test_job_still_running()
    test_job_failed()
    test_timeout_behavior()
    
    print("\n" + "=" * 80)
    print("✅ ALL TESTS PASSED!")
    print("=" * 80)
    print("""
SUMMARY OF FIX:
==============

Problem: Monitor loops infinitely when jobs complete quickly
Cause: _get_job_status() only checked squeue, returned 'RUNNING' fallback
Solution: Check squeue first, then sacct if job not found

New Behavior:
1. squeue shows RUNNING/PENDING → Return state
2. squeue empty → Check sacct for final state
3. sacct shows COMPLETED → Return 'COMPLETED'
4. sacct shows FAILED/etc → Return 'FAILED'
5. Loop exits when status is not RUNNING

Result: Parallel jobs complete properly without infinite loops!
""")
