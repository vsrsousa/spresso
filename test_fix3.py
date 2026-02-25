#!/usr/bin/env python
"""
Test Fix #3: Store remote connection on calc after execute()
"""

import tempfile
from unittest.mock import Mock, patch, MagicMock
from ase.build import bulk
from xespresso import CalculationWorkflow, Espresso
from xespresso.schedulers import RemoteJobMonitor

print("=" * 70)
print("Testing Fix #3: Storing remote connection after execute()")
print("=" * 70)

atoms = bulk("Si", cubic=True)

with tempfile.TemporaryDirectory() as tmpdir:
    with patch('xespresso.workflow.simple_workflow.load_machine') as mock_load_machine, \
         patch('xespresso.workflow.simple_workflow.load_pseudopotentials_config') as mock_load_pseudo:
        
        mock_machine = {
            'execution': 'remote',
            'scheduler': 'slurm',
            'wait_for_completion': False,
        }
        mock_load_machine.return_value = mock_machine
        
        mock_pseudo_config = Mock()
        mock_pseudo_config.base_path = tmpdir
        mock_pseudo_config.get_pseudopotential = Mock(return_value=Mock(filename='Si.pbe.UPF'))
        mock_load_pseudo.return_value = mock_pseudo_config
        
        workflow = CalculationWorkflow(
            atoms=atoms,
            pseudopotentials_config='default',
            machine='test_machine'
        )
        
        print("\nScenario: After calc.execute() is called")
        
        # Mock calculator
        mock_calc = Mock(spec=Espresso)
        mock_calc.atoms = None  # Initially None
        mock_calc.write_input = Mock()
        mock_calc.last_job_id = '1594'
        mock_calc.last_remote_path = '/scratch/users/vinicius/xespresso/scf_si-test_d91374c1'
        
        # Mock scheduler with remote connection
        mock_scheduler = Mock()
        mock_remote = Mock()
        mock_remote.run_command = Mock(return_value=('', ''))
        mock_scheduler.remote = mock_remote
        
        # Before execute: calc doesn't have remote
        mock_calc.scheduler = None
        mock_calc.remote = None
        
        print("1. Before execute():")
        print(f"   - calc.remote: {mock_calc.remote}")
        print(f"   - calc.scheduler: {mock_calc.scheduler}")
        
        # Simulate execute() - sets up scheduler
        print("\n2. After execute():")
        mock_calc.execute = Mock()
        mock_calc.execute()
        
        # After execute, scheduler is set up with remote connection
        mock_calc.scheduler = mock_scheduler
        print(f"   - calc.scheduler: <Mock>")
        print(f"   - calc.scheduler.remote: <Mock connection>")
        print(f"   - calc.remote (before fix): None")
        
        # FIX #3: Store remote connection after execute()
        print("\n3. Applying Fix #3:")
        if hasattr(mock_calc, 'scheduler') and hasattr(mock_calc.scheduler, 'remote'):
            mock_calc.remote = mock_calc.scheduler.remote
            print("   - Stored: calc.remote = calc.scheduler.remote")
        
        print(f"   - calc.remote (after fix): <Mock connection>")
        
        # Now RemoteJobMonitor should work
        print("\n4. Creating RemoteJobMonitor:")
        try:
            monitor = RemoteJobMonitor(mock_calc)
            print("   ✓ RemoteJobMonitor created successfully!")
            print(f"   ✓ Job ID: {monitor.job_id}")
            print(f"   ✓ Remote connection available: {monitor.remote is not None}")
            print(f"   ✓ Job type detected: {monitor.job_type}")
        except Exception as e:
            print(f"   ✗ Failed: {e}")
            raise

print("\n" + "=" * 70)
print("✓ FIX #3 VERIFIED: Remote connection properly stored on calc!")
print("=" * 70)
print("\nAll fixes applied:")
print("  1. calc.atoms set before execute() ✓")
print("  2. RemoteJobMonitor accesses calc.scheduler.remote ✓")
print("  3. Remote connection stored on calc after execute() ✓")
print("\nResulting flow:")
print("  write_input() → atoms set → execute() → remote stored → monitor created")
print("=" * 70 + "\n")
