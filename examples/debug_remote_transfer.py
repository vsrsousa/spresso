#!/usr/bin/env python
"""
Debug script to test remote file transfer in batch mode.

Usage:
    python examples/debug_remote_transfer.py

This script tests:
1. Machine configuration loading
2. Queue setup and remote connection
3. Input file generation
4. Remote file transfer
5. Batch job submission
"""

import os
import sys
from ase.build import bulk

# Enable verbose error reporting to see full stack traces
os.environ['XESPRESSO_VERBOSE_ERRORS'] = '1'

from xespresso.workflow.convergence_workflow import ConvergenceWorkflow
from xespresso.machines import load_machine


def test_machine_loading():
    """Test if machine configuration is loaded correctly."""
    print("\n" + "="*80)
    print("TEST 1: Machine Configuration Loading")
    print("="*80)
    
    try:
        machine = load_machine('medusa')
        print(f"✓ Machine loaded: {machine}")
        print(f"  - execution: {machine.get('execution')}")
        print(f"  - scheduler: {machine.get('scheduler')}")
        print(f"  - remote_host: {machine.get('remote_host')}")
        print(f"  - remote_user: {machine.get('remote_user')}")
        return machine
    except Exception as e:
        print(f"✗ Failed to load machine: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_convergence_workflow_setup():
    """Test ConvergenceWorkflow initialization with machine."""
    print("\n" + "="*80)
    print("TEST 2: ConvergenceWorkflow Initialization")
    print("="*80)
    
    try:
        atoms = bulk('Si', 'diamond', a=5.43)
        
        wf = ConvergenceWorkflow(
            atoms=atoms,
            pseudopotentials_config='SSSP_efficiency',
            precision='low',
            machine='medusa',
            code_version='7.4.1',
        )
        
        print(f"✓ ConvergenceWorkflow created")
        print(f"  - atoms: {wf.atoms.get_chemical_formula()}")
        print(f"  - machine: {wf.machine}")
        print(f"  - queue loaded: {wf.queue is not None}")
        
        if wf.queue:
            print(f"  - queue execution: {wf.queue.get('execution')}")
            print(f"  - queue scheduler: {wf.queue.get('scheduler')}")
        
        return wf
    except Exception as e:
        print(f"✗ Failed to create ConvergenceWorkflow: {e}")
        import traceback
        traceback.print_exc()
        return None


def test_batch_mode_detection(wf):
    """Test batch mode detection."""
    print("\n" + "="*80)
    print("TEST 3: Batch Mode Detection")
    print("="*80)
    
    if not wf:
        print("✗ ConvergenceWorkflow not available")
        return False
    
    try:
        is_remote = wf.queue and wf.queue.get('execution') == 'remote'
        print(f"✓ Batch mode detection:")
        print(f"  - is_remote: {is_remote}")
        print(f"  - will_use_batch: {is_remote and True}")
        
        return is_remote
    except Exception as e:
        print(f"✗ Batch mode detection failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_input_generation(wf):
    """Test input file generation."""
    print("\n" + "="*80)
    print("TEST 4: Input File Generation")
    print("="*80)
    
    if not wf:
        print("✗ ConvergenceWorkflow not available")
        return False
    
    try:
        from xespresso.workflow.calculation_workflow import CalculationWorkflow
        
        # Create CalculationWorkflow
        calc_wf = CalculationWorkflow(
            atoms=wf.atoms,
            pseudopotentials=wf.pseudopotentials,
            protocol=wf.protocol,
            queue=wf.queue,
        )
        
        print(f"✓ CalculationWorkflow created")
        print(f"  - atoms: {calc_wf.atoms.get_chemical_formula()}")
        print(f"  - queue: {calc_wf.queue is not None}")
        
        # Try to create an Espresso calculator (this writes input files)
        from xespresso import Espresso
        
        params = {
            'pseudopotentials': calc_wf.pseudopotentials,
            'label': 'test_scf',
            'calculation': 'scf',
            'input_data': calc_wf.input_data.copy(),
            'kpts': calc_wf._get_kpts(),
            'queue': calc_wf.queue,
        }
        
        calc = Espresso(**params)
        calc.write_input(calc_wf.atoms)
        
        print(f"✓ Input files generated in: {calc.directory}")
        
        # Check what files were created
        input_file = f"{calc.directory}/{calc.prefix}.{calc.package}i"
        if os.path.exists(input_file):
            print(f"  - Input file exists: {input_file}")
            with open(input_file, 'r') as f:
                lines = f.readlines()
                print(f"  - File size: {len(lines)} lines")
        else:
            print(f"  ✗ Input file not found: {input_file}")
        
        return True, calc, params
    except Exception as e:
        print(f"✗ Input generation failed: {e}")
        import traceback
        traceback.print_exc()
        return False, None, None


def test_remote_connection(wf):
    """Test remote connection setup."""
    print("\n" + "="*80)
    print("TEST 5: Remote Connection Setup")
    print("="*80)
    
    if not wf or not wf.queue:
        print("✗ ConvergenceWorkflow or queue not available")
        return False
    
    try:
        from xespresso.utils.auth import RemoteAuth
        
        remote = RemoteAuth(
            username=wf.queue['remote_user'],
            host=wf.queue['remote_host'],
            auth_config=wf.queue.get('remote_auth', {}),
        )
        
        remote.connect()
        
        print(f"✓ Remote connection established")
        print(f"  - host: {wf.queue['remote_host']}")
        print(f"  - user: {wf.queue['remote_user']}")
        
        # Test a simple command
        stdout, stderr = remote.run_command('pwd')
        print(f"  - remote working dir: {stdout.strip()}")
        
        return True, remote
    except Exception as e:
        print(f"✗ Remote connection failed: {e}")
        print(f"\n⚠️  This is likely the cause of file transfer failures!")
        print(f"    Check:")
        print(f"    1. SSH key is properly configured")
        print(f"    2. Remote host is reachable")
        print(f"    3. Remote user credentials are correct")
        import traceback
        traceback.print_exc()
        return False, None


def test_scheduler_execution(calc, wf):
    """Test scheduler setup and execution (without actual job submission)."""
    print("\n" + "="*80)
    print("TEST 6: Scheduler Setup")
    print("="*80)
    
    if not calc or not wf:
        print("✗ Calculator or ConvergenceWorkflow not available")
        return False
    
    try:
        from xespresso.scheduler import set_queue
        
        # Set up the scheduler (this initializes remote connection)
        set_queue(calc)
        
        print(f"✓ Scheduler configured")
        print(f"  - scheduler type: {calc.queue.get('scheduler')}")
        print(f"  - execution mode: {calc.queue.get('execution')}")
        print(f"  - job_file: {calc.scheduler.job_file}")
        
        return True
    except Exception as e:
        print(f"✗ Scheduler setup failed: {e}")
        print(f"\n⚠️  This is where file transfer would be initialized!")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all debug tests."""
    print("\n" + "="*80)
    print("REMOTE FILE TRANSFER DEBUG TEST")
    print("="*80)
    print(f"Python: {sys.executable}")
    print(f"PWD: {os.getcwd()}")
    
    # Test 1: Machine loading
    machine = test_machine_loading()
    if not machine:
        print("\n⚠️  Cannot continue without machine configuration")
        return
    
    # Test 2: ConvergenceWorkflow setup
    wf = test_convergence_workflow_setup()
    if not wf:
        print("\n⚠️  Cannot continue without ConvergenceWorkflow")
        return
    
    # Test 3: Batch mode detection
    is_remote = test_batch_mode_detection(wf)
    
    # Test 4: Input generation
    success, calc, params = test_input_generation(wf)
    if not success:
        print("\n⚠️  Cannot continue without generated inputs")
        return
    
    # Test 5: Remote connection
    success, remote = test_remote_connection(wf)
    if not success:
        print("\n⚠️  Remote connection failed - this would prevent file transfer")
        return
    
    # Test 6: Scheduler execution
    test_scheduler_execution(calc, wf)
    
    print("\n" + "="*80)
    print("DEBUG TEST COMPLETE")
    print("="*80)
    print("\nIf all tests passed, try running:")
    print("  ConvergenceWorkflow.optimize_parameters(atoms, machine='medusa', use_batch_mode=True)")


if __name__ == '__main__':
    main()
