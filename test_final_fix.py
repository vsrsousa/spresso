"""
Test the final fix: storing remote connection in calc.remote during job submission
"""
import os
os.environ['ASE_ESPRESSO_PSEUDO'] = '/home/vinicius/pseudos/SSSP_1.3.0_PBE_efficiency'

import sys
import logging

# Setup logging to see what's happening
logging.basicConfig(level=logging.INFO, format='%(levelname)-8s [%(name)s]: %(message)s')

from ase.build import bulk
from xespresso import Espresso
from xespresso.workflow.simple_workflow import CalculationWorkflow

print("\n" + "="*70)
print("TEST: Final Fix - Store remote connection in calc.remote")
print("="*70)

# Create a simple Silicon structure
atoms = bulk('Si', 'diamond', a=5.431)

# Create workflow with remote non-blocking execution
workflow = CalculationWorkflow(
    atoms,
    pseudopotentials={
        'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF'
    },
    machine='medusa',  # Load machine config
    # The non-blocking mode is already set in the machine config
)

print("\n1. Testing run_scf with remote non-blocking execution...")
print("-" * 70)

try:
    # This should now work with the fix:
    # 1. execute() sets calc.remote = self.remote in remote_mixin.py
    # 2. RemoteJobMonitor can then find calc.remote
    workflow.run_scf(label='scf/si-final-test')
    print("\n✅ SUCCESS: run_scf completed without errors!")
    print(f"   Last job ID: {workflow.last_calc.last_job_id}")
    print(f"   Remote path: {workflow.last_calc.last_remote_path}")
    print(f"   Remote connection exists: {hasattr(workflow.last_calc, 'remote') and workflow.last_calc.remote is not None}")
    
except ValueError as e:
    print(f"\n❌ FAILED with ValueError: {e}")
    sys.exit(1)
except AttributeError as e:
    print(f"\n❌ FAILED with AttributeError: {e}")
    sys.exit(1)
except Exception as e:
    print(f"\n❌ FAILED with unexpected error: {type(e).__name__}: {e}")
    sys.exit(1)

print("\n" + "="*70)
print("All tests passed! ✅")
print("="*70)
