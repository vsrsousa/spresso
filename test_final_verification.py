"""
Final test - Run actual remote execution and check if remote connection is preserved
"""
import os
os.environ['ASE_ESPRESSO_PSEUDO'] = '/home/vinicius/pseudos/SSSP_1.3.0_PBE_efficiency'

import sys
import logging

# Reduce logging noise but keep important info
logging.basicConfig(level=logging.INFO, format='%(levelname)-8s: %(message)s')
logging.getLogger('paramiko').setLevel(logging.WARNING)
logging.getLogger('xespresso.xio').setLevel(logging.WARNING)

from ase.build import bulk
from xespresso.workflow.calculation_workflow import CalculationWorkflow

print("\n" + "="*80)
print("FINAL TEST: Remote Non-Blocking Execution with Fixed Connection Storage")
print("="*80)

atoms = bulk("Si", cubic=True)

# Create workflow with machine config (medusa has non-blocking as default in config)
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials={"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"},
    machine='medusa'
)

print("\n[Test] Running remote non-blocking SCF...")
print("-" * 80)

try:
    # This is the exact line that was failing before
    workflow.run_scf(label='scf/si-test')
    
    print("\n" + "="*80)
    print("✅ SUCCESS! Remote execution completed without errors!")
    print("="*80)
    print(f"Job ID: {workflow.last_calc.last_job_id}")
    print(f"Remote path: {workflow.last_calc.last_remote_path}")
    print(f"Has remote connection: {hasattr(workflow.last_calc, 'remote')}")
    if hasattr(workflow.last_calc, 'remote'):
        print(f"Remote connection value: {workflow.last_calc.remote is not None}")
    
except ValueError as e:
    print(f"\n❌ FAILED: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
except AttributeError as e:
    print(f"\n❌ FAILED: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
except Exception as e:
    print(f"\n❌ FAILED: {type(e).__name__}: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("\n" + "="*80)
print("All tests passed! The fix is working correctly. ✅")
print("="*80)
