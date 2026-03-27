"""
Mock test: Complete workflow chain with Hubbard - does it appear correctly in input file?
"""

from ase.build import bulk
import numpy as np
from xespresso.workflow.convergence_workflow import ConvergenceWorkflow
import tempfile
import os
from pathlib import Path

print("="*80)
print("CONVERGENCE WORKFLOW - HUBBARD PARAMETER FLOW TEST")
print("="*80)

# Create Fe structure
atoms = bulk("Fe", cubic=True)

print("\n" + "="*80)
print("TEST 1: Old Format Flow (QE 6.8)")
print("="*80)

try:
    workflow_old = ConvergenceWorkflow(
        atoms=atoms,
        code_version='6.8',  # Old QE format
        pseudopotentials={'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'},
        precision='low',
        initial_kspacing=0.3,
        hubbard_config={'Fe': 4.3},  # Old format for QE 6.8
    )
    
    print("✅ ConvergenceWorkflow created")
    print(f"   code_version: 6.8")
    print(f"   hubbard_config: {{'Fe': 4.3}}")
    
    # Try phase 1 with mocked calculation
    # Check what's in the first calculation's input_data
    if hasattr(workflow_old, 'phase1_config'):
        print(f"\n   Phase 1 Config:")
        phase1_cfg = workflow_old.phase1_config
        print(f"   - input_data['qe_version']: {phase1_cfg.get('input_data', {}).get('qe_version', 'NOT SET')}")
        print(f"   - input_data['lda_plus_u']: {phase1_cfg.get('input_data', {}).get('lda_plus_u', 'NOT SET')}")
        if 'input_ntyp' in phase1_cfg.get('input_data', {}):
            hubbard_u = phase1_cfg['input_data']['input_ntyp'].get('Hubbard_U', {})
            print(f"   - input_data['input_ntyp']['Hubbard_U']: {hubbard_u}")
    
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*80)
print("TEST 2: New Format Flow (QE 7.2)")
print("="*80)

try:
    workflow_new = ConvergenceWorkflow(
        atoms=atoms,
        code_version='7.2',  # New QE format
        pseudopotentials={'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'},
        precision='low',
        initial_kspacing=0.3,
        hubbard_config={'u': {'Fe': 4.3}},  # New format for QE 7.2
    )
    
    print("✅ ConvergenceWorkflow created")
    print(f"   code_version: 7.2")
    print(f"   hubbard_config: {{'u': {{'Fe': 4.3}}}}")
    
    # Check phase 1 config
    if hasattr(workflow_new, 'phase1_config'):
        print(f"\n   Phase 1 Config:")
        phase1_cfg = workflow_new.phase1_config
        print(f"   - input_data['qe_version']: {phase1_cfg.get('input_data', {}).get('qe_version', 'NOT SET')}")
        print(f"   - input_data['lda_plus_u']: {phase1_cfg.get('input_data', {}).get('lda_plus_u', 'NOT SET')}")
        if 'hubbard' in phase1_cfg.get('input_data', {}):
            print(f"   - input_data['hubbard']: {phase1_cfg['input_data']['hubbard']}")
        else:
            print(f"   - input_data['hubbard']: NOT SET ❌")
    
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*80)
print("TEST 3: Check CalculationWorkflow receives hubbard_config correctly")
print("="*80)

try:
    from xespresso.workflow.calculation_workflow import CalculationWorkflow
    
    # Simulate what happened inside ConvergenceWorkflow
    calc_old = CalculationWorkflow(
        atoms=atoms,
        code_version='6.8',
        pseudopotentials={'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'},
        precision='low',
        hubbard_config={'Fe': 4.3},  # Old format
    )
    
    print("✅ CalculationWorkflow (QE 6.8) created directly")
    print(f"   lda_plus_u: {calc_old.input_data.get('lda_plus_u')}")
    print(f"   Hubbard_U: {calc_old.input_data.get('input_ntyp', {}).get('Hubbard_U', {})}")
    
    # New format
    calc_new = CalculationWorkflow(
        atoms=atoms,
        code_version='7.2',
        pseudopotentials={'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'},
        precision='low',
        hubbard_config={'u': {'Fe': 4.3}},  # New format
    )
    
    print("\n✅ CalculationWorkflow (QE 7.2) created directly")
    print(f"   qe_version: {calc_new.input_data.get('qe_version')}")
    print(f"   hubbard: {calc_new.input_data.get('hubbard')}")
    print(f"   lda_plus_u: {calc_new.input_data.get('lda_plus_u', 'NOT SET')}")
    
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print("""
The flow should be:

1. User calls ConvergenceWorkflow(..., code_version='7.2', hubbard_config={'u': {'Fe': 4.3}})

2. ConvergenceWorkflow stores these and creates CalculationWorkflow instances with:
   - code_version='7.2'
   - hubbard_config={'u': {'Fe': 4.3}}

3. CalculationWorkflow.__init__():
   - Sets input_data['qe_version']='7.2' from code_version
   - Detects hubbard_config has 'u' key → NEW FORMAT
   - Sets input_data['hubbard']={'u': {'Fe': 4.3}}
   - Does NOT set lda_plus_u (new format doesn't need it)

4. When Hubbard class reads input_data:
   - Sees qe_version='7.2' (>= 7.0) → NEW FORMAT
   - Sees input_data['hubbard'] present → Generates HUBBARD card
   
Expected in final input file:
   HUBBARD {atomic}
     U Fe 4.3
""")
