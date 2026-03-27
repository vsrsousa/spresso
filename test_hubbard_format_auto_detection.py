"""
Test: Auto-detection of Hubbard format based on code_version
This test shows how code_version is automatically converted to qe_version
and affects the Hubbard parameter format selection.
"""

import sys
from unittest.mock import patch, MagicMock
from ase.build import bulk
import numpy as np

print("="*80)
print("TEST: Auto-detection of Hubbard Format Based on code_version")
print("="*80)

# Mock CalculationWorkflow and ConvergenceWorkflow initialization to avoid actual calculations
with patch('xespresso.workflow.calculation_workflow.logger') as mock_calc_logger, \
     patch('xespresso.workflow.convergence_workflow.logger') as mock_conv_logger:
    
    from xespresso.workflow.calculation_workflow import CalculationWorkflow
    from xespresso.workflow.convergence_workflow import ConvergenceWorkflow
    
    atoms = bulk('Fe', cubic=True)
    
    print("\n" + "="*80)
    print("TEST 1: CalculationWorkflow with code_version='7.2'")
    print("="*80)
    
    try:
        # Create CalculationWorkflow with QE 7.2
        calc_wf_72 = CalculationWorkflow(
            atoms=atoms,
            pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},
            code_version='7.2',  # ← Should set qe_version='7.2'
        )
        
        print(f"✓ CalculationWorkflow initialized with code_version='7.2'")
        print(f"  input_data['qe_version'] = {calc_wf_72.input_data.get('qe_version')}")
        
        # Check if qe_version was auto-set
        if calc_wf_72.input_data.get('qe_version') == '7.2':
            print("  ✅ PASS: qe_version auto-set to '7.2'")
        else:
            print("  ❌ FAIL: qe_version was not set correctly")
    
    except Exception as e:
        print(f"  Error: {e}")
    
    print("\n" + "="*80)
    print("TEST 2: CalculationWorkflow with code_version='6.8'")
    print("="*80)
    
    try:
        # Create CalculationWorkflow with QE 6.8
        calc_wf_68 = CalculationWorkflow(
            atoms=atoms,
            pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},
            code_version='6.8',  # ← Should set qe_version='6.8'
        )
        
        print(f"✓ CalculationWorkflow initialized with code_version='6.8'")
        print(f"  input_data['qe_version'] = {calc_wf_68.input_data.get('qe_version')}")
        
        if calc_wf_68.input_data.get('qe_version') == '6.8':
            print("  ✅ PASS: qe_version auto-set to '6.8'")
        else:
            print("  ❌ FAIL: qe_version was not set correctly")
    
    except Exception as e:
        print(f"  Error: {e}")
    
    print("\n" + "="*80)
    print("TEST 3: Hubbard Format Selection Based on qe_version")
    print("="*80)
    
    from xespresso.hubbard import HubbardConfig
    
    # Test 3a: QE 7.2 should use new format
    print("\nTEST 3a: QE 7.2 - Should use NEW format (HUBBARD card)")
    input_data_72 = {
        'qe_version': '7.2',
        'hubbard': {
            'u': {'Fe-3d': 4.3}
        }
    }
    
    config_72 = HubbardConfig.from_input_data(input_data_72, qe_version='7.2')
    uses_new_72 = config_72.should_use_new_format()
    
    print(f"  HubbardConfig.should_use_new_format() = {uses_new_72}")
    if uses_new_72:
        print("  ✅ PASS: Correctly detected new format for QE 7.2")
        card_str = config_72.to_new_format_card()
        print(f"  Generated card:\n    {chr(10).join(card_str)}")
    else:
        print("  ❌ FAIL: Should use new format for QE 7.2")
    
    # Test 3b: QE 6.8 should use old format
    print("\nTEST 3b: QE 6.8 - Should use OLD format (SYSTEM namelist)")
    input_data_68 = {
        'qe_version': '6.8',
        'input_ntyp': {
            'Hubbard_U': {'Fe': 4.3}
        }
    }
    
    config_68 = HubbardConfig.from_input_data(input_data_68, qe_version='6.8')
    uses_new_68 = config_68.should_use_new_format()
    
    print(f"  HubbardConfig.should_use_new_format() = {uses_new_68}")
    if not uses_new_68:
        print("  ✅ PASS: Correctly detected old format for QE 6.8")
        print(f"  Old format data: {config_68.u_params}")
    else:
        print("  ❌ FAIL: Should use old format for QE 6.8")
    
    print("\n" + "="*80)
    print("TEST 4: Complete Workflow with code_version")
    print("="*80)
    
    try:
        # Test with explicit input_data that has hubbard config
        calc_with_hubbard = CalculationWorkflow(
            atoms=atoms,
            pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},
            code_version='7.2',
            input_data={
                'hubbard': {
                    'u': {'Fe-3d': 4.3}
                }
            }
        )
        
        print(f"✓ CalculationWorkflow with hubbard_config and code_version='7.2'")
        print(f"  qe_version in input_data: {calc_with_hubbard.input_data.get('qe_version')}")
        print(f"  hubbard in input_data: {calc_with_hubbard.input_data.get('hubbard')}")
        
        # Simulate building the Hubbard string
        config_from_workflow = HubbardConfig.from_input_data(
            calc_with_hubbard.input_data,
            qe_version=calc_with_hubbard.input_data.get('qe_version')
        )
        
        uses_new_workflow = config_from_workflow.should_use_new_format()
        print(f"  Would use new format: {uses_new_workflow}")
        
        if uses_new_workflow:
            print("  ✅ PASS: Workflow correctly uses new format due to qe_version")
        else:
            print("  ❌ FAIL: Should use new format")
            
    except Exception as e:
        print(f"  Error: {e}")

print("\n" + "="*80)
print("SUMMARY: code_version → qe_version → Hubbard Format Selection")
print("="*80)
print("""
Flow:
  code_version='7.2'
       ↓
  CalculationWorkflow.__init__()
       ↓
  Auto-set: input_data['qe_version'] = '7.2'
       ↓
  HubbardConfig.from_input_data(input_data, qe_version='7.2')
       ↓
  Detected QE >= 7.0 → Use NEW format (HUBBARD card)

---

  code_version='6.8'
       ↓
  CalculationWorkflow.__init__()
       ↓
  Auto-set: input_data['qe_version'] = '6.8'
       ↓
  HubbardConfig.from_input_data(input_data, qe_version='6.8')
       ↓
  Detected QE < 7.0 → Use OLD format (SYSTEM namelist)

Result: User just passes code_version and format is automatically selected!
""")
print("="*80)
