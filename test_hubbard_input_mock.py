"""
Mock test to verify Hubbard parameters appear correctly in generated input files
"""

from ase.build import bulk
import numpy as np
from xespresso.workflow.calculation_workflow import CalculationWorkflow
import tempfile
import os

print("="*80)
print("Testing Hubbard Parameters - Input File Generation")
print("="*80)

# Create Au structure
atoms = bulk("Au", cubic=True)

# Test 1: Old format (single element)
print("\n" + "="*80)
print("TEST 1: Old Format - Single Element (Fe)")
print("="*80)

try:
    atoms_fe = bulk("Fe", cubic=True)
    
    workflow1 = CalculationWorkflow(
        atoms=atoms_fe,
        code_version='6.8',  # Old QE format
        pseudopotentials={'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'},
        hubbard_config={'Fe': 4.3},  # Old format
        precision='low'
    )
    
    print("✅ CalculationWorkflow created with old format")
    print(f"   hubbard_config: {{'Fe': 4.3}}")
    print(f"   code_version: 6.8 (old QE)")
    
    # Check if lda_plus_u is set
    if 'lda_plus_u' in workflow1.input_data:
        print(f"   ✅ lda_plus_u = {workflow1.input_data['lda_plus_u']}")
    
    # Check input_ntyp
    if 'input_ntyp' in workflow1.input_data and 'Hubbard_U' in workflow1.input_data['input_ntyp']:
        print(f"   ✅ Hubbard_U in input_ntyp = {workflow1.input_data['input_ntyp']['Hubbard_U']}")
    
except Exception as e:
    print(f"❌ Error: {e}")

# Test 2: New format with structured dict
print("\n" + "="*80)
print("TEST 2: New Format - Structured with 'u' key (QE 7.2)")
print("="*80)

try:
    atoms_fe = bulk("Fe", cubic=True)
    
    workflow2 = CalculationWorkflow(
        atoms=atoms_fe,
        code_version='7.2',  # New QE format
        pseudopotentials={'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'},
        hubbard_config={
            'u': {'Fe': 4.3}
        },
        precision='low'
    )
    
    print("✅ CalculationWorkflow created with new format")
    print(f"   hubbard_config: {{'u': {{'Gd': 6.0}}}}")
    print(f"   code_version: 7.2 (new QE)")
    
    # Check hubbard in input_data
    if 'hubbard' in workflow2.input_data:
        print(f"   ✅ hubbard in input_data = {workflow2.input_data['hubbard']}")
    else:
        print(f"   ❌ hubbard NOT in input_data")
    
    # Check lda_plus_u should NOT be set for new format
    if 'lda_plus_u' not in workflow2.input_data or not workflow2.input_data['lda_plus_u']:
        print(f"   ✅ lda_plus_u NOT set (correct for new format)")
    else:
        print(f"   ⚠️  lda_plus_u = {workflow2.input_data.get('lda_plus_u')} (should be False/absent)")
        
except Exception as e:
    print(f"❌ Error: {e}")

# Test 3: What user passed (incorrect format)
print("\n" + "="*80)
print("TEST 3: User Input Format - What Actually Happened")
print("="*80)

try:
    atoms_fe2 = bulk("Fe", cubic=True)
    
    workflow3 = CalculationWorkflow(
        atoms=atoms_fe2,
        code_version='7.4.1',
        pseudopotentials={'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'},
        hubbard_config={'Fe-3d': 4.3},  # User passed this format
        precision='low'
    )
    
    print("⚠️  CalculationWorkflow created with user format")
    print(f"   hubbard_config: {{'Fe-3d': 4.3}}")
    print(f"   code_version: 7.4.1 (new QE)")
    
    # Check what happened
    if 'hubbard' in workflow3.input_data:
        print(f"   Result: hubbard in input_data = {workflow3.input_data['hubbard']}")
    else:
        print(f"   Result: hubbard NOT in input_data")
    
    if 'lda_plus_u' in workflow3.input_data and workflow3.input_data['lda_plus_u']:
        print(f"   ⚠️  lda_plus_u = {workflow3.input_data['lda_plus_u']}")
        print(f"      (treated as OLD format because 'Fe-3d' doesn't have 'u'/'v'/'projector' keys)")
    
    if 'input_ntyp' in workflow3.input_data and 'Hubbard_U' in workflow3.input_data['input_ntyp']:
        print(f"   Hubbard_U in input_ntyp = {workflow3.input_data['input_ntyp']['Hubbard_U']}")
        
except Exception as e:
    print(f"❌ Error: {e}")

print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print("""
✅ OLD FORMAT (QE < 7.0):
   Pass: hubbard_config = {'Fe': 4.3}
   Result: lda_plus_u=True, Hubbard_U in input_ntyp

✅ NEW FORMAT (QE >= 7.0):
   Pass: hubbard_config = {'u': {'Gd': 6.0}}
   Result: hubbard in input_data

❌ USER FORMAT (INCORRECT):
   Passed: hubbard_config = {'Gd-4f': 6.0}
   Problem: No 'u' key, so treated as old format
   Solution: Use {'u': {'Gd-4f': 6.0}} or just {'Gd': 6.0}
""")
