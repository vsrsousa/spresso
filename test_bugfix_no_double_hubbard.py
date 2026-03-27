"""
Test: Verify ONLY new format HUBBARD card appears, not old hubbard_u
"""

from ase.build import bulk
from xespresso.workflow.calculation_workflow import CalculationWorkflow
from xespresso.xio import write_espresso_in
import tempfile
import os

print("="*80)
print("BUG FIX TEST: Should NOT have both hubbard_u AND HUBBARD card")
print("="*80)

atoms = bulk("Fe", cubic=True)

with tempfile.TemporaryDirectory() as tmpdir:
    workflow = CalculationWorkflow(
        atoms=atoms,
        code_version='7.4.1',  # NEW format (v7.x)
        pseudopotentials={'Fe': 'Fe.pbe.UPF'},
        precision='low',
        hubbard_config={'Fe': 4.3},
        outdir=tmpdir,
    )
    
    print("\n✅ Created CalculationWorkflow")
    print(f"   code_version: 7.4.1")
    print(f"   qe_version in input_data: {workflow.input_data.get('qe_version')}")
    
    # Write input file
    input_file = os.path.join(tmpdir, 'test.pwi')
    write_espresso_in(
        input_file,
        atoms,
        input_data=workflow.input_data,
        pseudopotentials=workflow.pseudopotentials,
        kpts=(2, 2, 2)
    )
    
    with open(input_file, 'r') as f:
        content = f.read()
    
    print("\n📄 Checking generated input file:")
    print("-" * 80)
    
    # Check for hubbard_u in SYSTEM namelist (should NOT be there!)
    has_hubbard_u = 'hubbard_u(' in content.lower()
    has_lda_plus_u = 'lda_plus_u' in content.lower()
    has_hubbard_card = 'HUBBARD {' in content
    has_u_parameter = 'U Fe' in content
    
    print(f"   Has 'hubbard_u(': {has_hubbard_u} (should be FALSE ❌)")
    print(f"   Has 'lda_plus_u': {has_lda_plus_u} (should be FALSE ❌)")
    print(f"   Has 'HUBBARD {{': {has_hubbard_card} (should be TRUE ✅)")
    print(f"   Has 'U Fe-': {has_u_parameter} (should be TRUE ✅)")
    
    print("\n   SYSTEM namelist:")
    print("   " + "-" * 76)
    in_system = False
    for line in content.split('\n'):
        if '&SYSTEM' in line:
            in_system = True
        if in_system:
            print(f"   {line}")
        if '/' in line and in_system and line.strip() == '/':
            break
    
    print("\n   HUBBARD card:")
    print("   " + "-" * 76)
    in_hubbard = False
    for line in content.split('\n'):
        if 'HUBBARD' in line:
            in_hubbard = True
        if in_hubbard:
            print(f"   {line}")
            if line.strip() and not line.startswith(' ') and 'HUBBARD' not in line:
                break

print("\n" + "="*80)
print("RESULT")
print("="*80)

if has_hubbard_u or has_lda_plus_u:
    print("""
❌ PROBLEM: Old format parameters (hubbard_u or lda_plus_u) still present!
   QE 7.4.1 will reject this input.
""")
elif has_hubbard_card and has_u_parameter:
    print("""
✅ FIXED! 
   - Old format parameters removed ✅
   - HUBBARD card with new format present ✅
   - QE 7.4.1 will accept this input ✅
""")
else:
    print("""
❌ INCOMPLETE: Missing HUBBARD card!
""")
