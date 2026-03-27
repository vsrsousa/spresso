"""
Test: Does passing {'Fe': 4.3} with code_version='7.2' generate HUBBARD card?
"""

from ase.build import bulk
from xespresso.workflow.calculation_workflow import CalculationWorkflow
from xespresso.xio import write_espresso_in
import tempfile
import os

print("="*80)
print("TEST: Simple dict {'Fe': 4.3} with code_version='7.2'")
print("="*80)

atoms = bulk("Fe", cubic=True)

with tempfile.TemporaryDirectory() as tmpdir:
    try:
        # This is what user wants to pass
        workflow = CalculationWorkflow(
            atoms=atoms,
            code_version='7.2',  # NEW QE format
            pseudopotentials={'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'},
            precision='low',
            hubbard_config={'Fe': 4.3},  # Simple dict, no 'u' key
            outdir=tmpdir,
        )
        
        print("✅ CalculationWorkflow created")
        print(f"   code_version: 7.2")
        print(f"   hubbard_config: {{'Fe': 4.3}}")
        
        # Check what's in input_data
        print(f"\n   input_data['qe_version']: {workflow.input_data.get('qe_version')}")
        print(f"   input_data['hubbard']: {workflow.input_data.get('hubbard', 'NOT SET')}")
        print(f"   input_data['lda_plus_u']: {workflow.input_data.get('lda_plus_u', 'NOT SET')}")
        if 'input_ntyp' in workflow.input_data:
            print(f"   input_data['input_ntyp']['Hubbard_U']: {workflow.input_data['input_ntyp'].get('Hubbard_U', {})}")
        
        # Write input file
        input_file = os.path.join(tmpdir, 'test.pwi')
        write_espresso_in(
            input_file,
            atoms,
            input_data=workflow.input_data,
            pseudopotentials=workflow.pseudopotentials,
            kpts=(2, 2, 2)
        )
        
        print("\n📄 Generated Input File:")
        print("-" * 80)
        
        with open(input_file, 'r') as f:
            content = f.read()
        
        # Show SYSTEM namelist and last part of file
        lines = content.split('\n')
        
        print("\n   SYSTEM namelist:")
        in_system = False
        for i, line in enumerate(lines, 1):
            if '&SYSTEM' in line:
                in_system = True
            if in_system:
                print(f"   Line {i:3d}: {line}")
            if '/' in line and in_system:
                break
        
        print("\n   Last 20 lines (checking for HUBBARD card):")
        for i, line in enumerate(lines[-20:], len(lines)-19):
            print(f"   Line {i:3d}: {line}")
        
        # Verdict
        if 'HUBBARD' in content:
            print("\n✅ HUBBARD card FOUND - YES IT WORKS!")
        elif 'hubbard_u(1)' in content.lower():
            print("\n⚠️  Old format Hubbard_U FOUND - NOT new format!")
            print("   This is WRONG for QE 7.2!")
        else:
            print("\n❌ No Hubbard parameters found!")
            
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

print("\n" + "="*80)
print("ANSWER")
print("="*80)
print("""
If HUBBARD card appears:
✅ YES, it works - user can pass {'Fe': 4.3} with code_version='7.2'

If old format appears (hubbard_u(1)):
❌ NO, it doesn't work - treated as old format despite code_version='7.2'

If nothing appears:
❌ NO, parameters are lost
""")
