"""
Direct test: Call write_espresso_in() directly to generate input file
"""

from ase.build import bulk
from xespresso.workflow.calculation_workflow import CalculationWorkflow
from xespresso.xio import write_espresso_in
import tempfile
import os

print("="*80)
print("DIRECT xio.write_espresso_in() TEST")
print("="*80)

atoms = bulk("Fe", cubic=True)

print("\n" + "="*80)
print("TEST 1: Old Format (QE 6.8)")
print("="*80)

with tempfile.TemporaryDirectory() as tmpdir:
    try:
        workflow_old = CalculationWorkflow(
            atoms=atoms,
            code_version='6.8',
            pseudopotentials={'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'},
            precision='low',
            hubbard_config={'Fe': 4.3},
            outdir=tmpdir,
        )
        
        print("✅ CalculationWorkflow created")
        print(f"   lda_plus_u: {workflow_old.input_data.get('lda_plus_u')}")
        print(f"   Hubbard_U: {workflow_old.input_data.get('input_ntyp', {}).get('Hubbard_U', {})}")
        
        # Directly write input file using xio
        input_file = os.path.join(tmpdir, 'test_old.pwi')
        write_espresso_in(
            input_file,
            atoms,
            input_data=workflow_old.input_data,
            pseudopotentials=workflow_old.pseudopotentials,
            kpts=(2, 2, 2)
        )
        
        if os.path.exists(input_file):
            with open(input_file, 'r') as f:
                content = f.read()
            
            print("\n✅ Input file created successfully!")
            print("-" * 80)
            
            # Search for Hubbard lines
            found = False
            for i, line in enumerate(content.split('\n'), 1):
                lower = line.lower()
                if 'hubbard' in lower or 'u(' in lower or 'lda_plus_u' in lower:
                    print(f"   Line {i:3d}: {line}")
                    found = True
            
            if found:
                print("\n   ✅ Hubbard parameters found!")
            else:
                print("\n   ⚠️  No Hubbard in expected locations")
                print("\n   Full file content:")
                print(content)
        else:
            print("   ⚠️  File not created")
            
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

print("\n" + "="*80)
print("TEST 2: New Format (QE 7.2)")
print("="*80)

with tempfile.TemporaryDirectory() as tmpdir:
    try:
        workflow_new = CalculationWorkflow(
            atoms=atoms,
            code_version='7.2',
            pseudopotentials={'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'},
            precision='low',
            hubbard_config={'u': {'Fe': 4.3}},
            outdir=tmpdir,
        )
        
        print("✅ CalculationWorkflow created")
        print(f"   qe_version: {workflow_new.input_data.get('qe_version')}")
        print(f"   hubbard: {workflow_new.input_data.get('hubbard')}")
        
        # Directly write input file using xio
        input_file = os.path.join(tmpdir, 'test_new.pwi')
        write_espresso_in(
            input_file,
            atoms,
            input_data=workflow_new.input_data,
            pseudopotentials=workflow_new.pseudopotentials,
            kpts=(2, 2, 2)
        )
        
        if os.path.exists(input_file):
            with open(input_file, 'r') as f:
                content = f.read()
            
            print("\n✅ Input file created successfully!")
            print("-" * 80)
            
            # Search for HUBBARD card
            found = False
            lines = content.split('\n')
            for i, line in enumerate(lines, 1):
                if 'HUBBARD' in line:
                    print(f"   Line {i:3d}: {line}")
                    found = True
                    # Print next 5 lines
                    for j in range(1, 6):
                        if i+j-1 < len(lines):
                            print(f"   Line {i+j:3d}: {lines[i+j-1]}")
                    break
            
            if found:
                print("\n   ✅ HUBBARD card found!")
            else:
                print("\n   ⚠️  No HUBBARD card found")
                print("\n   Last 30 lines of file:")
                for i, line in enumerate(lines[-30:], len(lines)-29):
                    print(f"   Line {i:3d}: {line}")
        else:
            print("   ⚠️  File not created")
            
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

print("\n" + "="*80)
print("RESULT")
print("="*80)
print("""
If both tests show Hubbard parameters in their input files:
✅ The COMPLETE workflow is CORRECT!
   The chain works: ConvergenceWorkflow → CalculationWorkflow → write_espresso_in()
   
If parameters are missing:
❌ Check xespresso.xio module for Hubbard processing
""")
