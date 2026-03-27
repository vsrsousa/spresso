"""
Mock test: Verify Hubbard appears CORRECTLY in the actual generated input file
"""

from ase.build import bulk
from xespresso.workflow.calculation_workflow import CalculationWorkflow
import tempfile
import os

print("="*80)
print("INPUT FILE GENERATION TEST - Hubbard Parameters")
print("="*80)

# Create Fe structure
atoms = bulk("Fe", cubic=True)

print("\n" + "="*80)
print("TEST 1: Old Format (QE 6.8) - Check input file content")
print("="*80)

try:
    with tempfile.TemporaryDirectory() as tmpdir:
        calc_old = CalculationWorkflow(
            atoms=atoms,
            code_version='6.8',
            pseudopotentials={'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'},
            precision='low',
            hubbard_config={'Fe': 4.3},
            outdir=tmpdir,
        )
        
        print("✅ CalculationWorkflow (QE 6.8) created")
        print(f"   input_data['lda_plus_u']: {calc_old.input_data.get('lda_plus_u')}")
        print(f"   input_data['input_ntyp']['Hubbard_U']: {calc_old.input_data.get('input_ntyp', {}).get('Hubbard_U', {})}")
        
        # Try to generate input file
        try:
            # Check internal qe module to see if it writes correctly
            from xespresso.xio.qeinput import qeinput
            
            input_text = qeinput.qeinput(atoms, calc_old.input_data, calc_old.pseudopotentials)
            
            print("\n📄 Generated Input File Content:")
            print("-" * 80)
            
            # Show relevant parts
            for i, line in enumerate(input_text.split('\n')):
                if 'lda_plus_u' in line.lower() or 'hubbard' in line.lower() or 'u(' in line.lower():
                    print(f"   Line {i:3d}: {line}")
            
            if 'Hubbard_U' in input_text or 'U(' in input_text:
                print("\n   ✅ HUBBARD parameters found in input file")
            else:
                print("\n   ⚠️  No Hubbard parameters in input file")
                
        except Exception as e:
            print(f"   ℹ️  Could not generate input file: {type(e).__name__}: {e}")
        
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*80)
print("TEST 2: New Format (QE 7.2) - Check input file content")
print("="*80)

try:
    with tempfile.TemporaryDirectory() as tmpdir:
        calc_new = CalculationWorkflow(
            atoms=atoms,
            code_version='7.2',
            pseudopotentials={'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'},
            precision='low',
            hubbard_config={'u': {'Fe': 4.3}},
            outdir=tmpdir,
        )
        
        print("✅ CalculationWorkflow (QE 7.2) created")
        print(f"   input_data['qe_version']: {calc_new.input_data.get('qe_version')}")
        print(f"   input_data['hubbard']: {calc_new.input_data.get('hubbard')}")
        
        # Try to generate input file
        try:
            from xespresso.xio.qeinput import qeinput
            
            input_text = qeinput.qeinput(atoms, calc_new.input_data, calc_new.pseudopotentials)
            
            print("\n📄 Generated Input File Content:")
            print("-" * 80)
            
            # Show relevant parts
            for i, line in enumerate(input_text.split('\n')):
                if 'hubbard' in line.lower() or 'lda_plus_u' in line.lower():
                    print(f"   Line {i:3d}: {line}")
            
            if 'HUBBARD' in input_text:
                print("\n   ✅ HUBBARD card found in input file")
            else:
                print("\n   ⚠️  No HUBBARD card in input file")
                
        except Exception as e:
            print(f"   ℹ️  Could not generate input file: {type(e).__name__}: {e}")
        
except Exception as e:
    print(f"❌ Error: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "="*80)
print("CONCLUSION")
print("="*80)
print("""
If both tests show Hubbard parameters in the generated input files:
✅ The flow is CORRECT end-to-end
   - ConvergenceWorkflow passes hubbard_config
   - CalculationWorkflow processes it correctly
   - Input file is generated with Hubbard in correct format

If parameters are missing:
❌ Hubbard module (qeinput) might not recognize input_data['hubbard']
   - Check xespresso/xio/qeinput.py for Hubbard handling
""")
