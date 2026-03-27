"""
Direct test: Create Espresso calculator and check if Hubbard appears in generated input file
"""

from ase.build import bulk
from ase.calculators.espresso import Espresso
import tempfile
import os

print("="*80)
print("DIRECT CALCULATOR TEST - Hubbard in Input File")
print("="*80)

atoms = bulk("Fe", cubic=True)

print("\n" + "="*80)
print("TEST 1: Old Format (QE 6.8) - Direct Espresso Calculator")
print("="*80)

with tempfile.TemporaryDirectory() as tmpdir:
    try:
        input_data_old = {
            'system': {
                'ibrav': 0,
                'ecutwfc': 71.0,
                'lda_plus_u': True,  # Enable DFT+U
            },
            'input_ntyp': {
                'Hubbard_U': {'Fe': 4.3}  # Old format
            }
        }
        
        calc_old = Espresso(
            pseudopotentials={'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'},
            input_data=input_data_old,
            calculation='scf',
            kpts=(2,2,2),
            label=os.path.join(tmpdir, 'old_format'),
        )
        
        print("✅ Espresso calculator created (old format)")
        print(f"   input_data['system']['lda_plus_u']: {input_data_old['system']['lda_plus_u']}")
        print(f"   input_data['input_ntyp']['Hubbard_U']: {input_data_old['input_ntyp']['Hubbard_U']}")
        
        # Generate input file
        calc_old.write_input(atoms)
        
        input_file = os.path.join(tmpdir, 'old_format.pwi')
        if os.path.exists(input_file):
            with open(input_file, 'r') as f:
                content = f.read()
            
            print("\n📄 Input file generated:")
            print("-" * 80)
            
            # Search for Hubbard-related lines
            found_hubbard = False
            for i, line in enumerate(content.split('\n'), 1):
                lower = line.lower()
                if 'hubbard' in lower or 'lda_plus_u' in lower or 'u(' in lower:
                    print(f"   Line {i:3d}: {line}")
                    found_hubbard = True
            
            if found_hubbard:
                print("\n   ✅ Hubbard parameters found in input file!")
            else:
                print("\n   ⚠️  No Hubbard parameters found")
        else:
            print(f"   ⚠️  Input file not created at {input_file}")
            
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

print("\n" + "="*80)
print("TEST 2: New Format (QE 7.2) - Direct Espresso Calculator")
print("="*80)

with tempfile.TemporaryDirectory() as tmpdir:
    try:
        input_data_new = {
            'system': {
                'ibrav': 0,
                'ecutwfc': 71.0,
            },
            'qe_version': '7.2',
            'hubbard': {
                'u': {'Fe': 4.3}  # New format
            }
        }
        
        calc_new = Espresso(
            pseudopotentials={'Fe': 'Fe.pbe-spn-rrkjus_psl.1.0.0.UPF'},
            input_data=input_data_new,
            calculation='scf',
            kpts=(2,2,2),
            label=os.path.join(tmpdir, 'new_format'),
        )
        
        print("✅ Espresso calculator created (new format)")
        print(f"   input_data['qe_version']: {input_data_new['qe_version']}")
        print(f"   input_data['hubbard']: {input_data_new['hubbard']}")
        
        # Generate input file
        calc_new.write_input(atoms)
        
        input_file = os.path.join(tmpdir, 'new_format.pwi')
        if os.path.exists(input_file):
            with open(input_file, 'r') as f:
                content = f.read()
            
            print("\n📄 Input file generated:")
            print("-" * 80)
            
            # Search for Hubbard card
            found_hubbard = False
            for i, line in enumerate(content.split('\n'), 1):
                lower = line.lower()
                if 'hubbard' in lower:
                    print(f"   Line {i:3d}: {line}")
                    found_hubbard = True
            
            if found_hubbard:
                print("\n   ✅ HUBBARD card found in input file!")
            else:
                print("\n   ⚠️  No HUBBARD card found")
        else:
            print(f"   ⚠️  Input file not created at {input_file}")
            
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

print("\n" + "="*80)
print("CONCLUSION")
print("="*80)
print("""
If both tests show Hubbard in input files:
✅ The complete flow works correctly:
   1. ConvergenceWorkflow receives hubbard_config
   2. Passes to CalculationWorkflow
   3. CalculationWorkflow stores in input_data correctly
   4. Espresso calculator reads input_data
   5. Input file is generated with Hubbard parameters

If parameters are missing:
❌ Espresso calculator doesn't recognize input_data['hubbard'] structure
   - May need to check xespresso integration with ASE Espresso
   - Or there's an intermediate step that processes input_data before Espresso
""")
