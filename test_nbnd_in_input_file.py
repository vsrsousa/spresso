#!/usr/bin/env python3
"""
Test that nbnd appears in the &SYSTEM section of the input file when enhance_nbands=True.
"""

from ase.build import bulk
import tempfile
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)-8s [%(name)s]: %(message)s')
logger = logging.getLogger(__name__)

from xespresso.workflow import CalculationWorkflow

# Create dummy pseudopotentials
temp_dir = tempfile.mkdtemp()
pseudo_dir = Path(temp_dir) / "pseudos"
pseudo_dir.mkdir()

si_pseudo = pseudo_dir / "Si.pbe.UPF"
si_pseudo.write_text("""created by Si pseudopotential
z_valence = "4"
     4.00 (valence charge)
""", )

atoms = bulk('Si', 'diamond', a=5.43)
pseudos = {'Si': str(si_pseudo)}

print("="*70)
print("TEST: nbnd appears in &SYSTEM section of input file")
print("="*70)

# Create workflow with enhance_nbands=True
wf = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials=pseudos,
    protocol='fast',
    enhance_nbands=True
)

# Create a test directory
test_dir = Path(temp_dir) / "test_scf"
test_dir.mkdir()

# Dry run to generate input file without executing
print("\n→ Running SCF dry_run with enhance_nbands=True...")
calc = wf.run_scf(label='scf', dry_run=True)

# Read the input file and check for nbnd
input_file = Path('scf') / 'scf.pwi'
if input_file.exists():
    print(f"✓ Input file created: {input_file}")
    
    with open(input_file, 'r') as f:
        content = f.read()
    
    print("\nInput file &SYSTEM section:")
    print("-" * 70)
    
    # Extract &SYSTEM section
    if "&SYSTEM" in content:
        start = content.find("&SYSTEM")
        end = content.find("/", start)
        system_section = content[start:end+1]
        print(system_section)
        print("-" * 70)
        
        # Check if nbnd is present
        if "nbnd" in system_section:
            # Extract the nbnd value
            import re
            match = re.search(r'nbnd\s*=\s*(\d+)', system_section)
            if match:
                nbnd_value = int(match.group(1))
                expected_nbnd = len(atoms) * 4  # 2 Si atoms × 4 valence electrons
                
                print(f"\n✓ SUCCESS: nbnd found in &SYSTEM section!")
                print(f"  nbnd = {nbnd_value}")
                print(f"  Expected: {expected_nbnd} (2 atoms × 4 valence electrons)")
                
                if nbnd_value == expected_nbnd:
                    print(f"  ✓ PASS: Values match!")
                else:
                    print(f"  ✗ FAIL: Values don't match (got {nbnd_value}, expected {expected_nbnd})")
        else:
            print("\n✗ FAIL: nbnd NOT found in &SYSTEM section!")
    else:
        print("✗ ERROR: &SYSTEM section not found in input file")
else:
    print(f"✗ ERROR: Input file not created at {input_file}")

# Cleanup
import shutil
shutil.rmtree(temp_dir)

print("\n" + "="*70)
print("TEST COMPLETE")
print("="*70)
