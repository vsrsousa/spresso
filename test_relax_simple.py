#!/usr/bin/env python3
"""
Simple test of CalculationWorkflow.relax() for non-magnetic system
"""

from ase.build import molecule
from xespresso import CalculationWorkflow
import tempfile
import os

# Create a simple non-magnetic system: H2 molecule
atoms = molecule('H2', vacuum=5.0)

print("=" * 60)
print("Testing CalculationWorkflow.relax() - Non-magnetic H2")
print("=" * 60)

# Create workflow with fast protocol
wf = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials={"H": "H.pbe-kjpaw_psl.1.0.0.UPF"},
    protocol="fast",  # Use fast for quick test
)

print(f"\nInitial atomic positions:\n{atoms.positions}")
print(f"Initial cell:\n{atoms.cell}")

# Run relax with constraints
with tempfile.TemporaryDirectory() as tmpdir:
    label = os.path.join(tmpdir, "h2_relax")
    
    # Write input files
    wf.write_input(label=label)
    
    print(f"\n✓ Input files written to {label}/")
    print(f"  - {os.path.basename(label)}.pwi (input file)")
    print(f"  - {os.path.basename(label)}.asei (structure)")
    
    # Check if occupations is set correctly
    pwi_file = os.path.join(tmpdir, os.path.basename(label), f"{os.path.basename(label)}.pwi")
    if os.path.exists(pwi_file):
        with open(pwi_file, 'r') as f:
            content = f.read()
            # Check key parameters
            params = {
                'ecutwfc': [l for l in content.split('\n') if 'ecutwfc' in l.lower()],
                'occupations': [l for l in content.split('\n') if 'occupations' in l.lower() and '=' in l],
                'smearing': [l for l in content.split('\n') if 'smearing' in l.lower() and '=' in l],
                'degauss': [l for l in content.split('\n') if 'degauss' in l.lower() and '=' in l],
            }
            
            print(f"\n✓ Key parameters in PWI file:")
            for param, lines in params.items():
                if lines:
                    print(f"  {param}: {lines[0].strip()}")
                else:
                    print(f"  {param}: NOT FOUND ⚠️")

print("\n" + "=" * 60)
print("Setup complete! Now you can:")
print("  1. Modify the script to add remote execution")
print("  2. Check that structure relaxes (positions change)")
print("  3. Verify convergence with different protocols")
print("=" * 60)
