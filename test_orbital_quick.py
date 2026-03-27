"""
Quick test: Fe and Gd orbital auto-detection
"""

from ase.build import bulk
from xespresso.workflow.calculation_workflow import CalculationWorkflow
from xespresso.xio import write_espresso_in
import tempfile
import os

print("="*80)
print("ORBITAL AUTO-DETECTION: Fe and Gd")
print("="*80)

# Test Fe (3d)
print("\n TEST 1: Fe → should be Fe-3d")
print("-"*80)

with tempfile.TemporaryDirectory() as tmpdir:
    atoms = bulk("Fe", cubic=True)
    
    workflow = CalculationWorkflow(
        atoms=atoms,
        code_version='7.2',
        pseudopotentials={'Fe': 'Fe.pbe.UPF'},
        precision='low',
        hubbard_config={'Fe': 4.3},
        outdir=tmpdir,
    )
    
    input_file = os.path.join(tmpdir, 'test_fe.pwi')
    write_espresso_in(input_file, atoms, input_data=workflow.input_data, 
                     pseudopotentials=workflow.pseudopotentials, kpts=(2,2,2))
    
    with open(input_file, 'r') as f:
        for line in f:
            if 'HUBBARD' in line or ('U ' in line and 'Fe' in line):
                print(f"  {line.rstrip()}")

# Test Gd (4f lanthanide)
print("\n TEST 2: Gd → should be Gd-4f (lanthanide)")
print("-"*80)

with tempfile.TemporaryDirectory() as tmpdir:
    atoms = bulk("Gd", cubic=True)
    
    workflow = CalculationWorkflow(
        atoms=atoms,
        code_version='7.2',
        pseudopotentials={'Gd': 'Gd.pbe.UPF'},
        precision='low',
        hubbard_config={'Gd': 6.0},
        outdir=tmpdir,
    )
    
    input_file = os.path.join(tmpdir, 'test_gd.pwi')
    write_espresso_in(input_file, atoms, input_data=workflow.input_data,
                     pseudopotentials=workflow.pseudopotentials, kpts=(2,2,2))
    
    with open(input_file, 'r') as f:
        for line in f:
            if 'HUBBARD' in line or ('U ' in line and 'Gd' in line):
                print(f"  {line.rstrip()}")

print("\n" + "="*80)
print("✅ Auto-detection should work correctly!")
print("="*80)
