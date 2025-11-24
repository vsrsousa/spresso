"""
Example demonstrating the simplified workflow for Quantum ESPRESSO calculations.

This example shows how to:
1. Run SCF calculations with quality presets
2. Run relaxation calculations
3. Use k-spacing instead of explicit k-points
4. Work with CIF files
"""

from ase.build import bulk
from xespresso import CalculationWorkflow, quick_scf, quick_relax, PRESETS
import numpy as np

# Example 1: Simple SCF calculation with quality preset
print("="*60)
print("Example 1: Simple SCF calculation with quality preset")
print("="*60)

atoms = bulk("Si", cubic=True)
pseudopotentials = {"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"}

# Create a workflow with 'fast' quality preset
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials=pseudopotentials,
    protocol='fast'
)

print("\nPreset information:")
for key, value in workflow.get_preset_info().items():
    print(f"  {key}: {value}")

# Run SCF calculation (commented out for demonstration)
# calc = workflow.run_scf(label='scf/silicon-fast')
# print(f"Energy: {calc.results['energy']} eV")


# Example 2: SCF calculation with custom k-spacing
print("\n" + "="*60)
print("Example 2: SCF calculation with custom k-spacing")
print("="*60)

# Using k-spacing like in the problem statement
# from ase.io.espresso import kspacing_to_grid
# kpts = kspacing_to_grid(atoms, 0.20/(2*np.pi))

# With our workflow, you can specify k-spacing directly
workflow = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials=pseudopotentials,
    protocol='moderate',
    kspacing=0.3  # Angstrom^-1
)

print(f"K-spacing: {workflow.kspacing} Angstrom^-1")
print(f"Resulting k-points: {workflow._get_kpts()}")

# Run SCF (commented out)
# calc = workflow.run_scf(label='scf/silicon-kspacing')


# Example 3: Structure relaxation with different quality levels
print("\n" + "="*60)
print("Example 3: Structure relaxation with quality levels")
print("="*60)

# Fast relaxation for initial guess
workflow_fast = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials=pseudopotentials,
    protocol='fast'
)

# Accurate relaxation for final structure
workflow_accurate = CalculationWorkflow(
    atoms=atoms,
    pseudopotentials=pseudopotentials,
    protocol='accurate'
)

print("\nFast preset:")
for key, value in workflow_fast.get_preset_info().items():
    print(f"  {key}: {value}")

print("\nAccurate preset:")
for key, value in workflow_accurate.get_preset_info().items():
    print(f"  {key}: {value}")

# Run relaxation (commented out)
# calc = workflow_fast.run_relax(label='relax/silicon-fast')
# calc = workflow_accurate.run_relax(label='relax/silicon-accurate')


# Example 4: Using quick helper functions
print("\n" + "="*60)
print("Example 4: Using quick helper functions")
print("="*60)

# Quick SCF calculation
# calc = quick_scf(
#     atoms,
#     pseudopotentials,
#     label='scf/silicon-quick',
#     protocol='moderate'
# )

# Quick relaxation
# calc = quick_relax(
#     atoms,
#     pseudopotentials,
#     label='relax/silicon-quick',
#     protocol='moderate',
#     relax_type='vc-relax'  # relax cell and ions
# )

print("Quick functions allow one-line calculations")


# Example 5: Working with CIF files
print("\n" + "="*60)
print("Example 5: Working with CIF files")
print("="*60)

# If you have a CIF file, you can create a workflow directly from it
# workflow = CalculationWorkflow.from_cif(
#     'structure.cif',
#     pseudopotentials={'Fe': 'Fe.pbe-spn.UPF'},
#     protocol='moderate'
# )
# calc = workflow.run_scf(label='scf/from-cif')

# Or use the quick function
# calc = quick_scf(
#     'structure.cif',
#     {'Fe': 'Fe.pbe-spn.UPF'},
#     protocol='fast'
# )

print("Workflows can be created directly from CIF files")


# Example 6: View available presets
print("\n" + "="*60)
print("Example 6: Available quality presets")
print("="*60)

for quality_name, settings in PRESETS.items():
    print(f"\n{quality_name.upper()}:")
    for key, value in settings.items():
        print(f"  {key}: {value}")


print("\n" + "="*60)
print("Examples complete!")
print("="*60)
